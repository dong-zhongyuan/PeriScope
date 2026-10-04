"""CCWM — potential-form equilibrium simulator (pd_product track, ADR-0001 APPROVED).

Derived by edit from ad_omics_v2 unified_model_v3.py（来源基底目录已按 2026-09-18 用户指令删除，
派生关系以 git 历史与 ADR-0002 为准）: the free-form neural ODE (DriftNet) is
replaced by a potential gradient flow dz/dt = -grad U (structurally
non-divergent, dU/dt <= 0), and the disease-label random point pairing of the
AD engine is replaced by block-structured supervision:

- Block B (CITE-seq, same-cell true pairing): weighted denoising score matching
  on z_blood via -grad U (the equilibrium score) + protein head regression;
- Block A (PD brain+blood, disease-label only): blood coordinates fixed, brain
  coordinates minimized to equilibrium, decoded brain RNA matched to observed
  brain RNA at DISTRIBUTION level with stratified Sinkhorn (cond2 x compartment).

Claim semantics (docs/CCWM_PREREGISTRATION.md): model-space perturbation
response only; u is CITE-seq ADT, NOT plasma protein; causal grade
NOT-EVALUABLE (no held-out action data).
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F
from torch import nn


class MLP(nn.Module):
    def __init__(self, din: int, dout: int, hidden: int = 256, depth: int = 2, dropout: float = 0.1) -> None:
        super().__init__()
        layers: list[nn.Module] = []
        d = din
        for _ in range(depth):
            layers += [nn.Linear(d, hidden), nn.ReLU(inplace=True), nn.Dropout(dropout)]
            d = hidden
        layers.append(nn.Linear(d, dout))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class ProteinHead(MLP):
    """The measured ADT CLR vector lies in the zero-sum subspace."""
    def forward(self,x):
        y=super().forward(x)
        return y-y.mean(dim=-1,keepdim=True)


class Encoder(nn.Module):
    """RNA -> (mu, logvar) of compartment latent z."""

    def __init__(self, din: int, z_dim: int, hidden: int = 256, depth: int = 2) -> None:
        super().__init__()
        self.body = MLP(din, 2 * z_dim, hidden, depth)

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        out = self.body(x)
        mu, logvar = out.chunk(2, dim=-1)
        return mu, logvar.clamp(-8.0, 8.0)


def reparameterize(mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
    return mu + torch.exp(0.5 * logvar) * torch.randn_like(mu)


FREE_BITS = 0.05  # per-dim KL floor (nats) — v2.1


def kl_normal_free_bits(mu: torch.Tensor, logvar: torch.Tensor, lam: float = FREE_BITS) -> torch.Tensor:
    """Per-dimension KL floor in the penalty; does not guarantee latent usage."""
    kl_dim = 0.5 * (mu.pow(2) + logvar.exp() - 1.0 - logvar)  # (B, z_dim)
    return torch.clamp(kl_dim, min=lam).sum(dim=-1).mean()

def kl_normal_legacy(mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
    return 0.5 * (mu.pow(2) + logvar.exp() - 1.0 - logvar).sum(dim=-1).mean()


def ot_cost(x, y, eps=0.1, iters=50, scale=1.0):
    """Uniform entropic OT dual using one externally fitted cost scale."""
    c = torch.cdist(x, y).square() / scale
    n, m = c.shape
    f, g = c.new_zeros(n), c.new_zeros(m)
    loga, logb = -__import__('math').log(n), -__import__('math').log(m)
    for _ in range(iters):
        f = eps * loga - eps * torch.logsumexp((g[None,:]-c)/eps,1)
        g = eps * logb - eps * torch.logsumexp((f[:,None]-c)/eps,0)
    return f.mean() + g.mean()


def sinkhorn_divergence(x, y, eps=0.1, iters=50, scale=1.0):
    return ot_cost(x,y,eps,iters,scale) - .5*ot_cost(x,x,eps,iters,scale) - .5*ot_cost(y,y,eps,iters,scale)


@dataclass
class CCWMConfig:
    n_genes: int
    n_proteins: int
    cond_dim: int = 8          # onehot cond2(2) + compartment(4) + source(2)
    z_dim: int = 64            # per compartment; total latent state 2*z_dim = 128
    hidden: int = 256
    depth: int = 2
    k_steps: int = 4
    eta: float = 0.5
    alpha_l2: float = 0.01
    coupling: bool = True
    sigma_rel: tuple = (0.1, 0.3, 0.5)
    cond_dropout: float = 0.1
    sinkhorn_eps: float = 0.1
    sinkhorn_iters: int = 50
    beta_recon: float = 0.1
    beta_kl: float = 0.01
    beta_protein: float = 1.0
    beta_dsm: float = 1.0
    beta_sinkhorn: float = 1.0
    # v2.0 engine ports from ad_omics_v2 (prereg par.10)
    zero_weight: float = 20.0   # Zero-Masked MSE nonzero-gene weight (anti-collapse core)
    free_bits: float = 0.05      # per-dim KL floor (v2.1)
    kl_anneal_epochs: int = 50  # beta-annealing window (v2.1)
    esm_dim: int = 640
    potential_beta: float = 5.0
    ot_scale: float = 1.0
    sample_latent_a: bool = False


class CCWM(nn.Module):
    """dz/dt = -grad U; equilibrium z* read out by decoding; intervention =
    clamping u and re-solving the equilibrium (model-space only)."""

    def __init__(self, cfg: CCWMConfig) -> None:
        super().__init__()
        self.cfg = cfg
        z, g, p, cd = cfg.z_dim, cfg.n_genes, cfg.n_proteins, cfg.cond_dim
        # v2.2: ESM2 projection REMOVED (blueprint hard-constraint 2:
        # "跨组学不用 ESM2：RNA 直入、蛋白直入"). Encoder/decoder operate
        # directly in gene space. Anti-collapse via Zero-Masked MSE on
        # center-only data (which preserves sparse expression structure).
        enc_in = g
        dec_out_dim = g
        self.enc_blood = Encoder(enc_in, z, cfg.hidden, cfg.depth)
        self.enc_brain = Encoder(enc_in, z, cfg.hidden, cfg.depth)
        self.dec_blood = MLP(z, dec_out_dim, cfg.hidden, cfg.depth)
        self.dec_brain = MLP(z, dec_out_dim, cfg.hidden, cfg.depth)
        # v2.2: data arrives pre-centered (center_only at build); gene_mean = 0
        # (kept as buffer so zm_mse/centered_cos work unchanged; subtracting 0 is a no-op)
        self.register_buffer("gene_mean", torch.zeros(g))
        self.protein_head = ProteinHead(z, p, cfg.hidden, cfg.depth)
        self.u_blood = MLP(z + p + cd, 1, cfg.hidden, cfg.depth)
        self.u_brain = MLP(z + p + cd, 1, cfg.hidden, cfg.depth)
        self.u_coupling = MLP(2 * z + p + cd, 1, cfg.hidden, cfg.depth) if cfg.coupling else None
        # Mixed latent/protein derivatives must be nonzero during gradient-flow training.
        # Deterministic potentials make line-search comparisons meaningful.
        for potential in [self.u_blood, self.u_brain, self.u_coupling]:
            if potential is not None:
                for i,layer in enumerate(potential.net):
                    if isinstance(layer, nn.ReLU): potential.net[i] = nn.Softplus(beta=cfg.potential_beta)
                    elif isinstance(layer, nn.Dropout): potential.net[i] = nn.Identity()
        self.register_buffer('blood_center', torch.zeros(g))
        self.register_buffer('brain_center', torch.zeros(g))
        # v1.1: per-step learnable descent scale (softplus-parameterized; softplus(-2.2)
        # ~= 0.10). v1.0's fixed eta=0.5 overshot relative to ||z|| (~1.2) and the energy
        # ROSE during the solve (results/index/CCWM_V1_DIAGNOSTICS.md).
        self.step_sizes = nn.Parameter(torch.full((cfg.k_steps,), -2.2))

    # ---- v2.0 helpers: ESM2 projection, centering, zero-masked losses ----
    def zm_mse(self, pred_c: torch.Tensor, target_c: torch.Tensor, raw_target: torch.Tensor) -> torch.Tensor:
        """Zero-Masked MSE on CENTERED targets; mask from RAW expression (ad engine)."""
        nonzero = (raw_target.abs() > 1e-6).float()
        w = nonzero * self.cfg.zero_weight + (1.0 - nonzero)
        return ((pred_c - target_c) ** 2 * w).mean()

    def centered_cos(self, pred_c: torch.Tensor, target_c: torch.Tensor) -> torch.Tensor:
        return 1.0 - F.cosine_similarity(pred_c, target_c, dim=-1).mean()

    # ---- potential & equilibrium ----
    def potential(self, zb: torch.Tensor, zbr: torch.Tensor, u: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """U(z_b, z_B, u, c) = U_blood + U_brain + U_coupling + (alpha/2)||z||^2; (B,) energies."""
        e = self.u_blood(torch.cat([zb, u, c], dim=-1)).squeeze(-1)
        e = e + self.u_brain(torch.cat([zbr, u, c], dim=-1)).squeeze(-1)
        if self.u_coupling is not None:
            e = e + self.u_coupling(torch.cat([zb, zbr, u, c], dim=-1)).squeeze(-1)
        return e + 0.5 * self.cfg.alpha_l2 * (zb.pow(2).sum(dim=-1) + zbr.pow(2).sum(dim=-1))

    def solve_equilibrium(
        self, zb: torch.Tensor, zbr0: torch.Tensor, u: torch.Tensor, c: torch.Tensor,
        return_traj: bool = False,
    ):
        """K unrolled descent steps on U over brain coordinates with blood
        coordinates FIXED (Block A semantics), per-step learnable scale and
        Armijo backtracking (halve t until the energy decreases; v1.1 fix).
        Gradients flow through the ACCEPTED step (create_graph), so training
        still backpropagates through the unrolled solve. traj (optional) records
        (energy_before, energy_after, accepted_t) per step for gate-0 auditing."""
        zb = zb if zb.requires_grad else zb.detach().requires_grad_(True)
        zbr = zbr0.clone()
        if not zbr.requires_grad: zbr.requires_grad_(True)
        traj = []
        for k in range(self.cfg.k_steps):
            e = self.potential(zb, zbr, u, c)
            g = torch.autograd.grad(e.sum(), zbr, create_graph=True)[0]
            # learnable scale keeps its graph; backtracking multiplies a 2^-j guard
            t_grad = F.softplus(self.step_sizes[k])
            with torch.no_grad():
                e0 = e.detach()
                g2 = g.detach().square().sum(-1)
                guard = torch.ones_like(e0)
                for _ in range(8):
                    trial = self.potential(zb,zbr-(guard*float(t_grad))[:,None]*g,u,c)
                    ok = trial <= e0 - 1e-4*guard*float(t_grad)*g2 + 1e-8
                    if bool(ok.all()): break
                    guard = torch.where(ok,guard,guard*.5)
                trial = self.potential(zb,zbr-(guard*float(t_grad))[:,None]*g,u,c)
                guard = torch.where(trial<=e0+1e-8,guard,torch.zeros_like(guard))
            zbr = zbr - (guard*t_grad)[:,None]*g
            if return_traj:
                with torch.no_grad():
                    traj.append((float(e0.mean()),float(self.potential(zb,zbr,u,c).mean()),float((guard*t_grad).mean())))
        if return_traj:
            return zb, zbr, traj
        return zb, zbr

    def _drop(self, x: torch.Tensor) -> torch.Tensor:
        if not self.training or self.cfg.cond_dropout <= 0:
            return x
        keep = (torch.rand(x.shape[0], 1, device=x.device, dtype=x.dtype) >= self.cfg.cond_dropout).float()
        return x * keep

    # ---- losses ----
    def block_b(self, x: torch.Tensor, u: torch.Tensor, c: torch.Tensor) -> dict:
        """CITE-seq same-cell block: weighted DSM (score = -grad_b U at z_tilde),
        protein head, Zero-Masked reconstruction on centered targets (v2.0)."""
        mu, logvar = self.enc_blood(x)
        z = reparameterize(mu, logvar)
        sigma_rel = self.cfg.sigma_rel[torch.randint(len(self.cfg.sigma_rel), (1,)).item()]
        z_std = z.detach().std(dim=0, keepdim=True).clamp_min(1e-4)
        sig = sigma_rel * z_std
        noise = torch.randn_like(z)
        z_tilde = (z + sig * noise).detach().requires_grad_(True)
        e = self.potential(z_tilde, z_tilde, u, self._drop(c))
        g = torch.autograd.grad(e.sum(), z_tilde, create_graph=True)[0]
        # weighted DSM: || eps - sigma * grad U ||^2  (s_theta = -gradU matches -eps/sigma)
        dsm = (noise - sig * g).pow(2).mean()
        u_hat = self.protein_head(mu)
        prot = F.mse_loss(u_hat, u) + (1 - F.cosine_similarity(u_hat, u, dim=-1).mean())
        xr = self.dec_blood(z)
        rec = self.zm_mse(xr, x - self.gene_mean, x + self.blood_center) + self.centered_cos(xr, x - self.gene_mean)
        return {"dsm": dsm, "protein": prot, "recon": rec, "kl": kl_normal_free_bits(mu, logvar)}

    def block_a(
        self, x_blood: torch.Tensor, c_blood: torch.Tensor, x_brain_obs: torch.Tensor, c_brain: torch.Tensor
    ) -> dict:
        """Distribution-level disease block: encode blood -> fix zb -> solve brain
        equilibrium from zb init -> decode -> Sinkhorn divergence vs observed
        brain RNA of the same stratum. Protein input is inferred from blood RNA."""
        mu, lv = self.enc_blood(x_blood)
        n = x_brain_obs.shape[0]
        source = reparameterize(mu,lv) if self.cfg.sample_latent_a else mu
        zb = source[:n] if source.shape[0] >= n else source.repeat(2,1)[:n]
        # v2.0 û inferred control (blueprint §3.2 one-head-two-uses): the protein
        # head supplies the control input instead of zeros, so U's u-channel
        # receives training signal on disease data
        with torch.no_grad():
            u0 = self.protein_head(mu).detach()
        _, zbr_star = self.solve_equilibrium(zb, zb.clone(), u0, self._drop(c_brain[:n]))
        x_pred = self.dec_brain(zbr_star)
        obs_c = x_brain_obs - self.gene_mean
        sk = sinkhorn_divergence(x_pred - self.gene_mean, obs_c, self.cfg.sinkhorn_eps, self.cfg.sinkhorn_iters, self.cfg.ot_scale)
        # v1.3 variance-floor material: differentiable prediction-variance ratio
        var_ratio = x_pred.var(dim=0).mean() / x_brain_obs.var(dim=0).mean().clamp_min(1e-8).detach()
        return {"sinkhorn": sk, "x_pred": x_pred, "var_ratio": var_ratio}

    def recon(self, x_blood: torch.Tensor, x_brain: torch.Tensor) -> dict:
        mu_b, lv_b = self.enc_blood(x_blood)
        z_b = reparameterize(mu_b, lv_b)
        mu_B, lv_B = self.enc_brain(x_brain)
        z_B = reparameterize(mu_B, lv_B)
        xb_hat = self.dec_blood(z_b)
        xB_hat = self.dec_brain(z_B)
        rb = self.zm_mse(xb_hat, x_blood - self.gene_mean, x_blood + self.blood_center) + self.centered_cos(xb_hat, x_blood - self.gene_mean)
        rB = self.zm_mse(xB_hat, x_brain - self.gene_mean, x_brain + self.brain_center) + self.centered_cos(xB_hat, x_brain - self.gene_mean)
        return {"recon": 0.5 * (rb + rB), "kl": kl_normal_free_bits(mu_b, lv_b) + kl_normal_free_bits(mu_B, lv_B)}


@torch.no_grad()
def per_donor_spearman(u_hat: torch.Tensor, y_true: torch.Tensor, donor: torch.Tensor) -> dict:
    """G2-compatible metric: per donor, Spearman rho per protein across cells,
    averaged over proteins; then mean over donors. Ties broken by average rank
    via double argsort (adequate for continuous predictions)."""
    def rank_along_cells(t: torch.Tensor) -> torch.Tensor:
        from scipy.stats import rankdata
        return torch.as_tensor(rankdata(t.detach().cpu().numpy(),axis=0,method='average'),device=t.device,dtype=t.dtype)

    per_donor = {}
    for d in donor.unique().tolist():
        m = donor == d
        rp = rank_along_cells(u_hat[m])
        rt = rank_along_cells(y_true[m])
        rp = (rp - rp.mean(0)) / (rp.std(0,unbiased=False) + 1e-8)
        rt = (rt - rt.mean(0)) / (rt.std(0,unbiased=False) + 1e-8)
        rho = (rp * rt).mean(0)  # pearson of ranks == spearman (tie-free)
        per_donor[int(d)] = float(rho.mean())
    mean_rho = sum(per_donor.values()) / max(len(per_donor), 1)
    return {"per_donor": per_donor, "mean_rho": mean_rho}

