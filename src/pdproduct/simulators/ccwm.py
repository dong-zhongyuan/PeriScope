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
    """v2.1 free-bits KL: each latent dim must carry >= lam nats of information.
    sum(max(KL_dim, lam)) — decoder cannot ignore any dimension."""
    kl_dim = 0.5 * (mu.pow(2) + logvar.exp() - 1.0 - logvar)  # (B, z_dim)
    return torch.clamp(kl_dim, min=lam).sum(dim=-1).mean()

def kl_normal_legacy(mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
    return 0.5 * (mu.pow(2) + logvar.exp() - 1.0 - logvar).sum(dim=-1).mean()


def ot_cost(x: torch.Tensor, y: torch.Tensor, eps: float = 0.1, iters: int = 50) -> torch.Tensor:
    """Log-stabilized Sinkhorn OT cost, uniform weights, cost = squared euclidean
    normalized to mean 1 (scale-free across strata)."""
    c = torch.cdist(x, y) ** 2
    c = c / c.detach().mean().clamp_min(1e-8)
    n, m = c.shape
    f = torch.zeros(n, device=x.device, dtype=x.dtype)
    g = torch.zeros(m, device=x.device, dtype=x.dtype)
    for _ in range(iters):
        f = -eps * torch.logsumexp((g.unsqueeze(0) - c) / eps, dim=1)
        g = -eps * torch.logsumexp((f.unsqueeze(1) - c) / eps, dim=0)
    return f.mean() + g.mean()


def sinkhorn_divergence(x: torch.Tensor, y: torch.Tensor, eps: float = 0.1, iters: int = 50) -> torch.Tensor:
    """Debiased Sinkhorn divergence: OT(x,y) - 0.5 OT(x,x) - 0.5 OT(y,y)."""
    return ot_cost(x, y, eps, iters) - 0.5 * ot_cost(x, x, eps, iters) - 0.5 * ot_cost(y, y, eps, iters)


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
        self.protein_head = MLP(z, p, cfg.hidden, cfg.depth)
        self.u_blood = MLP(z + p + cd, 1, cfg.hidden, cfg.depth)
        self.u_brain = MLP(z + p + cd, 1, cfg.hidden, cfg.depth)
        self.u_coupling = MLP(2 * z + p + cd, 1, cfg.hidden, cfg.depth) if cfg.coupling else None
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
        zbr = zbr0.detach().clone().requires_grad_(True)
        traj = []
        for k in range(self.cfg.k_steps):
            e = self.potential(zb, zbr, u, c)
            g = torch.autograd.grad(e.sum(), zbr, create_graph=True)[0]
            # learnable scale keeps its graph; backtracking multiplies a 2^-j guard
            t_grad = F.softplus(self.step_sizes[k])
            with torch.no_grad():
                e0 = float(e.mean())
                g2 = float(g.detach().pow(2).sum(dim=-1).mean())
                guard = 1.0
                if g2 > 0.0:
                    for _ in range(6):
                        e_try = float(self.potential(zb, zbr - guard * float(t_grad) * g, u, c).mean())
                        if e_try <= e0 - 1e-4 * guard * float(t_grad) * g2:
                            break
                        guard *= 0.5
            zbr = zbr - guard * t_grad * g
            if return_traj:
                with torch.no_grad():
                    traj.append((e0, float(self.potential(zb, zbr, u, c).mean()), guard * float(t_grad)))
        if return_traj:
            return zb, zbr, traj
        return zb, zbr

    def _drop(self, x: torch.Tensor) -> torch.Tensor:
        if self.cfg.cond_dropout <= 0:
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
        u_hat = self.protein_head(z)
        prot = F.mse_loss(u_hat, u) + (1 - F.cosine_similarity(u_hat, u, dim=-1).mean())
        xr = self.dec_blood(z)
        rec = self.zm_mse(xr, x - self.gene_mean, x) + self.centered_cos(xr, x - self.gene_mean)
        return {"dsm": dsm, "protein": prot, "recon": rec, "kl": kl_normal_free_bits(mu, logvar)}

    def block_a(
        self, x_blood: torch.Tensor, c_blood: torch.Tensor, x_brain_obs: torch.Tensor, c_brain: torch.Tensor
    ) -> dict:
        """Distribution-level disease block: encode blood -> fix zb -> solve brain
        equilibrium from zb init -> decode -> Sinkhorn divergence vs observed
        brain RNA of the same stratum. u is zero on the brain side (no brain ADT)."""
        mu, _ = self.enc_blood(x_blood)
        n = x_brain_obs.shape[0]
        zb = mu[:n] if mu.shape[0] >= n else mu.repeat(2, 1)[:n]
        # v2.0 û inferred control (blueprint §3.2 one-head-two-uses): the protein
        # head supplies the control input instead of zeros, so U's u-channel
        # receives training signal on disease data
        with torch.no_grad():
            u0 = self.protein_head(mu).detach()
        _, zbr_star = self.solve_equilibrium(zb, zb.clone(), u0, self._drop(c_brain[:n]))
        x_pred = self.dec_brain(zbr_star)
        obs_c = x_brain_obs - self.gene_mean
        sk = sinkhorn_divergence(x_pred - self.gene_mean, obs_c, self.cfg.sinkhorn_eps, self.cfg.sinkhorn_iters)
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
        rb = self.zm_mse(xb_hat, x_blood - self.gene_mean, x_blood) + self.centered_cos(xb_hat, x_blood - self.gene_mean)
        rB = self.zm_mse(xB_hat, x_brain - self.gene_mean, x_brain) + self.centered_cos(xB_hat, x_brain - self.gene_mean)
        return {"recon": 0.5 * (rb + rB), "kl": kl_normal_free_bits(mu_b, lv_b) + kl_normal_free_bits(mu_B, lv_B)}


@torch.no_grad()
def per_donor_spearman(u_hat: torch.Tensor, y_true: torch.Tensor, donor: torch.Tensor) -> dict:
    """G2-compatible metric: per donor, Spearman rho per protein across cells,
    averaged over proteins; then mean over donors. Ties broken by average rank
    via double argsort (adequate for continuous predictions)."""
    def rank_along_cells(t: torch.Tensor) -> torch.Tensor:
        # average ranks are unnecessary for tie-free floats; use ordinal ranks
        return t.argsort(dim=0).argsort(dim=0).float()

    per_donor = {}
    for d in donor.unique().tolist():
        m = donor == d
        rp = rank_along_cells(u_hat[m])
        rt = rank_along_cells(y_true[m])
        rp = (rp - rp.mean(0)) / (rp.std(0) + 1e-8)
        rt = (rt - rt.mean(0)) / (rt.std(0) + 1e-8)
        rho = (rp * rt).mean(0)  # pearson of ranks == spearman (tie-free)
        per_donor[int(d)] = float(rho.mean())
    mean_rho = sum(per_donor.values()) / max(len(per_donor), 1)
    return {"per_donor": per_donor, "mean_rho": mean_rho}


def direction_concordance(
    model: "CCWM",
    blood_X: torch.Tensor, blood_cond2: torch.Tensor, blood_comp: torch.Tensor, blood_donor: torch.Tensor,
    brain_X: torch.Tensor, brain_cond2: torch.Tensor, brain_comp: torch.Tensor, brain_donor: torch.Tensor,
    comp: int = 0, n_per_donor: int = 256, n_perm: int = 2000, seed: int = 0,
) -> dict:
    """v1.1 readout on held-out donors: does the model's PREDICTED disease
    displacement (blood cells per condition -> encode -> solve -> decode ->
    donor-level pseudobulk -> condition contrast) align in GENE-LEVEL DIRECTION
    with the OBSERVED brain disease displacement (donor-level pseudobulk
    contrast)? This is the estimand unpaired data can identify (vs Sinkhorn,
    which v1.0 diagnostics showed is satisfiable by stratum-mean collapse).
    Returns spearman r, permutation p, and a specificity control vs the
    observed NEURONAL displacement (comp=3)."""
    import numpy as np
    from scipy.stats import spearmanr

    g = torch.Generator().manual_seed(seed)

    def predicted_means() -> torch.Tensor:
        """per-condition donor-level decoded means from blood side"""
        outs = {0: [], 1: []}
        m = blood_comp == comp
        for d in blood_donor[m].unique().tolist():
            idx = ((blood_donor == d) & m).nonzero(as_tuple=True)[0]
            idx = idx[torch.randperm(len(idx), generator=g)[:n_per_donor]]
            x = blood_X[idx]
            c = torch.zeros(len(idx), model.cfg.cond_dim, device=x.device)
            c[:, 0 if int(blood_cond2[idx[0]]) == 0 else 1] = 1.0
            c[:, 2 + comp] = 1.0  # source=blood block stays zero
            u0 = torch.zeros(len(idx), model.cfg.n_proteins, device=x.device)
            with torch.no_grad():
                mu, _ = model.enc_blood(x)
            with torch.enable_grad():
                _, z_star = model.solve_equilibrium(mu, mu.clone(), u0, c)
            with torch.no_grad():
                outs[1 if int(blood_cond2[idx[0]]) == 1 else 0].append(model.dec_brain(z_star.detach()).mean(dim=0))
        return torch.stack(outs[0]).mean(0), torch.stack(outs[1]).mean(0)

    with torch.no_grad():
        # donor -> condition map (donor-level, from the full arrays)
        d2c = {}
        m = blood_comp == comp
        for d in blood_donor[m].unique().tolist():
            d2c[int(d)] = int(blood_cond2[(blood_donor == d) & m][0])

    def observed_diff(comp_id: int) -> np.ndarray:
        mm = brain_comp == comp_id
        means, conds = [], []
        for d in brain_donor[mm].unique().tolist():
            sel = (brain_donor == d) & mm
            means.append(brain_X[sel].mean(dim=0))
            conds.append(int(brain_cond2[sel][0]))
        means = torch.stack(means)
        conds = torch.tensor(conds)
        return (means[conds == 1].mean(0) - means[conds == 0].mean(0)).cpu().numpy()

    pred0, pred1 = predicted_means()
    pred_diff = (pred1 - pred0).cpu().numpy()
    obs_my = observed_diff(comp)
    r = float(spearmanr(pred_diff, obs_my).statistic)
    rng = np.random.RandomState(seed)
    null = np.array([spearmanr(pred_diff, rng.permutation(obs_my)).statistic for _ in range(n_perm)])
    p = float((np.abs(null) >= abs(r)).mean())
    obs_neural = observed_diff(3)
    r_ctrl = float(spearmanr(pred_diff, obs_neural).statistic)
    return {"direction_r_myeloid": r, "permutation_p": p, "direction_r_neural_control": r_ctrl,
            "n_perm": n_perm}
