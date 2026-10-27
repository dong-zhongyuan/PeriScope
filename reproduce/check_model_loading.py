"""Load the registered weights and run a small synthetic CPU inference check.

This checks checkpoint/code compatibility, shapes and finite outputs. Synthetic
inputs do not reproduce a scientific benchmark or a target nomination.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "_work/optimization_20261003"))
from ccwm import CCWM, CCWMConfig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", required=True, type=Path)
    args = parser.parse_args()
    registry = json.loads((ROOT / "results/optimization_20261003/model_registry.json").read_text())
    with np.load(ROOT / "reproduce/inputs/brain_train_preprocessing.npz", allow_pickle=False) as z:
        states = list(z["state_vocab"].astype(str))
        n_genes = len(z["genes"])
    torch.set_num_threads(2)
    for seed, row in registry["models"].items():
        path = args.models / f"seed_{seed}/model.pt"
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            raise SystemExit(f"Checkpoint SHA-256 mismatch: seed {seed}")
        checkpoint = torch.load(path, map_location="cpu", weights_only=True)
        model = CCWM(CCWMConfig(**checkpoint["cfg"]))
        model.load_state_dict(checkpoint["state_dict"], strict=True)
        model.eval()
        assert model.cfg.n_genes == n_genes and model.cfg.cond_dim == len(states) + 4
        for parameter in model.parameters():
            parameter.requires_grad_(False)
        x = torch.zeros(2, n_genes)
        with torch.no_grad():
            zb, _ = model.enc_blood(x)
            protein = model.protein_head(zb)
        condition = torch.zeros(2, model.cfg.cond_dim)
        condition[:, 0] = 1
        condition[:, 2 + states.index("astro")] = 1
        condition[:, -1] = 1
        with torch.enable_grad():
            _, brain = model.solve_equilibrium(zb, zb.clone(), protein, condition)
            expression = model.dec_brain(brain)
        assert expression.shape == (2, n_genes) and torch.isfinite(expression).all()
        print(f"PASS: seed {seed}; strict load and synthetic CPU inference, {n_genes} genes", flush=True)
    print("PASS: all ten registered models")


if __name__ == "__main__":
    main()
