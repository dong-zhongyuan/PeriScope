"""All paths in one place. Scripts must not hard-code absolute paths.

(The sister project hard-coded its server paths into nearly every script, which is why
nobody could later tell which copy of a file was the live one.)
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
REFERENCE = ROOT / "reference"
PAPERS = ROOT / "papers"

#: interpreter that has numpy/scipy/pandas/statsmodels/matplotlib (shared with the hh project)
PYTHON = Path(os.environ.get("PDV3_PYTHON", ROOT.parent / "hh/env/bin/python"))


def paper(name):
    d = PAPERS / name
    if not d.is_dir():
        raise FileNotFoundError(f"no such paper directory: {d}")
    return d


def results_dir(paper_name, make=True):
    d = paper(paper_name) / "results"
    if make:
        d.mkdir(parents=True, exist_ok=True)
    return d


def figures_dir(paper_name, make=True):
    d = paper(paper_name) / "figures"
    if make:
        d.mkdir(parents=True, exist_ok=True)
    return d


def data_dir(source, sub=None, make=False):
    d = DATA / source / sub if sub else DATA / source
    if make:
        d.mkdir(parents=True, exist_ok=True)
    return d


def env_override(name, default):
    """Allow a path to be pointed elsewhere without editing code (e.g. on the server)."""
    return Path(os.environ.get(name, default))
