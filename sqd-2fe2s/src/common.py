"""Shared paths, constants and JSON helpers."""

from __future__ import annotations

import json
import pathlib

import numpy as np

HARTREE2CM = 219474.6313632
CHEMICAL_ACCURACY_HA = 1.6e-3

ROOT = pathlib.Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
GEOM = ROOT / "geometry"


def save_json(path: pathlib.Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    def enc(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.bool_):
            return bool(o)
        raise TypeError(type(o))

    path.write_text(json.dumps(obj, indent=2, default=enc))


def load_json(path: pathlib.Path) -> dict:
    return json.loads(path.read_text())
