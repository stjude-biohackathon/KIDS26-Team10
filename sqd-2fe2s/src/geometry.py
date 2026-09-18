"""Build an idealized [2Fe-2S(SMe)4]^2- model cluster.

The Fe2S2 rhombic core is placed in the xy-plane with the Fe-Fe vector along x.
Terminal thiolates complete a distorted tetrahedron at each Fe.  Bond lengths are
taken from the crystallographic range for oxidized (Fe(III)/Fe(III)) synthetic
[2Fe-2S] analogues:

    Fe-Fe        2.70 A
    Fe-S(mu2)    2.20 A
    Fe-S(term)   2.30 A
    S-C          1.81 A
    C-H          1.09 A
    Fe-S-C     105 deg

Only the Fe 3d shell enters the active space, so the methyl conformation is not
a sensitive parameter; the geometry is built to be chemically sane rather than
relaxed.
"""

from __future__ import annotations

import numpy as np

FE_FE = 2.70
FE_SB = 2.20
FE_ST = 2.30
S_C = 1.81
C_H = 1.09
FE_S_C = np.deg2rad(105.0)
S_C_H = np.deg2rad(109.5)

# half-angle of the terminal S-Fe-S wedge, measured from the outward bisector
ST_HALF_ANGLE = np.deg2rad(57.0)


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def _perp(v: np.ndarray) -> np.ndarray:
    """Some unit vector orthogonal to v."""
    trial = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(_unit(v), trial)) > 0.9:
        trial = np.array([1.0, 0.0, 0.0])
    return _unit(trial - np.dot(trial, _unit(v)) * _unit(v))


def _methyl(s_pos: np.ndarray, fe_pos: np.ndarray, twist: float) -> list[tuple[str, np.ndarray]]:
    """Place a -CH3 on a thiolate S, bent FE_S_C away from the S->Fe direction."""
    u = _unit(fe_pos - s_pos)          # S -> Fe
    p = _perp(u)
    q = np.cross(u, p)
    # rotate the in-plane reference by `twist` so the four methyls are not coplanar
    p = np.cos(twist) * p + np.sin(twist) * q
    d = np.cos(FE_S_C) * u + np.sin(FE_S_C) * p
    c_pos = s_pos + S_C * _unit(d)

    # three H staggered about the S-C axis
    axis = _unit(c_pos - s_pos)
    a1 = _perp(axis)
    a2 = np.cross(axis, a1)
    atoms = [("C", c_pos)]
    for k in range(3):
        phi = 2.0 * np.pi * k / 3.0 + twist
        radial = np.cos(phi) * a1 + np.sin(phi) * a2
        h_dir = np.cos(np.pi - S_C_H) * axis + np.sin(np.pi - S_C_H) * radial
        atoms.append(("H", c_pos + C_H * _unit(h_dir)))
    return atoms


def build_cluster() -> list[tuple[str, np.ndarray]]:
    """Return [(symbol, xyz_angstrom), ...] for [2Fe-2S(SMe)4]^2-."""
    x_fe = FE_FE / 2.0
    y_sb = np.sqrt(FE_SB**2 - x_fe**2)

    fe = [np.array([+x_fe, 0.0, 0.0]), np.array([-x_fe, 0.0, 0.0])]
    sb = [np.array([0.0, +y_sb, 0.0]), np.array([0.0, -y_sb, 0.0])]

    atoms: list[tuple[str, np.ndarray]] = []
    atoms += [("Fe", fe[0]), ("Fe", fe[1])]
    atoms += [("S", sb[0]), ("S", sb[1])]

    # terminal thiolates: outward bisector along +/-x, wedge opening along z
    dz = FE_ST * np.sin(ST_HALF_ANGLE)
    dx = FE_ST * np.cos(ST_HALF_ANGLE)
    terminal = []
    for sign, fe_pos in zip((+1.0, -1.0), fe):
        for zsign in (+1.0, -1.0):
            terminal.append((fe_pos, fe_pos + np.array([sign * dx, 0.0, zsign * dz])))

    for k, (fe_pos, s_pos) in enumerate(terminal):
        atoms.append(("S", s_pos))
        atoms += _methyl(s_pos, fe_pos, twist=k * np.pi / 2.0)

    return atoms


def to_xyz(atoms: list[tuple[str, np.ndarray]], comment: str = "") -> str:
    lines = [str(len(atoms)), comment]
    for sym, pos in atoms:
        lines.append(f"{sym:2s} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}")
    return "\n".join(lines) + "\n"


def to_pyscf_atom(atoms: list[tuple[str, np.ndarray]]) -> list[tuple[str, tuple[float, float, float]]]:
    return [(sym, tuple(float(c) for c in pos)) for sym, pos in atoms]


def report(atoms: list[tuple[str, np.ndarray]]) -> str:
    """Sanity-check distances for the core."""
    pos = {i: p for i, (_, p) in enumerate(atoms)}
    d = lambda i, j: float(np.linalg.norm(pos[i] - pos[j]))
    out = [
        f"atoms            : {len(atoms)}",
        f"Fe-Fe            : {d(0, 1):.3f} A",
        f"Fe-S(mu2)        : {d(0, 2):.3f} A",
        f"S(mu2)-S(mu2)    : {d(2, 3):.3f} A",
        f"Fe-S(term)       : {d(0, 4):.3f} A",
        f"S-C              : {d(4, 5):.3f} A",
    ]
    v1 = pos[2] - pos[0]
    v2 = pos[3] - pos[0]
    ang = np.degrees(np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
    out.append(f"S-Fe-S (bridge)  : {ang:.1f} deg")
    return "\n".join(out)


if __name__ == "__main__":
    import pathlib
    import sys

    atoms = build_cluster()
    print(report(atoms))
    out = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if out is not None:
        out.write_text(to_xyz(atoms, "[2Fe-2S(SMe)4]2- idealized model, charge -2"))
        print(f"\nwrote {out}")
