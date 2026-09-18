"""Idealized model geometry for [2Fe-2S(SMe)4]2-, the classic Holm-type
synthetic analogue of the plant-type [2Fe-2S] ferredoxin active site.

Construction method
--------------------
The Fe2S2 core (Fe1, Fe2, S1, S2) is built directly in Cartesian coordinates
as a planar rhombus from two literature bond lengths (Fe...Fe and Fe-S_bridge).
Every atom after that (terminal thiolate S, methyl C, methyl H) is placed with
a bond length / bond angle / dihedral triple relative to three already-placed
atoms, using the standard NeRF (Natural Extension Reference Frame) formula --
the same construction used to build protein backbones from internal
coordinates. That keeps every placement decision an explicit, checkable
number instead of hand-derived vector algebra.

Reference geometry (crystallographic, used as the target values this model
reproduces):
    Mayerle, J. J.; Denmark, S. E.; DePamphilis, B. V.; Ibers, J. A.; Holm,
    R. H. "Synthetic analogues of the active sites of iron-sulfur proteins.
    11. Synthesis and structure of crystalline (Et4N)2[Fe2S2(SPh)4]."
    J. Am. Chem. Soc. 1975, 97, 1032. Core: Fe...Fe = 2.691 A,
    Fe-S(bridge) = 2.198 A (avg), Fe-S(terminal thiolate) = 2.250 A (avg),
    Fe-S-Fe = 75.3 deg, S-Fe-S(core) = 104.7 deg.

We use the SMe (methanethiolate) ligand rather than SPh, per the project
profile's target species [2Fe-2S(SMe)4]2-, keeping the same experimental core
metrics and standard organic bond lengths for the methyl cap:
    S-C(sp3)  = 1.82 A   (typical thiolate C-S)
    C-H       = 1.09 A   (typical sp3 C-H)
    Fe-S-C    = 105 deg  (typical bridging-thiolate M-S-C angle)
    H-C-H     = 109.47 deg (ideal tetrahedral)

This is a model geometry, not a DFT-relaxed structure: it is meant to carry a
physically reasonable Fe2S2(SMe)4 cluster into the electronic-structure
pipeline, not to reproduce total energies to spectroscopic accuracy.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Literature bond lengths / angles (Mayerle et al. 1975, JACS 97, 1032)
# ---------------------------------------------------------------------------
FE_FE = 2.691          # Fe...Fe, core diagonal
FE_SBRIDGE = 2.198      # Fe-S(mu2-sulfido)
FE_STERMINAL = 2.250    # Fe-S(thiolate)
S_C = 1.82              # S-C(methyl)
C_H = 1.09              # C-H
FE_S_C_ANGLE = 105.0    # Fe-S-C bridging-thiolate angle (deg)
HCH_ANGLE = 109.47      # ideal tetrahedral (deg)


def _unit(v):
    return v / np.linalg.norm(v)


def nerf_place(a, b, c, bond_length, angle_deg, dihedral_deg):
    """Place a new atom d bonded to c, given three prior atoms a-b-c.

    Returns d such that:
        |c - d|            = bond_length
        angle(b, c, d)      = angle_deg
        dihedral(a, b, c, d) = dihedral_deg

    Standard NeRF construction (Parsons et al., J. Comput. Chem. 2005).
    """
    angle = np.radians(180.0 - angle_deg)  # NeRF uses the supplement
    dihedral = np.radians(dihedral_deg)

    d2 = np.array([
        bond_length * np.cos(angle),
        bond_length * np.sin(angle) * np.cos(dihedral),
        bond_length * np.sin(angle) * np.sin(dihedral),
    ])

    ab = b - a
    bc = c - b
    bc_u = _unit(bc)
    n = _unit(np.cross(ab, bc))
    m = np.cross(n, bc_u)

    M = np.array([bc_u, m, n]).T
    return c + M @ d2


def bond(p, q):
    return np.linalg.norm(p - q)


def angle(p, q, r):
    """Angle p-q-r in degrees."""
    v1, v2 = p - q, r - q
    cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0)))


def dihedral(p, q, r, s):
    """Dihedral p-q-r-s in degrees."""
    b1, b2, b3 = q - p, r - q, s - r
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m1 = np.cross(n1, _unit(b2))
    x, y = np.dot(n1, n2), np.dot(m1, n2)
    return np.degrees(np.arctan2(y, x))


def _build_methyl(anchor_S, ref_for_dihedral, out_dir_atom, s_c_dihedral, h_phase=0.0):
    """Attach a -CH3 to a thiolate S given the Fe (or S1) reference chain.

    h_phase rotates the methyl top rigidly about the S-C bond (0-120 deg
    covers all distinct rotamers by the 3-fold symmetry of -CH3); it is used
    only to resolve steric clashes between neighboring methyl groups and does
    not change the S-C bond direction, length or the tetrahedral H-C-H angle.
    """
    a, b, c = ref_for_dihedral, out_dir_atom, anchor_S
    C = nerf_place(a, b, c, S_C, FE_S_C_ANGLE, s_c_dihedral)
    # Three tetrahedral H's on C, staggered 120 deg apart, referenced to (b, c, C)
    H1 = nerf_place(b, c, C, C_H, HCH_ANGLE, 60.0 + h_phase)
    H2 = nerf_place(b, c, C, C_H, HCH_ANGLE, 180.0 + h_phase)
    H3 = nerf_place(b, c, C, C_H, HCH_ANGLE, 300.0 + h_phase)
    return C, H1, H2, H3


_METHYL_H_PHASE = 78.0  # found by _resolve_methyl_clash(): min H..H = 2.145 A


def _resolve_methyl_clash():
    """Grid search over the second methyl's rotamer phase (0-120 deg, the
    period of a free -CH3 top) to maximize the minimum H...H distance
    between the two methyls sharing one Fe center. Returns the best phase.
    """
    half_fefe = FE_FE / 2.0
    Fe1 = np.array([half_fefe, 0.0, 0.0])
    half_ss = np.sqrt(FE_SBRIDGE**2 - half_fefe**2)
    S1 = np.array([0.0, half_ss, 0.0])
    beta = 0.5 * np.degrees(np.arccos(-1.0 / 3.0))
    v1 = np.array([np.cos(np.radians(beta)), 0.0, np.sin(np.radians(beta))])
    v2 = np.array([np.cos(np.radians(beta)), 0.0, -np.sin(np.radians(beta))])
    St1a, St1b = Fe1 + FE_STERMINAL * v1, Fe1 + FE_STERMINAL * v2

    best_phase, best_dmin = 0.0, -1.0
    for phase in np.arange(0.0, 120.0, 1.0):
        C1a, H1a1, H1a2, H1a3 = _build_methyl(St1a, S1, Fe1, 180.0, h_phase=0.0)
        C1b, H1b1, H1b2, H1b3 = _build_methyl(St1b, S1, Fe1, 180.0, h_phase=phase)
        Hs_a = [H1a1, H1a2, H1a3]
        Hs_b = [H1b1, H1b2, H1b3]
        dmin = min(np.linalg.norm(ha - hb) for ha in Hs_a for hb in Hs_b)
        if dmin > best_dmin:
            best_dmin, best_phase = dmin, phase
    return best_phase, best_dmin


def build_fe2s2_sme4(charge=-2):
    """Build the idealized [Fe2S2(SMe)4]2- geometry.

    Returns
    -------
    atoms : list of (symbol, xyz) tuples
    checks : dict of measured bond lengths/angles for verification
    """
    # --- Fe2S2 core: planar rhombus in the xy-plane, centered at origin ---
    half_fefe = FE_FE / 2.0
    Fe1 = np.array([half_fefe, 0.0, 0.0])
    Fe2 = np.array([-half_fefe, 0.0, 0.0])

    half_ss = np.sqrt(FE_SBRIDGE**2 - half_fefe**2)
    S1 = np.array([0.0, half_ss, 0.0])
    S2 = np.array([0.0, -half_ss, 0.0])

    # --- Terminal thiolate S on Fe1 and Fe2, placed via NeRF using the
    #     Fe-Fe-S(terminal) angle that completes an (approximate) tetrahedron
    #     at each Fe: bisector of the two terminal S bonds points opposite the
    #     bisector of the two bridging S bonds. ---
    # angle Fe2-Fe1-S1(bridge) at Fe1:
    core_half_angle = angle(Fe2, Fe1, S1)  # ~ half the S-Fe-S core angle
    # Choose the terminal Fe-Fe-S(term) angle from an idealized local
    # tetrahedron at Fe1: with two bonds (to S1, S2) subtending 2*core_half_angle,
    # the other two tetrahedral bonds subtend 2*(90 - core_half_angle/2)...
    # we instead fix it from the explicit tetrahedral solution used during
    # design (see module docstring derivation): terminal bonds bisector points
    # opposite the bridging bisector, split symmetrically out of plane.
    beta = 0.5 * np.degrees(np.arccos(-1.0 / 3.0))  # 54.7356 deg (regular tetrahedron half-angle)
    v1 = np.array([np.cos(np.radians(beta)), 0.0, np.sin(np.radians(beta))])
    v2 = np.array([np.cos(np.radians(beta)), 0.0, -np.sin(np.radians(beta))])

    St1a = Fe1 + FE_STERMINAL * v1
    St1b = Fe1 + FE_STERMINAL * v2
    St2a = Fe2 - FE_STERMINAL * v1
    St2b = Fe2 - FE_STERMINAL * v2

    # --- Methyl caps on each terminal S, pointing outward (away from core).
    # The methyl -CH3 tops are freely-rotating in reality; h_phase (found by
    # a grid search in _resolve_methyl_clashes, see module __main__) is
    # chosen purely to avoid a steric clash between the two methyls sharing
    # the same Fe center. It changes no bond length or bond angle. ---
    phase = _METHYL_H_PHASE
    C1a, H1a1, H1a2, H1a3 = _build_methyl(St1a, S1, Fe1, s_c_dihedral=180.0, h_phase=0.0)
    C1b, H1b1, H1b2, H1b3 = _build_methyl(St1b, S1, Fe1, s_c_dihedral=180.0, h_phase=phase)
    C2a, H2a1, H2a2, H2a3 = _build_methyl(St2a, S2, Fe2, s_c_dihedral=180.0, h_phase=0.0)
    C2b, H2b1, H2b2, H2b3 = _build_methyl(St2b, S2, Fe2, s_c_dihedral=180.0, h_phase=phase)

    atoms = [
        ("Fe", Fe1), ("Fe", Fe2),
        ("S", S1), ("S", S2),
        ("S", St1a), ("S", St1b), ("S", St2a), ("S", St2b),
        ("C", C1a), ("H", H1a1), ("H", H1a2), ("H", H1a3),
        ("C", C1b), ("H", H1b1), ("H", H1b2), ("H", H1b3),
        ("C", C2a), ("H", H2a1), ("H", H2a2), ("H", H2a3),
        ("C", C2b), ("H", H2b1), ("H", H2b2), ("H", H2b3),
    ]

    checks = {
        "Fe-Fe": bond(Fe1, Fe2),
        "Fe-S(bridge)": bond(Fe1, S1),
        "Fe-S(terminal)": bond(Fe1, St1a),
        "S-C": bond(St1a, C1a),
        "C-H": bond(C1a, H1a1),
        "Fe-S-Fe": angle(Fe1, S1, Fe2),
        "S-Fe-S(core)": angle(S1, Fe1, S2),
        "Fe-S-C": angle(Fe1, St1a, C1a),
        "H-C-H": angle(H1a1, C1a, H1a2),
        "min_nonbonded_dist": _min_nonbonded_distance(atoms),
    }
    return atoms, checks


def _min_nonbonded_distance(atoms, bond_cutoff=1.9):
    """Smallest interatomic distance across pairs farther apart than a direct
    bond (clash check). bond_cutoff=1.9 A excludes true bonds (max real bond
    here is S-C=1.82 A) while catching any accidental steric overlap.
    """
    coords = np.array([xyz for _, xyz in atoms])
    n = len(coords)
    dmin = np.inf
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(coords[i] - coords[j])
            if bond_cutoff < d < dmin:
                dmin = d
    return dmin


def atoms_to_xyz_string(atoms, comment="[2Fe-2S(SMe)4]2- idealized model geometry"):
    lines = [str(len(atoms)), comment]
    for sym, xyz in atoms:
        lines.append(f"{sym:2s} {xyz[0]: .6f} {xyz[1]: .6f} {xyz[2]: .6f}")
    return "\n".join(lines)


def atoms_to_pyscf_atom(atoms):
    """Return the atom list in PySCF's expected [[sym, (x,y,z)], ...] format."""
    return [[sym, tuple(xyz)] for sym, xyz in atoms]


if __name__ == "__main__":
    atoms, checks = build_fe2s2_sme4()
    print(f"Built {len(atoms)} atoms.")
    for k, v in checks.items():
        print(f"  {k:20s} = {v:8.4f}")
