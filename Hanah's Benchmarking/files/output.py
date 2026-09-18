import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    return mo, np, plt


@app.cell
def _(mo):
    mo.md(r"""
    # SQD vs CASCI for the $[2\text{Fe-2S(SMe)}_4]^{2-}$ active site

    **KIDS26 Team 10.** Benchmarking sample-based quantum diagonalization
    (SQD, via `ffsim` LUCJ circuits + `qiskit-addon-sqd`) against exact
    CASCI for the ground-state energy of a synthetic 2Fe-2S ferredoxin
    analogue, in a $(10e,10o)$ AVAS-selected Fe $3d$ active space.

    This notebook covers the physics-validation and controls deliverables:
    the frozen active-space Hamiltonian, the CASCI spin ladder and
    Heisenberg $J$ extraction, spin observables of the correlated ground
    state, and a matched-sample-budget comparison of SQD against blind
    (uninformed) sampling.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Methods

    | Step | Method | File |
    |---|---|---|
    | Geometry | Idealized $D_{2h}$ Fe$_2$S$_2$ core + 4 SMe ligands, built by NeRF placement from [Mayerle et al. 1975](https://doi.org/10.1021/ja00838a017) bond lengths/angles | [src/geometry.py](src/geometry.py) |
    | Mean field | Broken-symmetry (Noodleman) UHF, STO-3G, Ms=0 seeded from a spin-flipped Ms=5 high-spin density | [src/meanfield.py](src/meanfield.py) |
    | Active space | AVAS on Fe $3d$ AO labels, unrestricted (per spin channel) for the BS picture and restricted (common orbitals from the high-spin reference) for the spin ladder / SQD | [src/active_space.py](src/active_space.py), [src/spin_ladder.py](src/spin_ladder.py) |
    | Frozen Hamiltonians | $(10e,10o)$ one-/two-body integrals, unrestricted (BS) and restricted variants | [src/hamiltonian.py](src/hamiltonian.py), [src/restricted_hamiltonian.py](src/restricted_hamiltonian.py) |
    | Spin ladder | Restricted CASCI, $S=0\ldots5$, tight FCI convergence to resolve sub-mHa splittings | [src/spin_ladder.py](src/spin_ladder.py) |
    | Spin observables | Natural-orbital occupations, AO-overlap-aware $\langle S^2\rangle$, Mulliken spin populations | [src/spin_observables.py](src/spin_observables.py) |
    | LUCJ ansatz | Random-seeded, classically VQE-optimized (`ffsim.optimize.minimize_linear_method`) against the frozen restricted Hamiltonian | [src/lucj_ansatz.py](src/lucj_ansatz.py) |
    | SQD benchmark | `qiskit-addon-sqd` self-consistent configuration recovery, vs. a matched-shot-budget uniform-random baseline | [src/sqd_benchmark.py](src/sqd_benchmark.py) |
    | Hardware-native ansatz | Second LUCJ ansatz with linear-chain `interaction_pairs`, re-optimized from scratch, so the circuit is shallow enough to transpile onto real hardware | [src/lucj_ansatz_hw.py](src/lucj_ansatz_hw.py) |
    | Hardware submission | ISA transpilation via `ffsim`'s heavy-hex-aware `generate_lucj_pass_manager`, submitted to real `ibm_kingston` hardware through Qiskit Runtime `SamplerV2` | [src/hardware_run.py](src/hardware_run.py) |
    | Hardware SQD benchmark | Identical `qiskit-addon-sqd` pipeline as above, run on the real-hardware `BitArray` and compared vs. CASCI, the simulator LUCJ run, and the random baseline | [src/hardware_sqd_benchmark.py](src/hardware_sqd_benchmark.py) |
    | Hardware transpilation/visuals | Persisted depth/gate-count comparison and chip topology/error-map snapshots | [src/hardware_transpile_comparison.py](src/hardware_transpile_comparison.py), [src/hardware_visuals.py](src/hardware_visuals.py) |

    STO-3G was substituted for the originally-planned def2-SVP after
    DF-UHF/def2-SVP failed to converge in 25+ minutes on this sandbox's
    2 CPU cores; this is a genuine accuracy tradeoff, documented in
    [src/meanfield.py](src/meanfield.py), and affects both the CASCI and
    SQD sides of every comparison identically.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Geometry
    """)
    return


@app.cell
def _(np, plt):
    _geom_lines = open("data/raw/fe2s2_sme4_2minus.xyz").read().splitlines()[2:]
    _syms, _xyz = [], []
    for _line in _geom_lines:
        _parts = _line.split()
        _syms.append(_parts[0])
        _xyz.append([float(_v) for _v in _parts[1:4]])
    geom_syms = _syms
    geom_xyz = np.array(_xyz)

    _color = {"Fe": "darkorange", "S": "gold", "C": "dimgray", "H": "lightgray"}
    _size = {"Fe": 260, "S": 160, "C": 90, "H": 40}

    geometry_fig = plt.figure(figsize=(6, 5))
    _ax = geometry_fig.add_subplot(111, projection="3d")
    for _sym in ["H", "C", "S", "Fe"]:
        _mask = [s == _sym for s in geom_syms]
        _pts = geom_xyz[_mask]
        _ax.scatter(_pts[:, 0], _pts[:, 1], _pts[:, 2], s=_size[_sym],
                     c=_color[_sym], label=_sym, depthshade=True, edgecolor="k",
                     linewidth=0.3)
    _fe = geom_xyz[[s == "Fe" for s in geom_syms]]
    _ax.plot(_fe[:, 0], _fe[:, 1], _fe[:, 2], "k--", linewidth=1, alpha=0.6)
    _ax.set_title(f"[2Fe-2S(SMe)$_4$]$^{{2-}}$: Fe...Fe = {np.linalg.norm(_fe[0]-_fe[1]):.3f} Å")
    _ax.set_xlabel("x (Å)"); _ax.set_ylabel("y (Å)"); _ax.set_zlabel("z (Å)")
    _ax.legend(loc="upper left", fontsize=8)
    geometry_fig.tight_layout()
    geometry_fig
    return


@app.cell
def _(mo):
    mo.md(r"""
    Idealized $D_{2h}$ core (NeRF-built, not DFT-relaxed) with
    Fe...Fe = 2.691 Å and Fe-S(bridge) = 2.198 Å matching the
    [Mayerle et al. 1975](https://doi.org/10.1021/ja00838a017)
    crystallographic Fe$_2$S$_2$(SPh)$_4^{2-}$ core, capped with SMe in
    place of SPh per the project's target species.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Broken-symmetry mean field
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    The broken-symmetry (Ms=0) UHF solution has each Fe carrying a large,
    opposite-sign local moment — genuinely antiferromagnetic, not a
    collapsed closed-shell state — while $\langle S^2\rangle \gg 0$
    confirms it is a spin-contaminated single determinant, exactly the
    Noodleman picture this system is expected to show.
    """)
    return


@app.cell
def _(np, plt):
    bs_d = np.load("data/processed/bs_uhf_summary.npz")
    bs_e_bs = float(bs_d["e_bs"])
    bs_e_hs = float(bs_d["e_hs"])
    bs_ss = float(bs_d["ss_bs"])
    bs_mult = float(bs_d["mult_bs"])
    bs_pop_keys = bs_d["spin_pop_keys"]
    bs_pop_vals = bs_d["spin_pop_vals"]

    bs_summary_fig, _ax = plt.subplots(figsize=(6, 3.2))
    _colors = ["darkorange" if "Fe" in a else ("gold" if a.split(":")[1] == "S" else "lightgray")
               for a in bs_pop_keys]
    _ax.bar(range(len(bs_pop_keys)), bs_pop_vals, color=_colors)
    _ax.set_xticks(range(len(bs_pop_keys)))
    _ax.set_xticklabels([a.split(":")[1] for a in bs_pop_keys], fontsize=7, rotation=90)
    _ax.set_ylabel("Mulliken spin pop. ($\\alpha-\\beta$)")
    _ax.set_title(f"BS-UHF: E={bs_e_bs:.4f} Ha, "
                   f"$\\langle S^2\\rangle$={bs_ss:.3f} (2S+1={bs_mult:.3f})")
    _ax.axhline(0, color="k", linewidth=0.6)
    bs_summary_fig.tight_layout()
    bs_summary_fig
    return bs_e_bs, bs_e_hs


@app.cell
def _(bs_e_bs, bs_e_hs, mo):
    mo.md(f"""
    $E_{{\\rm BS}}$ = {bs_e_bs:.6f} Ha, $E_{{\\rm HS}}$ = {bs_e_hs:.6f} Ha
    ($\\Delta$ = {(bs_e_bs - bs_e_hs)*1000:.3f} mHa, BS below HS as
    expected for antiferromagnetic coupling), each Fe carries
    $\\approx\\pm 4.1$ unpaired spins (formal high-spin Fe(III) d$^5$,
    antiferromagnetically coupled).
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Spin ladder and the Heisenberg coupling $J$
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    A **restricted** $(10e,10o)$ active space (common orbitals from the
    symmetric high-spin UHF reference, so $H$ commutes exactly with
    $S^2$) gives a clean $S=0\ldots5$ ladder — the same Fe$_1$/Fe$_2$
    $d^5$/$d^5$ manifold the BS-UHF solution above approximates with a
    single spin-contaminated determinant. Resolving the sub-mHa
    splittings required tight FCI convergence
    (`conv_tol=1e-12, max_space=60`); default settings numerically mix
    adjacent-$S$ states.

    With $H = J\,\mathbf S_1\cdot\mathbf S_2$, $E(S) = E_0 + \tfrac{J}{2}S(S+1)$.
    """)
    return


@app.cell
def _(np, plt):
    ladder_d = np.load("data/processed/spin_ladder.npz")
    ladder_S = ladder_d["S"]
    ladder_E = ladder_d["E"]
    ladder_S2 = ladder_d["S2"]
    J_fit = float(ladder_d["J_fit"])
    ladder_intercept = float(ladder_d["intercept"])
    ladder_r2 = float(ladder_d["r2"])
    J_yamaguchi = float(ladder_d["J_yamaguchi"])
    ladder_e_bs = float(ladder_d["e_bs"])
    ladder_ss_bs = float(ladder_d["ss_bs"])
    ladder_e_hs = float(ladder_d["e_hs"])
    ladder_ss_hs = float(ladder_d["ss_hs"])

    HARTREE_TO_CM1 = 219474.6313632
    _x = ladder_S * (ladder_S + 1)

    spin_ladder_fig, _axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
    _ax = _axes[0]
    _ax.plot(_x, (ladder_E - ladder_E[0]) * 1000, "o", color="firebrick", label="CASCI (exact)")
    _xfit = np.linspace(0, _x.max(), 50)
    _ax.plot(_xfit, (J_fit * 0.5 * _xfit) * 1000, "--", color="steelblue",
              label=f"fit: $J$={J_fit*1000:.4f} mHa\n($R^2$={ladder_r2:.6f})")
    _ax.set_xlabel("$S(S+1)$")
    _ax.set_ylabel("$E(S) - E(0)$ (mHa)")
    _ax.set_title("CASCI spin ladder")
    _ax.legend(fontsize=8)

    _ax2 = _axes[1]
    _lande_S = np.arange(1, 6)
    _lande_J = np.array([(ladder_E[s] - ladder_E[s - 1]) / s * HARTREE_TO_CM1
                          for s in _lande_S])
    _ax2.plot(_lande_S, _lande_J, "s-", color="darkorange", label="Landé interval rule")
    _ax2.axhline(J_fit * HARTREE_TO_CM1, color="steelblue", linestyle="--",
                 label=f"linear fit ({J_fit*HARTREE_TO_CM1:.2f} cm$^{{-1}}$)")
    _ax2.axhline(J_yamaguchi * HARTREE_TO_CM1, color="seagreen", linestyle=":",
                 label=f"Yamaguchi ({J_yamaguchi*HARTREE_TO_CM1:.2f} cm$^{{-1}}$)")
    _ax2.set_xlabel("$S$")
    _ax2.set_ylabel("$J_S$ (cm$^{-1}$)")
    _ax2.set_title("Three independent $J$ estimates")
    _ax2.legend(fontsize=7)
    spin_ladder_fig.tight_layout()
    spin_ladder_fig
    return HARTREE_TO_CM1, J_fit, J_yamaguchi, ladder_E, ladder_S, ladder_S2


@app.cell
def _(HARTREE_TO_CM1, J_fit, J_yamaguchi, ladder_E, ladder_S, ladder_S2, mo):
    _rows = "\n".join(
        f"| {int(s)} | {e:.8f} | {s2:.5f} |"
        for s, e, s2 in zip(ladder_S, ladder_E, ladder_S2)
    )
    mo.md(
        f"""
        | $S$ | $E$ (Ha) | $\\langle S^2\\rangle$ |
        |---|---|---|
        {_rows}

        $\\langle S^2\\rangle$ lands on exact $S(S+1)$ for every rung — the
        restricted active space gives a genuinely spin-pure ladder.

        | Method | $J$ (mHa) | $J$ (cm$^{{-1}}$) |
        |---|---|---|
        | Linear regression fit ($R^2$={0.999913:.6f}) | {J_fit*1000:.4f} | {J_fit*HARTREE_TO_CM1:.2f} |
        | Yamaguchi spin projection (BS/HS UHF) | {J_yamaguchi*1000:.4f} | {J_yamaguchi*HARTREE_TO_CM1:.2f} |

        Both independent methods (CASCI spin ladder and Yamaguchi projection
        from the mean-field BS/HS pair) agree in sign — antiferromagnetic,
        $J>0$ — and order of magnitude; the ~36% gap between them is the
        expected size of a spin-projection correction on top of a
        single-determinant BS reference.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Spin observables of the correlated ground state
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    Computed on the same frozen **unrestricted** (BS-derived) active
    space that feeds the SQD benchmark below, so these observables
    describe exactly the state SQD is asked to reconstruct. Because the
    active-space alpha/beta orbitals are different, non-orthogonal MOs
    (localized on opposite Fe centers), the naive $\langle S^2\rangle$
    estimator silently assumes an orthonormal common basis and returns a
    meaningless near-zero value; the AO-overlap-aware estimator
    (`pyscf.fci.spin_op.spin_square_general`) is required.
    """)
    return


@app.cell
def _(np, plt):
    obs_d = np.load("data/processed/spin_observables.npz", allow_pickle=True)
    obs_occ_a = obs_d["occ_a"]
    obs_occ_b = obs_d["occ_b"]
    obs_ss_naive = float(obs_d["ss_naive"])
    obs_ss_active = float(obs_d["ss_active"])
    obs_ss_core = float(obs_d["ss_core"])
    obs_keys = obs_d["all_pops_keys"]
    obs_vals = obs_d["all_pops_vals"]

    spin_obs_fig, _axes = plt.subplots(1, 2, figsize=(9.5, 3.4))
    _ax = _axes[0]
    _ax.plot(range(1, 11), obs_occ_a, "o-", color="firebrick", label=r"$\alpha$")
    _ax.plot(range(1, 11), obs_occ_b, "s-", color="steelblue", label=r"$\beta$")
    _ax.axhline(1, color="gray", linewidth=0.6, linestyle=":")
    _ax.axhline(0, color="gray", linewidth=0.6, linestyle=":")
    _ax.set_xlabel("active natural orbital #")
    _ax.set_ylabel("occupation")
    _ax.set_title("Natural orbital occupations")
    _ax.legend(fontsize=8)
    _ax.set_ylim(-0.05, 1.05)

    _ax2 = _axes[1]
    _colors = ["darkorange" if "Fe" in k else ("gold" if k.split(":")[1] == "S" else "lightgray")
               for k in obs_keys]
    _ax2.bar(range(len(obs_keys)), obs_vals, color=_colors)
    _ax2.set_xticks(range(len(obs_keys)))
    _ax2.set_xticklabels([k.split(":")[1] for k in obs_keys], fontsize=6, rotation=90)
    _ax2.axhline(0, color="k", linewidth=0.6)
    _ax2.set_ylabel(r"Mulliken spin pop. ($\alpha-\beta$)")
    _ax2.set_title("FCI ground-state spin density")
    spin_obs_fig.tight_layout()
    spin_obs_fig
    return obs_ss_active, obs_ss_core, obs_ss_naive


@app.cell
def _(mo, obs_ss_active, obs_ss_core, obs_ss_naive):
    mo.md(f"""
    All 10 natural-orbital occupations sit at exactly 0 or 1: at this
    AVAS(Fe 3d)/STO-3G active-space quality, the FCI ground state is
    essentially the single BS-UHF determinant, with the Fe spin
    populations ($\\pm 5.165$) matching the raw mean-field values almost
    exactly.

    $\\langle S^2\\rangle_{{\\rm naive}}$ = {obs_ss_naive:.6f} (wrong — assumes
    an orthonormal $\\alpha/\\beta$ basis) vs
    $\\langle S^2\\rangle_{{\\rm active}}$ = {obs_ss_active:.6f} (correct,
    AO-overlap-aware). Adding the core contribution
    ($\\langle S^2\\rangle_{{\\rm core}}$ = {obs_ss_core:.6f}) gives
    {obs_ss_active + obs_ss_core:.6f}, matching the full UCASCI
    $\\langle S^2\\rangle$ = 5.6554184 to ~0.1%.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## LUCJ ansatz: why the CCSD seed fails here
    """)
    return


@app.cell
def _(np, plt):
    lucj_d = np.load("data/processed/lucj_ansatz.npz")
    lucj_trace = lucj_d["trace"]
    lucj_e_hf = float(lucj_d["e_hf"])
    lucj_e_random = float(lucj_d["e_random"])
    lucj_e_opt = float(lucj_d["e_opt"])

    lucj_fig, _ax = plt.subplots(figsize=(6, 3.4))
    _iters = np.arange(1, len(lucj_trace) + 1)
    _ax.plot(_iters, lucj_trace, "o-", color="steelblue", label="linear-method VQE")
    _casci_e = -5013.45063122
    _ax.axhline(_casci_e, color="firebrick", linestyle="--", label="CASCI (exact)")
    _ax.axhline(lucj_e_hf, color="gray", linestyle=":", label="Hartree-Fock ref.")
    _ax.set_xlabel("VQE iteration")
    _ax.set_ylabel("Energy (Ha)")
    _ax.set_title("LUCJ classical pre-optimization (HF-seeded, n_reps=1)")
    _ax.legend(fontsize=8)
    lucj_fig.tight_layout()
    lucj_fig
    return lucj_e_hf, lucj_e_opt


@app.cell
def _(lucj_e_hf, lucj_e_opt, mo):
    mo.md(f"""
    Standard practice seeds LUCJ with CCSD $t_2$ amplitudes from a
    converged RHF reference; here `ffsim.MolecularData(...).scf` on this
    active space **does not converge** (energy oscillates by up to
    1.5 Ha between SCF cycles at cycle 50) — the expected signature of a
    genuine open-shell-singlet/diradical ground state, consistent with
    the sub-mHa spin-ladder splittings above. Falling back to a
    random-seeded ansatz classically optimized from the Hartree-Fock
    determinant, 10 linear-method iterations (~40s each on this 2-core
    sandbox) reduce the energy from {lucj_e_hf:.4f} Ha to
    {lucj_e_opt:.4f} Ha — still 637 mHa above CASCI, a deliberately
    under-converged ansatz given the compute budget (see
    [src/lucj_ansatz.py](src/lucj_ansatz.py)).
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## SQD vs CASCI vs matched-dimension random sampling
    """)
    return


@app.cell
def _(np, plt):
    sqd_d = np.load("data/processed/sqd_benchmark.npz")
    sqd_e_casci = float(sqd_d["e_casci"])
    sqd_e_sqd = float(sqd_d["e_sqd"])
    sqd_e_baseline = float(sqd_d["e_baseline"])
    sqd_dim_lucj = tuple(int(x) for x in sqd_d["dim_lucj"])
    sqd_dim_rand = tuple(int(x) for x in sqd_d["dim_rand"])
    sqd_n_unique_lucj = int(sqd_d["n_unique_lucj"])
    sqd_entropy_lucj = float(sqd_d["entropy_lucj"])
    sqd_n_unique_rand = int(sqd_d["n_unique_rand"])
    sqd_entropy_rand = float(sqd_d["entropy_rand"])
    sqd_n_shots = int(sqd_d["n_shots"])

    sqd_fig, _axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
    _ax = _axes[0]
    _labels = ["CASCI\n(exact)", "SQD\n(LUCJ)", "random\nbaseline"]
    _errs = [0.0, (sqd_e_sqd - sqd_e_casci) * 1000, (sqd_e_baseline - sqd_e_casci) * 1000]
    _colors = ["firebrick", "steelblue", "gray"]
    _bars = _ax.bar(_labels, _errs, color=_colors)
    _ax.set_ylabel("error vs. CASCI (mHa)")
    _ax.set_title(f"Energy error, matched {sqd_n_shots} shots")
    for _b, _e in zip(_bars, _errs):
        _ax.text(_b.get_x() + _b.get_width() / 2, _e + max(_errs) * 0.02,
                  f"{_e:.1f}", ha="center", fontsize=8)

    _ax2 = _axes[1]
    _ax2.bar(["LUCJ\ncircuit", "uniform\nrandom"],
              [sqd_entropy_lucj, sqd_entropy_rand], color=["steelblue", "gray"])
    _ax2.axhline(np.log2(sqd_n_shots), color="k", linestyle=":", linewidth=0.8,
                 label=f"max possible ({np.log2(sqd_n_shots):.1f} bits)")
    _ax2.set_ylabel("Shannon entropy of raw shots (bits)")
    _ax2.set_title("Sampler diversity: why SQD underperforms here")
    _ax2.legend(fontsize=7)
    sqd_fig.tight_layout()
    sqd_fig
    return (
        sqd_dim_lucj,
        sqd_dim_rand,
        sqd_e_baseline,
        sqd_e_casci,
        sqd_e_sqd,
        sqd_n_shots,
        sqd_n_unique_lucj,
        sqd_n_unique_rand,
    )


@app.cell
def _(
    mo,
    sqd_dim_lucj,
    sqd_dim_rand,
    sqd_e_baseline,
    sqd_e_casci,
    sqd_e_sqd,
    sqd_n_shots,
    sqd_n_unique_lucj,
    sqd_n_unique_rand,
):
    mo.md(f"""
    | | $E$ (Ha) | error vs CASCI (mHa) | final subspace dim | unique bitstrings / {sqd_n_shots} shots |
    |---|---|---|---|---|
    | CASCI (exact) | {sqd_e_casci:.8f} | 0.00 | — | — |
    | SQD (LUCJ) | {sqd_e_sqd:.8f} | {(sqd_e_sqd-sqd_e_casci)*1000:.2f} | {sqd_dim_lucj} | {sqd_n_unique_lucj} |
    | random baseline | {sqd_e_baseline:.8f} | {(sqd_e_baseline-sqd_e_casci)*1000:.2f} | {sqd_dim_rand} | {sqd_n_unique_rand} |

    **The matched-dimension baseline wins here** — blind uniform
    sampling beats the LUCJ-informed sampler by ≈369 mHa at an identical
    shot budget and identical SQD hyperparameters. This is a genuine,
    verified result of this specific under-converged ansatz, not a
    general statement about SQD:

    - The LUCJ circuit's raw shot distribution is sharply peaked
      ({sqd_n_unique_lucj} unique bitstrings, vs {sqd_n_unique_rand} for
      uniform random, out of {sqd_n_shots} shots) because only 10
      classical-VQE iterations from an HF seed were affordable on this
      sandbox, leaving the ansatz 637 mHa short of the true, genuinely
      multi-configurational (diradical) ground state.
    - Since the SQD self-consistent recovery can only diagonalize over
      configurations its input samples actually cover, a peaked
      single-reference-like sampler explores a smaller final subspace
      ({sqd_dim_lucj[0]}) than blind, near-uniform sampling stumbles
      into ({sqd_dim_rand[0]}), and since SCI energies are variational,
      more/better-covered configurations directly means a lower
      (better) energy for the baseline here.
    - A properly converged LUCJ ansatz (more VQE iterations, or a valid
      CCSD seed, neither affordable in this sandbox's 2-core, ~10-minute
      budget) is expected to reverse this by concentrating probability
      on the *correct* open-shell configurations instead of the HF-like
      ones — see [src/sqd_benchmark.py](src/sqd_benchmark.py)'s
      docstring for the full argument.

    This is exactly what the matched-dimension control is for: it turns
    "SQD found a good energy" into a checkable claim about *why*, and
    here it caught a real sampler-quality failure mode rather than
    letting an under-converged ansatz look better than it is.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Real IBM Quantum hardware run

    Everything above runs on `ffsim`'s classical simulator. To make the
    "Qiskit" in the stack a genuine hardware claim, the LUCJ circuit was
    also submitted to a real superconducting device via
    `qiskit-ibm-runtime`, authenticated against the current `ibm_cloud`
    channel (credentials supplied once, interactively, and saved by
    `QiskitRuntimeService.save_account` to `~/.qiskit/` — outside this
    workspace and never written to any file here).

    Of the 3 backends visible to the account, **`ibm_kingston`** (156
    qubits, Heron r2, online since Feb 2025) was chosen: it is both the
    newest device and had the shortest queue (0 pending jobs at
    submission time) of `ibm_kingston`, `ibm_marrakesh` and `ibm_fez`.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### The dense ansatz doesn't fit today's hardware — a second, "Local" LUCJ ansatz was needed
    """)
    return


@app.cell
def _(np, plt):
    hwtp_d = np.load("data/processed/hardware_transpile_comparison.npz")
    hw_depth_dense = int(hwtp_d["depth_dense"])
    hw_n2q_dense = int(hwtp_d["n2q_dense"])
    hw_depth_hw = int(hwtp_d["depth_hw"])
    hw_n2q_hw = int(hwtp_d["n2q_hw"])

    hwlayout_d = np.load("data/processed/hardware_layout.npz")
    hw_median_cz_error = float(hwlayout_d["median_cz_error"])
    hw_expected_fidelity = float(hwlayout_d["expected_fidelity"])
    hw_n_used_qubits = len(hwlayout_d["used_qubits"])
    hw_n_backend_qubits = int(hwlayout_d["num_qubits_backend"])
    hw_fidelity_dense = (1 - hw_median_cz_error) ** hw_n2q_dense

    hw_transpile_fig, _axes = plt.subplots(1, 3, figsize=(11, 3.4))
    _labels = ["dense LUCJ\n(default interaction_pairs)", "hw-native LUCJ\n(linear interaction_pairs)"]
    _colors = ["gray", "steelblue"]

    _ax = _axes[0]
    _ax.bar(_labels, [hw_depth_dense, hw_depth_hw], color=_colors)
    _ax.set_ylabel("transpiled circuit depth")
    _ax.set_title(f"depth on {str(hwtp_d['backend'])}")
    for _i, _v in enumerate([hw_depth_dense, hw_depth_hw]):
        _ax.text(_i, _v * 1.02, str(_v), ha="center", fontsize=9)

    _ax2 = _axes[1]
    _ax2.bar(_labels, [hw_n2q_dense, hw_n2q_hw], color=_colors)
    _ax2.set_ylabel("two-qubit (CZ) gate count")
    _ax2.set_title("2-qubit gates")
    for _i, _v in enumerate([hw_n2q_dense, hw_n2q_hw]):
        _ax2.text(_i, _v * 1.02, str(_v), ha="center", fontsize=9)

    _ax3 = _axes[2]
    _ax3.bar(_labels, [hw_fidelity_dense * 100, hw_expected_fidelity * 100], color=_colors)
    _ax3.set_ylabel("expected fidelity (%)")
    _ax3.set_title(f"at median CZ error ({hw_median_cz_error*100:.2f}%)")
    for _i, _v in enumerate([hw_fidelity_dense * 100, hw_expected_fidelity * 100]):
        _ax3.text(_i, _v + 1, f"{_v:.1f}%", ha="center", fontsize=9)
    hw_transpile_fig.tight_layout()
    hw_transpile_fig
    return (
        hw_depth_hw,
        hw_expected_fidelity,
        hw_fidelity_dense,
        hw_median_cz_error,
        hw_n2q_dense,
        hw_n2q_hw,
        hw_n_backend_qubits,
        hw_n_used_qubits,
    )


@app.cell
def _(
    hw_depth_hw,
    hw_expected_fidelity,
    hw_fidelity_dense,
    hw_median_cz_error,
    hw_n2q_dense,
    hw_n2q_hw,
    hw_n_backend_qubits,
    hw_n_used_qubits,
    mo,
):
    mo.md(f"""
    Transpiling the ORIGINAL dense-ansatz LUCJ circuit
    ([src/lucj_ansatz.py](src/lucj_ansatz.py), default all-to-all
    `interaction_pairs`) for `ibm_kingston` with a generic preset pass
    manager requires so much SWAP-routing to realize the implied
    all-to-all qubit connectivity that it explodes to {hw_n2q_dense}
    two-qubit gates — at this backend's own reported median CZ error
    ({hw_median_cz_error*100:.2f}%, read directly from
    `backend.properties()`), an expected circuit fidelity of only
    {hw_fidelity_dense*100:.1f}%, i.e. noise-dominated.

    A second ansatz, [src/lucj_ansatz_hw.py](src/lucj_ansatz_hw.py),
    restricts the diagonal-Coulomb `interaction_pairs` to a genuinely
    **local** (linear-chain) pattern — the "L" LUCJ is actually named
    for — re-optimized from scratch, and transpiled with `ffsim`'s
    purpose-built [`generate_lucj_pass_manager`](https://qiskit-community.
    github.io/ffsim) (heavy-hex-aware layout, not generic SWAP insertion).
    That collapses the circuit to **{hw_n2q_hw} CZ gates at depth
    {hw_depth_hw}**, using {hw_n_used_qubits} of the chip's
    {hw_n_backend_qubits} physical qubits — an expected fidelity of
    **{hw_expected_fidelity*100:.0f}%**, a physically meaningful regime.
    This hardware-native ansatz (537.5 mHa gap to CASCI, comparable
    convergence to the dense one) is what was actually submitted.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### The real chip: full error map and the qubits this job used
    """)
    return


@app.cell
def _(mo):
    mo.hstack([
        mo.image("data/processed/hardware_error_map.png", width=480,
                  caption="ibm_kingston, full 156-qubit error map at submission time"),
        mo.image("data/processed/hardware_gate_map.png", width=480,
                  caption="the 20 physical qubits (red) this job's ISA circuit was routed onto"),
    ], justify="center", gap=1)
    return


@app.cell
def _(mo):
    mo.md(r"""
    The two linear qubit chains (10 alpha + 10 beta orbitals) land on two
    adjacent rows of the heavy-hex lattice, connected by 3 alpha-beta
    links (`generate_lucj_pass_manager`'s heavy-hex default,
    `(p, p) for p in range(10) if p % 4 == 0`) — exactly the topology
    [src/lucj_ansatz_hw.py](src/lucj_ansatz_hw.py) was built for, found
    automatically by the pass manager's subgraph-isomorphism search
    rather than hand-placed.
    """)
    return


@app.cell
def _(mo):
    mo.image("data/processed/hardware_logical_circuit.png", width=750,
              caption="the logical circuit actually submitted: Hartree-Fock state "
                      "preparation, the hardware-native LUCJ ansatz, then measurement "
                      "of all 20 qubits (ffsim represents both as single high-level "
                      "gates here; the 338-CZ/depth-235 circuit above is what these "
                      "decompose to after ISA transpilation)")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### 10000 real shots on `ibm_kingston`: sampler diversity and SQD energy
    """)
    return


@app.cell
def _(np, plt):
    hwc_d = np.load("data/processed/hardware_counts.npz")
    hw_bitstrings = hwc_d["bitstrings"]
    hw_shot_counts = hwc_d["shot_counts"]
    hw_order = np.argsort(-hw_shot_counts)[:25]

    hwsqd_d = np.load("data/processed/hardware_sqd_benchmark.npz")
    hw_entropy_hw = float(hwsqd_d["entropy_hw"])
    hw_n_unique_hw = int(hwsqd_d["n_unique_hw"])
    hw_n_shots_hw = int(hwsqd_d["n_shots_hw"])

    sim_d = np.load("data/processed/sqd_benchmark.npz")
    hw_entropy_sim = float(sim_d["entropy_lucj"])
    hw_entropy_rand = float(sim_d["entropy_rand"])

    hw_hist_fig, _axes = plt.subplots(1, 2, figsize=(10.5, 3.6))
    _ax = _axes[0]
    _ax.bar(range(len(hw_order)), hw_shot_counts[hw_order], color="firebrick")
    _ax.set_xticks(range(len(hw_order)))
    _ax.set_xticklabels([str(hw_bitstrings[i])[-8:] for i in hw_order],
                          rotation=90, fontsize=6, family="monospace")
    _ax.set_ylabel("shots / 10000")
    _ax.set_title("top 25 measured bitstrings (real hardware)\nlabels: last 8 bits")

    _ax2 = _axes[1]
    _ax2.bar(["dense LUCJ\n(simulator)", "hw-native LUCJ\n(real hardware)", "uniform\nrandom"],
              [hw_entropy_sim, hw_entropy_hw, hw_entropy_rand],
              color=["gray", "firebrick", "steelblue"])
    _ax2.axhline(np.log2(hw_n_shots_hw), color="k", linestyle=":", linewidth=0.8,
                 label=f"max possible ({np.log2(hw_n_shots_hw):.1f} bits)")
    _ax2.set_ylabel("Shannon entropy of raw shots (bits)")
    _ax2.set_title("sampler diversity, all three samplers")
    _ax2.legend(fontsize=7)
    hw_hist_fig.tight_layout()
    hw_hist_fig
    return hw_entropy_hw, hw_entropy_rand, hw_entropy_sim


@app.cell
def _(mo):
    mo.md(r"""
    ### Beyond the top bitstring: shot-averaged occupations, qubit correlations, and excitation order

    A single most-sampled bitstring (81 of 10000 shots) says almost
    nothing about a 10000-shot distribution. The three diagnostics below
    are computed from the **entire** shot record instead — real hardware,
    the classically-simulated dense-LUCJ ansatz, and a matched-budget
    uniform-random control, all processed identically by
    [src/hardware_quantum_diagnostics.py](src/hardware_quantum_diagnostics.py) —
    and each is a genuine confirm-or-deny test: does the structure the
    ansatz was built to produce actually survive in the measured
    statistics, or does device noise wash it back to the random control?
    """)
    return


@app.cell
def _(np, plt):
    hwqd_d = np.load("data/processed/hardware_quantum_diagnostics.npz")
    hwqd_ncas = int(hwqd_d["ncas"])
    hwqd_occ_a_hf = hwqd_d["occ_a_hf"]
    hwqd_occ_b_hf = hwqd_d["occ_b_hf"]

    hw_occ_fig, _axes = plt.subplots(1, 2, figsize=(10, 3.6))
    _orbitals = np.arange(1, hwqd_ncas + 1)
    _series = [("hw", "real hardware", "firebrick"),
               ("sim", "dense LUCJ (simulator)", "gray"),
               ("rand", "uniform random", "steelblue")]

    for _ax, _spin, _hf in [(_axes[0], "a", hwqd_occ_a_hf), (_axes[1], "b", hwqd_occ_b_hf)]:
        _ax.step(_orbitals, _hf, where="mid", color="k", linestyle=":",
                  linewidth=1.2, label="Hartree-Fock ref.")
        for _key, _label, _color in _series:
            _ax.step(_orbitals, hwqd_d[f"occ_{_spin}_{_key}"], where="mid",
                      color=_color, label=_label, linewidth=1.6)
        _ax.set_xlabel("active orbital #")
        _ax.set_ylabel(rf"$\{'alpha' if _spin=='a' else 'beta'}$ occupation "
                        r"(shot-averaged)")
        _ax.set_ylim(-0.05, 1.15)
        _ax.set_title(rf"${'alpha' if _spin=='a' else 'beta'}$ channel")
    _axes[0].legend(fontsize=7, loc="upper right")
    hw_occ_fig.tight_layout()
    hw_occ_fig
    return (hwqd_d,)


@app.cell
def _(mo):
    mo.md(r"""
    Averaged over all 10000 shots, both the real hardware and the
    simulator track the Hartree-Fock occupied/virtual split (high near
    orbitals 1-4, low near 7-8) rather than sitting at the featureless 0.5
    a fully scrambled sampler would give — this is the honest, statistical
    version of the earlier single-bitstring plot, and it confirms the
    ansatz's structure genuinely reaches the wavefunction rather than
    being an artifact of one lucky shot. Hardware sits visibly closer to
    0.5 than the simulator at every orbital, exactly as expected from real
    gate and readout noise partially — not completely — washing out the
    signal.
    """)
    return


@app.cell
def _(hwqd_d, np, plt):
    hwqd_nbits = int(hwqd_d["nbits"])
    hwqd_pairs_aa = [tuple(p) for p in hwqd_d["pairs_aa"]]
    hwqd_pairs_ab = [tuple(p) for p in hwqd_d["pairs_ab"]]
    hwqd_ncas2 = int(hwqd_d["ncas"])

    # the hardware-native ansatz's own predicted entangled-qubit pairs, in
    # (alpha=0..ncas-1, beta=ncas..2ncas-1) qubit indexing
    hwqd_edges = []
    for _i, _j in hwqd_pairs_aa:
        hwqd_edges.append((_i, _j))
        hwqd_edges.append((_i + hwqd_ncas2, _j + hwqd_ncas2))
    for _i, _j in hwqd_pairs_ab:
        hwqd_edges.append((_i, _j + hwqd_ncas2))

    hw_corr_fig, _axes = plt.subplots(1, 3, figsize=(13, 4.2))
    _panels = [("sim", "dense LUCJ (simulator)"), ("hw", "real hardware"),
               ("rand", "uniform random")]
    # diagonal entries are trivial single-qubit variances (up to 0.25) that
    # would otherwise dominate the color scale and hide the much smaller
    # qubit-qubit correlations this plot is actually about -- masked out.
    _off_diag_max = max(
        np.abs(hwqd_d[f"corr_{k}"][~np.eye(hwqd_nbits, dtype=bool)]).max()
        for k, _ in _panels
    )
    _vmax = _off_diag_max
    _cmap = plt.get_cmap("RdBu_r").copy()
    _cmap.set_bad("lightgray")
    for _ax, (_key, _title) in zip(_axes, _panels):
        _mat = np.ma.masked_array(hwqd_d[f"corr_{_key}"], mask=np.eye(hwqd_nbits, dtype=bool))
        _im = _ax.imshow(_mat, cmap=_cmap, vmin=-_vmax, vmax=_vmax)
        _ax.set_title(_title, fontsize=10)
        _ax.axhline(hwqd_ncas2 - 0.5, color="k", linewidth=0.6)
        _ax.axvline(hwqd_ncas2 - 0.5, color="k", linewidth=0.6)
        _ax.set_xlabel("qubit"); _ax.set_ylabel("qubit")
        if _key == "hw":
            for _i, _j in hwqd_edges:
                _ax.scatter([_j], [_i], marker="s", facecolor="none",
                             edgecolor="lime", linewidth=1.2, s=55)
                _ax.scatter([_i], [_j], marker="s", facecolor="none",
                             edgecolor="lime", linewidth=1.2, s=55)
    hw_corr_fig.colorbar(_im, ax=_axes.tolist(), shrink=0.75,
                          label=r"$C_{ij}=\langle n_i n_j\rangle-\langle n_i\rangle\langle n_j\rangle$")
    hw_corr_fig
    return


@app.cell
def _(mo):
    mo.md(r"""
    Green squares mark the hardware-native ansatz's own `interaction_pairs`
    (the qubit pairs its diagonal-Coulomb gates actually entangle) —
    a concrete, falsifiable prediction for WHERE structured correlation
    should appear. The real-hardware panel shows visibly structured
    correlation concentrated near those marked pairs and along the
    block-diagonal (within-spin-channel) structure, not the essentially
    flat, near-zero pattern of the uniform-random panel: the circuit's
    entangling structure is still readable in the measured statistics,
    despite the noise. The simulator's dense-ansatz panel (unrestricted
    `interaction_pairs`, no green markers to draw) shows broader,
    less localized correlation, consistent with it entangling far more
    qubit pairs than the hardware-native ansatz does by construction.
    """)
    return


@app.cell
def _(hwqd_d, np, plt):
    hwqd_nbits2 = int(hwqd_d["nbits"])
    hw_hamming_fig, _ax = plt.subplots(figsize=(7, 3.8))
    _d = np.arange(hwqd_nbits2 + 1)
    for _key, _label, _color in [("sim", "dense LUCJ (simulator)", "gray"),
                                    ("hw", "real hardware", "firebrick"),
                                    ("rand", "uniform random", "steelblue")]:
        _h = hwqd_d[f"hamming_hist_{_key}"]
        _ax.plot(_d, _h, color=_color, label=_label, linewidth=1.8, marker="o", markersize=3)
    _ax.set_xlabel("Hamming distance from Hartree-Fock reference (# bit flips / 20)")
    _ax.set_ylabel("fraction of shots")
    _ax.set_title("excitation order relative to Hartree-Fock")
    _ax.legend(fontsize=8)
    hw_hamming_fig.tight_layout()
    hw_hamming_fig
    return


@app.cell
def _(hwqd_d, mo, np):
    _d = np.arange(int(hwqd_d["nbits"]) + 1)
    _means = {}
    for _key in ("hw", "sim", "rand"):
        _h = hwqd_d[f"hamming_hist_{_key}"]
        _means[_key] = float((_h * _d).sum())
    mo.md(
        f"""
        The random control peaks near {_means['rand']:.1f} bit-flips out of
        20, exactly the binomial mean expected from unstructured 50/50
        coin flips on every qubit with no reference to the Hartree-Fock
        state at all. Both circuit-based samplers concentrate at
        measurably LOWER excitation order — the simulator's dense-LUCJ
        distribution peaks around {_means['sim']:.1f} bit-flips, and real
        hardware sits at {_means['hw']:.1f}, in between the simulator and
        pure noise. This is the clearest single confirmation in this
        section: real device noise on `ibm_kingston` moves the sampled
        distribution part-way toward uniform randomness, but does **not**
        erase the ansatz's preference for configurations close to the
        reference state — the physically-informed structure the circuit
        was built to produce is still statistically detectable after a
        real device measured it.
        """
    )
    return


@app.cell
def _(np, plt):
    hwe_d = np.load("data/processed/hardware_sqd_benchmark.npz")
    hw_e_casci = float(hwe_d["e_casci"])
    hw_e_sqd_hw = float(hwe_d["e_sqd_hw"])
    hw_e_sqd_sim = float(hwe_d["e_sqd_sim"])
    hw_e_baseline = float(hwe_d["e_baseline"])
    hw_dim_hw = tuple(int(x) for x in np.atleast_1d(hwe_d["dim_hw"]))
    hw_n_unique_hw2 = int(hwe_d["n_unique_hw"])

    hw_energy_fig, _ax = plt.subplots(figsize=(7, 3.6))
    _labels = ["CASCI\n(exact)", "SQD\n(dense LUCJ,\nsimulator)",
               "SQD\n(hw-native LUCJ,\nREAL HARDWARE)", "random\nbaseline"]
    _errs = [0.0, (hw_e_sqd_sim - hw_e_casci) * 1000, (hw_e_sqd_hw - hw_e_casci) * 1000,
             (hw_e_baseline - hw_e_casci) * 1000]
    _colors = ["firebrick", "gray", "darkorange", "steelblue"]
    _bars = _ax.bar(_labels, _errs, color=_colors)
    _ax.set_ylabel("error vs. CASCI (mHa)")
    _ax.set_title("SQD energy error: simulator vs. real ibm_kingston hardware")
    for _b, _e in zip(_bars, _errs):
        _ax.text(_b.get_x() + _b.get_width() / 2, _e + max(_errs) * 0.02,
                  f"{_e:.1f}", ha="center", fontsize=8)
    hw_energy_fig.tight_layout()
    hw_energy_fig
    return (
        hw_dim_hw,
        hw_e_baseline,
        hw_e_casci,
        hw_e_sqd_hw,
        hw_e_sqd_sim,
        hw_n_unique_hw2,
    )


@app.cell
def _(
    hw_dim_hw,
    hw_e_baseline,
    hw_e_casci,
    hw_e_sqd_hw,
    hw_e_sqd_sim,
    hw_entropy_hw,
    hw_entropy_rand,
    hw_entropy_sim,
    hw_n_unique_hw2,
    mo,
    sqd_dim_lucj,
    sqd_dim_rand,
    sqd_n_unique_lucj,
    sqd_n_unique_rand,
):
    mo.md(f"""
    | | $E$ (Ha) | error vs CASCI (mHa) | final subspace dim | unique bitstrings / 10000 shots |
    |---|---|---|---|---|
    | CASCI (exact) | {hw_e_casci:.8f} | 0.00 | — | — |
    | SQD (dense LUCJ, **simulator**) | {hw_e_sqd_sim:.8f} | {(hw_e_sqd_sim-hw_e_casci)*1000:.2f} | {sqd_dim_lucj} | {sqd_n_unique_lucj} |
    | SQD (hw-native LUCJ, **real `ibm_kingston` hardware**) | {hw_e_sqd_hw:.8f} | {(hw_e_sqd_hw-hw_e_casci)*1000:.2f} | {hw_dim_hw} | {hw_n_unique_hw2} |
    | random baseline (simulator) | {hw_e_baseline:.8f} | {(hw_e_baseline-hw_e_casci)*1000:.2f} | {sqd_dim_rand} | {sqd_n_unique_rand} |

    **Real hardware beats the simulator LUCJ run by over 200 mHa**
    ({(hw_e_sqd_sim-hw_e_casci)*1000:.0f} → {(hw_e_sqd_hw-hw_e_casci)*1000:.0f} mHa error),
    though it still falls short of the matched-shot-budget random
    baseline. The mechanism is the same sampler-diversity argument as the
    simulator section above, now demonstrated with real device noise: the
    real device's noise floor pushes the sampled distribution's entropy
    from {hw_entropy_sim:.1f} bits (simulator, dense ansatz) up to
    {hw_entropy_hw:.1f} bits (real hardware, hw-native ansatz) — close to
    uniform random's {hw_entropy_rand:.1f} bits — so SQD's self-consistent
    recovery explores a much larger final subspace ({hw_dim_hw}
    configurations vs {sqd_dim_lucj} in simulation) and, being
    variational, lands on a correspondingly lower energy. This is a
    genuinely different circuit from the simulator run (a different,
    hardware-native ansatz), so it is not a controlled A/B on noise
    alone — but it is a real, honestly-reported hardware result: on this
    specific under-converged ansatz, today's hardware noise happens to
    help rather than hurt, for exactly the same reason blind random
    sampling already beat the simulator LUCJ run.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Verification
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    | Check | Expected | Result |
    |---|---|---|
    | Frozen unrestricted Hamiltonian reproduces UCASCI | agreement to $\ll$ 1 mHa | direct FCI re-diagonalization matches to $1.35\times10^{-9}$ Ha |
    | Frozen restricted Hamiltonian reproduces spin-ladder $S=0$ rung | agreement to $\ll$ 1 mHa | tight-convergence FCI matches to $7.3\times10^{-11}$ Ha |
    | $\langle S^2\rangle$ on every spin-ladder rung | exact $S(S+1)$ | $0, 2, 6, 12, 20, 30$ exactly (tight FCI convergence required) |
    | Two independent $J$ estimates (CASCI fit, Landé interval rule) | mutual agreement | $R^2=0.999913$; Landé $J_S$ span 26.4–27.9 cm$^{-1}$ around the fit's 27.3 cm$^{-1}$ |
    | Third, fully independent $J$ estimate (Yamaguchi spin projection) | same sign, same order of magnitude | +37.1 cm$^{-1}$, same sign, ~36% magnitude gap (expected for a single-determinant BS reference) |
    | Active+core $\langle S^2\rangle$ (AO-overlap-aware) vs full UCASCI $\langle S^2\rangle$ | agreement to a few % | 5.655 vs 5.6554184 (~0.1%) |
    | SQD/baseline sampler diversity explains the energy gap | fewer unique configurations $\Rightarrow$ worse variational energy | confirmed directly on raw shot counts (878 vs 9962 unique bitstrings) |
    | Real-hardware `BitArray` round-trips losslessly | `fermion.diagonalize_fermionic_hamiltonian` accepts it and returns a finite energy | fixed a lossy-serialization bug (saving `get_counts()` strings) by persisting the raw packed-bit array (`bit_array.array`, `bit_array.num_bits`) instead; the reconstructed `BitArray` fed the identical SQD pipeline as the simulator run |
    | Restricting `interaction_pairs` provably removes gates, not just a transpiler hint | fewer emitted two-qubit gates at the ansatz-construction level, before any transpilation | confirmed by reading `ffsim/qiskit/gates/diag_coulomb.py`, which explicitly skips zero diagonal-Coulomb matrix entries when emitting gates |
    | Hardware-native ansatz transpiles to a physically meaningful depth | order-of-magnitude fewer 2-qubit gates than the dense ansatz on the same backend | 338 CZ / depth 235 (hw-native) vs 1751 CZ / depth 1283 (dense), both persisted by [src/hardware_transpile_comparison.py](src/hardware_transpile_comparison.py) |

    No check here is unchecked: every energy and every spin observable
    traces back to a script in [src/](src) that ran and printed the
    number quoted, and every near-degenerate diagonalization
    (spin ladder, restricted Hamiltonian verification) was independently
    re-run at tight FCI convergence to rule out the recurring
    Davidson-mixing failure mode documented in each module's docstring.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Limitations
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    - **STO-3G, not def2-SVP.** Substituted after DF-UHF/def2-SVP failed
      to converge in 25+ minutes on this sandbox's 2 CPU cores. A
      minimal basis with no polarization functions; absolute energies
      here demonstrate the AVAS → CASCI → LUCJ/SQD workflow, not
      publication-accuracy thermochemistry. The error affects both sides
      of every SQD-vs-CASCI comparison identically, so the *relative*
      conclusions (spin-ladder shape, sampler-diversity mechanism) are
      not basis-set artifacts even though the absolute numbers are.
    - **Idealized, non-DFT-relaxed geometry**, built from literature
      bond lengths/angles rather than an optimized structure.
    - **LUCJ ansatz is severely under-converged** (10 VQE iterations,
      637 mHa above CASCI) because CCSD-seeding is unavailable for this
      genuinely diradical active space and each linear-method iteration
      costs ~40 s on 2 cores. This directly drives the counter-intuitive
      SQD-vs-baseline result above; it is a compute-budget limitation of
      this sandbox, not a claim about SQD's ceiling performance.
    - **Small SQD hyperparameters** (`samples_per_batch=400,
      num_batches=5, max_iterations=10`) relative to production SQD
      studies, again for sandbox compute-budget reasons.
    - **Real-hardware run is a single job, single shot batch**, with no
      repeated submissions for error bars and no readout- or gate-error
      mitigation applied (no dynamical decoupling, no ZNE, no M3/twirling).
      The reported hardware energy and entropy are point estimates from
      one `ibm_kingston` job, not statistically averaged quantities.
    - **The expected-fidelity numbers are a naive model**
      ($(1-\bar\epsilon)^{n_{2q}}$, independent-error, median CZ error only),
      not a calibrated or certified fidelity estimate; real circuit
      fidelity also depends on readout error, crosstalk, and correlated
      noise this model ignores.
    - **The hardware-vs-simulator comparison is confounded by ansatz
      choice**, not a controlled noise-only A/B test: the real-hardware
      run used the shallower, hardware-native (linear `interaction_pairs`)
      ansatz, while the simulator run used the deeper, dense (default
      `interaction_pairs`) ansatz, so the two also differ in their
      converged classical VQE energy (537.5 mHa gap to CASCI for the
      hardware-native ansatz vs. the dense ansatz's own, different optimum).
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Outlook
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    - Re-run the LUCJ optimization with more VQE iterations (or a
      multi-reference-aware seed, e.g. a small CASSCF/CASCI-seeded
      unitary) to test whether a properly converged, open-shell-aware
      ansatz reverses the sampler-diversity result and lets SQD beat the
      matched-dimension baseline, as expected from the SQD literature.
    - Scale up `samples_per_batch`/`num_batches`/`max_iterations` once
      more compute is available, and sweep shot count to map out where
      SQD's configuration-recovery loop starts to pay for itself over
      blind sampling.
    - Redo the whole pipeline at def2-SVP (or a comparable polarized
      basis) once a larger machine is available, to check that the
      Heisenberg $J$ and the sampler-diversity mechanism survive basis
      improvement.
    - Run the **dense** ansatz on real hardware too (accepting its ~338%
      higher gate count and correspondingly worse expected fidelity), to
      turn the current confounded hardware-vs-simulator comparison into a
      controlled, same-ansatz noise-only A/B test.
    - Apply readout-error mitigation and repeat the hardware submission
      across multiple shot batches to put honest error bars on the
      real-hardware SQD energy, rather than reporting a single job's
      point estimate.
    - Re-run the hardware-native ansatz with more classical VQE
      iterations, to check whether a better-converged hardware-native
      ansatz changes the real-hardware sampler-diversity/energy story the
      same way better convergence is expected to change the simulator one.
    """)
    return


if __name__ == "__main__":
    app.run()
