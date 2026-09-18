# Sample-based Quantum Diagonalization (SQD) for Iron-Sulfur Clusters 
**KIDS26 · Team 10 | Day 2: Final results and demo**

We built an SQD workflow for [Fe₂S₂(SMe)₄]²⁻, tested it on simulators and IBM quantum hardware, and compared its energies with classical references. The smallest active space matched CASCI at the reported precision. Larger calculations remained less accurate, and the measured 32-qubit comparison favored CASCI in both accuracy and runtime.

**[Project report](sqd-2fe2s/REPORT.md) · [Final presentation](presentation/BH26_Team10_Presentation.pptx) · [Code and detailed instructions](sqd-2fe2s/README.md)**

## Project Profile

- **Project name:** SQD for Iron-Sulfur Clusters
- **Question, problem, or opportunity:** How accurately can sample-based quantum diagonalization (SQD) reproduce the ground-state energy of [2Fe-2S(SMe)₄]²⁻?
- **Data, inputs, or evidence:** Model cluster geometry, STO-3G basis, PySCF Hamiltonian, AVAS-selected (10-electron, 10-orbital) Fe 3d active space, and simulated or hardware-generated ffsim/LUCJ samples.
- **Output:** An executable workflow, CASCI benchmarks, convergence and sampling analysis, orbital occupancies, hardware results, and a final report.
- **Tools and stack:** Python, PySCF, AVAS, ffsim/LUCJ, Qiskit, `qiskit-addon-sqd`, and IBM Quantum Runtime; GPU calculations used `gpu4pyscf`.
- **Team lead:** Amandeep Singh Bhatia [deepquantum88](https://github.com/deepquantum88)
- **Team members and roles:** Emily Clifton (documentation), Soham Bopardikar (code and simulation), Amandeep Singh Bhatia (analysis and presentation). See [team roles](project-management/team.md).
- **Communication:** [Slack: #team10](https://stjudebiohackathon.slack.com/archives/C0BSA3JQGS1)

## Vision and Mission

- **Vision:** Evaluate quantum-assisted approaches to iron–sulfur electronic structure.
- **Mission:** Benchmark SQD against classical references and identify how sampling, subspace size, and hardware noise affect accuracy and computational cost.

## About

Iron–sulfur clusters are found in enzyme active sites. Their magnetic coupling and Fe–S covalency make their electronic structure challenging for single-reference methods.

Our workflow connects PySCF → AVAS → ffsim/LUCJ → qiskit-addon-sqd. We compare SQD with CASCI using the same Hamiltonian and fixed orbitals. CASCI provides an exact reference within the selected active space.

## Final Results

Errors below are relative to the reference for each active space. The hardware column lists runs without added error suppression.

| Active space | Qubits | Simulator SQD error | Hardware SQD error | Reference |
| --- | --- | --- | --- | --- |
| 10 electrons, 10 orbitals | 20 | 0.0000 mHa at reported precision | 0.0000 mHa at reported precision | Exact CASCI |
| 22 electrons, 16 orbitals | 32 | +34.5 mHa | +529.7 mHa | Exact CASCI |
| 26 electrons, 18 orbitals | 36 | Not run | +653.4 mHa | Extrapolated singlet from four computed spin sectors |

- **Smallest case:** Hardware post-processing recovered the full 63,504-determinant space. A uniform-random sampling control also matched the reference, so this result does not demonstrate quantum advantage or subspace compression.
- **Runtime:** In the reported 32-qubit comparison, CASCI took 0.96 hours and SQD took 4.54 hours on an otherwise idle machine using one core. SQD remained 34.5 mHa above CASCI after ten recovery iterations.
- **Hardware sampling:** The fraction of shots with the target electron count and spin projection fell from 5.90% at 20 qubits to 0.47% at 32 and 0.10% at 36.
- **Magnetic properties:** The tested SQD spin-sector calculations did not reliably reproduce exchange coupling. Reported coupling estimates use classical spin-sector energies.

The 1.6 mHa accuracy target was met for the smallest case only. All calculations used a minimal basis, the active spaces are not strictly nested, and the 36-qubit singlet reference is extrapolated. Conclusions apply to the settings tested. See the [report](sqd-2fe2s/REPORT.md) for methods and limitations and [results record](sqd-2fe2s/docs/RESULTS_MASTER.md) for numerical provenance.

## Run the Demo

From the repository root, with Python and the required dependencies available:

```bash
cd sqd-2fe2s
pip install -r requirements.txt
python3 run_all.py --quick
```

The quick run is reported to take about two minutes. Install dependencies before the presentation. The full initial (10e,10o) workflow uses `python3 run_all.py` and is reported to take about twelve minutes; larger calculations and QPU submissions are separate.

To replay the saved 20-qubit hardware samples without submitting a new quantum job:

```bash
python3 src/stage7_hardware.py --tag sto-3g_fe3d --replay
```

See [hardware instructions](sqd-2fe2s/docs/HARDWARE.md) for details. For the short presentation, existing [figures](sqd-2fe2s/results/figures/) and the report provide a backup if the live run exceeds the available time.

## Roadmap and Handoff

| When | Focus |
| --- | --- |
| Day 1 | Initial (10e,10o) CASCI reference and SQD implementation |
| Day 2 | Benchmarking, convergence analysis, and expanded calculations |
| Day 3 | Final results, documentation, presentation, and demo |

Possible follow-up work includes completing the missing 36-qubit classical spin sectors, testing larger SQD subspaces and shallower circuits, and examining basis-set effects.

[References](sqd-2fe2s/docs/REFERENCES.md) · [Run summary](sqd-2fe2s/docs/RUN_SUMMARY.md) · [License](LICENSE.md)


