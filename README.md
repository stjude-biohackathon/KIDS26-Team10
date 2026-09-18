# Sample-based Quantum Diagonalization (SQD) for Iron-Sulfer Clusters 
**KIDS26 · Team 10 | Day 2: Benchmarking**

## Project Profile

- **Project name:** SQD for Iron-Sulfer Clusters
- **Question, problem, or opportunity:** How accurately can sample-based quantum diagonalization (SQD) reproduce the ground-state energy of [2Fe-2S(SMe)₄]²⁻?
- **Data, inputs, or evidence:** Model cluster geometry, a PySCF Hamiltonian, an AVAS-selected (10-electron, 10-orbital) Fe 3d active space, and ffsim/LUCJ samples.
- **Expected output:** A reproducible SQD workflow benchmarked against CASCI, with energy convergence, sampling efficiency, and orbital occupancies.
- **Tools and stack:** Python, PySCF, AVAS, ffsim/LUCJ, Qiskit, and qiskit-addon-sqd.
- **Team lead:** Amandeep Singh Bhatia [deepquantum88](https://github.com/deepquantum88)
- **Team members and roles:** [`project-management/team.md`](https://github.com/stjude-biohackathon/KIDS26-Team10/blob/main/project-management/team.md)
- **Communication:** [Slack: #team10](https://stjudebiohackathon.slack.com/archives/C0BSA3JQGS1)

## Vision and Mission

- **Vision:** Evaluate quantum-assisted approaches to iron–sulfur electronic structure.
- **Mission:** Run and benchmark an SQD pipeline against an exact classical reference in the same active space, targeting an energy error within 1.6 mHa.

## About

Iron–sulfur clusters are found in enzyme active sites. Their magnetic coupling and Fe–S covalency make their electronic structure challenging for single-reference methods.

Our workflow connects PySCF → AVAS → ffsim/LUCJ → qiskit-addon-sqd. We compare SQD with CASCI using the same Hamiltonian and fixed orbitals. CASCI provides an exact reference within the selected active space.

## Roadmap and Milestones

| When | Focus | Expected outcome |
| --- | --- | --- |
| Day 1 | Use completed (10-electron, 10-orbital) CASCI results as the reference and develop the SQD implementation | Classical reference available; SQD implementation underway |
| Day 2 | Compare SQD with CASCI; analyze convergence, valid-sample fractions, and orbital occupancies | Quantified energy error and benchmark plots |
| Day 3 | Finalize runs, document methods and limitations, and prepare the presentation | Reproducible workflow and results summary |




