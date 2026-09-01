# Simulation inputs and analysis for GPR15L–polySia and GPR15L–GPR15 molecular dynamics

Scripts, input files, and analysis notebooks accompanying:

> **Polysialylation of GPR15 regulates T cell trafficking to inflamed mucosa in ulcerative colitis**
> Makhsous S., et al. — see the manuscript's *Materials and Methods*,
> sections "Molecular dynamics simulations of GPR15L–polysialic acid interactions",
> "Dissociation free energy (DFE) calculations", and
> "Molecular simulations of GPR15L in the presence and absence of polysialylation".

This repository contains **exactly the input files and shell scripts used to run the
simulations** (recovered verbatim from the original run directories on the NTNU
triumvirate cluster and NMRbox), together with the Jupyter notebooks used for the
manuscript figures. No trajectories or other large binary outputs are included.

---

## Repository layout

```
GPR15_polySia_simulations/
├── aqueous_DP5_DP10_DP20/        GPR15L + polySia chains of increasing length
│   ├── DP5/   complex.pdb, POL_48SA_0SA.pdb, start.top/topol.top, index.ndx,
│   │          mdp/, top/, plumed.dat, run.sh, unbiased_run.sh, continue.sh
│   ├── DP10/  ... + POL_98SA_0SA.pdb, pp.sh
│   └── DP20/  ... + POL_198SA_0SA.pdb, pp.sh
├── bilayer_receptor/             GPR15±glycan embedded in POPC + GPR15L
│   ├── RL_complex/               unglycosylated receptor  (tleap.in, mdout.mdp, mdp/)
│   ├── core2_RL_complex/         T32 core-2 O-glycan      (tleap.in, mdout.mdp, mdp/)
│   ├── sia10_core2_RL_complex/   T32 core-2 glycan + 10 sialic acids (D10)
│   │                             (tleap.in, mdp/, prot.amb2gmx/, clean structure)
│   ├── free_ligand/              GPR15L alone in solution (reference, 1 µs)
│   └── clean_bilayer_RL_OLT.pdb  shared membrane+protein starting structure
├── system_preparation/
│   ├── doGlycans_polySia/        polySia chain construction (doGlycans)
│   ├── glycam06h_ff/             GLYCAM06(h) force-field files used with GROMACS
│   ├── prepare.sh                protein topology preparation (GROMACS)
│   ├── packmol.inp + memgen.sh   POPC bilayer packing (Packmol-MemGen route)
│   └── run_tleap_acpype.sh       Amber→GROMACS topology conversion
└── analysis/
    ├── gpr15lg_contacts.ipynb            Fig. 3B — GPR15L–polySia contact maps
    ├── analyze-fes.ipynb                 funnel-restrained FES inspection (development)
    ├── receptor_and_complexes_analysis.ipynb  Figs. 3C–J — DFE, conformational
    │                                     space (z, φ), Arg34–Asp148 distances,
    │                                     C-terminal RMSD vs cryo-EM 9WXM
    └── freq_/contacts_*                  saved contact outputs (DP5/DP10/DP20)
```

## System descriptions

| System | Contents | Production |
|---|---|---|
| `aqueous_DP5_DP10_DP20` | GPR15L (from Boltz/AlphaFold model) + α2,8-polySia DP5/DP10/DP20, rhombic dodecahedron, NaCl | Unbiased, 1.0 µs (625 ns DP20), Nosé–Hoover/Parrinello–Rahman @ 300 K |
| `aqueous_*/replicas` | Same complexes, metadynamics on COM distance | 13/15/22 replicas × 10 ns (DP5/DP10/DP20), σ=0.005 nm, W=0.42 kJ/mol, pace 45 steps, wall 7.3 nm |
| `bilayer_receptor/RL_complex` | Cryo-EM GPR15–GPR15L (PDB 9WXM) + AF N-terminus in POPC | 3 replicas × ≥300 ns analyzed |
| `bilayer_receptor/core2_RL_complex` | As above + T32 core-2 O-glycan | 3 replicas × ≥300 ns analyzed |
| `bilayer_receptor/sia10_core2_RL_complex` | As above + 10 sialic acids (D10) on the core-2 glycan | 3 replicas × ≥300 ns analyzed |
| `bilayer_receptor/free_ligand` | GPR15L alone | 1.0 µs |

## Software

- GROMACS 2022.3 + PLUMED 2.8.3 (`gmx_mpi`, CUDA builds on both clusters)
- AmberTools (tleap; ff14SB, GLYCAM06j-1/h, Lipid21, TIP3P)
- doGlycans (polySia GLYCAM06h ITP generation)
- Packmol-MemGen / Packmol (membrane assembly)
- ACPYPE + ParmEd (Amber↔GROMACS conversion)
- MDAnalysis (trajectory analysis; notebooks)


## Reproducing the simulations

1. **Aqueous GPR15L–polySia complexes**

   ```bash
   cd aqueous_DP5_DP10_DP20/DP5
   # 1) build solvated, neutralized system (see unbiased_run.sh header):
   #    editconf -> solvate -> genion -> emin.mdp
   # 2) one long unbiased reference trajectory:
   ./unbiased_run.sh          # slurm job; grompp mdp/{nvt,npt,md}.mdp then mdrun
   # 3) metadynamics replica ensemble:
   ./run.sh                   # regenerates plumed.dat, loops rep_i with seeded velocities;
                              # grompp nvt_seeded.mdp -> npt.mdp -> md_short.mdp (+plumed)
   ```

   `run.sh` writes `plumed.dat` itself and is self-contained apart from cluster
   environment variables (`$gmx`, `$plumed`, `PLUMED_KERNEL`) — adjust the paths at
   the top for your installation (the committed copies point at NTNU triumvirate /
   NMRbox installations). Seeds are derived from `BASE_SEED` to make the ensemble
   reproducible.

2. **Membrane systems** — rebuild path (as performed):

   ```
   system_preparation/packmol.inp  (packmol-memgen equivalent packing of
                                    POPC/water/ions around the protein)
   bilayer_receptor/<system>/tleap.in   (tleap: build glycan with χ/φ imposing,
                                        bond to Thr32, save prmtop/inpcrd)
   prot.amb2gmx/                        (acpype conversion; rungmx.sh notes)
   bilayer_receptor/<system>/mdp/       (ions -> minim -> nvt -> npt -> md)
   ```

   The committed `clean_bilayer_RL_OLT.pdb` is the solvated/ion-free lipid+protein
   structure; water and ions were added with `gmx solvate`/`genion` using
   `mdp/ions.mdp`.

3. **Free ligand in solution**

   ```bash
   cd bilayer_receptor/free_ligand
   gmx pdb2gmx -f lig.pdb      # as recorded in the typescript log
   # then mdp/{emin,nvt,npt,md}.mdp via freelig/run_nmrbox.sh (NMRbox SLURM)
   ```

## Reproducing the figures

Notebooks are executed top-to-bottom with the trajectories in place:

| Manuscript item | Notebook |
|---|---|
| Fig. 3B (contact frequency maps DP5/DP10/DP20) | `analysis/gpr15lg_contacts.ipynb` |
| Fig. 3C & S3A–C (DFE profiles, convergence, bootstrap) | `analysis/receptor_and_complexes_analysis.ipynb` ("Complexes" section) |
| Fig. 3F,G (conformational space z/φ) | same notebook, "Receptor simulations" section |
| Fig. 3H (representative bound pose; rendering) | PyMOL/VMD session (not scripted here) |
| Fig. 3I,J (C-terminal RMSD vs 9WXM) | same notebook |

Summary contact-frequency files (`analysis/freq_complex_polySia_DP*.txt`) record
residues with contact frequency > 0.85, e.g.
`select binding_site, resi 2+4+5+6+7+8+11+15+34+35+36` (DP5).


## Clusters

Simulations were executed on NTNU "triumvirate" (GROMACS 2022.3 + PLUMED 2.8.3,
CUDA) and on NMRbox (nmrbox.org). SBATCH headers in `run.sh` /
`run_nmrbox.sh` reflect those queues and serve as documentation rather than
portable submission scripts.

