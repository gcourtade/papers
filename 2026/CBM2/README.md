# Decoding CBM2 binding modes: how similar folds achieve distinct substrate preferences
Data and scripts used to analyse data for the manuscript.

Authors: Susanne Hansen Troøyen, Andrew Philip Rennison, Tina Dahlgren, Peter Fojan, Marie Sofie Møller, Gaston Courtade


Citation:
```


``` 


# Repository layout

```text
├── 2026
│   ├── CBM2
│   │   ├── AFM                                                AFM data/scripts
│   │   │   ├── AFM_raw_data                                   Raw data
│   │   │   ├── Plot_AFM.ipynb                                 Plotting AFM height profiles
│   │   │   └── ...                                           Gwyddion projects (.gwy), images, height profiles (.txt)
│   │   ├── Alignment                                          Sequence alignment of ScCBM2 and BaCBM2 (Clustal Omega)
│   │   │   ├── alignment.png, alignment.svg                   Alignment from Jalview
│   │   │   ├── CBM_conservation.ipynb                         Analysis of sequence alignment
│   │   │   └── ...                                           Clustal Omega output files
│   │   ├── Binding_assays                                     Pull-down assay data/scripts
│   │   │   ├── Cellulose_chitin_assays
│   │   │   │   ├── data
│   │   │   │   ├── plots
│   │   │   │   └── Cellulose_chitin_assays.ipynb               Analysis of cellulose/chitin pull-down assays
│   │   │   ├── PET_assays
│   │   │   │   ├── data
│   │   │   │   ├── plots
│   │   │   │   └── PET_assays.ipynb                            Analysis of PET pull-down assays
│   │   │   └── pulldown_assays_raw_data.xlsx                  Raw plate reader data
│   │   ├── Bioinformatic_analysis                             MSA analysis of specific/promiscuous CBM2s (see *)
│   │   │   ├── Bioinformatics_fasta
│   │   │   │   ├── per_substrate_files                         FASTA sequences for CBM2s binding each substrate
│   │   │   │   │   └── Clustal_omega_alignments                Aligned sequences and sequence logos
│   │   │   │   └── specific_promiscuous_files                  FASTA sequences for substrate-specific CBM2s
│   │   │   │       └── Clustal_omega_alignments                Aligned sequences and sequence logos
│   │   │   ├── CBM2_bioinformatics.ipynb                       Filtering binders and generating FASTA files
│   │   │   └── Supplementarytable1_Rennison_highthroughput.xlsx  HTP data from Rennison et al.
│   │   ├── Figures                                            All figures
│   │   │   └── SI_Figures                                     Supplementary figures
│   │   ├── NMR
│   │   │   ├── CCPN_CSP
│   │   │   │   ├── Data                                       Chemical shift data exported from CCPN
│   │   │   │   └── CCPN_project.zip                           CCPN project and NMR spectra for BaCBM2 + PET/BHET
│   │   │   ├── DynamicsCenter_data                            Relaxation data
│   │   │   │   ├── 15N
│   │   │   │   ├── BaCBM2
│   │   │   │   └── ScCBM2_free
│   │   │   ├── files                                          Output files from NMR_plots_relax_PRE.ipynb
│   │   │   │   ├── BaCBM2 BHET.txt
│   │   │   │   └── BaCBM2 PET.txt
│   │   │   ├── Plots
│   │   │   ├── pymol_changes                                  PyMOL .pml coloring scripts from NMR_plots_relax_PRE.ipynb
│   │   │   ├── pymol_csp_i_lw                                 PyMOL .pml coloring scripts from NMR_plots_csp.ipynb
│   │   │   ├── pymol_relax                                    PyMOL .pml coloring scripts from NMR_plots_relax_PRE.ipynb
│   │   │   ├── NMR_plots_csp.ipynb                             Analysis of CSP data
│   │   │   └── NMR_plots_relax_PRE.ipynb                       Analysis of relaxation/PRE data
│   │   ├── Particle_size                                      Particle size determination
│   │   │   ├── Particle_size_substrates.xlsx                  Raw laser diffraction data
│   │   │   └── particle_size.ipynb                            Analysis of particle size
│   │   ├── Simulations                                        MD simulation data
│   │   │   └── Simulations                                    RDFs from five blocks for each variant
│   │   │       ├── BaCBM2
│   │   │       ├── ScCBM2
│   │   │       ├── ScCBM2_double                              ScCBM2-E10T-H69W
│   │   │       ├── ScCBM2_E272T                               ScCBM2-E10T
│   │   │       ├── ScCBM2_H331W                               ScCBM2-H69W
│   │   │       ├── analyze_rdf.ipynb                          Analysis of RDFs
│   │   │       └── ...                                       Plots
│   │   ├── XRD                                                XRD data/scripts
│   │   │   ├── TOPAS                                          TOPAS project files (.pro)
│   │   │   ├── Visualization_crystals                         PyMOL sessions for cellulose crystals
│   │   │   ├── XRD_processed_data                             XRD profiles (.txt), coordinates (.str), plots (.svg)
│   │   │   │   ├── Rietveldt_method.ipynb                     Plotting XRD profiles
│   │   │   │   └── Segal_method.ipynb                         Plotting Segal method XRD profiles
│   │   │   └── XRD_raw_data                                   Raw XRD data (.raw, .brml)
│   │   ├── PyMOL_sessions.zip                                 PyMOL sessions for figures (see **)
│   │   └── README.md
```

```text
* MSA analysis CBM2 sequences:
Binders to each substrate was selected from high-throughput data from
Rennison et al. https://doi.org/10.1002/cssc.202500468,
and filtered using CBM2_bioinformatics.ipynb.

** PyMOL_sessions.zip contains: 
CBM2_structure_align.pse          Structural alignment of ScCBM2, BaCBM2, BaCBM2-t
Electrostatic_potential.pse       Surface electrostatics
fold_bacbm2_model_0.cif           BaCBM2 AlphaFold model
fold_sccbm2_model_0.cif           ScCBM2 AlphaFold model
NMR_1H-R1R2.pse                   1H relaxation rates
NMR_15N-R1.pse                    15N relaxation rates
NMR_BaCBM2_BHET_changes.pse       Relaxation and PRE changes for BaCBM2 + BHET
NMR_BaCBM2_csp_i_lw.pse           CSP, intensity and linewidth changes for BaCBM2 + PET/BHET
NMR_BaCBM2_PET_changes.pse        Relaxation and PRE changes for BaCBM2 + PET
NMR_PRE.pse                       sPRE data
S_binding_plane.py                PyMOL script for calculating binding surface area
S_protein_footprint.py            PyMOL script for calculating full protein footprint
SASA_footprints.pse               Visualization of surface areas and protein footprints
Surface_area_calc_results.txt     Surface area calculation results
Water_density_BaCBM2.pse         Water density visualization for BaCBM2
Water_density_ScCBM2.pse         Water density visualization for ScCBM2
```
