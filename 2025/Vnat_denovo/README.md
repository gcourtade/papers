# Development of a *Vibrio natriegens*-based plate-clearing assay for rapid screening of PET-hydrolyzing enzymes  
This GitHub repository contains the raw data and scripts for the following paper. Please cite accordingly if you use any of the code or data in this repository.

Bolstad I. M., Troøyen S. H., Luciano D., Petersen E., Courtade G. (2025) Development of a *Vibrio natriegens*-based plate-clearing assay for rapid screening of PET-hydrolyzing enzymes. TBD

## Contents
[**All_proteins**](https://github.com/gcourtade/masters/tree/main/2024/IBM_denovo_PET/All_proteins) contains .pdb files of all _de novo_ proteins in this project. These are output-files from [RFdiffusion](https://doi.org/10.1038/s41586-023-06415-8) + [ProteinMPNN](https://www.science.org/doi/10.1126/science.add2187).

[**CD**](https://github.com/gcourtade/masters/tree/main/2024/IBM_denovo_PET/CD) contains code and raw data from the circular dichroism analysis. The notebook _CD\_secondarystructure.ipynb_ plots the CD spectrum of B22 and additionally makes a pie chart of the secondary structure deconvolution of the spectrum performed by [CDNN 2.1](https://pubmed.ncbi.nlm.nih.gov/1409538/). The notebook _CD\_thermostability.ipynb_ plots the melt curve of B22 and finds an estimate of the apparent melting temperature.

[**NMR**](https://github.com/gcourtade/masters/tree/main/2024/IBM_denovo_PET/NMR) contains raw data and code from the NMR-analyses. The jupyter notebook _Plot\_1H\_NMRspectra.ipynb_ plots the 1H-NMR spectrum of B22 with and without BHET. The folders _Timeresolved\_BHET_ and _Timeresolved\_PET_ contain raw data and code to plot the time-resolved NMR with BHET and PET as substrate, respectively. The code in _Timeresolved\_PET.ipynb_ is modified from Hellesnes et al. [(1)](https://pubs.acs.org/doi/epdf/10.1021/acs.biochem.2c00619), and so is the raw data for FsC and FsC-L182A (_137\_pH6p5\_TSP.txt_ and _148\_pH6p5\_L182A\_TSP.txt_).

[**Plate\_clearance**](https://github.com/gcourtade/masters/tree/main/2024/IBM_denovo_PET/Plate_clearance) contains data for the quantitative plate-clearing assay. The notebook _plate\_clearance.ipynb_ plots the measurements from _All\_plate\_diameters.csv_ as ratios of clearance zone diameter to colony diameters for each timepoint, as well as a barplot for the final timepoint (72h).

[**Thermofluor**](https://github.com/gcourtade/masters/tree/main/2024/IBM_denovo_PET/Thermofluor) contains data for the thermofluor assay. The notebook _Tm\_thermofluor.ipynb_ plots the raw data from _Measured\_temp\_fluorescence.xlsx_ and fits a Boltzmann sigmoidal curve to the data, and finds Tm as the inflection point.
