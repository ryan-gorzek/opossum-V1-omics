
# Comparative transcriptomics of opossum and mouse V1

Code repository accompanying
**Gorzek & Trachtenberg (2025)**
*Comparative transcriptomics reveals differences in cortical cell type organization between metatherian and eutherian mammals*
**bioRxiv**: [https://doi.org/10.1101/2025.10.03.680397](https://doi.org/10.1101/2025.10.03.680397)

This repository contains the analysis code used to generate all main and supplementary figures in the manuscript, spanning single-nucleus RNA-seq, spatial transcriptomics, and immunofluorescence data.

## Status

Code upload is in progress.

All analysis notebooks for main and supplementary figures will be deposited in this repository. Core utilities used across figures are already available in the [comparatome](https://github.com/ryan-gorzek/comparatome) R package.

## Data Availability

All raw and processed data will be made publicly available upon publication:

* **snRNA-seq**: GEO **GSE299387**
* **Stereo-seq**: GEO **GSE299386**
* **Immunofluorescence images**: Figshare (link will be provided here)

## Repository Structure

```
preprocessing/      snRNA-seq and Stereo-seq preprocessing pipelines
figure 1/           Atlas overview
figure 2/           Cross-species IT mapping and archetypes
figure 3/           Spatial transcriptomics
figure 4/           Immunofluorescence
figure S1-S4/       Supplementary figures
comparatome/        R package used throughout the project
```

Each `figure *` directory contains self-contained notebooks or scripts corresponding directly to panel(s) in the manuscript.

## Software

Analyses were performed in the following environments:

* **R** ≥ 4.3.0

  * Seurat 4.3.0
* **Python** ≥ 3.9

  * used for spatial region selection (preprocessing) and spatial archetype analysis (figure 3)
* **MATLAB** r2023a

  * used for ParTI archetype analysis

Some analyses were performed with custom R tools:

* **Package**: [comparatome](https://github.com/ryan-gorzek/comparatome)

---

## Citation

If you use this code or data, please cite the preprint listed above.
