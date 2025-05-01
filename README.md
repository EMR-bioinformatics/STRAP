# STRAP Molecular Analysis

This repository contains scripts to support the publication: 'Deep molecular profiling of synovial biopsies in the STRAP trial identifies signatures predictive of treatment response to biologic therapies in Rheumatoid Arthritis'


## Repository Structure

```
.
├── 0_DataPreprocessing
    ├── 1_salmon_mapping.sh
    ├── 2_prepare_counts.R
    └── 3_exploratory.R
├── 1_DEGanalysis
    └── 1_DEGanalysis.R
├── 2_modularanalysis
    ├── 1_load_mods_qmod.R
    ├── 2_modular_analysis.R
    └── 3_qmod.R
├── 3_Machine_Learning_Models
    ├── STRAP_models_FINAL.R
    ├── STRAP_models_check_in_R4RA_final.R
    ├── STRAP_models_tables.R
    ├── clinical_models.R
    ├── comBat_norm.R
    ├── compare_mod.R
    ├── pooled_models.R
    └── create_pseudoRNAseq_coefficients_merged_github.R
├── 4_Subset_Module_Scores
    ├── 1_AMP_module_scores_github.R
    └── 2_AMP_module_plots_github.R
├── 5_common and 3-way DEGs
    ├── Figure2a_STRAP.DEGs.R
    ├── Figure2b_pathways.R
    ├── Figure2e_forestPlot.R
    └── Figures2cd_polarPlots.R
├── 6_molecular_Clusters
    ├── Figure4a_STRAP.clusters.R
    ├── Figure4b_R4RA.clusters.R
    ├── Figure4c_STRAP.clusters.pathways.R
    ├── Figure4d_VennDiagram.R
    ├── Figure4e_STRAP.PCAs.R
    └── SupplementaryTable2_STRAP.clusters.table.Rmd
```

## Data exploration website

The supporting website is available at: [https://strap.hpc.qmul.ac.uk/](https://strap.hpc.qmul.ac.uk/)

Raw RNAseq files: [E-MTAB-13733](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13733)
