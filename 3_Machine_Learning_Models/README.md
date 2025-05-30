**compare_mod.R** is the pipeline function which compares nestedcv models using repeated nested CV

**STRAP_models_FINAL.R** is used to generate the machine learning models predictive of drug response to etanercept, rituximab and tocilizumab in STRAP

**STRAP_models_tables.R** is to retrieve metrics table of the models

**comBat_norm.R** is the merging of R4RA and STRAP RNA-seq for validation

**STRAP_models_check_in_R4RA_final.R** is the validation of STRAP models in R4RA

**clinical_models.R** is used to generate models using clinical parameters alone

**pooled_models.R** is used to generate models using pooled STRAP and R4RA data

**create_pseudoRNAseq_coefficients_merged_github.R** is used to calculate predefined linear regression components (intercept and slope coefficient)for each gene. Using these models, nCounter assay data was converted to RNA-Seq scale (referred to as pseudo-RNA-Seq) based on stored regression models for each gene.
