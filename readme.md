# Molecule-dynamic-based Aging Clock and Aging Roadmap Forecast with ***Sundial***

This repository accompanies our manuscript "Molecule-dynamic-based Aging Clock and Aging Roadmap Forecast with ***Sundial***".

## Data preparation

If you want to reproduce the full analysis based on individual-level data, you will also need to apply for access to the following dataset:

- The UK Biobank (Sudlow *et al.* (2015), apply from https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access)

Data preprocess in the file folder: datapre. It contains all `R` scripts to extract the necessary data from UK Biobank.

## Supervisedclock

We calculated protein-predicted BA for the filtered cohort using five-fold cross-validation. Within each fold, a Light-GBM model was trained using the tuned hyperparameters and predicted age values were generated for the test set of that fold. We then combined the predicted age values from each of the folds to create a measure of BA for the entire sample. (Ludger *et al.* (2024), *Argentieri et al.* (2024))

## feature selection

Feature selection combined Boruta and recursive feature elimination (RFE) to identify proteins predictive of chronological age. Boruta used shadow features as noise benchmarks, iteratively removing features with lower SHAP values until only relevant ones remained. RFE then refined the selection using five-fold cross-validation on the UK Biobank dataset, removing the least important proteins based on SHAP values until five remained. The minimal subset for accurate prediction was determined to be 20 proteins, as reducing this number significantly impaired performance.

## sundialclock

step 1: knn to construct the graph

step 2: transition matirx to diffusion fied

step3: pseudo-age to pedicted biological age by distribution mapping

step4: visualization

## references

1. Sudlow, C., Gallacher, J., Allen, N., Beral, V., Burton, P., Danesh, J., Downey, P., Elliott, P., Green, J., Landray, M., et al. (2015). UK Biobank: An Open Access Resource for Identifying the Causes of a Wide Range of Complex Diseases of Middle and Old Age. PLoS Med 12, e1001779. https://doi.org/10.1371/JOURNAL.PMED.1001779
2. Ludger J.E. Goeminne, Anastasiya Vladimirova, Alec Eames, Alexander Tyshkovskiy, M. Austin Argentieri, Kejun Ying, Mahdi Moqri, Vadim N. Gladyshev (2024). Plasma protein-based organ-specific aging and mortality models unveil diseases as accelerated aging of organismal systems, *Cell Metabolism*, doi: https://doi.org/10.1016/j.cmet.2024.10.005.
3. Argentieri, M.A., Xiao, S., Bennett, D. *et al.* Proteomic aging clock predicts mortality and risk of common age-related diseases in diverse populations. *Nat Med*  **30** , 2450–2460 (2024). https://doi.org/10.1038/s41591-024-03164-7
