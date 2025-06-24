# Multiple Imputation and ANCOVA Analysis in R

This repository provides R functions and scripts for performing multiple imputation (MI) and ANCOVA with Rubin's Rule pooling in longitudinal clinical trial data. The codebase supports the simulation of dummy data with user-specified missingness, flexible implementation of MI in two steps (MCMC and monotone regression), and downstream pooled analysis.

## Repository Structure

- `function/`  
  Contains core R functions for the multiple imputation workflow.

- `00generate_dummy.R`  
  Script to generate dummy longitudinal trial data, define missing data mechanisms, and prepare the gold standard complete dataset.

- `01imputation.R`  
  Script to perform multiple imputation on the dummy dataset using the custom functions in the `function/` folder.

- `02ANCOVA.R`  
  Script to conduct ANCOVA analysis on imputed datasets and use Rubin's Rule pooling to obtain final inference.

## Workflow

1. **Data Simulation**  
   Run `00generate_dummy.R` to create a fully observed dummy dataset, and introduce user-defined or randomized missing values to simulate realistic trial data scenarios.

2. **Multiple Imputation**  
   Run `01imputation.R` to apply the two-step imputation process (MCMC and monotone regression), using the core functions from `function/`.

3. **Analysis**  
   Run `02ANCOVA.R` to perform ANCOVA for each imputed dataset, pool results using Rubin’s Rule, and generate tables or plots for inference and method comparison.

## Suggested Requirements

- R (>= 4.4)
- Suggested packages: `renv`

## Getting Started

Clone the repository and run the code to install the needed packages:

```r
renv::restore()
```

And then you can try to run the R scripts by order.

## Citation
If you use or adapt this code, please cite the corresponding paper or source.

## License
This project is licensed under the MIT License.

## Acknowledgments
Special thanks to all contributors and reviewers involved in the development and validation of these methods.