# Mathematical modelling of cytomegalovirus vaccination scenarios in high seroprevalence populations highlights the importance of immunizing infants
This repository accompanies the article “Mathematical modeling of cytomegalovirus vaccination scenarios in high-seroprevalence populations highlights the importance of immunizing infants.”

**R Version**

The original simulations were run with R 4.2.2. To avoid compatibility issues, we recommend installing R 4.2.2.

**Reproducing the Parameter Fit**

Run cedar_lenormand.r, which contains three additional scripts:

1. contact_matrix_use_continuous_progress.r
Required R packages: Hmisc, parallel, MASS, plyr, dplyr, ggplot2, pomp, reshape2, tidyverse, tidyr, matrixStats, extdplyr

Required data files:

ave_mom_child.txt — mean number of mother–infant contacts

breast_feeding_exponential.txt — contains best-fit parameters resulting from fitting an exponential curve to the breastfeeding data

contact_edited3BRA.txt — time-dependent contact intensities

toilet_training_gompertz.txt — contains best-fit parameters resulting from fitting a Gompertz curve to the potty-training data

NHANES_CMV_antibodies_combined_weighted.csv — seroprevalence data (scenario-specific: All-of-Brazil, South Brazil, North Brazil)

2. EasyABC_internal_modified.r
Additional package: EasyABC

3. ABC_sequential.r
   
**Visualization**

Use readin_cmv_results.r to visualize the results (e.g., parameter distributions).

**Vaccine Scenario Simulations**

natural_vaccine.r — simulates a vaccine mimicking natural immunity

sterilizing_vaccine.r — simulates a lifelong sterilizing vaccine

**Visualization - vaccine scenarios**

To reproduce the graphs for the vaccine simulations, use readin_vaccine_new.r script.

**Running on HPCs**

Use cmv_cluster.sh to run cedar_lenormand script.

For vaccine scenario simulations, use vaccine_cluster.sh and specify which script to run (natural_vaccine.r or sterilizing_vaccine.r).

