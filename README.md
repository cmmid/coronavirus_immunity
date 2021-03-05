# coronavirus_immunity
 code for paper - fitting seasonal coronavirus and simulating covid


All the fitting, analysis and figure generation can be run from "master_script.R"


Package versions are: 

coda_0.19-3
doParallel_1.0.15
iterators_1.0.12
foreach_1.5.0
here_0.1
gridExtra_2.3
ggplot2_3.3.2
tmvtnorm_1.4-10
gmm_1.6-5
sandwich_2.5-1 
Matrix_1.2-18
mvtnorm_1.1-1 
Rcpp_1.0.6
rootSolve_1.8.2.1
deSolve_1.28
MASS_7.3-51.6
data.table_1.13.4


All analysis except for the fitting was done in R version 4.0.0 on macOS Catalina 10.15.7. The fitting was done on R 3.4 on AWS ec2 machines, running a linux AMI. On AWS the package mvtnorm was version 1.0.8.

