# A Goodness-of-fit Test for Log-Linearity in Cox Proportional Hazard Model under Monotonic Covariate Effects #

This repository contains  ```R``` codes for the article “A Goodness-of-fit Test for Log-Linearity in Cox Proportional Hazard Model under Monotonic Covariate Effects.” 
We proposed a goodness-of-fit (GOF) for the linearity versus isotonic effects in the Cox PH Model. 
Deviation of the partial likelihood ratios suggests a non-linear effect. 
The critical values were determined by Monte Carlo random samples from the ordinary Cox Model with linear interpolated Breslow estimation of the baseline hazard function. We implemented two ```R``` functions, [StepFun_uniquify.R](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/StepFun_uniquify.R) and [Linear Inverse.R](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Linear_Inverse.R), to interpolate between consecutive jumps in the cumulated baseline and apply the inverse methods when generating Monte Carlo samples). 
In the end, some extensions and discussions beyond the paper were provided. 
<!--Lastly, we discuss the proposed GOF tests beyond the GOF tests. -->
<!-- This article has been submitted for publication. -->

<!-- Prior to using R programs on this repository, please download the main R program [EGJ_USO_Library.R](https://raw.githubusercontent.com/cftang9/MSUSO/master/EGJ_USO_Library.r).  -->

## Part 1. Reproducing size and power comparison results with a residual-based GOF test in the manuscript. 
[Directly from Baseline](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline Exp(1)/Sim_Uni_0Z_n50.R)


## Part 2. Reproducing the assessments of the robustness of the GOF tests with baseline hazards. 


<!--
## Beyond this work: 

### Further discussion of censoring times: ties and high-censoring rates

### Further discussion of further discussion of the partial linear models
-->

## Reference: 
1. 
2. 
3. 

