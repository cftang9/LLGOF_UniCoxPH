# A Goodness-of-fit Test for Log-Linearity in Cox Proportional Hazard Model under Monotonic Covariate Effects #

This repository contains  ```R``` codes for the article “A Goodness-of-fit Test for Log-Linearity in Cox Proportional Hazard Model under Monotonic Covariate Effects.” 
We proposed a goodness-of-fit (GOF) test for the log-linearity versus monotonic effects on the hazard function in the Cox Proportional Hazard Model (PH) Model: 

$$
  H_0: \lambda_i(t) = \lambda_0(t) \exp(Z_i \beta) ~~ \mbox{versus} ~~ H_1: \lambda_i(t) = \lambda_0(t) \exp (g(Z_i)), 
$$

where $g$ is a monotonic function. 
Deviation of the partial likelihood ratios suggests exploring non-log-linear but monotonic effects, such as log-concave and log-convex, on the hazard rate. 
The critical values were determined by Monte Carlo random samples from the ordinary Cox Model with linear interpolated Breslow estimation of the baseline hazard function. We implemented two ```R``` functions, [StepFun_uniquify.R](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/StepFun_uniquify.R) and [Linear Inverse.R](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Linear_Inverse.R), to interpolate between consecutive jumps in the Breslow cumulated baseline estimation and then apply the inverse methods when generating Monte Carlo samples). 
In the end, some extensions and discussions beyond the paper were provided. 
<!--Lastly, we discuss the proposed GOF tests beyond the GOF tests. -->
<!-- This article has been submitted for publication. -->

<!-- Prior to using R programs on this repository, please download the main R program [EGJ_USO_Library.R](https://raw.githubusercontent.com/cftang9/MSUSO/master/EGJ_USO_Library.r).  -->

## Part 1. Reproducing size and power comparison results with a residual-based GOF test in the manuscript. 

Under the Cox PH model, 
We generate random samples from the Cox PH models, using the baseline hazard from the Exponential distribution with a mean of $1$, denoted by Exp(1). 
Under $H_0$, we consider $\beta = $
[0](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline%20Exp(1)/Sim_Uni_0Z_n50.R), 
[1](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline%20Exp(1)/Sim_Uni_1Z_n50.R), and 
[5](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline%20Exp(1)/Sim_Uni_5Z_n50.R).
Under $H_1$, we consider $g(Z) = $ 
[Z^2](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline%20Exp(1)/Sim_Uni_Z2_n50.R), 
[exp(Z)](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline%20Exp(1)/Sim_Uni_eZ_n50.R), 
[log(Z)](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline%20Exp(1)/Sim_Uni_lZ_n50.R), 
[6(Z)^(1/2)](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline%20Exp(1)/Sim_Uni_lZ_n50.R), and 
[sin(Z pi)/4](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline%20Exp(1)/Sim_Uni_lZ_n50.R), 

with sample size $n=50, 100, 200, 500$ and $1000$. 



## Part 2. Reproducing the assessments of the robustness of the GOF tests with baseline hazards. 

[<img src="https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/BHF.png" width="600" />](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/BHF.png)


[Baseline from Gompertz(1,2) distribution](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline%20G(1%2C2)/Sim_Uni_0Z_n50_Gompertz.R)

[Baseline from Gompertz(2,0.5) distribution](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/Baseline%20G(2%2C0.5)/Sim_Uni_0Z_n50_Gompertz.R)


<!-- ![BHF.png](https://github.com/cftang9/LLGOF_UniCoxPH/blob/main/BHF.png) -->


<!--
## Beyond this work: 

### Further discussion of censoring times: ties and high-censoring rates

### Further discussion of further discussion of the partial linear models
-->

## Reference: 
1. 
2. 
3. 

