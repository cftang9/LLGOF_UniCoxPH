rm(list=ls(all=TRUE))
library(survival)
library(isoSurv)
library(timereg)
library(parallel)
source("Linear_Inverse.R")
source("StepFun_uniquify.R")

n = 50; 

gfunction = function(z){
  gz = log(z);
  return(gz)}
lower_unif_censoring = 0; 
upper_unif_censoring = 1.5; 

Rn_power = 0; 

B = 500; 
Tn = array(,B); 
start.time = Sys.time()
for(b in 1:B){
  set.seed(02022025000+b)
  Z = array(runif(n,min=0.1,max=2),c(n,1)); 
  C = array(runif(n,lower_unif_censoring,upper_unif_censoring),c(n,1)); 
  gZ = gfunction(Z); 
  Exp_i = exp(gZ); 
  
  U = runif(n); 
  time_0 = log(1+2*(-log(U)/Exp_i))/2; 
  time = pmin(time_0,C); 
  status = array(1,c(n,1)); status[C<time_0] = 0; 
  Data = data.frame(cbind(time,status,Z)); 
  colnames(Data) = c("time","status","Z"); 
  Data = Data[order(Data$Z),]
  
  fit <- cox.aalen(Surv(time, status)~prop(Z), data = Data, residuals = 1, rate.sim=0, n.sim=100)
  mt_resids = cum.residuals(fit, data=Data, cum.resid = 1, n.sim = 100)
  Rn_power = as.numeric(mt_resids$pval.test<0.05)/B + Rn_power
  
  interval.time = (Sys.time()-start.time)/b
  if(b%%1==0){
    print(b)
    print(interval.time)
    print(Sys.time()+interval.time*(B-b))
    print(Rn_power/b*B)
  }
}

Rn_power 

# 1000 # 0.974
#  500 # 0.676
#  200 # 0.312
#  100 # 0.128
#   50 # 0.060