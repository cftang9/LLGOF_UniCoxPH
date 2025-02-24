rm(list=ls(all=TRUE))
library(survival)
library(isoSurv)
library(parallel)
source("Linear_Inverse.R")
source("StepFun_uniquify.R")
n = 50; 

gfunction = function(z){
  gz = 1*z;
  return(gz)}
lower_unif_censoring = 0; 
upper_unif_censoring = 0.82; 

PLRT.power = 0; 

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
  
  iso.fit = isoph(Surv(time, status) ~ iso(Z,shape="inc"), data=Data, maxiter=50, eps = 10^(-3))
  Data$psi.hat = iso.fit$iso.cov$psi.hat; 
  Data$Iso.Exp_i = exp(Data$psi.hat)
  Data$Iso.Exp_i[Data$Iso.Exp_i==0] = 10^(-8)
  Data = Data[order(Data$time,decreasing=T),]
  Data$Iso.cum.Exp_i = cumsum(Data$Iso.Exp_i)
  
  cox.fit = coxph(Surv(time, status) ~ Z, timefix=F, data=Data);
  Data$Cox.Exp_i = exp(Data$Z*max(cox.fit$coefficients,0))
  Data = Data[order(Data$time,decreasing=T),]
  Data$Cox.cum.Exp_i = cumsum(Data$Cox.Exp_i)
  Data = Data[order(Data$time,decreasing=F),]
  Data$Cox.Brewlow = cumsum(Data$status/Data$Cox.cum.Exp_i)
  Cum_BL_Haz = StepFun_uniquify(Data$time,Data$Cox.Brewlow)
  Cen_KM = survfit(Surv(time, 1-status) ~ 1, data = Data)
  Cum_Cen_Haz = StepFun_uniquify(Cen_KM$time,1-Cen_KM$surv)
  
  NLL.Cox = sum(Data$status*(log(Data$Cox.Exp_i)-log(Data$Cox.cum.Exp_i)))
  NLL.Iso = sum(Data$status*(log(Data$Iso.Exp_i)-log(Data$Iso.cum.Exp_i)))
  Tn[b] = NLL.Iso - NLL.Cox; 
  
  ###bootstrapping
  BB = 500;
  Tn_b = array(,BB);
  
  do_bb = function(bb){
    set.seed(02082025000+bb)
    time_0b = Linear_Inverse(y=-log(runif(n))/Data$Cox.Exp_i,x=Cum_BL_Haz$uni_x,fx=Cum_BL_Haz$uni_y)
    C_b = Linear_Inverse(y=runif(n),x=Cum_Cen_Haz$uni_x,fx=Cum_Cen_Haz$uni_y)
    C_b = C_b*Data$status + Data$time*(1-Data$status);
    time_b = pmin(time_0b,C_b);
    status_b = array(1,c(n,1)); status_b[C_b<time_0b] = 0;
    
    Data_b = data.frame(cbind(time_b,status_b,Data$Z));
    colnames(Data_b) = c("time_b","status_b","Z");
    Data_b = Data_b[order(Data_b$Z),]
    
    iso.fit_b = isoph(Surv(time_b, status_b) ~ iso(Z,shape="inc"), data=Data_b, maxiter=50, eps = 10^(-3))
    Data_b$psi.hat = iso.fit_b$iso.cov$psi.hat;
    Data_b$Iso.Exp_i = exp(Data_b$psi.hat)
    Data_b$Iso.Exp_i[Data_b$Iso.Exp_i==0] = 10^(-8)
    Data_b = Data_b[order(Data_b$time,decreasing=T),]
    Data_b$Iso.cum.Exp_i = cumsum(Data_b$Iso.Exp_i)
    
    cox.fit_b = coxph(Surv(time_b, status_b) ~ Z, timefix=F, data=Data_b);
    Data_b$Cox.Exp_i = exp(Data_b$Z*max(cox.fit_b$coefficients,0))
    Data_b = Data_b[order(Data_b$time,decreasing=T),]
    Data_b$Cox.cum.Exp_i = cumsum(Data_b$Cox.Exp_i)
    
    NLL.Cox_b = sum(Data_b$status_b*(log(Data_b$Cox.Exp_i)-log(Data_b$Cox.cum.Exp_i)))
    NLL.Iso_b = sum(Data_b$status_b*(log(Data_b$Iso.Exp_i)-log(Data_b$Iso.cum.Exp_i)))
    Tn_b = NLL.Iso_b - NLL.Cox_b; Tn_b;
    
    return(Tn_b)
  }
  Test_bb = mclapply(1:BB, do_bb, mc.cores = 16)
  Tn_b = unlist(Test_bb)
  
  PLRT.power = as.numeric(Tn[b] > quantile(p=0.95,Tn_b))/B + PLRT.power;
  
  interval.time = (Sys.time()-start.time)/b
  if(b%%1==0){
    print(b)
    print(interval.time)
    print(Sys.time()+interval.time*(B-b))
    print(c(Tn[b],quantile(prob=0.95,Tn_b)))
    print(PLRT.power/b*B)
  }
}

PLRT.power #

write.csv(PLRT.power, "Sim_Uni_1Z_n50_Gompertz.csv")