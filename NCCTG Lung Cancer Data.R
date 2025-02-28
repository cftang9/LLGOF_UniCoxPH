rm(list=ls(all=TRUE))
library(survival)
library(isoSurv)
library(timereg)
library(parallel)
source("Linear_Inverse.R")
source("StepFun_uniquify.R")

lung_univ = lung[,c(2,3,8)]
lung_univ = na.omit(lung_univ)
Data = lung_univ
set.seed(02172025)
#lung_univ$pat.karno = lung_univ$pat.karno + rnorm(length(lung_univ$pat.karno))/100
n = length(lung_univ$time)

hfit = coxph(Surv(time, status) ~ pat.karno, data= lung_univ)
hfit$coefficients 
ifit = isoph(Surv(time, status) ~ iso(pat.karno,shape="dec"), data= lung_univ)

plot((ifit$iso.cov)$pat.karno, (ifit$iso.cov)$psi.hat,type="s",ylim=c(-0.8,1.8))
abline(a=-hfit$coefficients*ifit$Zk ,b=hfit$coefficients, col="red")

lung_univ$status = lung_univ$status-1; 
Data = lung_univ
colnames(Data) = c("time","status","Z"); 
Data = Data[order(Data$Z),]
#Data$Z = Data$Z + rnorm(length(Data$Z))/100

iso.fit = isoph(Surv(time, status) ~ iso(Z,shape="dec"), data=Data, maxiter=50, eps = 10^(-3))
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

Tn = NLL.Iso - NLL.Cox; 

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
  
  iso.fit_b = isoph(Surv(time_b, status_b) ~ iso(Z,shape="dec"), data=Data_b, maxiter=50, eps = 10^(-3))
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

Test_bb = sapply(1:BB, do_bb)
Tn_b = unlist(Test_bb); quantile(Tn_b,probs=0.95) # 3.670195 

hist(Tn_b); abline(v=Tn); # 7.481391

fit <- cox.aalen(Surv(time, status)~prop(Z), data = Data, residuals = 1, rate.sim=0, n.sim=100)
mt_resids = cum.residuals(fit, data=Data, cum.resid = 1, n.sim = 100)
as.numeric(mt_resids$pval.test<0.05) #0.697>0.05

Data$Z2 = log(Data$Z)
hist(Data$Z2)
cox.fit2 = coxph(Surv(time, status) ~ Z2, timefix=F, data=Data);
summary(cox.fit2)

Data$Z3 = exp(Data$Z)
hist(Data$Z3)
cox.fit3 = coxph(Surv(time, status) ~ Z3, timefix=F, data=Data);
summary(cox.fit3)

