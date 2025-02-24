rm(list=ls(all=TRUE))
library(VGAM)
u = seq(0,2,by=0.01); 

dW = dweibull(u,shape=3,scale=1)/(1-pweibull(u,shape=3,scale=1))
dG1 = dgompertz(u, scale = 2, shape = 0.5, log = FALSE)/(1-pgompertz(u, scale = 2, shape = 0.5, log = FALSE))
dG2 = dgompertz(u, scale = 1, shape = 2, log = FALSE)/(1-pgompertz(u, scale = 1, shape = 2, log = FALSE))

par(mar=c(4.2,4.2,0.1,0.1))
plot(u,dG1,type="l",ylim=c(0,27),xlab="t",ylab=expression(lambda[0](t)),
     col="red",lwd=1)
lines(c(0,2),c(1,1),lty=1,lwd=1)
#lines(c(0,2),c(.5,.5),lty=1,lwd=1,col="green")
lines(u,dG2,col="blue",lty=1,lwd=1)
legend("topleft", 
       legend = c("Exp(1)","Gom(1,2)","Gom(2,.5)"),
       lty=c(1,1,1),
       col=c("black","blue","red")
)