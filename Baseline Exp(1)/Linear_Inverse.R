Linear_Inverse = function(y,x,fx){
  ### x and fx must be strictly positive and increasing 
  # x = seq(1,10); 
  # fx = x+runif(10); #fx = sort(fx); 
  # y = x + rexp(10); 
  # 
  ind = sapply(y, function(y){sum(y>fx)})
  x0 = c(0,x); fx0 = c(0,fx); 
  ind0 = ind + 1; 
  
  m = ((fx0[ind0+1] - fx0[ind0])/(x0[ind0+1] - x0[ind0]))
  m[is.na(m)] = Inf; 
  
  nx0 = length(x0); 
  Linear_Inverse = ((y-fx0[ind0])/m + x0[ind0]) + (ind0==nx0)*runif(length(ind0),min=0,max=0.0001)
  # Linear_Inverse
  # diff(Linear_Inverse)
  # plot(c(0,x),c(0,fx),type="l",xlim=c(0,11),ylim=c(0,11))
  # points(c(0,x),c(0,fx),col="blue")
  # points(Linear_Inverse,y)
  return(Linear_Inverse)
}