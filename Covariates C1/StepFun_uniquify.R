StepFun_uniquify = function(x,y){
  # This uniquification is based on positive values with smallest value zero. 
  Ind = (abs(diff(c(0,y))<=10^(-9)))
  uni_x = x[Ind==0]; uni_y = y[Ind==0]
  uni_x = c(0,uni_x); uni_y = c(0,uni_y)
  return(list(uni_x=uni_x, uni_y=uni_y))
}

# 
# x = 1:13
# y = c(1,2,2,2,2,2,3,4,5,5,7,7,7); 
# W = StepFun_uniquify(x,y)
# plot(x,y)
# points(W$uni_x,W$uni_y,col="red")