# Part a
# My function is a uniform distribution from 0 to 1 scaled by a factor of 2

# Part b
xval = runif(1000000)
yval = dbeta(xval,3,2)
u = runif(1000000,0,2)
sum(ifelse(u < yval & xval < 0.25, 1, 0))/sum(ifelse(yval[u < yval],1,0))

# Part c
plot(xval[u < yval],u[u < yval],col="green",pch=20)
points(xval[u > yval],u[u > yval],col="blue",pch=20)
  sum(ifelse(u < yval,1,0))/1000000

quantile(xval[u < yval],c(0.025,0.975))

# Part d
pbeta(0.25,3,2)
qbeta(c(0.025,0.975),3,2)

