library(Bolstad)
library(Bolstad2)

# Linear Regression
bayes.lin.reg(dat3$y, dat3$x, "normal", "normal", 0 , 3, -20,3, 1, pred.x = 10)

# Acceptance/Rejection
p = runif(100000)
theta = ifelse(p < 1/2, log(2*p), -log(2*(1-p)))
yval1 = exp(-abs(theta))
yval2 = exp(-abs(theta)-(theta-0.5)^2/2)
u = c()
for(val in yval1)
  u = c(u,runif(1,0,val))
plot(theta[u < yval2],u[u < yval2],col="green",pch=20,ylim = c(0,1))
points(theta[u > yval2],u[u > yval2],col="blue",pch=20)
theta_keep = theta[u < yval2]
hist(theta_keep)

# MH
f <- function(x,y,z){
  (22.5-z)^(x-1)*(28-z)^(x-1)*exp(-(22.5+28-2*z)/y-x-y/2-z^2/8)/(gamma(x)*y^x)^2
} 

#how many replications we want
nrep <-20000
#If we graph the target density we will notice the variance is about 5 or so.
stdev <- 0.15
#this gives me empty vectors of size nrep
x <- rep(NA, nrep)
y <- rep(NA, nrep)
z <- rep(NA, nrep)

#start point (we should pick it within the range)
x[1] <- 5
y[1] <- 6
z[1] <- 5
for(i in 1:(nrep - 1)){
  #Pick my new point from a bivariate normal distribution. If we pick from two normals then (x,y) will be
  #bivariate normal if they are independent.
  # the number 1 indicates we want one point
  # we make sure the new distribution is centered at the old x
  print(i)
  x_prime <- rnorm(1, x[i], stdev)
  y_prime <- rnorm(1, y[i], stdev)
  z_prime <- rnorm(1, z[i], stdev)
  #pick a uniform random variable from U(0,1) (that is the default)
  u <- runif(1)
  
  # Define the ratio
  ratio <- f(x_prime, y_prime, z_prime)/f(x[i],y[i],z[i])
  
  #ifelse returns the first value if the condition is true and the second if it isn't
  x[i+1] <- ifelse(u > ratio, x[i], x_prime)
  y[i+1] <- ifelse(u > ratio, y[i], y_prime)
  z[i+1] <- ifelse(u > ratio, z[i], z_prime)
}

# Gibbs Sampling
#how many replications we want
nrep <-20000
#If we graph the target density we will notice the variance is about 5 or so.
stdev <- 0.15
#this gives me empty vectors of size nrep
x <- rep(NA, nrep)
y <- rep(NA, nrep)

#start point (we should pick it within the range)
x[1] <- 1
y[1] <- 1

x_vals <- c(0.0041,1.077,0.2393,0.0556,.5093,0.1252,0.5179,0.4699,1.1619,0.9864,1.025,1.8864,1.0819,0.0196,0.5420,0.0090,1.1341,0.3209,1.2337,0.4501)
for(i in 1:(nrep - 1)){
  #Pick my new point from a normal distribution centered at the old x
  # the number 1 indicates we want one point
  x[i+1] <- rgamma(1,shape=21,scale=1/(1+y[i]*sum(x_vals)))
  y[i+1] <- rgamma(1,shape=21,scale=1/(1+x[i+1]*sum(x_vals)))
}

mean(x[1001:nrep])
var(x[1001:nrep])
mean(y[1001:nrep])
var(y[1001:nrep])

#95% credible interval from resulting data
quantile(x, 0.025)
quantile(x, 0.975)
quantile(y, 0.025)
quantile(y, 0.975)

plot(x,y)
plot(x[1001:nrep],y[1001:nrep])
sum(ifelse(x[1001:nrep]<2 & y[1001:nrep] < 1.5,1,0))/(nrep-1000)
