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
