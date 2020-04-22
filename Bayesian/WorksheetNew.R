#how many replications we want
nrep <-100000
#If we graph the target density we will notice the variance is about 5 or so.
stdev <- 0.15
#this gives me empty vectors of size nrep
x <- rep(NA, nrep)
y <- rep(NA, nrep)

#start point (we should pick it within the range)
x[1] <- 10
y[1] <- 10
mux <- 6
muy <- 36
sigx <- 3
sigy <- 1.6
rho <- -0.6
for(i in 1:(nrep - 1)){
  #Pick my new point from a normal distribution centered at the old x
  # the number 1 indicates we want one point
  x[i+1] <- rnorm(1, mux+(sigx/sigy)*rho*(y[i]-muy), sqrt((1-rho^2)*sigx^2))
  y[i+1] <- rnorm(1, muy+(sigy/sigx)*rho*(x[i+1]-mux), sqrt((1-rho^2)*sigy^2))
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
sum(ifelse(x[1001:nrep]<5 & y[1001:nrep] < 40,1,0))/(nrep-1000)
