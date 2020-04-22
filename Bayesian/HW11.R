#Problem 1

# Sigma_0 < Sigma_1
# Sigma_0 < Sigma_2

install.packages("Bolstad2")
#Problem 2
library(Bolstad2)

##The below are arguments for norMixMH() function:
##normMixMH(theta0, theta1, p, candidate, steps = 100000, type = "rw")

#theta0 contains the mean and st deviation of the first component of the normal mixture N(1,1) here
theta0 <- c(1,sqrt(1))

#theta1 contains the mean and st deviation of the second component of the normal mixture N(6,9) here
theta1 <- c(6,3)

#p a value between 0 and 1 representing the mixture proportion, so that the true density is:
# p  f(mu1, sigma1) + (1-p) f(mu2, sigma2)
p <- 0.75

#candidate: a victor of length 2 containing the mean and st deviation of the candidate density
candidate <- c(0, 0.5)

#steps number of steps

# type = 'ind' or 'rw' depending on whether an independent candidate density or a random
#walk candidate density is to be used. We will only use random walk density here.

MCMCsampleRW <- normMixMH(theta0, theta1, p, candidate, steps = 50000, type = "rw")

summary(MCMCsampleRW)
var(MCMCsampleRW)
#Problem 3
library(Bolstad2)

##The below are arguments for norMixMH() function:
##normMixMH(theta0, theta1, p, candidate, steps = 100000, type = "rw")

#theta0 contains the mean and st deviation of the first component of the normal mixture N(1,1) here
theta0 <- c(-1,2)

#theta1 contains the mean and st deviation of the second component of the normal mixture N(6,9) here
theta1 <- c(4,1)

#p a value between 0 and 1 representing the mixture proportion, so that the true density is:
# p  f(mu1, sigma1) + (1-p) f(mu2, sigma2)
p <- 0.6

#candidate: a victor of length 2 containing the mean and st deviation of the candidate density
candidate <- c(0, 1.5)

#steps number of steps

# type = 'ind' or 'rw' depending on whether an independent candidate density or a random
#walk candidate density is to be used. We will only use random walk density here.

MCMCsampleRW <- normMixMH(theta0, theta1, p, candidate, steps = 50000, type = "rw")

#Problem 4
?bivnormMH
BVMH <- bivnormMH(0.3, steps = 50000)
summary(BVMH$targetSample)
cov(BVMH$targetSample)
length(BVMH$targetSample[length(BVMH$targetSample)-1000:length(BVMH$targetSample)])
plot(BVMH$targetSample$x[length(BVMH$targetSample$x)-1000:length(BVMH$targetSample$x)],BVMH$targetSample$y[length(BVMH$targetSample$y)-1000:length(BVMH$targetSample$y)])
