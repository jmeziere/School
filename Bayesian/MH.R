## Example done in class
library(ggplot2)
library(tidyverse); library(viridis); library(viridisLite); 
library(ggExtra)
# Question 1

# Part a
# x - alpha
# y- beta
f <- function(x,y){
  max(0,(cos(x)+sin(y))/exp((x/10)^2+(y/10)^2))
} 

#how many replications we want
nrep <-50000
#If we graph the target density we will notice the variance is about 5 or so.
stdev <- 2
#this gives me empty vectors of size nrep
x <- rep(NA, nrep)
y <- rep(NA, nrep)

#start point (we should pick it within the range)
x[1] <- 5
y[1] <- 6

for(i in 1:(nrep - 1)){
  #Pick my new point from a bivariate normal distribution. If we pick from two normals then (x,y) will be
  #bivariate normal if they are independent.
  # the number 1 indicates we want one point
  # we make sure the new distribution is centered at the old x
  x_prime <- rnorm(1, x[i], stdev)
  y_prime <- rnorm(1, y[i], stdev)
  
  #pick a uniform random variable from U(0,1) (that is the default)
  u <- runif(1)
  
  # Define the ratio
  ratio <- f(x_prime, y_prime)/f(x[i],y[i])
  
  #ifelse returns the first value if the condition is true and the second if it isn't
  x[i+1] <- ifelse(u > ratio, x[i], x_prime)
  y[i+1] <- ifelse(u > ratio, y[i], y_prime)
}

##Simple plot
plot(y[1001:nrep]~x[1001:nrep])

##More advanced plots:
df <- data.frame(x[1001:nrep] , y[1001:nrep])
# Scatterplot with Marginal Histograms
p <- ggplot(df, aes(x[1001:nrep], y[1001:nrep])) +
  geom_point(size = .1) +
  theme_bw()
ggMarginal(p, type = "histogram", color = "white", bins = 50)

## Does a heat map of the distribution showing high and low points
ggplot(df, aes(x[1001:nrep], y[1001:nrep])) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  labs(title = "Density Contour Plot") +
  scale_fill_viridis() +
  theme_bw()

# Part b
x[1001:nrep]
sum(ifelse(x[1001:nrep] < 2 & y[1001:nrep] < 5,1,0))/(nrep-1000)
# I found this value by counting all of the points with x < 2 and y < 5 and then dividing by the total number

# Question 2  

# Part b
# x - alpha
# y - beta
# z - theta
# Part a
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

# Store the data in a tibble
sim.dat <- tibble(x[1001:nrep],y[1001:nrep],z[1001:nrep])
# Create the plot
qplot(x[1001:nrep], y[1001:nrep], colour=z[1001:nrep], data=sim.dat) +
  scale_colour_gradient(low="blue", high="red") +
  geom_point(size = 0.1, stroke = 0, shape = 16)

# Store the data in a tibble
sim.dat <- tibble(y[1001:nrep],z[1001:nrep],x[1001:nrep])
# Create the plot
qplot(y[1001:nrep], z[1001:nrep], colour=x[1001:nrep], data=sim.dat) +
  scale_colour_gradient(low="blue", high="red") +
  geom_point(size = 0.1, stroke = 0, shape = 16)

# Store the data in a tibble
sim.dat <- tibble(z[1001:nrep],x[1001:nrep],y[1001:nrep])
# Create the plot
qplot(z[1001:nrep], x[1001:nrep], colour=y[1001:nrep], data=sim.dat) +
  scale_colour_gradient(low="blue", high="red") +
  geom_point(size = 0.1, stroke = 0, shape = 16)

# Part c
hist(x[1001:nrep])
hist(y[1001:nrep])
hist(z[1001:nrep])

# Part d
sum(ifelse(x[1001:nrep]<1.293,1,0))/(nrep-1000)
sum(ifelse(y[1001:nrep]<7.526,1,0))/(nrep-1000)
sum(ifelse(z[1001:nrep]<-2.111,1,0))/(nrep-1000)

# Part e
# It seems unwise to try to use less observations than parameters to estimate those parameters.