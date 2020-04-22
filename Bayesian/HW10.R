# Problem 1

# Part a

x=runif(1000000,0,14)
theta=sin(runif(1000000,0,pi/2))*14
2/(sum(ifelse(x < theta,1,0))/length(x))

# Problem 2

# Part d

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

# Part e

mean(theta_keep)
sd(theta_keep)^2

# Part f

sum(ifelse(u < yval2,1,0))/100000


# Problem 3

# Part b

x = rnorm(10000,4,3)
yval1 = dnorm(x,4,3)*6.8
yval2 = 0.75*exp(-(x-3)^2/2)+ 0.25*exp(-(x-6)^2/2)
u = c()
for(val in yval1)
  u = c(u,runif(1,0,val))
plot(x[u < yval2],u[u < yval2],col="green",pch=20,ylim = c(0,1))
points(x[u > yval2],u[u > yval2],col="blue",pch=20)
x_keep = x[u < yval2]
hist(x_keep)

# Part c

mean(x_keep)
sd(x_keep)^2

# Part d

sum(ifelse(u < yval2,1,0))/10000

# Problem 4

# Part c 

x = rt(10000,1)+0.4
yval1 = dt(x,1)*pi*exp(-0.4^2/2)
yval2 = exp(-abs(x))*exp(-(x-0.4)^2/2)
u = c()
for(val in yval1)
  u = c(u,runif(1,0,val))
plot(x[u < yval2],u[u < yval2],col="green",pch=20,ylim = c(0,1))
points(x[u > yval2],u[u > yval2],col="blue",pch=20)
x_keep = x[u < yval2]
hist(x_keep)

# Part d

mean(x_keep)
sd(x_keep)^2

# Part e

sum(ifelse(u < yval2,1,0))/10000
