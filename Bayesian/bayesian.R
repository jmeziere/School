library(Bolstad)
binodp(3, 8, seq(0, 1, .2), rep(1/6, 6))
# a. The conditional probabilities can be found in the n = 8 y = 3 section
# b. the likelihoods are under the column Likelihood
# c. The joint distribution is found by taking all possibilities of the observed quantity and the probability and calculating the binomial from that.
# d. The marginal distribution is found by taking all of the conditionals possible, multiplying by the marginal of x, and then adding them all up.
# e. The posteriors are found by normalizing the likelihood times the prior

# Question 6.2a
binodp(2, 7, seq(0, 1, .2),
       c(0.0000000, 0.2628337,0.4989733, 0.2217659,0.0164271,0.0000000))
# Question 6.2b
binodp(5, 15, seq(0, 1, .2), rep(1/6, 6))

# c. This shows that it doesn't matter if you do it in two separate steps. You'll get the same answer.