source("multihypothesis.R")

  # Bayesian setting
H     = 500  # horizon
th    = c(0.3, 0.5, 0.7)  # hypotheses
thgam = th  # Bayesian problem
gam   = rep(1, 3) / 3  # least informative a priori distribution
lam   = c(100,100,100)  # Lagrange multipliers
l     = array(lam, dim = c(3, 3))  # convert in a matrix
test  = OptTest(H, l, th, gam, thgam)  # design optimal Bayesian test
list(setting = "Bayes", ErrorProbabilities = sapply(1:3, function(i) 1 - PAccept(test, th[i], i)),  # calculate error probabilities
     ESS=sapply(th, function(x) ESS(test, x)), MaxSteps=maxNumber(test))  # calculate ESS
     
  # modified Kiefer-Weiss setting
H     = 500  # horizon
th    = c(0.3, 0.5, 0.7)  # hypotheses
thgam = c(0.4,0.6)  # modified Kiefer-Weiss problem
gam   = rep(1, 2) / 2  # equal weights
lam   = c(100,100,100)  # Lagrange multipliers
l     = array(lam, dim = c(3, 3))  # convert in a matrix
test  = OptTest(H, l, th, gam, thgam)  # use modified Kiefer-Weiss setting
list(setting = "modified Kiefer-Weiss", ErrorProbabilities = sapply(1:3, function(i) 1 - PAccept(test, th[i], i)),  # calculate error probabilities
     ESS=sapply(th, function(x) ESS(test, x)), MaxSteps=maxNumber(test))  # calculate ESS

  # details in https://doi.org/10.48550/arXiv.2212.05151
