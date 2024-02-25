

pmf <- function(n, th){
  # joint probability of Bernoully trials as a function of the sufficient statistic s=\sum_{i=1}^n x_i
  s = seq(0, n)
  return(dbinom(s, n, th))
}


llRatio <- function(n, th1, th2){
  # log likelihood ratio
  log(pmf(n, th1) / pmf(n, th2))
}


OptTest <- function(H, lam, th, gam, thgam){
  # computation of the optimal test in the Bayes/modified Kiefer-Weiss problem
  # lam matrix of Lagrange Multipliers 
  # gam vector of weights
  # thgam vector of \vartheta the weights are applied to
  # th hypothesized values
  # H horizon

  k = length(th)
  m = length(thgam)
  lagr = list()
  accept = list()
  cont = list()
  z = sapply(th, function(x) pmf(H, x))
  for(i in 1:k){
    s0 = rep(0, H + 1)
    for(j in 1:k){
      if(j == i)
        next
      s0 = s0 + lam[j, i] * z[,j]
    }
    if(i == 1)
      mn = s0
    else
      mn = pmin(mn, s0)
    accept[[H]] = ifelse(mn >= s0, rep(i, H + 1), accept[[H]])
  }

  lagr[[H]] = mn
  cont[[H]] = rep(FALSE, H + 1)
  if(H > 1)
    for(n in (H-1):1){
      z = sapply(th, function(x) pmf(n, x))
      for(i in 1:k){
        s0 = rep(0, n + 1)
        for(j in 1:k){
          if(j == i)
            next
          s0 = s0 + lam[j, i] * z[,j]
        }
        if(i == 1)
          mn = s0
        else
          mn = pmin(mn, s0)

        accept[[n]] = ifelse(mn >= s0, rep(i, n + 1), accept[[n]])
      }
      x = seq(0, n + 1)
      z = sapply(thgam, function(x) pmf(n, x))
      s1 = rep(0, n + 1)

      for(i in 1:m){
        s1 = s1 + gam[i] * z[,i]
      }

      lagr[[n]] = pmin(
        s1 + head(lagr[[n + 1]] / (n + 1) * (n + 1 - x), n + 1) + tail(lagr[[n + 1]]/(n + 1) * (x),n + 1),
        mn
      )
      cont[[n]] = lagr[[n]] < mn
    }

  # returns the optimal test 
  return(rbind(cont, accept))
}


DBCTest<-function(H, lam, th, gam, thgam){
  fn0 <- function(l) gam[l] * pmf(n,thgam[l])
  fn1 <- function(j, l) lam[l, j] * pmf(n,th[l])

  k = length(th)
  m = length(gam)
  cont = list()
  accept = list()
  n = 1
  repeat{
    cont[[n]] = rep(TRUE, n + 1)
    accept[[n]] = rep(0, n + 1)
    s0 = fn0(1)
    if (m > 1)
      for(l in 2:m) s0 = s0 + fn0(l)
    for(i in 1:k){    
      s1 = rep(0, n + 1)
      for (j in 1:k){
        if (j == i) 
          next
        s1 = s1 + fn1(i, j)
      }
      cont[[n]] = cont[[n]] & s1 > s0
      if (i == 1){ 
        mn = s1 
        accept[[n]] = rep(1, n + 1)
      } 
      else{
         accept[[n]] = ifelse(mn > s1, rep(i, n + 1), accept[[n]])
         mn = pmin(mn, s1)
      }
    }
      if( n>= H)break
  if (all(!cont[[n]])) 
    break
  n=n+1
  }
  cont[[n]]=rep(FALSE,n+1)
  test=rbind(cont,accept)
  return(test)
}


MSPRT <- function(H, lA, th){

  cont = list()
  accept = list()
  N = length(th)
  for(n in 1:H){
    a = rep(0, n + 1)

    for(i in 1:N){
      b = rep(TRUE, n + 1)
      for(j in 1:N){
        if (j == i) next
        b = b & llRatio(n, th[i], th[j]) > lA[j]
      }
      a = a + ifelse(b, i, 0)
    }
    accept[[n]] = a
    cont[[n]] = accept[[n]] == 0
  }

  return(rbind(cont, accept))
}


maxNumber <- function(test){

  H = length(test[1,])
  n = 1
  repeat{
    if(all(!test[1,][[n]]))
      break
    if(n == H)
      break
    n = n + 1
  }

  return(n)
}


PAccept <- function(test, th, i){
  # probability to accept hypothesis i given th
  # probability to accept nothing if i==0

  cont = test[1,]
  accept = test[2,]

  H = length(cont)
  a = list()
  a[[H]] = ifelse(accept[[H]] == i, pmf(H, th), 0)
  if(H > 1)
  for(n in (H - 1):1){
    x = seq(0, n + 1)
    h = head(a[[n + 1]]/(n + 1) * (n + 1 - x), n + 1)
    t = tail(a[[n + 1]]/(n + 1) * (x), n + 1)

    a[[n]] = ifelse(cont[[n]], h + t, ifelse(accept[[n]] == i, pmf(n, th), 0))
  }

  return(head(a[[1]], 1) + tail(a[[1]], 1))
}


ESS <- function(test, th){
  # ESS of a test at point th

  cont = test[1,]
  H = length(cont)
  ess = list()
  ess[[H]] = rep(0, H + 1)
  if(H == 1)
    return(1)
  for(n in (H-1):1){
    x = seq(0, n + 1)
    ess[[n]] = ifelse(
      cont[[n]],
      (pmf(n, th) + head(ess[[n + 1]] / (n + 1) * (n + 1 - x), n + 1) + tail(ess[[n + 1]] / (n + 1) * (x), n + 1)),
      0
    )
  }

  return(1 + head(ess[[1]], 1) + tail(ess[[1]], 1))
}


prob_to_stop_after <- function(test, th,k){
  # probability to stop after stage k
  # given the true value of parameter th

  cont = test[1,]
  H = length(cont)
  if(k>=H)return(0)
  ess = list()
  ess[[k]] =ifelse(cont[[k]], pmf(k, th), 0)
  if(k <= 0) 
    return(1)
  if(k > 1)
  for(n in (k-1):1){
    x = seq(0, n + 1)
    ess[[n]] = ifelse(
      cont[[n]],
      ( head(ess[[n + 1]] / (n + 1) * (n + 1 - x), n + 1) + tail(ess[[n + 1]] / (n + 1) * (x), n + 1)),
      0
    )
  }

  return( head(ess[[1]], 1) + tail(ess[[1]], 1))
}


monte_carlo_simulation <- function(K, test, hyp, nMC) {
  # test simulation
  # hyp = true success probability
  # returns rates of acceptations
  # standard error of the rates
  # the ESS and its standard errors
  # K number of hypotheses
  # nMC number of replications for Monte Carlo

  cont = test[1,]
  accept = test[2,]
  H = length(cont)
  ss = 0
  totaccepted = rep(0,K)
  totn = 0

  for(i in 1 : nMC) {
    s = 0 # accumulated sum for test run
    n = 0 # number of steps in current run
    accepted = rep(0, K) # number of accepting in the run (0 or 1)

    for (stage in 1:H){
      # run starts
      # generate
      summand = rbinom(1, 1, hyp) # 1 bernoulli
      s = s + summand # accumulated sum
      n = n + 1   # one step more

      if(stage == H){
        # the last stage of the run
        accepted[accept[[stage]][s+1 ]] = accepted[accept[[stage]][s+1 ]] + 1
      }
      else{
        if(cont[[stage]][s +1] == FALSE){
          # stop by the optimal stopping rule
          accepted[accept[[stage]][s + 1]] = accepted[accept[[stage]][s + 1]] + 1
          break # accepted or rejected, stop the run; no more stages
        }
      }

    }
    totaccepted = totaccepted + accepted
    totn = totn + n
    ss = ss + n^2

  }
  nrep = as.double(nMC)
  OC = (totaccepted) / nrep
  seOC = sqrt(OC * (1 - OC) / nrep)
  ESS = totn / nrep
  sdESS = sqrt((ss - totn^2 / nrep) / nrep)

  return(list(OC=OC, SEOC=seOC, ESS=ESS, SEESS=sdESS/sqrt(nrep)))
}
