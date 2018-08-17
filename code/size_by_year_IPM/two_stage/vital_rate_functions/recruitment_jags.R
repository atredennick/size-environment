model{
  ##  LIKELIHOOD AND PROCESS MODELS
  for(i in 1:N){
    y[i] ~ dnegbin(q[i],theta)
    ynew[i] ~ dnegbin(q[i],theta) # simulate a new dataset for posterior checks
    q[i] <- theta/(theta+lambda[i])
    lambda[i] <- trueP1[i]*fecundity[i]
    log(fecundity[i]) <- intcpt.yr[year[i]] + intcpt.gr[Group[i]] + dd*trueP2[i] 
    trueP1[i] <- parents1[i]*u+parents2[i]*(1-u)
    trueP2[i] <- sqrt(trueP1[i])
  }
  
  ##  PRIOR MODELS
  # Intercepts
  u ~ dunif(0, 1)
  theta ~ dgamma(0.001, 0.001)
  g.tau ~ dgamma(2, 0.5)
  for(g in 1:Ngroups){
    intcpt.gr[g] ~ dnorm(0, g.tau)
  }
  intcpt.mu ~ dnorm(0, 0.001)
  intcpt.tau ~ dgamma(0.001, 0.001)
  for(m in 1:Nyrs){
    intcpt.yr[m] ~ dnorm(intcpt.mu, intcpt.tau)  
  }
  
  # Density-dependence
  dd ~ dnorm(0, 0.001)
  
  # Posterior predictive checks
  sd_y <- sd(y[])
  sd_ynew <- sd(ynew[]) 
  pvalue_sd <- step(sd_ynew - sd_y)
  
  mean_y <- mean(y[])
  mean_ynew <- mean(ynew[])
  pvalue_mean <- step(mean_ynew - mean_y)
}