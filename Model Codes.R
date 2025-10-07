# File includes functions relating to the nimble model 

#################### dCT ######################################
# dCt is a nimble function returning the density function of
# the uncensored CT value
# rCT is a nimble function returning a randomly generated 
# uncensored CT value

# x is the uncensored CT value
# mu is the expected value of CT
# sd is the variance of CT
# sd_cont is sigma_c from Model
# type is an indicator (1 = contaminated replicate,
# 2 = not contaminated or inhibited, 3 = inhibited)

dCt <- nimbleFunction(run = function(x = double(0), 
                                     mu = double(0),
                                     sd = double(0),
                                     sd_cont = double(0),
                                     type = double(0),
                                     log = integer(0, default = 0)) {
  returnType(double(0))
  
  if (type == 1){ # if contaminated
    
    if (x > mu){ # expect x to be less than mu
      dct <- 0
    } else {
      xt <- abs(x - mu)
      dct <- sqrt(2) * exp(-0.5 * (xt/sd_cont)^2) / (sd_cont*sqrt(pi))
    }
    
  } else if (type == 2){ # if good replicate
    
    xt <- x - mu
    dct <- exp(-0.5 * (xt/sd)^2) / (sd*sqrt(2*pi))
    
  } else { # if inhibited
    
    if (x < mu){ # expect x to be greater than mu
      dct <- 0
    } else {
      xt <- abs(x - mu)
      dct <- sqrt(2) * exp(-0.5 * (xt/sd_cont)^2) / (sd_cont*sqrt(pi))
    }
    
  }
  
  if (log){
    dct <- log(dct)
  }
  
  return(dct)
  
})

rCt <- nimbleFunction(run = function(n = integer(0), 
                                     mu = double(0),
                                     sd = double(0),
                                     sd_cont = double(0),
                                     type = double(0)) {
  returnType(double(0))
  
  if(n != 1) print("rCt only allows n = 1; using n = 1.")
  
  if (type == 1){ # if contaminated
    
    cont.amount <- rnorm(1, mean = 0, sd = sd_cont)
    x <- mu - abs(cont.amount)
    
  } else if (type == 2){ # if good replicate
    
    x <- rnorm(1, mean = mu, sd = sd)
    
  } else { # if inhibited
    
    cont.amount <- rnorm(1, mean = 0, sd = sd_cont)
    x <- mu + abs(cont.amount)
    
  }
  
  return(x)
  
})

##################### d_exp_uniform ###########################
# dbetab0 is a nimble function returning the density function of
# the exponential-uniform density function for betab0
# rbetab0 is a nimble function returning a randomly generated value
# with density dbetab0

# x is the log-DNA concentration
# b is the maximum DNA concentration for the uniform distribution

dbetab0 <- nimbleFunction(run = function(x = double(0),
                                         b = double(0),
                                         log = integer(0, default = 0)){
  
  # exp(betab0) should be uniform(0, b), so d(betab0) <- 1/b * exp(betab0) 
  returnType(double(0))
  
  dbetab0 <- exp(x)/b
  
  return(dbetab0)
})

rbetab0 <- nimbleFunction(run = function(n = integer(0),
                                         b = double(0)){
  # random generating function for dbetab0
  
  if(n!=1){
    stop('rbetab0 is only set up for n = 1')
  }
  
  returnType(double(0))
  
  rbetab0 <- log(runif(1, 0, b))
  
  return(rbetab0)
  
})

##################### nimble models ###########################
# Full_Code is Model 1
# Sigmay_Code is Model 2 (constant plate variance)
# Nocont_Code is Model 3 (ignore contamination and inhibition)

Full_Code <- nimbleCode({
  
  ############### Priors:
  tau2 ~ dinvgamma(a_sigma.tau, scale = b_sigma.tau) # default is shape, rate
  tau2.1 ~ dinvgamma(a_sigma.tau1, scale = b_sigma.tau1)
  sigma2 ~ dinvgamma(a_sigma, scale = b_sigma)
  
  sd_rho2 ~ dinvgamma(a_sigma.rho, scale = b_sigma.rho)
  rho0 ~ dnorm(1,1)
  
  if (nT != 1){
    if (num_sites == 1){
      rho ~ dnorm(rho0, sd_rho)
    } else {
      for (i in 1:num_sites){
        rho[i] ~ dnorm(rho0, sd_rho)
      }
    }
  }
  
  tau <- sqrt(tau2)
  tau.1 <- sqrt(tau2.1)
  sigma <- sqrt(sigma2)
  sd_rho <- sqrt(sd_rho2)
  
  a ~ dnorm(a0, sd_a)
  b ~ dnorm(b0, sd_b)
  
  #betab0 ~ dnorm(mean_betab0, sd_betab0)
  betab0 ~ dbetab0(b.max)
  
  if (numP == 1){
    alpha1 ~ dnorm(alpha10, sd = sigma_alpha)
    alpha2 ~ dnorm(alpha20, sd = sigma_alpha)
  } else {
    for (i in 1:numP){
      alpha1[i] ~ dnorm(alpha10, sd = sigma_alpha)
      alpha2[i] ~ dnorm(alpha20, sd = sigma_alpha)
    }
  }
  
  if (ncovb != 0){
    if (ncovb == 1){
      betab ~ dnorm(0,1)
    }
    else {
      for(i in 1:ncovb){
        betab[i] ~ dnorm(0,1)
      }
    }
  }
  
  if (ncovw != 0){
    if (ncovw == 1){
      betaw ~ dnorm(0,1)
    }
    else {
      for(i in 1:ncovw){
        betaw[i] ~ dnorm(0,1)
      }
    }
  }
  
  ################# Latent Parameters:
  if (nT != 1){
    if (ncovb == 0){
      if (num_sites == 1){
        l[1, 1] ~ dnorm(betab0, tau.1)
        for (j in 2:nT){
          l[1,j] ~ dnorm(rho*l[1,j-1], var = tau2)
        }
      } else {
        for (i in 1:num_sites){
          l[i,1] ~ dnorm(betab0, tau.1)
          for (j in 2:nT){
            l[i,j] ~ dnorm(rho[i]*l[i,j-1], var = tau2)
          }
        }
      }
    } else {
      if (num_sites == 1){
        l1.mu[1] <- inprod(Xb[1,1,], betab[])
        l[1,1] ~ dnorm(betab0 + l1.mu[1], tau.1)
        for (j in 2:nT){
          l.mu[1,j-1] <- inprod(Xb[1, j,] - rho*Xb[1, j-1,], betab[])
          l[1,j] ~ dnorm(rho*l[1,j-1] + l.mu[1,j-1], var = tau2)
        }
      } else {
        for (i in 1:num_sites){
          l1.mu[i] <- inprod(Xb[i,1,], betab[])
          l[i,1] ~ dnorm(betab0 + l1.mu[i], tau.1)
          for (j in 2:nT){
            l.mu[i,j-1] <- inprod(Xb[i, j,] - rho[i]*Xb[i, j-1,], betab[])
            l[i,j] ~ dnorm(rho[i]*l[i,j-1] + l.mu[i,j-1], var = tau2)
          }
        }
      }
    }
  } else {
    if (ncovb == 0){
      if (num_sites == 1){
        l[1, 1] ~ dnorm(betab0, tau.1)
      } else {
        for (i in 1:num_sites){
          l[i,1] ~ dnorm(betab0, tau.1)
        }
      }
    } else {
      if (num_sites == 1){
        l1.mu[1] <- inprod(Xb[1,1,], betab[])
        l[1,1] ~ dnorm(betab0 + l1.mu[1], tau.1)
      } else {
        for (i in 1:num_sites){
          l1.mu[i] <- inprod(Xb[i,1,], betab[])
          l[i,1] ~ dnorm(betab0 + l1.mu[i], tau.1)
        }
      }
    }
  }
  
  if (ncovw == 0){
    for (i in 1:num_samples){
      v[i] ~ dnorm(l[id_site_l[id_site[i]], id_site_time[id_site[i]]], var = sigma2)
    }
  } else {
    for (i in 1:num_samples){
      Xw.betaw[i] <- inprod(Xw[i,], betaw[])
      v[i] ~ dnorm(l[id_site_l[id_site[i]], id_site_time[id_site[i]]] + Xw.betaw[i], var = sigma2)
    }
  }
  
  for (i in 1:num_samples){
    w[i] <- exp(v[i]) 
  }
  
  # nimble chapter 5 on censorship
  # what dinterval distribution does
  
  # Lab contamination per replciate
  pi.type[1:3] ~ ddirch(alpha = alpha[1:3])
  
  constraint_data ~ dconstraint(pi.type[2] > pi.type[1] + pi.type[3])
  
  for (i in 1:num_replicates){
    type[i] ~ dcat(pi.type[1:3])
  }
  for (i in 1:num_replicates_star){
    type_star[i] ~ dcat(pi.type[1:3])
  }
  
  #################### Observations
  if (numP == 1){
    for (i in 1:num_replicates){
      mu[i] <- alpha1 + alpha2*log(w[id_sample[i]])
      sd[i] <- sqrt(exp(a + b*log(w[id_sample[i]])))
      Ct[i] ~ dCt(mu[i], sd[i], sd_cont, type[i])
      delta_inv[i] ~ dinterval(Ct[i], CT.max) 
    }
    for (i in 1:num_replicates_star){
      mu_star[i] <- alpha1 + alpha2*log(w_star[i])
      sd_star[i] <- sqrt(exp(a + b*log(w_star[i])))
      Ct_star[i] ~ dCt(mu_star[i], sd_star[i], sd_cont, type_star[i])
      delta_star_inv[i] ~ dinterval(Ct_star[i], CT.max)
    }
  } else {
    for (i in 1:num_replicates){
      mu[i] <- alpha1[P[i]] + alpha2[P[i]]*log(w[id_sample[i]])
      sd[i] <- sqrt(exp(a + b*log(w[id_sample[i]])))
      Ct[i] ~ dCt(mu[i], sd[i], sd_cont, type[i])
      delta_inv[i] ~ dinterval(Ct[i], CT.max)
    }
    for (i in 1:num_replicates_star){
      mu_star[i] <- alpha1[P_standards[i]] + alpha2[P_standards[i]]*log(w_star[i])
      sd_star[i] <- sqrt(exp(a + b*log(w_star[i])))
      Ct_star[i] ~ dCt(mu_star[i], sd_star[i], sd_cont, type_star[i])
      delta_star_inv[i] ~ dinterval(Ct_star[i], CT.max)
    }
  }
  
})

Sigmay_Code <- nimbleCode({
  
  ############### Priors:
  tau2 ~ dinvgamma(a_sigma.tau, scale = b_sigma.tau) # default is shape, rate
  tau2.1 ~ dinvgamma(a_sigma.tau1, scale = b_sigma.tau1)
  sigma2 ~ dinvgamma(a_sigma, scale = b_sigma)
  
  sd_rho2 ~ dinvgamma(a_sigma.rho, scale = b_sigma.rho)
  rho0 ~ dnorm(1,1)
  
  if (nT != 1){
    if (num_sites == 1){
      rho ~ dnorm(rho0, sd_rho)
    } else {
      for (i in 1:num_sites){
        rho[i] ~ dnorm(rho0, sd_rho)
      }
    }
  }
  
  tau <- sqrt(tau2)
  tau.1 <- sqrt(tau2.1)
  sigma <- sqrt(sigma2)
  sd_rho <- sqrt(sd_rho2)
  
  #betab0 ~ dnorm(mean_betab0, sd_betab0)
  betab0 ~ dbetab0(b.max)
  
  if (numP == 1){
    alpha1 ~ dnorm(alpha10, sd = sigma_alpha)
    alpha2 ~ dnorm(alpha20, sd = sigma_alpha)
    sigma.y2 ~ dinvgamma(a_sigma_y, scale = b_sigma_y)
    sigma.y <- sqrt(sigma.y2)
  } else {
    for (i in 1:numP){
      alpha1[i] ~ dnorm(alpha10, sd = sigma_alpha)
      alpha2[i] ~ dnorm(alpha20, sd = sigma_alpha)
      sigma.y2[i] ~ dinvgamma(a_sigma_y, scale = b_sigma_y)
      sigma.y[i] <- sqrt(sigma.y2[i])
    }
  }
  
  if (ncovb != 0){
    if (ncovb == 1){
      betab ~ dnorm(0,1)
    }
    else {
      for(i in 1:ncovb){
        betab[i] ~ dnorm(0,1)
      }
    }
  }
  
  if (ncovw != 0){
    if (ncovw == 1){
      betaw ~ dnorm(0,1)
    }
    else {
      for(i in 1:ncovw){
        betaw[i] ~ dnorm(0,1)
      }
    }
  }
  
  ################# Latent Parameters:
  if (nT != 1){
    if (ncovb == 0){
      if (num_sites == 1){
        l[1, 1] ~ dnorm(betab0, tau.1)
        for (j in 2:nT){
          l[1,j] ~ dnorm(rho*l[1,j-1], var = tau2)
        }
      } else {
        for (i in 1:num_sites){
          l[i,1] ~ dnorm(betab0, tau.1)
          for (j in 2:nT){
            l[i,j] ~ dnorm(rho[i]*l[i,j-1], var = tau2)
          }
        }
      }
    } else {
      if (num_sites == 1){
        l1.mu[1] <- inprod(Xb[1,1,], betab[])
        l[1,1] ~ dnorm(betab0 + l1.mu[1], tau.1)
        for (j in 2:nT){
          l.mu[1,j-1] <- inprod(Xb[1, j,] - rho*Xb[1, j-1,], betab[])
          l[1,j] ~ dnorm(rho*l[1,j-1] + l.mu[1,j-1], var = tau2)
        }
      } else {
        for (i in 1:num_sites){
          l1.mu[i] <- inprod(Xb[i,1,], betab[])
          l[i,1] ~ dnorm(betab0 + l1.mu[i], tau.1)
          for (j in 2:nT){
            l.mu[i,j-1] <- inprod(Xb[i, j,] - rho[i]*Xb[i, j-1,], betab[])
            l[i,j] ~ dnorm(rho[i]*l[i,j-1] + l.mu[i,j-1], var = tau2)
          }
        }
      }
    }
  } else {
    if (ncovb == 0){
      if (num_sites == 1){
        l[1, 1] ~ dnorm(betab0, tau.1)
      } else {
        for (i in 1:num_sites){
          l[i,1] ~ dnorm(betab0, tau.1)
        }
      }
    } else {
      if (num_sites == 1){
        l1.mu[1] <- inprod(Xb[1,1,], betab[])
        l[1,1] ~ dnorm(betab0 + l1.mu[1], tau.1)
      } else {
        for (i in 1:num_sites){
          l1.mu[i] <- inprod(Xb[i,1,], betab[])
          l[i,1] ~ dnorm(betab0 + l1.mu[i], tau.1)
        }
      }
    }
  }
  
  if (ncovw == 0){
    for (i in 1:num_samples){
      v[i] ~ dnorm(l[id_site_l[id_site[i]], id_site_time[id_site[i]]], var = sigma2)
    }
  } else {
    for (i in 1:num_samples){
      Xw.betaw[i] <- inprod(Xw[i,], betaw[])
      v[i] ~ dnorm(l[id_site_l[id_site[i]], id_site_time[id_site[i]]] + Xw.betaw[i], var = sigma2)
    }
  }
  
  for (i in 1:num_samples){
    w[i] <- exp(v[i]) 
  }
  
  # nimble chapter 5 on censorship
  # what dinterval distribution does
  
  # Lab contamination per replciate
  pi.type[1:3] ~ ddirch(alpha = alpha[1:3])
  
  constraint_data ~ dconstraint(pi.type[2] > pi.type[1] + pi.type[3])
  
  for (i in 1:num_replicates){
    type[i] ~ dcat(pi.type[1:3])
  }
  for (i in 1:num_replicates_star){
    type_star[i] ~ dcat(pi.type[1:3])
  }
  
  #################### Observations
  if (numP == 1){
    for (i in 1:num_replicates){
      mu[i] <- alpha1 + alpha2*log(w[id_sample[i]])
      sd[i] <- sigma.y
      Ct[i] ~ dCt(mu[i], sd[i], sd_cont, type[i])
      delta_inv[i] ~ dinterval(Ct[i], CT.max) 
    }
    for (i in 1:num_replicates_star){
      mu_star[i] <- alpha1 + alpha2*log(w_star[i])
      sd_star[i] <- sigma.y
      Ct_star[i] ~ dCt(mu_star[i], sd_star[i], sd_cont, type_star[i])
      delta_star_inv[i] ~ dinterval(Ct_star[i], CT.max)
    }
  } else {
    for (i in 1:num_replicates){
      mu[i] <- alpha1[P[i]] + alpha2[P[i]]*log(w[id_sample[i]])
      sd[i] <- sigma.y[P[i]]
      Ct[i] ~ dCt(mu[i], sd[i], sd_cont, type[i])
      delta_inv[i] ~ dinterval(Ct[i], CT.max)
    }
    for (i in 1:num_replicates_star){
      mu_star[i] <- alpha1[P_standards[i]] + alpha2[P_standards[i]]*log(w_star[i])
      sd_star[i] <- sigma.y[P_standards[i]]
      Ct_star[i] ~ dCt(mu_star[i], sd_star[i], sd_cont, type_star[i])
      delta_star_inv[i] ~ dinterval(Ct_star[i], CT.max)
    }
  }
  
})

Nocont_Code <- nimbleCode({
  
  ############### Priors:
  tau2 ~ dinvgamma(a_sigma.tau, scale = b_sigma.tau) # default is shape, rate
  tau2.1 ~ dinvgamma(a_sigma.tau1, scale = b_sigma.tau1)
  sigma2 ~ dinvgamma(a_sigma, scale = b_sigma)
  
  sd_rho2 ~ dinvgamma(a_sigma.rho, scale = b_sigma.rho)
  rho0 ~ dnorm(1,1)
  
  if (nT != 1){
    if (num_sites == 1){
      rho ~ dnorm(rho0, sd_rho)
    } else {
      for (i in 1:num_sites){
        rho[i] ~ dnorm(rho0, sd_rho)
      }
    }
  }
  
  tau <- sqrt(tau2)
  tau.1 <- sqrt(tau2.1)
  sigma <- sqrt(sigma2)
  sd_rho <- sqrt(sd_rho2)
  
  a ~ dnorm(a0, sd_a)
  b ~ dnorm(b0, sd_b)
  
  #betab0 ~ dnorm(mean_betab0, sd_betab0)
  betab0 ~ dbetab0(b.max)
  
  if (numP == 1){
    alpha1 ~ dnorm(alpha10, sd = sigma_alpha)
    alpha2 ~ dnorm(alpha20, sd = sigma_alpha)
  } else {
    for (i in 1:numP){
      alpha1[i] ~ dnorm(alpha10, sd = sigma_alpha)
      alpha2[i] ~ dnorm(alpha20, sd = sigma_alpha)
    }
  }
  
  if (ncovb != 0){
    if (ncovb == 1){
      betab ~ dnorm(0,1)
    }
    else {
      for(i in 1:ncovb){
        betab[i] ~ dnorm(0,1)
      }
    }
  }
  
  if (ncovw != 0){
    if (ncovw == 1){
      betaw ~ dnorm(0,1)
    }
    else {
      for(i in 1:ncovw){
        betaw[i] ~ dnorm(0,1)
      }
    }
  }
  
  ################# Latent Parameters:
  if (nT != 1){
    if (ncovb == 0){
      if (num_sites == 1){
        l[1, 1] ~ dnorm(betab0, tau.1)
        for (j in 2:nT){
          l[1,j] ~ dnorm(rho*l[1,j-1], var = tau2)
        }
      } else {
        for (i in 1:num_sites){
          l[i,1] ~ dnorm(betab0, tau.1)
          for (j in 2:nT){
            l[i,j] ~ dnorm(rho[i]*l[i,j-1], var = tau2)
          }
        }
      }
    } else {
      if (num_sites == 1){
        l1.mu[1] <- inprod(Xb[1,1,], betab[])
        l[1,1] ~ dnorm(betab0 + l1.mu[1], tau.1)
        for (j in 2:nT){
          l.mu[1,j-1] <- inprod(Xb[1, j,] - rho*Xb[1, j-1,], betab[])
          l[1,j] ~ dnorm(rho*l[1,j-1] + l.mu[1,j-1], var = tau2)
        }
      } else {
        for (i in 1:num_sites){
          l1.mu[i] <- inprod(Xb[i,1,], betab[])
          l[i,1] ~ dnorm(betab0 + l1.mu[i], tau.1)
          for (j in 2:nT){
            l.mu[i,j-1] <- inprod(Xb[i, j,] - rho[i]*Xb[i, j-1,], betab[])
            l[i,j] ~ dnorm(rho[i]*l[i,j-1] + l.mu[i,j-1], var = tau2)
          }
        }
      }
    }
  } else {
    if (ncovb == 0){
      if (num_sites == 1){
        l[1, 1] ~ dnorm(betab0, tau.1)
      } else {
        for (i in 1:num_sites){
          l[i,1] ~ dnorm(betab0, tau.1)
        }
      }
    } else {
      if (num_sites == 1){
        l1.mu[1] <- inprod(Xb[1,1,], betab[])
        l[1,1] ~ dnorm(betab0 + l1.mu[1], tau.1)
      } else {
        for (i in 1:num_sites){
          l1.mu[i] <- inprod(Xb[i,1,], betab[])
          l[i,1] ~ dnorm(betab0 + l1.mu[i], tau.1)
        }
      }
    }
  }
  
  if (ncovw == 0){
    for (i in 1:num_samples){
      v[i] ~ dnorm(l[id_site_l[id_site[i]], id_site_time[id_site[i]]], var = sigma2)
    }
  } else {
    for (i in 1:num_samples){
      Xw.betaw[i] <- inprod(Xw[i,], betaw[])
      v[i] ~ dnorm(l[id_site_l[id_site[i]], id_site_time[id_site[i]]] + Xw.betaw[i], var = sigma2)
    }
  }
  
  for (i in 1:num_samples){
    w[i] <- exp(v[i]) 
  }
  
  #################### Observations
  if (numP == 1){
    for (i in 1:num_replicates){
      mu[i] <- alpha1 + alpha2*log(w[id_sample[i]])
      sd[i] <- sqrt(exp(a + b*log(w[id_sample[i]])))
      Ct[i] ~ dnorm(mu[i], sd[i])
      delta_inv[i] ~ dinterval(Ct[i], CT.max) 
    }
    for (i in 1:num_replicates_star){
      mu_star[i] <- alpha1 + alpha2*log(w_star[i])
      sd_star[i] <- sqrt(exp(a + b*log(w_star[i])))
      Ct_star[i] ~ dnorm(mu_star[i], sd_star[i])
      delta_star_inv[i] ~ dinterval(Ct_star[i], CT.max)
    }
  } else {
    for (i in 1:num_replicates){
      mu[i] <- alpha1[P[i]] + alpha2[P[i]]*log(w[id_sample[i]])
      sd[i] <- sqrt(exp(a + b*log(w[id_sample[i]]))) 
      Ct[i] ~ dnorm(mu[i], sd[i])
      delta_inv[i] ~ dinterval(Ct[i], CT.max)
    }
    for (i in 1:num_replicates_star){
      mu_star[i] <- alpha1[P_standards[i]] + alpha2[P_standards[i]]*log(w_star[i])
      sd_star[i] <- sqrt(exp(a + b*log(w_star[i])))
      Ct_star[i] ~ dnorm(mu_star[i], sd_star[i])
      delta_star_inv[i] ~ dinterval(Ct_star[i], CT.max)
    }
  }
  
})
