#################### Example set up for simulation ###################

{
  ncovb = 2 # num site level covariates
  ncovw = 2 # num sample level covariates
  n = 10 # num sites
  nT = 20 # num time-points
  M = rep(5, n*nT) # num samples per site
  K = rep(5, sum(M)) # num replicates per sample
}

dataParams <- list(n = n,
                   nT = nT,
                   M = M,
                   K = K,
                   w_standards = c(3e+07, 3e+06, 3e+05, 3e+04, 3e+03, 3e+02, 3e+01), # standard concentrations
                   K_standards = 3, # num replicates per standard concentration
                   CT.max = 40 # Censoring limit
)

hyperParams <- list(tau = 1, # sd time series
                    sigma = 1, # sd samples
                    tau2.1 = 1, # variance across sites when t=1
                    betab = c(1, -1), # site covariate coeffs
                    betaw = c(1, -1), # sample covariate coeffs
                    betab0 = 6, # intercept across sites when t=1
                    alpha1.0 = 44, # mean plate intercept
                    alpha2.0 = -1.7, # mean plate slope
                    sigma_alpha1 = .1, # sd plate intercept
                    sigma_alpha2 = .01, # sd plate slope
                    rho = rep(1, n), # time series coeff
                    a = 0.2, # plate variance intercept
                    b = -0.25, # plate variance slope
                    lambda0 = 3e+3, # mean contamination concentration
                    sd_lambda = 100, # sd contamination concentration
                    p0 = c(0.05, 0.1), # prob contamination and inhibition
                    multiplier = -9/10 # inhibition effect
)


#################### simulate data ###########################
# Function that returns an example data set 

simulateData <- function(dataParams,
                         hyperParams){
  
  list2env(dataParams,environment()) 
  list2env(hyperParams,environment()) 
  
  Xb = array(NA, dim = c(n, nT, ncovb))
  Xw = array(NA, dim = c(sum(M), ncovw))
  
  # half continuous, half categorical covariates
  Xb[,,1:floor(ncovb/2)] = array(rnorm(n * nT * floor(ncovb/2)), dim = c(n, nT, floor(ncovb/2)))
  Xb[,,(floor(ncovb/2)+1):ncovb] = array(sample(0:1, size = n*nT*(ncovb - floor(ncovb/2)), replace = T), 
                                         dim = c(n, nT, (ncovb - floor(ncovb/2))))
  
  # half continuous, half categorical covariates
  Xw[,1:floor(ncovw/2)] = rnorm(sum(M) * floor(ncovw/2))
  Xw[,(floor(ncovw/2)+1):ncovw] = sample(0:1, size = sum(M)*(ncovw - floor(ncovw/2)),
                                         replace = T)
  
  # Create the L's:
  l = matrix(NA, ncol = nT, nrow = n)
  l[,1] = rnorm(n, betab0 + Xb[,1,]%*%betab, sd = sqrt(tau2.1))
  for (i in 1:n){
    for (j in 2:nT){
      l[i, j] = rnorm(1, rho[i]*l[i,j-1] + (Xb[i, j, ] - rho[i]*Xb[i, j-1, ])%*%betab, sd = tau)
    }
  }
  
  l = as.vector(t(l)) # collapse into vector
  id_site_l = rep(1:n, each = nT) # link l_it to site
  id_site_time = rep(1:nT, n) # link l_it to time
  
  id_site = rep(1:length(l), M) # link sample to l_it
  Xbetaw = Xw %*% betaw
  v = rnorm(sum(M), rep(l, M) + Xbetaw, sd = sigma) # v per sample
  
  w = exp(v) 
  
  # observations
  # assume that each site and time uses its own plate
  alpha1 = rnorm(n*nT, alpha1.0, sd = sigma_alpha1)
  alpha2 = rnorm(n*nT, alpha2.0, sd = sigma_alpha2)
  
  id_sample = rep(1:sum(M), K) # links replicate to sample
  P = rep(id_site, K) # assume plate for each site and time
  
  Ct = rep(NA, sum(K))
  lambda = rep(0, sum(K)) # indicator for contamination
  for (i in 1:sum(K)){ # loop through replicates
    u = runif(1)
    if (u <= p0[1]){# if outlier is contamination
      lambda[i] = rnorm(1, lambda0, sd_lambda)
    } else if (u >= (1 - p0[2])) {# if outlier is inhibition
      lambda[i] <- multiplier * w[id_sample[i]]
    }
    sd = sqrt(exp(a + b*log(w[id_sample[i]] + lambda[i]))) # heterogeneous variance
    mu = alpha1[P[i]] + alpha2[P[i]]*log(w[id_sample[i]] + lambda[i]) # add contamination to replicate
    Ct[i] = rnorm(1, mu, sd = sd)
    if (Ct[i] > CT.max){
      Ct[i] = NA
    }
  }
  
  delta = 1*!is.na(Ct) # returns 1 if amplified, 0 if not
  
  data_real <- data.frame(Site = rep(rep(id_site_l, M), K),
                          Time = rep(rep(id_site_time, M), K),
                          Sample = id_sample,
                          Plate = P,
                          Ct = Ct)
  
  # standards:
  # there will be a set of standards for every plate that we are using
  w_star = rep(w_standards, each = K_standards)
  w_star = rep(w_star, n*nT) # make replicates per plate
  P_standards = rep(1:(n*nT), each = length(w_standards)*K_standards)
  
  Ct_star = rep(NA, length(w_star))
  lambda_star = rep(0, length(w_star))
  for (i in 1:length(w_star)){ # loop through replicates
    u = runif(1)
    if (u <= p0[1]){# if outlier is contamination
      lambda_star[i] = rnorm(1, lambda0, sd_lambda)
    } else if (u >= (1 - p0[2])) {# if outlier is inhibition
      lambda_star[i] <- multiplier * w_star[i]
    }
    mu = alpha1[P_standards[i]] + alpha2[P_standards[i]]*log(w_star[i] + lambda_star[i])
    sd = sqrt(exp(a + b*log(w_star[i] + lambda_star[i])))
    Ct_star[i] = rnorm(1, mu, sd = sd)
    if (Ct_star[i] > CT.max){
      Ct_star[i] = NA
    }
  }
  
  delta_star = 1*!is.na(Ct_star)
  
  data_standard <- data.frame(Plate = P_standards,
                              Quantity = w_star,
                              Ct = Ct_star)
  
  # save true values -----------
  
  trueParams <- list(
    l_true = l,
    v_true = v,
    rho_true = rho,
    tau_true = tau,
    tau1_true = sqrt(tau2.1),
    alpha1_true = alpha1,
    alpha2_true = alpha2,
    a_true = a,
    b_true = b,
    beta_b_true = betab,
    beta_w_true = betaw,
    betab0_true = betab0,
    delta_true = delta,
    delta_star_true = delta_star,
    sigma_true = sigma,
    lambda_true = lambda,
    lambda_star_true = lambda_star,
    gamma_true = 1*(lambda != 0),
    gamma_star_true = 1*(lambda_star != 0),
    p0_true = p0,
    lambda0_true = lambda0,
    sd_lambda_true = sd_lambda
  )
  
  dataParams <- list("data_real" = data_real,
                     "data_standard" = data_standard,
                     "X_b" = Xb,
                     "X_w" = Xw,
                     "ncovb" = dim(Xb)[3],
                     "ncovw" = dim(Xw)[2],
                     "id_site" = id_site,
                     "id_sample" = id_sample,
                     "id_site_l" = id_site_l,
                     "id_site_time" = id_site_time,
                     "P" = data_real$Plate,
                     "P_star" = data_standard$Plate,
                     "numP" = max(P),
                     "num_sites" = n,
                     "nT" = nT,
                     "num_samples" = sum(M),
                     "num_replicates" = sum(K),
                     "num_replicates_star" = length(w_star),
                     "CT.max" = CT.max
                     
  )
  
  list("trueParams" = trueParams, 
       "dataParams" = dataParams)
}
