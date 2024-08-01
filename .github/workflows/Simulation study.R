# File contains functions relating to running the simulation
# study in the paper

#################### Set model Parameters #####################
# Naming convention: set_data1, set_data2, set_data3 relate to
# Models 1, 2, and 3 respectively.

set_data1 = function(dataParams, trueParams){
  # Data - with contamination and variance
  
  delta_inv = 1*(trueParams$delta_true==0)
  delta_star_inv = 1*(trueParams$delta_star_true==0)
  
  DEdata = list(Ct = dataParams$data_real$Ct, 
                Ct_star = dataParams$data_standard$Ct,
                delta_inv = delta_inv, 
                delta_star_inv = delta_star_inv,
                w_star = dataParams$data_standard$Quantity, 
                Xb = dataParams$X_b, 
                Xw = dataParams$X_w,
                constraint_data = 1)
  
  return(DEdata)
}

set_data2 = function(dataParams, trueParams){
  # Data - with contamination and no variance
  
  delta_inv = 1*(trueParams$delta_true==0)
  delta_star_inv = 1*(trueParams$delta_star_true==0)
  
  DEdata = list(Ct = dataParams$data_real$Ct, 
                Ct_star = dataParams$data_standard$Ct,
                delta_inv = delta_inv, 
                delta_star_inv = delta_star_inv,
                w_star = dataParams$data_standard$Quantity,
                Xb = dataParams$X_b, 
                Xw = dataParams$X_w,
                constraint_data = 1)
  
  return(DEdata)
}

set_data3 = function(dataParams, trueParams){
  # Data - with no contamination and variance
  
  delta_inv = 1*(trueParams$delta_true==0)
  delta_star_inv = 1*(trueParams$delta_star_true==0)
  
  DEdata = list(Ct = dataParams$data_real$Ct, 
                Ct_star = dataParams$data_standard$Ct,
                delta_inv = delta_inv, 
                delta_star_inv = delta_star_inv,
                w_star = dataParams$data_standard$Quantity, 
                Xb = dataParams$X_b, 
                Xw = dataParams$X_w)
  
  return(DEdata)
}

set_inits1 = function(dataParams, trueParams){
  # Inits - with contamination and variance
  
  Ct_inits = 1*(is.na(dataParams$data_real$Ct)) # gives 0 for amplified values
  Ct_inits[Ct_inits == 0] = NA
  Ct_inits[Ct_inits == 1] = dataParams$CT.max + 2 # random initial value for non amplified values
  
  Ct_star_inits = 1*(is.na(dataParams$data_standard$Ct)) # gives 0 for amplified values
  Ct_star_inits[Ct_star_inits == 0] = NA
  Ct_star_inits[Ct_star_inits == 1] = dataParams$CT.max + 2 # random initial value for non amplified values
  
  l.init = t(matrix(trueParams$l_true, nrow = nT))
  
  DEinits = function(){list(betab = rep(0, ncovb),
                            betaw = rep(0, ncovw),
                            betab0 = trueParams$betab0_true,
                            tau2 = trueParams$tau_true**2, 
                            tau2.1 = trueParams$tau1_true**2,
                            sd_rho2 = .1,
                            rho = rep(1, dataParams$num_sites),
                            sigma2 = trueParams$sigma_true**2, 
                            alpha1 = trueParams$alpha1_true,
                            alpha2 = trueParams$alpha2_true,
                            a = trueParams$a_true, b = trueParams$b_true,
                            l = l.init, v = trueParams$v_true,
                            Ct = Ct_inits, Ct_star = Ct_star_inits, 
                            type = rep(2, dataParams$num_replicates), 
                            type_star = rep(2, dataParams$num_replicates_star),
                            pi.type = c(0.01, 0.98, 0.01),
                            rho0 = 1
  )}
  
  return(DEinits)
  
}

set_inits2 = function(dataParams, trueParams){
  # Inits - with contamination and no variance
  
  Ct_inits = 1*(is.na(dataParams$data_real$Ct)) # gives 0 for amplified values
  Ct_inits[Ct_inits == 0] = NA
  Ct_inits[Ct_inits == 1] = dataParams$CT.max + 2 # random initial value for non amplified values
  
  Ct_star_inits = 1*(is.na(dataParams$data_standard$Ct)) # gives 0 for amplified values
  Ct_star_inits[Ct_star_inits == 0] = NA
  Ct_star_inits[Ct_star_inits == 1] = dataParams$CT.max + 2 # random initial value for non amplified values
  
  l.init = t(matrix(trueParams$l_true, nrow = nT))
  
  DEinits = function(){list(betab = rep(0, ncovb),
                            betaw = rep(0, ncovw),
                            betab0 = trueParams$betab0_true,
                            tau2 = trueParams$tau_true**2, 
                            tau2.1 = trueParams$tau1_true**2,
                            sd_rho2 = .1,
                            rho = rep(1, dataParams$num_sites),
                            sigma2 = trueParams$sigma_true**2, 
                            alpha1 = trueParams$alpha1_true,
                            alpha2 = trueParams$alpha2_true,
                            sigma.y2 = rep(.5, dataParams$numP),
                            l = l.init, v = trueParams$v_true,
                            Ct = Ct_inits, Ct_star = Ct_star_inits, 
                            type = rep(2, dataParams$num_replicates), 
                            type_star = rep(2, dataParams$num_replicates_star),
                            pi.type = c(0.01, 0.98, 0.01),
                            rho0 = 1
  )}
  
  return(DEinits)
  
} 

set_inits3 = function(dataParams, trueParams){
  # Inits - with no contamination and variance
  
  Ct_inits = 1*(is.na(dataParams$data_real$Ct)) # gives 0 for amplified values
  Ct_inits[Ct_inits == 0] = NA
  Ct_inits[Ct_inits == 1] = dataParams$CT.max + 2 # random initial value for non amplified values
  
  Ct_star_inits = 1*(is.na(dataParams$data_standard$Ct)) # gives 0 for amplified values
  Ct_star_inits[Ct_star_inits == 0] = NA
  Ct_star_inits[Ct_star_inits == 1] = dataParams$CT.max + 2 # random initial value for non amplified values
  
  l.init = t(matrix(trueParams$l_true, nrow = nT))
  
  DEinits = function(){list(betab = rep(0, ncovb),
                            betaw = rep(0, ncovw),
                            betab0 = trueParams$betab0_true,
                            tau2 = trueParams$tau_true**2, 
                            tau2.1 = trueParams$tau1_true**2,
                            sd_rho2 = .1,
                            rho = rep(1, dataParams$num_sites),
                            sigma2 = trueParams$sigma_true**2, 
                            alpha1 = trueParams$alpha1_true,
                            alpha2 = trueParams$alpha2_true,
                            a = trueParams$a_true, b =trueParams$ b_true,
                            l = l.init, v = trueParams$v_true,
                            Ct = Ct_inits, Ct_star = Ct_star_inits, 
                            rho0 = 1
  )}
  
  return(DEinits)
  
}

set_constants = function(dataParams, trueParams){
  
  DEconstants = list(numP = dataParams$numP, nT = dataParams$nT, 
                     ncovb = ncovb, ncovw = ncovw,
                     num_sites = dataParams$num_sites, 
                     num_samples = dataParams$num_samples,
                     id_site = dataParams$id_site, 
                     id_sample = dataParams$id_sample,
                     id_site_time = dataParams$id_site_time, 
                     id_site_l = dataParams$id_site_l,
                     num_replicates = dataParams$num_replicates, 
                     num_replicates_star = dataParams$num_replicates_star, 
                     P = dataParams$P, P_standards = dataParams$P_star,
                     a_sigma.tau = 2, b_sigma.tau = 1,
                     a_sigma = 2, b_sigma = 1,
                     a_sigma.tau1 = 2, b_sigma.tau1 = 1,
                     a_sigma.rho = 2, b_sigma.rho = .1,
                     a_sigma_y = 2, b_sigma_y = .5,
                     alpha10 = mean(trueParams$alpha1_true), 
                     alpha20 = mean(trueParams$alpha2_true),
                     sigma_alpha = 1,
                     CT.max = dataParams$CT.max, sd_cont = 30,
                     mean_betab0 = trueParams$betab0_true, sd_betab0 = 1,
                     a0 = trueParams$a_true, b0 = trueParams$b_true, 
                     sd_a = 1, sd_b = 1,
                     alpha = c(0.01, 0.98, 0.01))
  
  return(DEconstants)
  
}

#################### Run Simulation Comparing models #####################
# Sim_name: character string for name of saved output
# nchains = chains to run each simulation
# niter = number of iterations to run MCMC
# nburn = burn-in
# Nsims = how many simulations to run
# save.every = how many simulations to run before saving output
# (saves current progress)

run_test = function(Sim_name = 'Simulation', 
                    nchains = 1, niter = 11000, 
                    nburn = 1000, Nsims = 100,
                    save.every = 20){
  # Sim_name (str): the file name the simulation results will be saved under
  # nchains (int): number of chains to run
  # niters (int): number of iterations per chain
  # nburn (int): number of burn in iterations
  # N (int): number of simulations to run
  # save.every (int): save all chains after save.every simulations
  
  data <- simulateData(setupParams, hyperParams)
  
  trueParams <- data$trueParams
  dataParams <- data$dataParams
  
  DEconstants <- set_constants(dataParams, trueParams)
  
  DEdata1 = set_data1(dataParams, trueParams)
  DEdata2 = set_data2(dataParams, trueParams)
  DEdata3 = set_data3(dataParams, trueParams)
  
  DEinits1 = set_inits1(dataParams, trueParams)
  DEinits2 = set_inits2(dataParams, trueParams)
  DEinits3 = set_inits3(dataParams, trueParams)
  
  DEmodel1 = nimbleModel(Full_Code, 
                         constants = DEconstants, 
                         data = DEdata1, 
                         inits = DEinits1())
  
  DEmodel2 = nimbleModel(Sigmay_Code, 
                         constants = DEconstants, 
                         data = DEdata2, 
                         inits = DEinits2())
  
  DEmodel3 = nimbleModel(Nocont_Code, 
                         constants = DEconstants, 
                         data = DEdata3, 
                         inits = DEinits3())
  
  data_list = list()
  data_list[[1]] = data
  
  for (n_iter in 1:Nsims){
    print(n_iter)
    
    if (n_iter != 1){ # re-set data and inits
      data <- simulateData(setupParams, hyperParams)
      
      trueParams <- data$trueParams
      dataParams <- data$dataParams
      
      DEdata1 = set_data1(dataParams, trueParams)
      DEmodel1$resetData()
      DEmodel1$setData(DEdata1)
      
      DEdata2 = set_data2(dataParams, trueParams)
      DEmodel2$resetData()
      DEmodel2$setData(DEdata2)
      
      DEdata3 = set_data3(dataParams, trueParams)
      DEmodel3$resetData()
      DEmodel3$setData(DEdata3)
      
      DEinits1 = set_inits1(dataParams, trueParams)
      DEinits2 = set_inits2(dataParams, trueParams)
      DEinits3 = set_inits3(dataParams, trueParams)
      
      DEmodel1$setInits(DEinits1())
      DEmodel2$setInits(DEinits2())
      DEmodel3$setInits(DEinits3())
      
      data_list[[n_iter]] = data
      
    }
    
    # Run the model
    mcmcConf1 = configureMCMC(DEmodel1, monitors = c('rho', 'a', 'b', 
                                                     'betaw', 'betab', 'betab0',
                                                     'tau', 'sigma', 'l', 'tau.1', 'pi.type'
    ), print = FALSE, enableWAIC = TRUE)
    
    mcmcConf2 = configureMCMC(DEmodel2, monitors = c('rho', 'sigma.y',
                                                     'betaw', 'betab', 'betab0',
                                                     'tau', 'sigma', 'l', 'tau.1', 'pi.type'
    ), print = FALSE, enableWAIC = TRUE)
    
    mcmcConf3 = configureMCMC(DEmodel3, monitors = c('rho', 'a', 'b', 
                                                     'betaw', 'betab', 'betab0',
                                                     'tau', 'sigma', 'l', 'tau.1'
    ), print = FALSE, enableWAIC = TRUE)
    
    
    
    DEmcmc1 <- buildMCMC(mcmcConf1, print = FALSE) 
    DEmcmc2 <- buildMCMC(mcmcConf2, print = FALSE) 
    DEmcmc3 <- buildMCMC(mcmcConf3, print = FALSE) 
    
    cDEmodel1 <- compileNimble(DEmodel1)
    cDEmodel2 <- compileNimble(DEmodel2)
    cDEmodel3 <- compileNimble(DEmodel3)
    
    cDEmcmc1 <- compileNimble(DEmcmc1, project = cDEmodel1)
    cDEmcmc2 <- compileNimble(DEmcmc2, project = cDEmodel2)
    cDEmcmc3 <- compileNimble(DEmcmc3, project = cDEmodel3)
    
    DEresults1 <- runMCMC(cDEmcmc1, niter = niter, nburnin = nburn, nchains = nchains, WAIC = TRUE)
    DEresults2 <- runMCMC(cDEmcmc2, niter = niter, nburnin = nburn, nchains = nchains, WAIC = TRUE)
    DEresults3 <- runMCMC(cDEmcmc3, niter = niter, nburnin = nburn, nchains = nchains, WAIC = TRUE)
    
    # store results
    if (n_iter == 1){
      
      par_chain1 = list()
      par_chain2 = list()
      par_chain3 = list()
      
      par_names1 = colnames(DEresults1$samples)
      par_names2 = colnames(DEresults2$samples)
      par_names3 = colnames(DEresults3$samples)
      
      for (par in par_names1){
        par_chain1[[par]] = matrix(numeric(0), nrow = niter - nburn, ncol = Nsims)
      }
      for (par in par_names2){
        par_chain2[[par]] = matrix(numeric(0), nrow = niter - nburn, ncol = Nsims)
      }
      for (par in par_names3){
        par_chain3[[par]] = matrix(numeric(0), nrow = niter - nburn, ncol = Nsims)
      }
    }
    
    for (par in par_names1){
      par_chain1[[par]][,n_iter] = DEresults1$samples[,par]
    }
    for (par in par_names2){
      par_chain2[[par]][,n_iter] = DEresults2$samples[,par]
    }
    for (par in par_names3){
      par_chain3[[par]][,n_iter] = DEresults3$samples[,par]
    }
    
    if (n_iter%%save.every == 0){
      filename <- paste0(Sim_name, '_', as.character(n_iter%/%save.every), '.RData')
      save(par_chain1, par_chain2, par_chain3, data_list, file = filename)
    } # iteratively save the data (just in case something goes wrong half way through)
    
  }
  
}

