################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################



# Parallel tempering. 
# Metropolis Coupled Markov Chain Monte Carlo

# Args:
#   init_theta: either a vector of intial paramter values
#               or a list of initial parameter values for each chaing
#   data_to_fit: the data!
#   pop_size: population sizes by age group (vector)
#   MaxTime: length of time to run the chain
#   indep: period of time for which chains should wander independently
#   covmat: covariance matrix for the parameter proposals
#   n_chains: number of chains to run. Must equal length of temperature vector
#   adapt_rate: at what rate to adapt the temperatures. Default 1.
#   adapt_start: at what point to start adapting. Default infinite
#   virus: which virus to fit, seasonal or covid. In current aws format only seasonal works
#   model: SEIR, # fututre options to include SEIRR or SEIRRR
#   bb: whether betabinomial or binomial. In current format only binomial


# Returns:
#   chains: list containing matrix for each chain, first col is loglik + log prior prob,
#           remaining columns are fn parameters in order given in the pars[[i]]
#   Temperatures: list of temperature values every indep step

mcmcmc_fn <- function(init_theta, MaxTime = 1e3,
                      indep = 100, covmat, n_chains,
                      adapt_rate=1, adapt_start = Inf,
                      prev_swap_proposed, prev_swap_accepted, pt, Temperatures, 
                      previous_trace = 0, virus , model_type,bb,
                      ...){
  
  theta_names <- names(proposal_sd)
  
  Temperatures_store <- matrix(nrow=MaxTime/indep, ncol=n_chains)
  adapt_switch  <- T
  #  covmat <- covmat + diag(ncol(covmat))*0.01
  pars <- init_theta
  n_pars <- length(theta_names)
  # list to store likelihood of each chain
  if(n_chains == 1){
    likelihood <-data.frame()
  } else{
    likelihood <- list(length = n_chains) }
  # create the empty traces for each chain
  chains <- lapply(1:n_chains, function(i)  matrix(NA, nrow  = MaxTime,
                                                   ncol = (4 + n_pars)) )
  # set the first row of each trace to 0
  for (i in n_chains){ chains[[i]][,n_pars + 3] <- 0 }
  
  # The independent intervals, run model in parallel for each chain
  Interval <- matrix(1:MaxTime, nrow = indep)
  
  print(Sys.time())

  # the outer time loop over intervals
  for(s in 1:(MaxTime / indep)){
    # print(paste0("s is ", s))
    # Evolve chains independently for "indep" time steps
    out <- lapply( 1:n_chains,
                   function(i){
                     # matrix to store output
                     out <- matrix(NA, ncol = length(theta_names) + 4,
                                   nrow = indep)
                     # For each time step in indep
                     #    print(paste0("i is ", i))
                     for(t in 1:indep){
                      #  print(paste0("t is ", t))
                       
                       #if first step, calculate the likelhiood of the parameters
                       if (t == 1){
                         # get the log posterior of the initial values
                         if(virus == "seasonal"){
                           liklihood_init <- get_log_posterior_no_covid(theta = pars[[i]], bb,
                                                                        model_type = model_type)
                           } 
                         #  if(!is.finite(liklihood_init)){browser()}
                        
                         Pi <- liklihood_init
                       }  else {
                         Pi <- out[t-1,1]
                       }
                       
                       #propose new parameters
                    
                       proposed <- step_fn(pars[[i]], covmat,
                                           lower_bounds, upper_bounds,
                                           theta_names)
                       # calculate the temperature weighting for the chain
                       beta_to_use <- temp_func(Temperatures, i)
                       
                       # calculate the temperature weighting for the chain
                       if(virus == "seasonal"){
                         likelihood_proposed <- get_log_posterior_no_covid(theta = proposed, bb,
                                                                           model_type = model_type)
                         } 
                     #  print(paste0("ll proposed", likelihood_proposed, " old ", Pi))
                       #  if(!is.finite(likelihood_proposed)){browser()}
                       # Calculate acceptance ratio
                       log_acceptance <- likelihood_proposed - Pi
                       #    if(!is.finite(log_acceptance)){browser()}
                       #adjst the acceptance for the unsymmetrical proposal distribution
                       log_acceptance <- adjust_unsymmetrical(log_acceptance = log_acceptance,
                                                              parameters_old = pars[[i]],
                                                              parameters_new = proposed,
                                                              covmat = covmat,
                                                              lower_bounds, upper_bounds,
                                                              theta_names)
                       
                       #alpha is the weighted liklihood of the proposal by the temperature
                       alpha <- exp(beta_to_use * log_acceptance )
                       #check alpha is a number, if not make -inf
                       if(is.nan(alpha)){alpha <- -Inf}
                       # if new parameters accepted
                       if ( alpha  > runif(1) ){
                         #print("step accepted")
                         #update current parameters and likelihood
                         pars[[i]] <- proposed
                         Liklihood_stored <- (likelihood_proposed)
                         #   if(!is.finite(season_out_proposed)){browser()}
                         #  if(!is.finite(Prior_proposed)){browser()}
                         Accepted <- 1
                         #else store the likelihoood of the previous parameters
                       }else {
                        Liklihood_stored <- Pi
                         Accepted <- 0} 
                       #update the trace
                       out[t,] <- c(Liklihood_stored, unlist(pars[[i]][theta_names]),
                                    Accepted, 0,0)
               
                     }
                     # print(paste0("t is ", t))
                     out
                   })
    #print(paste0("s is ", s))
    # save the output from each indep into the main chains
    for(i in 1:n_chains){
      
      #copy the trace for the interval into the overall trace (chain)
      chains[[i]][Interval[,s],] <- out[[i]]
      # copy parameters from last step of indep trace

      pars[[i]] <- as.list(out[[i]][indep,][-c(1,n_pars+2, n_pars+3, n_pars+4)])
      names(pars[[i]]) <- theta_names
      # store the likelihood value for each so don't have to recaluclate later
      
      if(n_chains ==1){
        likelihood <- out[[i]][indep,][1]
      } else{
        likelihood[[i]]<- out[[i]][indep,][1]
      }
    }
    # Propose a chain swaps every "indep" time steps. num swaps - num chains

      if(n_chains > 1){
      for (swap_test in 1:n_chains){
        # Choose two out of total number of chains and store their numbers in 'pick'
        pick <- sample(1:(n_chains-1), 1)
        i <- pick; j <- pick+1 # for convience
        # calculate the swap potential

        if(is.finite(likelihood[[i]]) & is.finite(likelihood[[j]])){
          R <- exp((likelihood[[i]] - likelihood[[j]]) /
                     ( temp_func(Temperatures , j) - temp_func(Temperatures , i)))
        } else { R = -Inf}
        
        #store swap attempted
        chains[[i]][Interval[dim(Interval)[1],s], n_pars + 4] <- 1
        #paste0("R is ", R, ". ll i is ", likelihood[[i]], ". ll j is ", likelihood[[j]]))
        # accept or reject swap
        if(R > runif(1)){

          pars_temp <- pars[[i]]
          pars[[i]] <- pars[[j]]
          pars[[j]] <- pars_temp
          
          #swap the likrliooh. stops a swap back to the worse prameters
          likelihood_temp <- likelihood[[i]]
          likelihood[[i]] <- likelihood[[j]]
          likelihood[[j]] <- likelihood_temp

#store the number it swapped to
          chains[[i]][Interval[dim(Interval)[1],s],n_pars+3]<-j
          
        } } }
    if(Interval[1,s] > adapt_start){
      if(adapt_switch == T){
        print("Starting to adapt at this point.......")
        adapt_switch <- F
      }
      # calculate swap rate for each temperatre
      # print(prev_swap_proposed)
      # print(prev_swap_accepted)
      swap_rates <- sapply(1:n_chains, function(i)
        swapping_rates(chains[[i]], n_pars = n_pars,
                       prev_swap_proposed = prev_swap_proposed[i],
                       prev_swap_accepted = prev_swap_accepted[i]))
      # swap_rates[1] is between 1 and 2, swap_rates[2] is between 2 and 3 etc.
      # calculate kappa
      if(previous_trace==0){
        time_adapting<- Interval[dim(Interval)[1],s] + (pt-1)*MaxTime - adapt_start
      } else{
        time_adapting <- (Interval[dim(Interval)[1],s] + 
                            (pt-1)*MaxTime - adapt_start) + 
          (previous_trace - adapt_start)
      }
      
      kappa <- 1/(1+time_adapting)^(adapt_rate)
      # store current temperatures
      temps_store <- Temperatures
      for(tem in 2:(n_chains-1)){
        
        # calculate s = log(Ti - Ti+1) for first chain
        S2 <- log(temps_store[tem]- temps_store[tem-1])
        # calculate change in Schange = k[Ai - Ai+1] + S for first chain
        s2new <- kappa*(swap_rates[tem] - swap_rates[tem+1])+S2
        #calculate new temperature, based on new temperature for previous one
        Tnew <- exp(s2new)+Temperatures[tem-1]
        Temperatures[tem] <- Tnew
      }}
    #print at right point in time
    percent_points <- seq(from = 0, to = MaxTime / indep, by = (MaxTime / indep) /10 )
    if(s %in% percent_points){
      print(paste0("pt complete: ", s/(MaxTime/indep)*100, "% - ",
                   " time is ", Sys.time()))}
    
    Temperatures_store[s,] <- Temperatures
  }

  # Returns the full history of all chains
  swaps_proposed <- sapply(1:n_chains, function(i) swapping_props(chains[[i]], n_pars=n_pars,
                                                                  prev_swap_proposed = prev_swap_proposed[i]))
  swaps_accepted <- sapply(1:n_chains, function(i) swapping_acceps(chains[[i]], n_pars=n_pars,
                                                                   prev_swap_accepted = prev_swap_accepted[i]))
  #TODO return the swap attempts and swaps accepted for each chain
  chains[["temperatures"]] <- Temperatures_store
  chains[["swaps_proposed"]] <- swaps_proposed
  chains[["swaps_accepted"]] <- swaps_accepted
  return(chains)
}

trace_wrapper <- function(length_run ,
                          each_run ,
                          n_chains,
                          init_theta,
                          indep,
                          covmat,
                          adapt_rate,
                          adaption_starter,
                          proposal_sd,
                          Temperatures,
                          name, 
                          lower_bounds = lower_bounds, 
                          upper_bounds = upper_bounds, 
                          virus, 
                          model_type, 
                          bb, 
                          cont = F
){
  
  
  #calculate in which pt to start the adaption
  # which one should it start in and the remainder
  pt_start <- ceiling(adaption_starter/each_run)
  pt_rem <- adaption_starter %% each_run
  
  # wrapper around traces
  for (pt in 1:(length_run/each_run)) {
    
    # specify the adapt start for each run.
    if(pt == pt_start){
      if( pt_rem == 0 ) { adapt_start <- each_run } else{
        adapt_start <- pt_rem
      }
    } else if (pt < pt_start) {adapt_start <- Inf
    } else if (pt > pt_start) {adapt_start <- 0}
    # specify the init_theta as the main one if the first one
    if( pt == 1 ){
      if(cont == F){init_theta <- lapply(1:n_chains, function(x) as.list(init_theta))}
      prev_swap_proposed <- rep(0, n_chains)
      prev_swap_accepted <- rep(0, n_chains)
    }
    # Run the mcmc for the pt section
    trace <- mcmcmc_fn(init_theta = init_theta,
                       MaxTime = each_run,
                       indep = indep,
                       covmat = covmat,
                       n_chains = n_chains,
                       adapt_rate = adapt_rate,
                       adapt_start = adapt_start,
                       prev_swap_proposed = prev_swap_proposed,
                       prev_swap_accepted = prev_swap_accepted,
                       pt = pt,
                       Temperatures = Temperatures,
                       lower_bounds = lower_bounds,
                       upper_bounds = upper_bounds, 
                       virus = virus, 
                       model_type = model_type, bb = bb)

    # save the trace into total_trace
    if (pt ==1){total_trace <- trace} else {
      total_trace <- lapply(1:length(trace), function(x)
        add_on(total_trace[[x]], trace[[x]]))
    }

    
    # update the init_theta to the new one

    
    end_trace <- lapply(1:n_chains, function (x) 
      as.list(tail(total_trace[[x]],1)[,2:(length(proposal_sd)+1)]))
    for ( section in 1:n_chains){
      names(end_trace[[section]]) <- names(proposal_sd)}
    init_theta <- end_trace
    # update the temperatures to the correct ones
    
    Temperatures <- c(tail(trace[["temperatures"]],1))
    prev_swap_proposed <- trace[["swaps_proposed"]]
    prev_swap_accepted <- trace[["swaps_accepted"]]

    
    print(paste0("Overall: ", pt/(length_run/each_run)*100, "%"))
    save(total_trace, file=paste0(name, "_", array_num, "_", Sys.Date(),
                                  ".Rdata"))

    #print(total_trace)
  }
  return(total_trace)
}