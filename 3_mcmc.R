{###############################################################################
  #loading data and function files
  ###############################################################################
  source("1_dataloading_eda.R")
  source("2_functions.R")
  load("parms_ga_sigma.RData")
  ###############################################################################
}


webshot::install_phantomjs()
#font_import()
pacman::p_load(webshot,orca)
# install.packages("plotly")
library(plotly)
library(gridExtra)
library(ggplot2)
library(deSolve)
{
  # Global variable for text size
  TITLE_TEXT_SIZE <- 9
  X_TEXT_SIZE  <- 9
  x_tf = "Lastik Test SemiBold" #xaxis_title_font
  x_ts = 50 #xaxis_title_size
  x_vf = "Lastik Test Regular" #xaxis_value_font\
  x_vs = 50 #xaxis_value_size
  y_tf = "Lastik Test SemiBold" #yaxis_title_font
  y_ts = 50 #yaxis_title_size
  y_vf = "Lastik Test Regular" #yaxis_value_font
  y_vs = 50 #yaxis_value_size 
  lg_f = "Lastik Test Bold" #legend_font
  lg_s = 40 #legend_size
  
}

# standard deviations
sd_vec = c(
  betal_sd = 0.5,
  deltal_sd = 0.5,
  gammal_sd = 0.5,
  S_sd = 1000,
  El_sd = 1000,
  Il_sd = 1000,
  R_sd = 3000,
  deltal_mean_sd = 0.05,
  gammal_mean_sd = 0.05,
  betal_sd_prior = 0.0015,
  S_sd = 1500, #11
  El_sd = 1500, #12
  Il_sd = 1500, #13
  R_sd = 1500, #14
  tau_sd1 = 0.01,
  tau_sd2 = 0.01
)
{
tau1_init = log(1/10000)
tau2_init = log(1/10000)
# Setting up the inverse gamma prior parameters for sigma_rw
invg_shape <- 1
invg_rate <- 6.5 # rate parameter
hist(1/rgamma(100000, shape = invg_shape, rate = invg_rate), breaks = 1000, xlim = c(0,0.1))
tau1_init = parms$tau1
tau2_init = parms$tau2
sd_vec[10] = parms$sigma_rw
betal_eda = parms$beta_x
deltal_eda = parms$delta_x
gammal_eda = parms$gamma_x
El_eda[1] = parms$E0
Il_eda[1] = parms$I0

###############################################################################
###############################################################################
###############################################################################

# Setting up the gamma prior hyper parameters for tau parameters
gamma_tau1_shape <- 50  # Shape parameter close to zero
gamma_tau1_rate <- 1000    # rate parameter
gamma_tau2_shape <- 20  # Shape parameter close to zero
gamma_tau2_rate <- 1000    # rate parameter
hist(rgamma(100000, shape = gamma_tau1_shape, rate = gamma_tau1_rate), breaks = 1000)
hist(rgamma(100000, shape = gamma_tau2_shape, rate = gamma_tau2_rate), breaks = 1000)
}
###############################################################################
###############################################################################
# sub population size for computational efficiency
pop0 = 50000
###############################################################################

# MCMC parameters
num_iter_end <- 300000  # Number of MCMC iteration
T = as.numeric(length(betal_eda))
#timepoint <- 2:T  
# Start from the second timepoint because of the random walk structure
ve = 0.75
{  
  ###############################################################################
  # Initialize arrays to store MCMC samples
  betal_samples <- vector("list", num_iter_end)
  deltal_samples <- vector("list", num_iter_end)
  gammal_samples <- vector("list", num_iter_end)
  deltal_mean_samples <- vector("numeric", num_iter_end)
  gammal_mean_samples <- vector("numeric", num_iter_end)
  El_samples <- vector("list", num_iter_end)
  Il_samples <- vector("list", num_iter_end)
  S_samples = vector("list", num_iter_end);
  R_samples = vector("list", num_iter_end);
  El_sol_samples = vector("list", num_iter_end);
  Il_sol_samples = vector("list", num_iter_end);
  S_sol_samples = vector("list", num_iter_end);
  R_sol_samples = vector("list", num_iter_end);
  tau1_samples = vector("numeric", num_iter_end);
  tau2_samples = vector("numeric", num_iter_end);
  betal_rw_var_samples = vector("numeric", num_iter_end);
  
  
  # Initialize the parameters with zero index values 
  betal_samples[[1]] = betal_eda[1:(T-1)];
  deltal_samples[[1]] = deltal_eda[1:(T-1)];
  deltal_mean_samples[1] = mean(deltal_eda);
  gammal_samples[[1]] = gammal_eda[1:(T-1)];
  gammal_mean_samples[1] = mean(gammal_eda);
  El_samples[[1]] = El_eda[1:(T)];
  Il_samples[[1]] = Il_eda[1:(T)];
  S_samples[[1]] = S_eda[1:(T)];
  R_samples[[1]] = R_eda[1:(T)];
  El_sol_samples[[1]] = El_eda[1:(T)];
  Il_sol_samples[[1]] = Il_eda[1:(T)];
  R_sol_samples[[1]] = R_eda[1:(T)];
  tau1_samples[1] = tau1_init;
  tau2_samples[1] = tau2_init;
  betal_rw_var_samples[1] = parms$sigma_rw;
  
  
  v_x = c(0, (cases_vaccinated))[1:(T)] # vaccination out of pop0
}
prop_sd1 = 0.3
prop_sd2 = 0.3
prop_sd3 = 1

#######################################################################################################################
# Metropolis-Hastings within Gibbs
#######################################################################################################################
{
  save_interval <- 50000
  num_iter_start = 2
  pb <- txtProgressBar(min = 2, max = num_iter_end, style = 3)
  set.seed(12345)
  for (i in num_iter_start:num_iter_end) {
    ###################################################################################################################
    ###################################################################################################################
    ######################################### 1. update betal_t #######################################################
    ###################################################################################################################
    ###################################################################################################################
    #Updating betal[0] with time index [1] 
    #Propose new value from normal distribution with mean at previous iteration
    t = 1;
    betal_t_prop <- rnorm(1, mean = betal_eda[t], sd = prop_sd1)
    
    E0 = El_samples[[i-1]][t] * pop0/popsize # Initial number of exposed out of pop0
    I0 = Il_samples[[i-1]][t] * pop0/popsize # Initial number of infected out of pop0
    R0 = 0
    V0 = sum(v_x[1:t]) * vax_eff * pop0/popsize # Initial number of vaccinated out of pop0
    W =  pop0 - (E0 + I0 + R0 + V0)
    N = W + E0 + I0 + R0 + V0
    state_values = c(S = W, E = E0, I = I0, R = R0, V = V0)
    
    ode_state_parms_prop <- c(beta = exp(betal_t_prop), 
                              delta = exp(deltal_samples[[i-1]][t]),
                              gamma = exp(gammal_samples[[i-1]][t]),
                              v_t = v_x[t] * pop0/popsize
    )  
    
    sol_prop = ode_solution(state_values, parms = ode_state_parms_prop)
    #S_sol is real solution of ODE at (t+1).
    S_sol_prop = (sol_prop[2,2]) * popsize/pop0 ; 
    El_sol_prop = (sol_prop[2,3]) * popsize/pop0 ; 
    Il_sol_prop = (sol_prop[2,4]) * popsize/pop0 ; 
    R_sol_prop = (sol_prop[2,5]) * popsize/pop0;
    
    ode_state_parms_curr <- c(beta = exp(betal_samples[[i-1]][t]), 
                              delta = exp(deltal_samples[[i-1]][t]),
                              gamma = exp(gammal_samples[[i-1]][t]),
                              v_t = v_x[t] * pop0/popsize
    )  
    
    sol_curr = ode_solution(state_values, parms = ode_state_parms_curr)
    #S_sol is real solution of ODE at (t+1).
    S_sol_curr = (sol_curr[2,2]) * popsize/pop0 ; 
    El_sol_curr = (sol_curr[2,3]) * popsize/pop0 ; 
    Il_sol_curr = (sol_curr[2,4]) * popsize/pop0 ; 
    R_sol_curr = (sol_curr[2,5]) * popsize/pop0;
    
    
    
    # Calculate the log priors and log likelihoods
    # For current value
    betal_t_curr_log_postd <- betal_log_postd_0(betal = betal_samples[[i-1]][t], 
                                                betal_next = betal_samples[[i-1]][t+1],
                                                El_t1 = El_samples[[i-1]][t+1],
                                                El_t1sol = El_sol_curr,
                                                Il_t1 = Il_samples[[i-1]][t+1],
                                                Il_t1sol = Il_sol_curr,
                                                R_t1 = R_samples[[i-1]][t+1],
                                                R_t1sol = R_sol_curr,
                                                deltal = deltal_samples[[i-1]][t],
                                                gammal = gammal_samples[[i-1]][t],
                                                t,
                                                obs1 = cases_new[t+1],
                                                obs2 = cases_recovered[t+1],
                                                sd_vec = sd_vec, 
                                                tau1 = tau1_samples[i-1],
                                                tau2 = tau2_samples[i-1]
    )
    
    # Calculate the log priors and log likelihoods
    # For proposed value
    betal_t_prop_log_postd <- betal_log_postd_0(betal = betal_t_prop, 
                                                betal_next = betal_samples[[i-1]][t+1],
                                                El_t1 = El_samples[[i-1]][t+1],
                                                El_t1sol = El_sol_prop,
                                                Il_t1 = Il_samples[[i-1]][t+1],
                                                Il_t1sol = Il_sol_prop,
                                                R_t1 = R_samples[[i-1]][t+1],
                                                R_t1sol = R_sol_prop,
                                                deltal = deltal_samples[[i-1]][t],
                                                gammal = gammal_samples[[i-1]][t],
                                                t,
                                                obs1 = cases_new[t+1],
                                                obs2 = cases_recovered[t+1],
                                                sd_vec = sd_vec, 
                                                tau1 = tau1_samples[i-1],
                                                tau2 = tau2_samples[i-1]
    )
    
    log_q_beta_t_to_prop <- dnorm(betal_t_prop, 
                                  mean = betal_eda[t], 
                                  sd = sd_vec[1], log = TRUE)
    log_q_beta_prop_to_t <- dnorm(betal_samples[[i-1]][t],
                                  mean = betal_eda[t], 
                                  sd = sd_vec[1], log = TRUE)
    
    betal_t_prop_log_postd2 = betal_t_prop_log_postd + log_q_beta_prop_to_t
    betal_t_curr_log_postd2 = betal_t_curr_log_postd + log_q_beta_t_to_prop
    
    # Calculate the log acceptance probabilities
    accep_prob <- min(1, exp(betal_t_prop_log_postd2 - betal_t_curr_log_postd2))
    
    # Acceptance/rejection step for log_beta
    rand_u = runif(1) # Generate a random number from a uniform distribution
    
    if (!is.na(accep_prob) && rand_u < accep_prob) {
      betal_samples[[i]][t] <- betal_t_prop
      El_sol_samples[[i]][t+1] = El_sol_prop
      Il_sol_samples[[i]][t+1] = Il_sol_prop
      R_sol_samples[[i]][t+1] = R_sol_prop
    } else {
      betal_samples[[i]][t] <- betal_samples[[i-1]][t]
      El_sol_samples[[i]][t+1] = El_sol_curr
      Il_sol_samples[[i]][t+1] = Il_sol_curr
      R_sol_samples[[i]][t+1] = R_sol_curr
    }
    #############################################################################################
    #Updating betal[t] with time index [2:365] T=366
    
    for (t in 2:(T-2)) {
      
      #Propose new value from normal distribution with mean at previous iteration
      betal_t_prop <- rnorm(1, mean = betal_eda[t], sd = prop_sd1)
      
      
      E0 = El_samples[[i-1]][t] * pop0/popsize # Initial number of exposed out of pop0
      I0 = Il_samples[[i-1]][t] * pop0/popsize # Initial number of infected out of pop0
      R0 = R_samples[[i-1]][t] * pop0/popsize # Initial number of recovered out of pop0
      V0 = sum(v_x[1:t]) * vax_eff * pop0/popsize # Initial number of vaccinated out of pop0
      W =  pop0 - (E0 + I0 + R0 + V0)
      N = W + E0 + I0 + R0 + V0
      state_values = c(S = W, E = E0, I = I0, R = R0, V = V0)
      
      #parameters of ODE at current betal_t at i-1
      ode_state_parms_curr <- c(beta = exp(betal_samples[[i-1]][t]), 
                                delta = exp(deltal_samples[[i-1]][t]),
                                gamma = exp(gammal_samples[[i-1]][t]),
                                v_t = v_x[t] * pop0/popsize
      )
      
      #parameters of ODE at proposed value betal_t_prop 
      ode_state_parms_prop <- c(beta = exp(betal_t_prop), 
                                delta = exp(deltal_samples[[i-1]][t]),
                                gamma = exp(gammal_samples[[i-1]][t]),
                                v_t = v_x[t] * pop0/popsize
      ) 
      
      
      #solution of ODE at current value betal_t at i-1.
      sol_curr = ode_solution(state_values, parms = ode_state_parms_curr)
      S_sol_curr = (sol_curr[2,2]) * popsize/pop0 ; 
      El_sol_curr = (sol_curr[2,3]) * popsize/pop0 ; 
      Il_sol_curr = (sol_curr[2,4]) * popsize/pop0 ; 
      R_sol_curr = (sol_curr[2,5]) * popsize/pop0;
      
      #solution of ODE at proposed value betal_t at i-1.
      sol_prop = ode_solution(state_values, parms = ode_state_parms_prop)
      #S_sol is real solution of ODE at (t+1).
      S_sol_prop = (sol_prop[2,2]) * popsize/pop0 ; 
      El_sol_prop = (sol_prop[2,3]) * popsize/pop0 ; 
      Il_sol_prop = (sol_prop[2,4]) * popsize/pop0 ; 
      R_sol_prop = (sol_prop[2,5]) * popsize/pop0;
      
      
      # Calculate the log priors and log likelihoods
      # For current value
      betal_t_curr_log_postd <- betal_log_postd(betal = betal_samples[[i-1]][t], 
                                                betal_prev =  betal_samples[[i]][t-1], #NOTE: $$$$$$$$$$
                                                betal_next = betal_samples[[i-1]][t+1],
                                                El_t1 = El_samples[[i-1]][t+1],
                                                El_t1sol = El_sol_curr,
                                                Il_t1 = Il_samples[[i-1]][t+1],
                                                Il_t1sol = Il_sol_curr,
                                                R_t1 = R_samples[[i-1]][t+1],
                                                R_t1sol = R_sol_curr,
                                                deltal = deltal_samples[[i-1]][t],
                                                gammal = gammal_samples[[i-1]][t],
                                                t,
                                                obs1 = cases_new[t+1],
                                                obs2 = cases_recovered[t+1],
                                                sd_vec = sd_vec, 
                                                tau1 = tau1_samples[i-1],
                                                tau2 = tau2_samples[i-1]
      )
      
      # Calculate the log priors and log likelihoods
      # For proposed value
      betal_t_prop_log_postd <- betal_log_postd(betal = betal_t_prop, 
                                                betal_prev =  betal_samples[[i]][t-1],
                                                betal_next = betal_samples[[i-1]][t+1],
                                                El_t1 = El_samples[[i-1]][t+1],
                                                El_t1sol =  El_sol_prop,
                                                Il_t1 = Il_samples[[i-1]][t+1],
                                                Il_t1sol = Il_sol_prop,
                                                R_t1 = R_samples[[i-1]][t+1],
                                                R_t1sol = R_sol_prop,
                                                deltal = deltal_samples[[i-1]][t],
                                                gammal = gammal_samples[[i-1]][t],
                                                t,
                                                obs1 = cases_new[t+1],
                                                obs2 = cases_recovered[t+1],
                                                sd_vec = sd_vec, 
                                                tau1 = tau1_samples[i-1],
                                                tau2 = tau2_samples[i-1]
      )
      
      log_q_beta_t_to_prop <- dnorm(betal_t_prop, 
                                    mean = betal_eda[t], 
                                    sd = sd_vec[1], log = TRUE)
      log_q_beta_prop_to_t <- dnorm(betal_samples[[i-1]][t],
                                    mean = betal_eda[t], 
                                    sd = sd_vec[1], log = TRUE)
      
      betal_t_prop_log_postd2 = betal_t_prop_log_postd + log_q_beta_prop_to_t
      betal_t_curr_log_postd2 = betal_t_curr_log_postd + log_q_beta_t_to_prop
      
      # Calculate the log acceptance probabilities
      accep_prob <- min(1, exp(betal_t_prop_log_postd2 - betal_t_curr_log_postd2))
      
      # Acceptance/rejection step for log_beta
      rand_u = runif(1) # Generate a random number from a uniform distribution
      
      
      
      if (!is.na(accep_prob) && rand_u < accep_prob) {
        betal_samples[[i]][t] <- betal_t_prop
        El_sol_samples[[i]][t+1] = El_sol_prop
        Il_sol_samples[[i]][t+1] = Il_sol_prop
        R_sol_samples[[i]][t+1] = R_sol_prop
      } else {
        betal_samples[[i]][t] <- betal_samples[[i-1]][t]
        El_sol_samples[[i]][t+1] = El_sol_curr
        Il_sol_samples[[i]][t+1] = Il_sol_curr
        R_sol_samples[[i]][t+1] = R_sol_curr
      }
    }
    #############################################################################################
    #Updating betal[t] with time index [2:365] T=366
    
    t=T-1
    
    #Propose new value from normal distribution with mean at previous iteration
    betal_t_prop <- rnorm(1, mean = betal_eda[t], sd = prop_sd1)
    
    
    E0 = El_samples[[i-1]][t] * pop0/popsize # Initial number of exposed out of pop0
    I0 = Il_samples[[i-1]][t] * pop0/popsize # Initial number of infected out of pop0
    R0 = R_samples[[i-1]][t] * pop0/popsize # Initial number of recovered out of pop0
    V0 = sum(v_x[1:t]) * vax_eff * pop0/popsize # Initial number of vaccinated out of pop0
    W =  pop0 - (E0 + I0 + R0 + V0)
    N = W + E0 + I0 + R0 + V0
    state_values = c(S = W, E = E0, I = I0, R = R0, V = V0)
    
    #parameters of ODE at current betal_t at i-1
    ode_state_parms_curr <- c(beta = exp(betal_samples[[i-1]][t]), 
                              delta = exp(deltal_samples[[i-1]][t]),
                              gamma = exp(gammal_samples[[i-1]][t]),
                              v_t = v_x[t] * pop0/popsize
    )
    
    #parameters of ODE at proposed value betal_t_prop 
    ode_state_parms_prop <- c(beta = exp(betal_t_prop), 
                              delta = exp(deltal_samples[[i-1]][t]),
                              gamma = exp(gammal_samples[[i-1]][t]),
                              v_t = v_x[t] * pop0/popsize
    ) 
    
    
    #solution of ODE at current value betal_t at i-1.
    sol_curr = ode_solution(state_values, parms = ode_state_parms_curr)
    S_sol_curr = (sol_curr[2,2]) * popsize/pop0 ; 
    El_sol_curr = (sol_curr[2,3]) * popsize/pop0 ; 
    Il_sol_curr = (sol_curr[2,4]) * popsize/pop0 ; 
    R_sol_curr = (sol_curr[2,5]) * popsize/pop0;
    
    #solution of ODE at proposed value betal_t at i-1.
    sol_prop = ode_solution(state_values, parms = ode_state_parms_prop)
    #S_sol is real solution of ODE at (t+1).
    S_sol_prop = (sol_prop[2,2]) * popsize/pop0 ; 
    El_sol_prop = (sol_prop[2,3]) * popsize/pop0 ; 
    Il_sol_prop = (sol_prop[2,4]) * popsize/pop0 ; 
    R_sol_prop = (sol_prop[2,5]) * popsize/pop0;
    
    
    # Calculate the log priors and log likelihoods
    # For current value
    betal_t_curr_log_postd <- betal_log_postd_t(betal = betal_samples[[i-1]][t], 
                                                betal_prev =  betal_samples[[i]][t-1], #NOTE: $$$$$$$$$$
                                                El_t1 = El_samples[[i-1]][t+1],
                                                El_t1sol = El_sol_curr,
                                                Il_t1 = Il_samples[[i-1]][t+1],
                                                Il_t1sol = Il_sol_curr,
                                                R_t1 = R_samples[[i-1]][t+1],
                                                R_t1sol = R_sol_curr,
                                                deltal = deltal_samples[[i-1]][t],
                                                gammal = gammal_samples[[i-1]][t],
                                                t,
                                                obs1 = cases_new[t+1],
                                                obs2 = cases_recovered[t+1],
                                                sd_vec = sd_vec, 
                                                tau1 = tau1_samples[i-1],
                                                tau2 = tau2_samples[i-1]
    )
    
    # Calculate the log priors and log likelihoods
    # For proposed value
    betal_t_prop_log_postd <- betal_log_postd_t(betal = betal_t_prop, 
                                                betal_prev =  betal_samples[[i]][t-1],
                                                El_t1 = El_samples[[i-1]][t+1],
                                                El_t1sol =  El_sol_prop,
                                                Il_t1 = Il_samples[[i-1]][t+1],
                                                Il_t1sol = Il_sol_prop,
                                                R_t1 = R_samples[[i-1]][t+1],
                                                R_t1sol = R_sol_prop,
                                                deltal = deltal_samples[[i-1]][t],
                                                gammal = gammal_samples[[i-1]][t],
                                                t,
                                                obs1 = cases_new[t+1],
                                                obs2 = cases_recovered[t+1],
                                                sd_vec = sd_vec, 
                                                tau1 = tau1_samples[i-1],
                                                tau2 = tau2_samples[i-1]
    )
    
    log_q_beta_t_to_prop <- dnorm(betal_t_prop, 
                                  mean = betal_eda[t], 
                                  sd = sd_vec[1], log = TRUE)
    log_q_beta_prop_to_t <- dnorm(betal_samples[[i-1]][t],
                                  mean = betal_eda[t], 
                                  sd = sd_vec[1], log = TRUE)
    
    betal_t_prop_log_postd2 = betal_t_prop_log_postd + log_q_beta_prop_to_t
    betal_t_curr_log_postd2 = betal_t_curr_log_postd + log_q_beta_t_to_prop
    
    # Calculate the log acceptance probabilities
    accep_prob <- min(1, exp(betal_t_prop_log_postd2 - betal_t_curr_log_postd2))
    
    # Acceptance/rejection step for log_beta
    rand_u = runif(1) # Generate a random number from a uniform distribution
    
    
    
    if (!is.na(accep_prob) && rand_u < accep_prob) {
      betal_samples[[i]][t] <- betal_t_prop
      El_sol_samples[[i]][t+1] = El_sol_prop
      Il_sol_samples[[i]][t+1] = Il_sol_prop
      R_sol_samples[[i]][t+1] = R_sol_prop
    } else {
      betal_samples[[i]][t] <- betal_samples[[i-1]][t]
      El_sol_samples[[i]][t+1] = El_sol_curr
      Il_sol_samples[[i]][t+1] = Il_sol_curr
      R_sol_samples[[i]][t+1] = R_sol_curr
    }
    
    ###################################################################################################################
    ###################################################################################################################
    ######################################### 2. update deltal_t #######################################################
    ###################################################################################################################
    ###################################################################################################################
    #Updating deltal[0] with time index [1] 
    #Propose new value from normal distribution with mean at previous iteration
    
    for (t in 1:(T-1)) {
      
      #Propose new value from normal distribution with mean at previous iteration
      deltal_t_prop <- rnorm(1, mean = deltal_eda[t], sd = prop_sd2)
      
      E0 = El_samples[[i-1]][t] * pop0/popsize # Initial number of exposed out of pop0
      I0 = Il_samples[[i-1]][t] * pop0/popsize # Initial number of infected out of pop0
      R0 = R_samples[[i-1]][t] * pop0/popsize # Initial number of recovered out of pop0
      V0 = sum(v_x[1:t]) * vax_eff * pop0/popsize # Initial number of vaccinated out of pop0
      W =  pop0 - (E0 + I0 + R0 + V0)
      N = W + E0 + I0 + R0 + V0
      state_values = c(S = W, E = E0, I = I0, R = R0, V = V0)
      
      #parameters of ODE at current deltal_t at i-1
      ode_state_parms_curr <- c(beta = exp(betal_samples[[i]][t]), 
                                delta = exp(deltal_samples[[i-1]][t]),
                                gamma = exp(gammal_samples[[i-1]][t]),
                                v_t = v_x[t] * pop0/popsize
      )
      
      #parameters of ODE at proposed value deltal_t_prop 
      ode_state_parms_prop <- c(beta = exp(betal_samples[[i]][t]), 
                                delta = exp(deltal_t_prop),
                                gamma = exp(gammal_samples[[i-1]][t]),
                                v_t = v_x[t] * pop0/popsize
      ) 
      
      
      #solution of ODE at current value deltal_t at i-1.
      sol_curr = ode_solution(state_values, parms = ode_state_parms_curr)
      S_sol_curr = (sol_curr[2,2]) * popsize/pop0 ; 
      El_sol_curr = (sol_curr[2,3]) * popsize/pop0 ; 
      Il_sol_curr = (sol_curr[2,4]) * popsize/pop0 ; 
      R_sol_curr = (sol_curr[2,5]) * popsize/pop0;
      
      #solution of ODE at proposed value betal_t at i-1.
      sol_prop = ode_solution(state_values, parms = ode_state_parms_prop)
      #S_sol is real solution of ODE at (t+1).
      S_sol_prop = (sol_prop[2,2]) * popsize/pop0 ; 
      El_sol_prop = (sol_prop[2,3]) * popsize/pop0 ; 
      Il_sol_prop = (sol_prop[2,4]) * popsize/pop0 ; 
      R_sol_prop = (sol_prop[2,5]) * popsize/pop0;
      
      
      # Calculate the log priors and log likelihoods
      # For current value
      deltal_t_curr_log_postd <- deltal_log_postd(deltal = deltal_samples[[i-1]][t],
                                                  deltal_mean = deltal_mean_samples[i-1],
                                                  gammal = gammal_samples[[i-1]][t],
                                                  El_t1 = El_samples[[i-1]][t+1],
                                                  El_t1sol = El_sol_curr,
                                                  Il_t1 = Il_samples[[i-1]][t+1],
                                                  Il_t1sol = Il_sol_curr,
                                                  R_t1 = R_samples[[i-1]][t+1],
                                                  R_t1sol = R_sol_curr,
                                                  t,
                                                  obs1 = cases_new[t+1],
                                                  obs2 = cases_recovered[t+1],
                                                  sd_vec = sd_vec, 
                                                  tau1 = tau1_samples[i-1],
                                                  tau2 = tau2_samples[i-1]
      )
      
      # Calculate the log priors and log likelihoods
      # For proposed value
      deltal_t_prop_log_postd <- deltal_log_postd(deltal = deltal_t_prop,
                                                  gammal = gammal_samples[[i-1]][t],
                                                  deltal_mean = deltal_mean_samples[i-1],
                                                  El_t1 = El_samples[[i-1]][t+1],
                                                  El_t1sol =  El_sol_prop,
                                                  Il_t1 = Il_samples[[i-1]][t+1],
                                                  Il_t1sol = Il_sol_prop,
                                                  R_t1 = R_samples[[i-1]][t+1],
                                                  R_t1sol = R_sol_prop,
                                                  t,
                                                  obs1 = cases_new[t+1],
                                                  obs2 = cases_recovered[t+1],
                                                  sd_vec = sd_vec, 
                                                  tau1 = tau1_samples[i-1],
                                                  tau2 = tau2_samples[i-1]
      )
      
      log_q_delta_t_to_prop <- dnorm(deltal_t_prop, 
                                     mean = deltal_eda[t], 
                                     sd = sd_vec[2], log = TRUE)
      log_q_delta_prop_to_t <- dnorm(deltal_samples[[i-1]][t],
                                     mean = deltal_eda[t], 
                                     sd = sd_vec[2], log = TRUE)
      
      deltal_t_prop_log_postd2 = deltal_t_prop_log_postd + log_q_delta_prop_to_t
      deltal_t_curr_log_postd2 = deltal_t_curr_log_postd + log_q_delta_t_to_prop
      
      # Calculate the log acceptance probabilities
      accep_prob <- min(1, exp(deltal_t_prop_log_postd2 - deltal_t_curr_log_postd2))
      
      # Acceptance/rejection step for log_beta
      rand_u = runif(1) # Generate a random number from a uniform distribution
      
      
      
      if (!is.na(accep_prob) && rand_u < accep_prob) {
        deltal_samples[[i]][t] <- deltal_t_prop
        El_sol_samples[[i]][t+1] = El_sol_prop
        Il_sol_samples[[i]][t+1] = Il_sol_prop
        R_sol_samples[[i]][t+1] = R_sol_prop
      } else {
        deltal_samples[[i]][t] <- deltal_samples[[i-1]][t]
        El_sol_samples[[i]][t+1] = El_sol_curr
        Il_sol_samples[[i]][t+1] = Il_sol_curr
        R_sol_samples[[i]][t+1] = R_sol_curr
      }
    }
    
    ###################################################################################################################
    ###################################################################################################################
    ######################################### 3. update gammal_t #######################################################
    ###################################################################################################################
    ###################################################################################################################
    #Updating deltal[0] with time index [1] 
    #Propose new value from normal distribution with mean at previous iteration
    
    for (t in 1:(T-1)) {
      
      #Propose new value from normal distribution with mean at previous iteration
      gammal_t_prop <- rnorm(1, mean = gammal_eda[t], sd = prop_sd3)
      
      E0 = El_samples[[i-1]][t] * pop0/popsize # Initial number of exposed out of pop0
      I0 = Il_samples[[i-1]][t] * pop0/popsize # Initial number of infected out of pop0
      R0 = R_samples[[i-1]][t] * pop0/popsize # Initial number of recovered out of pop0
      V0 = sum(v_x[1:t]) * vax_eff * pop0/popsize # Initial number of vaccinated out of pop0
      W =  pop0 - (E0 + I0 + R0 + V0)
      N = W + E0 + I0 + R0 + V0
      state_values = c(S = W, E = E0, I = I0, R = R0, V = V0)
      
      #parameters of ODE at current deltal_t at i-1
      ode_state_parms_curr <- c(beta = exp(betal_samples[[i]][t]), 
                                delta = exp(deltal_samples[[i]][t]),
                                gamma = exp(gammal_samples[[i-1]][t]),
                                v_t = v_x[t] * pop0/popsize
      )
      
      #parameters of ODE at proposed value deltal_t_prop 
      ode_state_parms_prop <- c(beta = exp(betal_samples[[i]][t]), 
                                delta = exp(deltal_samples[[i]][t]),
                                gamma = exp(gammal_t_prop),
                                v_t = v_x[t] * pop0/popsize
      ) 
      
      
      #solution of ODE at current value deltal_t at i-1.
      sol_curr = ode_solution(state_values, parms = ode_state_parms_curr)
      S_sol_curr = (sol_curr[2,2]) * popsize/pop0 ; 
      El_sol_curr = (sol_curr[2,3]) * popsize/pop0 ; 
      Il_sol_curr = (sol_curr[2,4]) * popsize/pop0 ; 
      R_sol_curr = (sol_curr[2,5]) * popsize/pop0;
      
      #solution of ODE at proposed value betal_t at i-1.
      sol_prop = ode_solution(state_values, parms = ode_state_parms_prop)
      #S_sol is real solution of ODE at (t+1).
      S_sol_prop = (sol_prop[2,2]) * popsize/pop0 ; 
      El_sol_prop = (sol_prop[2,3]) * popsize/pop0 ; 
      Il_sol_prop = (sol_prop[2,4]) * popsize/pop0 ; 
      R_sol_prop = (sol_prop[2,5]) * popsize/pop0;
      
      
      # Calculate the log priors and log likelihoods
      # For current value
      gammal_t_curr_log_postd <- gammal_log_postd(deltal = deltal_samples[[i]][t],
                                                  gammal_mean = gammal_mean_samples[i-1],
                                                  gammal = gammal_samples[[i-1]][t],
                                                  El_t1 = El_samples[[i-1]][t+1],
                                                  El_t1sol = El_sol_curr,
                                                  Il_t1 = Il_samples[[i-1]][t+1],
                                                  Il_t1sol = Il_sol_curr,
                                                  R_t1 = R_samples[[i-1]][t+1],
                                                  R_t1sol = R_sol_curr,
                                                  t,
                                                  obs1 = cases_new[t+1],
                                                  obs2 = cases_recovered[t+1],
                                                  sd_vec = sd_vec, 
                                                  tau1 = tau1_samples[i-1],
                                                  tau2 = tau2_samples[i-1]
      )
      
      # Calculate the log priors and log likelihoods
      # For proposed value
      gammal_t_prop_log_postd <- gammal_log_postd(deltal = deltal_samples[[i]][t],
                                                  gammal = gammal_t_prop,
                                                  gammal_mean = gammal_mean_samples[i-1],
                                                  El_t1 = El_samples[[i-1]][t+1],
                                                  El_t1sol =  El_sol_prop,
                                                  Il_t1 = Il_samples[[i-1]][t+1],
                                                  Il_t1sol = Il_sol_prop,
                                                  R_t1 = R_samples[[i-1]][t+1],
                                                  R_t1sol = R_sol_prop,
                                                  t,
                                                  obs1 = cases_new[t+1],
                                                  obs2 = cases_recovered[t+1],
                                                  sd_vec = sd_vec, 
                                                  tau1 = tau1_samples[i-1],
                                                  tau2 = tau2_samples[i-1]
      )
      
      log_q_gamma_t_to_prop <- dnorm(gammal_t_prop, 
                                     mean = gammal_eda[t], 
                                     sd = sd_vec[3], log = TRUE)
      log_q_gamma_prop_to_t <- dnorm(gammal_samples[[i-1]][t],
                                     mean = gammal_eda[t], 
                                     sd = sd_vec[3], log = TRUE)
      
      gammal_t_prop_log_postd2 = gammal_t_prop_log_postd + log_q_gamma_prop_to_t
      gammal_t_curr_log_postd2 = gammal_t_curr_log_postd + log_q_gamma_t_to_prop
      
      # Calculate the log acceptance probabilities
      accep_prob <- min(1, exp(gammal_t_prop_log_postd2 - gammal_t_curr_log_postd2))
      
      # Acceptance/rejection step for log_beta
      rand_u = runif(1) # Generate a random number from a uniform distribution
      
      
      
      if (!is.na(accep_prob) && rand_u < accep_prob) {
        gammal_samples[[i]][t] <- gammal_t_prop
        El_sol_samples[[i]][t+1] = El_sol_prop
        Il_sol_samples[[i]][t+1] = Il_sol_prop
        R_sol_samples[[i]][t+1] = R_sol_prop
      } else {
        gammal_samples[[i]][t] <- gammal_samples[[i-1]][t]
        El_sol_samples[[i]][t+1] = El_sol_curr
        Il_sol_samples[[i]][t+1] = Il_sol_curr
        R_sol_samples[[i]][t+1] = R_sol_curr
      }
    }
    
    ####################################################################
    
    deltal_mean_samples[i] = deltal_mean_samples[i-1]
    gammal_mean_samples[i] = gammal_mean_samples[i-1]
    
    ####################################################################
    
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    ######################################### 6. update El_t #######################################################
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    
    { t=1
    # Propose new value from normal distribution with mean at previous iteration
    El_t_prop <- rnorm(1, mean = El_samples[[i-1]][t], sd = sd_vec[5])
    
    
    ode_parms <- c(beta = exp(betal_samples[[i]][t]),
                   delta = exp(deltal_samples[[i]][t]),
                   gamma = exp(gammal_samples[[i]][t]),
                   v_t = v_x[t] * pop0/popsize
    )
    
    # parameters of ODE at current El_t at i-1
    ode_state_parms_curr <- c(El_samples[[i-1]][t],
                              Il_samples[[i-1]][t],
                              R_samples[[i-1]][t],
                              sum(v_x[1:t])
    )
    
    # parameters of ODE at proposed value El_t_prop 
    ode_state_parms_prop <- c(El_t_prop,
                              Il_samples[[i-1]][t],
                              R_samples[[i-1]][t],
                              sum(v_x[1:t])
    )
    
    
    #solution of ODE at current value El_t at i-1.
    sol_curr = ode_solution_El_Il(parms = ode_parms,
                                  parms2 = ode_state_parms_curr
    )
    S_sol_curr = (sol_curr[2,2]) * popsize/pop0 ; 
    El_sol_curr = (sol_curr[2,3]) * popsize/pop0 ; 
    Il_sol_curr = (sol_curr[2,4]) * popsize/pop0 ; 
    R_sol_curr = (sol_curr[2,5]) * popsize/pop0;
    
    #solution of ODE at proposed value 
    sol_prop = ode_solution_El_Il(parms = ode_parms,
                                  parms2 = ode_state_parms_prop
    )
    #S_sol is real solution of ODE at (t+1).
    S_sol_prop = (sol_prop[2,2]) * popsize/pop0 ; 
    El_sol_prop = (sol_prop[2,3]) * popsize/pop0 ; 
    Il_sol_prop = (sol_prop[2,4]) * popsize/pop0 ; 
    R_sol_prop = (sol_prop[2,5]) * popsize/pop0;
    
    
    # Calculate the log priors and log likelihoods
    # For current value
    El_t_curr_log_postd <- El_log_postd(deltal = deltal_samples[[i]][t],
                                        gammal = gammal_samples[[i]][t],
                                        El_t = El_samples[[i-1]][t],
                                        El_tsol = El_eda[t], ###### NOTEEE !!!!!!
                                        El_t1 = El_samples[[i-1]][t+1],
                                        El_t1sol = El_sol_curr, 
                                        Il_t1 = Il_samples[[i-1]][t+1],
                                        Il_t1sol = Il_sol_curr, 
                                        R_t1 = R_samples[[i-1]][t+1],
                                        R_t1sol = R_sol_curr,  
                                        t,
                                        obs1 = cases_new[t+1],
                                        obs2 = cases_recovered[t+1],
                                        sd_vec = sd_vec, 
                                        tau1 = tau1_samples[i-1],
                                        tau2 = tau2_samples[i-1]
    )
    
    # Calculate the log priors and log likelihoods
    # For proposed value
    El_t_prop_log_postd <- El_log_postd(deltal = deltal_samples[[i]][t],
                                        gammal = gammal_samples[[i]][t],
                                        El_t = El_t_prop,
                                        El_tsol = El_eda[t], ###### NOTEEE !!!!!!
                                        El_t1 = El_samples[[i-1]][t+1],
                                        El_t1sol = El_sol_prop, 
                                        Il_t1 = Il_samples[[i-1]][t+1],
                                        Il_t1sol = Il_sol_prop, 
                                        R_t1 = R_samples[[i-1]][t+1],
                                        R_t1sol = R_sol_prop,
                                        t,
                                        obs1 = cases_new[t+1],
                                        obs2 = cases_recovered[t+1],
                                        sd_vec = sd_vec, 
                                        tau1 = tau1_samples[i-1],
                                        tau2 = tau2_samples[i-1]
    )
    # Calculate the log acceptance probabilities
    accep_prob <- min(1, exp(El_t_prop_log_postd - El_t_curr_log_postd))
    
    # Acceptance/rejection step for log_beta
    rand_u = runif(1) # Generate a random number from a uniform distribution
    
    
    
    if (!is.na(accep_prob) && rand_u < accep_prob) {
      El_samples[[i]][t] = El_t_prop
      El_sol_samples[[i]][t+1] = El_sol_prop
      Il_sol_samples[[i]][t+1] = Il_sol_prop
      R_sol_samples[[i]][t+1] = R_sol_prop
    } else {
      El_samples[[i]][t] = El_samples[[i-1]][t]
      El_sol_samples[[i]][t+1] = El_sol_curr
      Il_sol_samples[[i]][t+1] = Il_sol_curr
      R_sol_samples[[i]][t+1] = R_sol_curr
    }
    }
    
    ###################################################################################################################
    
    
    for (t in 2:(T-1)) {
      # Propose new value from normal distribution with mean at previous iteration
      El_t_prop <- rnorm(1, mean = El_samples[[i-1]][t], sd = sd_vec[5])
      
      
      ode_parms <- c(beta = exp(betal_samples[[i]][t]),
                     delta = exp(deltal_samples[[i]][t]),
                     gamma = exp(gammal_samples[[i]][t]),
                     v_t = v_x[t] * pop0/popsize
      )
      
      # parameters of ODE at current El_t at i-1
      ode_state_parms_curr <- c(El_samples[[i-1]][t],
                                Il_samples[[i-1]][t],
                                R_samples[[i-1]][t],
                                sum(v_x[1:t])
      )
      
      # parameters of ODE at proposed value El_t_prop 
      ode_state_parms_prop <- c(El_t_prop,
                                Il_samples[[i-1]][t],
                                R_samples[[i-1]][t],
                                sum(v_x[1:t])
      )
      
      
      #solution of ODE at current value El_t at i-1.
      sol_curr = ode_solution_El_Il(parms = ode_parms,
                                    parms2 = ode_state_parms_curr
      )
      S_sol_curr = (sol_curr[2,2]) * popsize/pop0 ; 
      El_sol_curr = (sol_curr[2,3]) * popsize/pop0 ; 
      Il_sol_curr = (sol_curr[2,4]) * popsize/pop0 ; 
      R_sol_curr = (sol_curr[2,5]) * popsize/pop0;
      
      #solution of ODE at proposed value 
      sol_prop = ode_solution_El_Il(parms = ode_parms,
                                    parms2 = ode_state_parms_prop
      )
      #S_sol is real solution of ODE at (t+1).
      S_sol_prop = (sol_prop[2,2]) * popsize/pop0 ; 
      El_sol_prop = (sol_prop[2,3]) * popsize/pop0 ; 
      Il_sol_prop = (sol_prop[2,4]) * popsize/pop0 ; 
      R_sol_prop = (sol_prop[2,5]) * popsize/pop0;
      
      
      # Calculate the log priors and log likelihoods
      # For current value
      El_t_curr_log_postd <- El_log_postd(deltal = deltal_samples[[i]][t],
                                          gammal = gammal_samples[[i]][t],
                                          El_t = El_samples[[i-1]][t],
                                          El_tsol = El_sol_samples[[i]][t], ###### NOTEEE !!!!!!
                                          El_t1 = El_samples[[i-1]][t+1],
                                          El_t1sol = El_sol_curr, 
                                          Il_t1 = Il_samples[[i-1]][t+1],
                                          Il_t1sol = Il_sol_curr, 
                                          R_t1 = R_samples[[i-1]][t+1],
                                          R_t1sol = R_sol_curr,
                                          t,
                                          obs1 = cases_new[t+1],
                                          obs2 = cases_recovered[t+1],
                                          sd_vec = sd_vec, 
                                          tau1 = tau1_samples[i-1],
                                          tau2 = tau2_samples[i-1]
      )
      
      # Calculate the log priors and log likelihoods
      # For proposed value
      El_t_prop_log_postd <- El_log_postd(deltal = deltal_samples[[i]][t],
                                          gammal = gammal_samples[[i]][t],
                                          El_t = El_t_prop,
                                          El_tsol = El_sol_samples[[i]][t], ###### NOTEEE !!!!!!
                                          El_t1 = El_samples[[i-1]][t+1],
                                          El_t1sol = El_sol_prop, 
                                          Il_t1 = Il_samples[[i-1]][t+1],
                                          Il_t1sol = Il_sol_prop, 
                                          R_t1 = R_samples[[i-1]][t+1],
                                          R_t1sol = R_sol_prop,
                                          t,
                                          obs1 = cases_new[t+1],
                                          obs2 = cases_recovered[t+1],
                                          sd_vec = sd_vec, 
                                          tau1 = tau1_samples[i-1],
                                          tau2 = tau2_samples[i-1]
      )
      # Calculate the log acceptance probabilities
      accep_prob <- min(1, exp(El_t_prop_log_postd - El_t_curr_log_postd))
      
      # Acceptance/rejection step for log_beta
      rand_u = runif(1) # Generate a random number from a uniform distribution
      
      
      if (!is.na(accep_prob) && rand_u < accep_prob) {
        El_samples[[i]][t] = El_t_prop
        El_sol_samples[[i]][t+1] = El_sol_prop
        Il_sol_samples[[i]][t+1] = Il_sol_prop
        R_sol_samples[[i]][t+1] = R_sol_prop
      } else {
        El_samples[[i]][t] = El_samples[[i-1]][t]
        El_sol_samples[[i]][t+1] = El_sol_curr
        Il_sol_samples[[i]][t+1] = Il_sol_curr
        R_sol_samples[[i]][t+1] = R_sol_curr
      }
    }
    
    ############################################################################################
    t=T
    # Propose new value from normal distribution with mean at previous iteration
    El_t_prop <- rnorm(1, mean = El_samples[[i-1]][t], sd = sd_vec[5])
    
    
    # Calculate the log priors and log likelihoods
    # For current value
    El_t_curr_log_postd <- El_log_postd_T(El_t = El_samples[[i-1]][t],
                                          El_tsol = El_sol_samples[[i]][t],
                                          sd_vec = sd_vec
    )
    
    # Calculate the log priors and log likelihoods
    # For proposed value
    El_t_prop_log_postd <- El_log_postd_T(El_t = El_t_prop,
                                          El_tsol = El_sol_samples[[i]][t],
                                          sd_vec = sd_vec
    )
    # Calculate the log acceptance probabilities
    accep_prob <- min(1, exp(El_t_prop_log_postd - El_t_curr_log_postd))
    
    # Acceptance/rejection step for log_beta
    rand_u = runif(1) # Generate a random number from a uniform distribution
    
    
    if (!is.na(accep_prob) && rand_u < accep_prob) {
      El_samples[[i]][t] = El_t_prop
    } else {
      El_samples[[i]][t] = El_samples[[i-1]][t]
    }
    
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    ######################################### 7. update Il_t #######################################################
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    
    {t=1
    # Propose new value from normal distribution with mean at previous iteration
    Il_t_prop <- rnorm(1, mean = Il_samples[[i-1]][t], sd = sd_vec[6])
    
    
    ode_parms <- c(beta = exp(betal_samples[[i]][t]), 
                   delta = exp(deltal_samples[[i]][t]),
                   gamma = exp(gammal_samples[[i]][t]),
                   v_t = v_x[t] * pop0/popsize
    )
    
    # parameters of ODE at current Il_t at i-1
    ode_state_parms_curr <- c(El_samples[[i]][t],
                              Il_samples[[i-1]][t],
                              R_samples[[i-1]][t],
                              sum(v_x[1:t])
    )
    
    # parameters of ODE at proposed value  
    ode_state_parms_prop <- c(El_samples[[i]][t],
                              Il_t_prop,
                              R_samples[[i-1]][t],
                              sum(v_x[1:t])
    )
    
    
    #solution of ODE at current value .
    sol_curr = ode_solution_El_Il(parms = ode_parms,
                                  parms2 = ode_state_parms_curr
    )
    S_sol_curr = (sol_curr[2,2]) * popsize/pop0 ; 
    El_sol_curr = (sol_curr[2,3]) * popsize/pop0 ; 
    Il_sol_curr = (sol_curr[2,4]) * popsize/pop0 ; 
    R_sol_curr = (sol_curr[2,5]) * popsize/pop0;
    
    #solution of ODE at proposed value .
    sol_prop = ode_solution_El_Il(parms = ode_parms,
                                  parms2 = ode_state_parms_prop
    )
    #S_sol is real solution of ODE at (t+1).
    S_sol_prop = (sol_prop[2,2]) * popsize/pop0 ; 
    El_sol_prop = (sol_prop[2,3]) * popsize/pop0 ; 
    Il_sol_prop = (sol_prop[2,4]) * popsize/pop0 ; 
    R_sol_prop = (sol_prop[2,5]) * popsize/pop0;
    
    
    # Calculate the log priors and log likelihoods
    # For current value
    Il_t_curr_log_postd <- Il_log_postd(deltal = deltal_samples[[i]][t],
                                        gammal = gammal_samples[[i]][t],
                                        Il_t = Il_samples[[i-1]][t],
                                        Il_tsol = Il_eda[t], ######## NOTE !!!!!!!!!!!!!
                                        Il_t1 = Il_samples[[i-1]][t+1],
                                        Il_t1sol = Il_sol_curr, 
                                        El_t1 = El_samples[[i-1]][t+1],
                                        El_t1sol = El_sol_curr, 
                                        R_t1 = R_samples[[i-1]][t+1],
                                        R_t1sol = R_sol_curr,
                                        t,
                                        obs1 = cases_new[t+1],
                                        obs2 = cases_recovered[t+1],
                                        sd_vec = sd_vec, 
                                        tau1 = tau1_samples[i-1],
                                        tau2 = tau2_samples[i-1]
    )
    
    
    # Calculate the log priors and log likelihoods
    # For proposed value
    Il_t_prop_log_postd <- Il_log_postd(deltal = deltal_samples[[i]][t],
                                        gammal = gammal_samples[[i]][t],
                                        Il_t = Il_t_prop,
                                        Il_tsol = Il_eda[t], ######## NOTE !!!!!!!!!!!!!
                                        Il_t1 = Il_samples[[i-1]][t+1],
                                        Il_t1sol = Il_sol_prop, 
                                        El_t1 = El_samples[[i-1]][t+1],
                                        El_t1sol = El_sol_prop, 
                                        R_t1 = R_samples[[i-1]][t+1],
                                        R_t1sol = R_sol_prop,
                                        t,
                                        obs1 = cases_new[t+1],
                                        obs2 = cases_recovered[t+1],
                                        sd_vec = sd_vec, 
                                        tau1 = tau1_samples[i-1],
                                        tau2 = tau2_samples[i-1]
    )
    # Calculate the log acceptance probabilities
    accep_prob <- min(1, exp(Il_t_prop_log_postd - Il_t_curr_log_postd))
    
    # Acceptance/rejection step for log_beta
    rand_u = runif(1) # Generate a random number from a uniform distribution
    
    if (!is.na(accep_prob) && rand_u < accep_prob) {
      Il_samples[[i]][t] = Il_t_prop
      El_sol_samples[[i]][t+1] = El_sol_prop
      Il_sol_samples[[i]][t+1] = Il_sol_prop
      R_sol_samples[[i]][t+1] = R_sol_prop
    } else {
      Il_samples[[i]][t] = Il_samples[[i-1]][t]
      El_sol_samples[[i]][t+1] = El_sol_curr
      Il_sol_samples[[i]][t+1] = Il_sol_curr
      R_sol_samples[[i]][t+1] = R_sol_curr
    }
    }
    
    ###################################################################################################################
    for (t in 2:(T-1)) {
      # Propose new value from normal distribution with mean at previous iteration
      Il_t_prop <- rnorm(1, mean = Il_samples[[i-1]][t], sd = sd_vec[6])
      
      
      ode_parms <- c(beta = exp(betal_samples[[i]][t]), 
                     delta = exp(deltal_samples[[i]][t]),
                     gamma = exp(gammal_samples[[i]][t]),
                     v_t = v_x[t] * pop0/popsize
      )
      
      # parameters of ODE at current Il_t at i-1
      ode_state_parms_curr <- c(El_samples[[i]][t],
                                Il_samples[[i-1]][t],
                                R_samples[[i-1]][t],
                                sum(v_x[1:t])
      )
      
      # parameters of ODE at proposed value  
      ode_state_parms_prop <- c(El_samples[[i]][t],
                                Il_t_prop,
                                R_samples[[i-1]][t],
                                sum(v_x[1:t])
      )
      
      
      #solution of ODE at current value .
      sol_curr = ode_solution_El_Il(parms = ode_parms,
                                    parms2 = ode_state_parms_curr
      )
      S_sol_curr = (sol_curr[2,2]) * popsize/pop0 ; 
      El_sol_curr = (sol_curr[2,3]) * popsize/pop0 ; 
      Il_sol_curr = (sol_curr[2,4]) * popsize/pop0 ; 
      R_sol_curr = (sol_curr[2,5]) * popsize/pop0;
      
      #solution of ODE at proposed value .
      sol_prop = ode_solution_El_Il(parms = ode_parms,
                                    parms2 = ode_state_parms_prop
      )
      #S_sol is real solution of ODE at (t+1).
      S_sol_prop = (sol_prop[2,2]) * popsize/pop0 ; 
      El_sol_prop = (sol_prop[2,3]) * popsize/pop0 ; 
      Il_sol_prop = (sol_prop[2,4]) * popsize/pop0 ; 
      R_sol_prop = (sol_prop[2,5]) * popsize/pop0;
      
      
      # Calculate the log priors and log likelihoods
      # For current value
      Il_t_curr_log_postd <- Il_log_postd(deltal = deltal_samples[[i]][t],
                                          gammal = gammal_samples[[i]][t],
                                          Il_t = Il_samples[[i-1]][t],
                                          Il_tsol = Il_sol_samples[[i]][t], ######## NOTE !!!!!!!!!!!!!
                                          Il_t1 = Il_samples[[i-1]][t+1],
                                          Il_t1sol = Il_sol_curr, 
                                          El_t1 = El_samples[[i-1]][t+1],
                                          El_t1sol = El_sol_curr, 
                                          R_t1 = R_samples[[i-1]][t+1],
                                          R_t1sol = R_sol_curr,
                                          t,
                                          obs1 = cases_new[t+1],
                                          obs2 = cases_recovered[t+1],
                                          sd_vec = sd_vec, 
                                          tau1 = tau1_samples[i-1],
                                          tau2 = tau2_samples[i-1]
      )
      
      
      # Calculate the log priors and log likelihoods
      # For proposed value
      Il_t_prop_log_postd <- Il_log_postd(deltal = deltal_samples[[i]][t],
                                          gammal = gammal_samples[[i]][t],
                                          Il_t = Il_t_prop,
                                          Il_tsol = Il_sol_samples[[i]][t], ######## NOTE !!!!!!!!!!!!!
                                          Il_t1 = Il_samples[[i-1]][t+1],
                                          Il_t1sol = Il_sol_prop, 
                                          El_t1 = El_samples[[i-1]][t+1],
                                          El_t1sol = El_sol_prop, 
                                          R_t1 = R_samples[[i-1]][t+1],
                                          R_t1sol = R_sol_prop,
                                          t,
                                          obs1 = cases_new[t+1],
                                          obs2 = cases_recovered[t+1],
                                          sd_vec = sd_vec, 
                                          tau1 = tau1_samples[i-1],
                                          tau2 = tau2_samples[i-1]
      )
      # Calculate the log acceptance probabilities
      accep_prob <- min(1, exp(Il_t_prop_log_postd - Il_t_curr_log_postd))
      
      # Acceptance/rejection step for log_beta
      rand_u = runif(1) # Generate a random number from a uniform distribution
      
      
      
      if (!is.na(accep_prob) && rand_u < accep_prob) {
        Il_samples[[i]][t] = Il_t_prop
        El_sol_samples[[i]][t+1] = El_sol_prop
        Il_sol_samples[[i]][t+1] = Il_sol_prop
        R_sol_samples[[i]][t+1] = R_sol_prop
      } else {
        Il_samples[[i]][t] = Il_samples[[i-1]][t]
        El_sol_samples[[i]][t+1] = El_sol_curr
        Il_sol_samples[[i]][t+1] = Il_sol_curr
        R_sol_samples[[i]][t+1] = R_sol_curr
      }
    }
    
    ############################################################################################
    t=T
    # Propose new value from normal distribution with mean at previous iteration
    Il_t_prop <- rnorm(1, mean = Il_samples[[i-1]][t], sd = sd_vec[5])
    
    
    # Calculate the log priors and log likelihoods
    # For current value
    Il_t_curr_log_postd <- Il_log_postd_T(Il_t = Il_samples[[i-1]][t],
                                          Il_tsol = Il_sol_samples[[i]][t],
                                          sd_vec = sd_vec
    )
    
    # Calculate the log priors and log likelihoods
    # For proposed value
    Il_t_prop_log_postd <- Il_log_postd_T(Il_t = Il_t_prop,
                                          Il_tsol = Il_sol_samples[[i]][t],
                                          sd_vec = sd_vec
    )
    # Calculate the log acceptance probabilities
    accep_prob <- min(1, exp(Il_t_prop_log_postd - Il_t_curr_log_postd))
    
    # Acceptance/rejection step for log_beta
    rand_u = runif(1) # Generate a random number from a uniform distribution
    
    Il_samples[[i]][t] = Il_samples[[i-1]][t]
    
    if (!is.na(accep_prob) && rand_u < accep_prob) {
      Il_samples[[i]][t] = Il_t_prop
    } else {
      Il_samples[[i]][t] = Il_samples[[i-1]][t]
    }
    
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    ######################################### 8. update Rl_t #######################################################
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    t= 1
    R_samples[[i]][t] = 0
    ###################################################################################################################
    for (t in 2:(T-1)) {
      # Propose new value from normal distribution with mean at previous iteration
      # R_t_prop <- rnorm(1, mean = R_samples[[i-1]][t], sd = sd_vec[7])
      R_t_prop <- rnorm(1, mean = R_samples[[i-1]][t], sd = sd_vec[7])
      
      if (R_t_prop < R_samples[[i]][t-1] || R_t_prop > R_samples[[i-1]][t+1]) {
        R_samples[[i]][t] = R_samples[[i-1]][t]  # Keep current value if restriction not met
      } else {
      
      ode_parms <- c(beta = exp(betal_samples[[i]][t]), 
                     delta = exp(deltal_samples[[i]][t]),
                     gamma = exp(gammal_samples[[i]][t]),
                     v_t = v_x[t] * pop0/popsize
      )
      
      # parameters of ODE at current Rl_t at i-1
      ode_state_parms_curr <- c(El_samples[[i]][t],
                                Il_samples[[i]][t],
                                R_samples[[i-1]][t],
                                sum(v_x[1:t])
      )
      
      # parameters of ODE at proposed value  
      ode_state_parms_prop <- c(El_samples[[i]][t],
                                Il_samples[[i]][t],
                                R_t_prop,
                                sum(v_x[1:t])
      )
      
      
      #solution of ODE at current value .
      sol_curr = ode_solution_El_Il(parms = ode_parms,
                                    parms2 = ode_state_parms_curr
      )
      S_sol_curr = (sol_curr[2,2]) * popsize/pop0 ; 
      El_sol_curr = (sol_curr[2,3]) * popsize/pop0 ; 
      Il_sol_curr = (sol_curr[2,4]) * popsize/pop0 ; 
      R_sol_curr = (sol_curr[2,5]) * popsize/pop0;
      
      #solution of ODE at proposed value .
      sol_prop = ode_solution_El_Il(parms = ode_parms,
                                    parms2 = ode_state_parms_prop
      )
      #S_sol is real solution of ODE at (t+1).
      S_sol_prop = (sol_prop[2,2]) * popsize/pop0 ; 
      El_sol_prop = (sol_prop[2,3]) * popsize/pop0 ; 
      Il_sol_prop = (sol_prop[2,4]) * popsize/pop0 ; 
      R_sol_prop = (sol_prop[2,5]) * popsize/pop0;
      
      
      # Calculate the log priors and log likelihoods
      # For current value
      R_t_curr_log_postd <- R_log_postd(deltal = deltal_samples[[i]][t],
                                        gammal = gammal_samples[[i]][t],
                                        R_t = R_samples[[i-1]][t],
                                        R_tsol = R_sol_samples[[i]][t], ######## NOTE !!!!!!!!!!!!!
                                        El_t1 = El_samples[[i-1]][t+1],
                                        El_t1sol = El_sol_curr, 
                                        Il_t1 = Il_samples[[i-1]][t+1],
                                        Il_t1sol = Il_sol_curr, 
                                        R_t1 = R_samples[[i-1]][t+1],
                                        R_t1sol = R_sol_curr,
                                        t,
                                        obs1 = cases_new[t+1],
                                        obs2 = cases_recovered[t+1],
                                        sd_vec = sd_vec, 
                                        tau1 = tau1_samples[i-1],
                                        tau2 = tau2_samples[i-1]
      )
      
      
      # Calculate the log priors and log likelihoods
      # For proposed value
      R_t_prop_log_postd <- R_log_postd(deltal = deltal_samples[[i]][t],
                                        gammal = gammal_samples[[i]][t],
                                        R_t = R_t_prop,
                                        R_tsol = R_sol_samples[[i]][t], ######## NOTE !!!!!!!!!!!!!
                                        El_t1 = El_samples[[i-1]][t+1],
                                        El_t1sol = El_sol_prop, 
                                        Il_t1 = Il_samples[[i-1]][t+1],
                                        Il_t1sol = Il_sol_prop, 
                                        R_t1 = R_samples[[i-1]][t+1],
                                        R_t1sol = R_sol_prop,
                                        t,
                                        obs1 = cases_new[t+1],
                                        obs2 = cases_recovered[t+1],
                                        sd_vec = sd_vec, 
                                        tau1 = tau1_samples[i-1],
                                        tau2 = tau2_samples[i-1]
      )
      

      # Calculate the log acceptance probabilities
      accep_prob <- min(1, exp(R_t_prop_log_postd - R_t_curr_log_postd))
      
      # Acceptance/rejection step for log_beta
      rand_u = runif(1) # Generate a random number from a uniform distribution
      
      
      
      if (!is.na(accep_prob) && rand_u < accep_prob) {
        
        R_samples[[i]][t] = R_t_prop
        El_sol_samples[[i]][t+1] = El_sol_prop
        Il_sol_samples[[i]][t+1] = Il_sol_prop
        R_sol_samples[[i]][t+1] = R_sol_prop
      } else {
        R_samples[[i]][t] = R_samples[[i-1]][t]
        El_sol_samples[[i]][t+1] = El_sol_curr
        Il_sol_samples[[i]][t+1] = Il_sol_curr
        R_sol_samples[[i]][t+1] = R_sol_curr
      }
      }
    }
    ############################################################################################
    t=T
    # Propose new value from normal distribution with mean at previous iteration
    R_t_prop <- rnorm(1, mean = R_samples[[i-1]][t], sd = sd_vec[7])
    
    if (R_t_prop < R_samples[[i]][t-1]) {
      R_samples[[i]][t] = R_samples[[i-1]][t]  # Keep current value if restriction not met
    } else {
      
    # Calculate the log priors and log likelihoods
    # For current value
    R_t_curr_log_postd <- R_log_postd_T(R_t = R_samples[[i-1]][t],
                                        R_tsol = R_sol_samples[[i]][t],
                                        sd_vec = sd_vec
    )
    
    
    # Calculate the log priors and log likelihoods
    # For proposed value
    R_t_prop_log_postd <- R_log_postd_T(R_t = R_t_prop,
                                        R_tsol = R_sol_samples[[i]][t],
                                        sd_vec = sd_vec
    )
    # Calculate the log acceptance probabilities
    accep_prob <- min(1, exp(R_t_prop_log_postd - R_t_curr_log_postd))
    
    # Acceptance/rejection step for log_beta
    rand_u = runif(1) # Generate a random number from a uniform distribution
    
    R_samples[[i]][t] = R_samples[[i-1]][t]
    
    if (!is.na(accep_prob) && rand_u < accep_prob) {
      
      R_samples[[i]][t] = R_t_prop
    }
    
    }
    tau1_samples[i] = tau1_samples[i-1]
    tau2_samples[i] = tau2_samples[i-1]
    # ###################################################################################################################
    # ###################################################################################################################
    # ########################################   update tau1/tau2   #####################################################
    # ###################################################################################################################
    # ###################################################################################################################
    # E0 = El_samples[[i]][1] * pop0/popsize # Initial number of exposed out of pop0
    # I0 = Il_samples[[i]][1] * pop0/popsize # Initial number of infected out of pop0
    # R0 = 0
    # V0 = 0
    # W =  pop0 - (E0 + I0 + R0 + V0)
    # N = W + E0 + I0 + R0 + V0
    # state_values = c(S = W, E = E0, I = I0, R = R0, V = V0)
    # 
    # parms_vec = list(delta_x = exp(deltal_samples[[i]]),
    #                  beta_x = exp(betal_samples[[i]]),
    #                  gamma_x = exp(gammal_samples[[i]]),
    #                  v_e = vax_eff,
    #                  v_t = v_x[1:T] * pop0/popsize
    # )
    # 
    # sol = ode_solution_for_vec(
    #   state_values = state_values,
    #   parms_vec = parms_vec,
    #   time_vec = seq(0, (T-1), by = 1)
    # )
    # 
    # #S_sol is real solution of ODE at (t+1).
    # S_sol = (sol[,2]) * popsize/pop0 ;
    # El_sol = (sol[,3]) * popsize/pop0 ;
    # Il_sol = (sol[,4]) * popsize/pop0 ;
    # R_sol = (sol[,5]) * popsize/pop0;
    # ###################################################################################################################
    # ###################################################################################################################
    # ########################################   update tau1   ##########################################################
    # ###################################################################################################################
    # ###################################################################################################################
    # 
    # shape_prop = (exp(tau1_samples[1])^2)/(sd_vec['tau_sd1']^2)
    # rate_prop = (exp(tau1_samples[1]))/(sd_vec['tau_sd1']^2)
    # 
    # tau1_prop <- rgamma(1, shape = shape_prop , rate = rate_prop)
    # # hist(rgamma(100000, shape = shape_prop , rate = rate_prop))
    # 
    # # Calculate the log priors and log likelihoods
    # # For current value
    # tau1_curr_log_postd <- tau1_log_postd(tau1_x = exp(tau1_samples[i-1]),
    #                                       shape_p = gamma_tau1_shape,
    #                                       rate_p = gamma_tau1_rate,
    #                                       obs1 = cases_new[2:T],
    #                                       deltal = deltal_samples[[i]][1:T-1],
    #                                       gammal = gammal_samples[[i]][1:T-1],
    #                                       El_sol = El_sol[2:T]
    # )
    # 
    # # Calculate the log priors and log likelihoods
    # # For proposed value
    # tau1_prop_log_postd <- tau1_log_postd(tau1_x = tau1_prop,
    #                                       shape_p = gamma_tau1_shape,
    #                                       rate_p = gamma_tau1_rate,
    #                                       obs1 = cases_new[2:T],
    #                                       deltal = deltal_samples[[i]][1:T-1],
    #                                       gammal = gammal_samples[[i]][1:T-1],
    #                                       El_sol = El_sol[2:T]
    # )
    # 
    # shape_curr = (tau1_prop^2)/(sd_vec['tau_sd1']^2)
    # rate_curr = (tau1_prop)/(sd_vec['tau_sd1']^2)
    # 
    # log_q_tau1_curr_to_prop <- dgamma(tau1_prop,
    #                                   shape = shape_prop,
    #                                   rate = rate_prop, log = TRUE)
    # log_q_tau1_prop_to_curr <- dgamma(exp(tau1_samples[i-1]),
    #                                   shape = shape_curr,
    #                                   rate = rate_curr, log = TRUE)
    # 
    # tau1_prop_log_postd2 = tau1_prop_log_postd + log_q_tau1_prop_to_curr
    # tau1_curr_log_postd2 = tau1_curr_log_postd + log_q_tau1_curr_to_prop
    # 
    # # Calculate the log acceptance probabilities
    # tau1_accep_prob <- min(1, exp(tau1_prop_log_postd2 - tau1_curr_log_postd2))
    # 
    # # Acceptance/rejection step for log_beta
    # rand_u = runif(1)
    # 
    # 
    # 
    # if (!is.na(tau1_accep_prob) && rand_u < tau1_accep_prob) {
    #   tau1_samples[i] <- log(tau1_prop)
    # } else {
    #   tau1_samples[i] <- tau1_samples[i-1]
    # }
    # 
    # 
    # 
    # ###################################################################################################################
    # ###################################################################################################################
    # ########################################   update tau2   ##########################################################
    # ###################################################################################################################
    # ###################################################################################################################
    # shape_prop = (exp(tau2_samples[1])^2)/(sd_vec['tau_sd2']^2)
    # rate_prop = (exp(tau2_samples[1]))/(sd_vec['tau_sd2']^2)
    # 
    # tau2_prop <- rgamma(1, shape = shape_prop , rate = rate_prop)
    # # hist(rgamma(100000, shape = shape_prop , rate = rate_prop))
    # 
    # # Calculate the log priors and log likelihoods
    # # For current value
    # tau2_curr_log_postd <- tau2_log_postd(tau2_x = exp(tau2_samples[i-1]),
    #                                       shape_p = gamma_tau2_shape,
    #                                       rate_p = gamma_tau2_rate,
    #                                       obs2 = cases_recovered[2:T],
    #                                       deltal = deltal_samples[[i]][1:T-1],
    #                                       gammal = gammal_samples[[i]][1:T-1],
    #                                       Il_sol = El_sol[2:T]
    # )
    # 
    # # Calculate the log priors and log likelihoods
    # # For proposed value
    # tau2_prop_log_postd <- tau2_log_postd(tau2_x = tau2_prop,
    #                                       shape_p = gamma_tau2_shape,
    #                                       rate_p = gamma_tau2_rate,
    #                                       obs2 = cases_recovered[2:T],
    #                                       deltal = deltal_samples[[i]][1:T-1],
    #                                       gammal = gammal_samples[[i]][1:T-1],
    #                                       Il_sol = El_sol[2:T]
    # )
    # 
    # shape_curr = (tau2_prop^2)/(sd_vec['tau_sd2']^2)
    # rate_curr = (tau2_prop)/(sd_vec['tau_sd2']^2)
    # 
    # log_q_tau2_curr_to_prop <- dgamma(tau2_prop,
    #                                   shape = shape_prop,
    #                                   rate = rate_prop, log = TRUE)
    # log_q_tau2_prop_to_curr <- dgamma(exp(tau2_samples[i-1]),
    #                                   shape = shape_curr,
    #                                   rate = rate_curr, log = TRUE)
    # 
    # tau2_prop_log_postd2 = tau2_prop_log_postd + log_q_tau2_prop_to_curr
    # tau2_curr_log_postd2 = tau2_curr_log_postd + log_q_tau2_curr_to_prop
    # 
    # # Calculate the log acceptance probabilities
    # tau2_accep_prob <- min(1, exp(tau2_prop_log_postd2 - tau2_curr_log_postd2))
    # 
    # # Acceptance/rejection step for log_beta
    # rand_u = runif(1)
    # 
    # 
    # 
    # if (!is.na(tau2_accep_prob) && rand_u < tau2_accep_prob) {
    #   tau2_samples[i] <- log(tau2_prop)
    # } else {
    #   tau2_samples[i] <- tau2_samples[i-1]
    # }
    
    ####################################################################
    ####################################################################
    ####################################################################
    ###################################################################################################################
    ###################################################################################################################
    ######################################## update sigma^2  ##########################################################
    ###################################################################################################################
    ###################################################################################################################
    
    sum_squared_diffs <- sum((diff(betal_samples[[i]]))^2)
    # Update the parameters for the inverse gamma posterior
    invg_shape_post <- invg_shape + length(betal_samples[[i]])/2
    invg_rate_post <- invg_rate + 0.5 * sum_squared_diffs
    
    # Sample a new value for sigma^2 from the inverse gamma posterior
    betal_rw_var_new <- 1 / rgamma(1, shape = invg_shape_post, rate = invg_rate_post)
    
    sd_vec[10] = betal_rw_var_new
    betal_rw_var_samples[i] = betal_rw_var_new
    
    
    
    setTxtProgressBar(pb, i)
    
    # Save entire environment every 10,000 iterations
    if (i > num_iter_start && i %% save_interval == 0) {
      save_file <- paste0("eda_saved_environment_", i, ".RData")
      save(list = ls(all.names = TRUE), file = save_file)
      cat("Saved entire environment at iteration", i, "\n")
    }
    
    if (i %% 10000 == 0) {
      filename <- paste0("1_beta_plot_iteration_", i, ".png")
      png(filename)
      plot(betal_eda[1:250], col = "red", type = "l", lwd=3, ylim = c(-3,-1))
      lines(betal_samples[[i]][1:250], lwd=3, col = "blue", type = "l")
      dev.off()
      filename2 <- paste0("2_delta_plot_iteration_", i, ".png")
      png(filename2)
      plot(deltal_eda[1:250], col = "red", type = "l", lwd=3, ylim = c(-3,-1))
      lines(deltal_samples[[i]][1:250], lwd=3 , col = "blue", type = "l")
      dev.off()
      filename3 <- paste0("3_gamma_plot_iteration_", i, ".png")
      png(filename3)
      plot(gammal_eda[1:250], col = "red", type = "l", lwd=3, ylim = c(-3,-1))
      lines(gammal_samples[[i]][1:250], lwd=3 , col = "blue", type = "l")
      dev.off()
      ##################################################
      post_samp <- exp(tau1_samples[1:i])
      filename4 <- paste0("4_tau1_", i, ".png")
      png(filename4)
      hist(post_samp)
      dev.off()
      ##################################################
      post_samp <- exp(tau2_samples[1:i])
      filename5 <- paste0("5_tau2_", i, ".png")
      png(filename5)
      hist(post_samp)
      dev.off()
      # ##################################################
      post_samp <- betal_rw_var_samples[1:i]
      filename6 <- paste0("6_betal_rw_var_", i, ".png")
      png(filename6)
      hist(post_samp)
      dev.off()
    }
    
  }    
  # Close the progress bar
  close(pb)
  
  # save the entire environment
  save(list = ls(all.names = TRUE), file=paste0("eda_K_", num_iter_end, ".RData"))
}

