###########################################################################################
###########################################################################################
############################# SEIRV MODEL AND SOLUTION FUNCTIONS ##########################
###########################################################################################
###########################################################################################

#The seirv model with variable 
seirv_1 <- function(timepoint, state_values, parms) {
  # create state variables (local variables)
  S = state_values[1]  # Susceptible
  E = state_values[2]  # Exposed
  I = state_values[3]  # Infectious asymptomatic
  R = state_values[4]  # Recovered symptomatic confirmed cases
  V = state_values[5]  # Vaccinated
  N = S + E + I + R + V
  with(as.list(parms), {  # variable names within parameters can be used 
    incidents = beta * S * I/N
    dS = -incidents - v_t  * (0.75)
    dE = incidents - delta * E
    dI = delta * E - gamma * I
    dR = gamma * I
    dV = v_t * (0.75)
    results = c(dS, dE, dI, dR, dV)
    return(list(results))
  })
}


# solution of ODE. The state values indices are same for all parameters.
ode_solution <- function(state_values,parms){
  #initial values for ODE
  
  times = c(0,1)
  results = lsoda(state_values, times, seirv_1, parms = parms)
  return(results)
}


ode_solution_El_Il <- function(parms,parms2){
  #initial values for ODE
  E0 = parms2[1] * pop0/popsize # Initial number of exposed out of pop0
  I0 = parms2[2] * pop0/popsize # Initial number of infected out of pop0
  R0 = parms2[3] * pop0/popsize # Initial number of recovered out of pop0
  V0 = parms2[4] * vax_eff * pop0/popsize # Initial number of vaccinated out of pop0
  W =  pop0 - (E0 + I0 + R0 + V0)
  N = W + E0 + I0 + R0 + V0
  state_values = c(S = W, E = E0, I = I0, R = R0, V = V0)
  
  times = c(0,1)
  results = lsoda(state_values, times, seirv_1, parms = parms)
  return(results)
}
###########################################################################################
###########################################################################################
############################# LOG LIKELIHOOD FUNCTION #####################################
###########################################################################################
###########################################################################################

# negative binomial log likelihood density
# El_sol, and Il_sol are re scaled back to original population. 
log_lik_dens <- function(obs1,
                         obs2,
                         deltal,
                         gammal,
                         El_t1sol,
                         Il_t1sol, 
                         tau1,
                         tau2
) 
{
  mu_obs1 = exp(deltal) * (El_t1sol)
  mu_obs2 = exp(gammal) * (Il_t1sol)
  
  # Calculate the log likelihood
  output = sum(dnbinom(obs1, size = 1/exp(tau1), 
                       mu = mu_obs1, log = TRUE),
               dnbinom(obs2, size = 1/exp(tau2), 
                       mu = mu_obs2, log = TRUE)
  )
  return(output)
}

###########################################################################################
###########################################################################################
#################### LOG PRIOR, LOG POSTERIOR FUNCTIONS FOR BETAL #########################
###########################################################################################
###########################################################################################
# Random walk log prior density with other parameters fixed at t 
# sum log liklihood and log prior from previous function 
# and showing the result for reducing code length
betal_log_postd_0 <- function(betal, 
                            betal_next, 
                            El_t1, 
                            El_t1sol, 
                            Il_t1, 
                            Il_t1sol, 
                            R_t1, 
                            R_t1sol,
                            deltal, 
                            gammal, 
                            t,
                            obs1,
                            obs2,
                            sd_vec, 
                            tau1,
                            tau2
) 
{
  
  # calculate the log likelihood
  log_lik = log_lik_dens(obs1 = obs1, # global variable
                         obs2 = obs2, # global variable
                         deltal, 
                         gammal, 
                         El_t1sol,
                         Il_t1sol,
                         tau1,  # global variable
                         tau2   # global variable
  )
  # calculate the prior
  betal_log_prior = sum(
    dnorm(betal_next, mean = betal, sd = sd_vec[10], log = TRUE),
    dnorm(El_t1, mean = El_t1sol, sd = sd_vec[12], log = TRUE),
    dnorm(Il_t1, mean = Il_t1sol, sd = sd_vec[13], log = TRUE),
    dnorm(R_t1, mean = R_t1sol, sd = sd_vec[14], log = TRUE)
  )
  
  # calculate the log posterior
  log_postd = log_lik + betal_log_prior
  
  return(log_postd)
}
#############################################
# sum log liklihood and log prior from previous function 
# and showing the result for reducing code length
betal_log_postd <- function(betal, 
                            betal_prev, 
                            betal_next, 
                            El_t1, 
                            El_t1sol, 
                            Il_t1, 
                            Il_t1sol, 
                            R_t1, 
                            R_t1sol, 
                            deltal, 
                            gammal, 
                            t,
                            obs1,
                            obs2,
                            sd_vec, 
                            tau1,
                            tau2
) 
{
  
  # calculate the log likelihood
  log_lik = log_lik_dens(obs1 = obs1, # global variable
                         obs2 = obs2, # global variable
                         deltal, 
                         gammal, 
                         El_t1sol,
                         Il_t1sol,
                         tau1,  # global variable
                         tau2   # global variable
  )
  # calculate the prior
  betal_log_prior = sum(
    dnorm(betal, mean = betal_prev, sd = sd_vec[10], log = TRUE),
    dnorm(betal_next, mean = betal, sd = sd_vec[10], log = TRUE),
    dnorm(El_t1, mean = El_t1sol, sd = sd_vec[12], log = TRUE),
    dnorm(Il_t1, mean = Il_t1sol, sd = sd_vec[13], log = TRUE),
    dnorm(R_t1, mean = R_t1sol, sd = sd_vec[14], log = TRUE)
  )
  
  # calculate the log posterior
  log_postd = log_lik + betal_log_prior
  
  return(log_postd)
}

#############################################
# sum log liklihood and log prior from previous function 
# and showing the result for reducing code length
betal_log_postd_t <- function(betal, 
                            betal_prev, 
                            El_t1, 
                            El_t1sol, 
                            Il_t1, 
                            Il_t1sol, 
                            R_t1, 
                            R_t1sol, 
                            deltal, 
                            gammal, 
                            t,
                            obs1,
                            obs2,
                            sd_vec, 
                            tau1,
                            tau2
) 
{
  
  # calculate the log likelihood
  log_lik = log_lik_dens(obs1 = obs1, # global variable
                         obs2 = obs2, # global variable
                         deltal, 
                         gammal, 
                         El_t1sol,
                         Il_t1sol,
                         tau1,  # global variable
                         tau2   # global variable
  )
  # calculate the prior
  betal_log_prior = sum(
    dnorm(betal, mean = betal_prev, sd = sd_vec[10], log = TRUE),
    dnorm(El_t1, mean = El_t1sol, sd = sd_vec[12], log = TRUE),
    dnorm(Il_t1, mean = Il_t1sol, sd = sd_vec[13], log = TRUE),
    dnorm(R_t1, mean = R_t1sol, sd = sd_vec[14], log = TRUE)
  )
  
  # calculate the log posterior
  log_postd = log_lik + betal_log_prior
  
  return(log_postd)
}
###########################################################################################
###########################################################################################
#################### LOG PRIOR, LOG POSTERIOR FUNCTIONS FOR DELTAL ########################
###########################################################################################
###########################################################################################
# sum log liklihood and log prior from previous function 
# and showing the result for reducing code length
deltal_log_postd <- function(deltal,
                             gammal,
                             deltal_mean,
                             El_t1,
                             El_t1sol, 
                             Il_t1, 
                             Il_t1sol, 
                             R_t1, 
                             R_t1sol,
                             t,
                             obs1,
                             obs2,
                             sd_vec,
                             tau1,
                             tau2
) 
{
  # calculate the log likelihood
  log_lik = log_lik_dens(obs1 = obs1, # global variable
                         obs2 = obs2, # global variable
                         deltal, 
                         gammal, 
                         El_t1sol,
                         Il_t1sol,
                         tau1,  # global variable
                         tau2   # global variable
  )
  
  # calculate the prior
  deltal_log_prior = sum(dnorm(deltal, mean = deltal_mean, 
                               sd = sd_vec[2], log = TRUE),
                         dnorm((El_t1), mean = El_t1sol, 
                               sd = sd_vec[12], log = TRUE),
                         dnorm((Il_t1), mean = Il_t1sol, 
                               sd = sd_vec[13], log = TRUE),
                         dnorm((R_t1), mean = R_t1sol, 
                               sd = sd_vec[14], log = TRUE)
  )
  
  # calculate the log posterior
  log_postd = log_lik + deltal_log_prior
  
  return(log_postd)
}

###########################################################################################
###########################################################################################
#################### LOG PRIOR, LOG POSTERIOR FUNCTIONS FOR GAMMAL ########################
###########################################################################################
###########################################################################################

# sum log liklihood and log prior from previous function 
# and showing the result for reducing code length
gammal_log_postd <- function(deltal,
                             gammal,
                             gammal_mean,
                             El_t1,
                             El_t1sol, 
                             Il_t1, 
                             Il_t1sol, 
                             R_t1, 
                             R_t1sol,
                             t,
                             obs1,
                             obs2,
                             sd_vec,
                             tau1,
                             tau2
) 
{
  # calculate the log likelihood
  log_lik = log_lik_dens(obs1 = obs1, # global variable
                         obs2 = obs2, # global variable
                         deltal, 
                         gammal, 
                         El_t1sol,
                         Il_t1sol,
                         tau1,  # global variable
                         tau2   # global variable
  )
  # calculate the prior
  
  gammal_log_prior = sum(dnorm(gammal, mean = gammal_mean, 
                               sd = sd_vec[3], log = TRUE),
                         dnorm((El_t1), mean = El_t1sol, 
                               sd = sd_vec[12], log = TRUE),
                         dnorm((Il_t1), mean = Il_t1sol, 
                               sd = sd_vec[13], log = TRUE),
                         dnorm((R_t1), mean = R_t1sol, 
                               sd = sd_vec[14], log = TRUE)
  )

  # calculate the log posterior
  log_postd = log_lik + gammal_log_prior
  
  return(log_postd)
}


###########################################################################################
###########################################################################################
################ LOG PRIOR, LOG POSTERIOR FUNCTIONS FOR DELTAL_MEAN #######################
###########################################################################################
###########################################################################################


deltal_mean_log_postd <- function(deltal, 
                                  deltal_mean,
                                  deltal_eda,
                                  sd_vec
)
{
  output = sum(dnorm(deltal_mean, mean = deltal, 
                     sd = sd_vec[8], log = TRUE),
               dnorm(deltal_mean, mean = deltal_eda,
                     sd = sd_vec[8], log = TRUE)
  )
  return(output)
}

###########################################################################################
###########################################################################################
################ LOG PRIOR, LOG POSTERIOR FUNCTIONS FOR GAMMAL_MEAN #######################
###########################################################################################
###########################################################################################

gammal_mean_log_postd <- function(gammal, 
                                  gammal_mean,
                                  gammal_eda,
                                  sd_vec
)
{
  output = sum(dnorm(gammal, mean = gammal_mean,
                     sd = sd_vec[9], log = TRUE),
               dnorm(gammal_mean, mean = gammal_eda, 
                     sd = sd_vec[9], log = TRUE)
  )
  return(output)
}


###########################################################################################
###########################################################################################
################ LOG PRIOR, LOG POSTERIOR FUNCTIONS FOR El_t #######################
###########################################################################################
###########################################################################################
###########################################################################################
# log prior plus log likelihood 
El_log_postd <- function(deltal,
                         gammal,
                         El_t,
                         El_tsol, 
                         El_t1, 
                         El_t1sol, 
                         Il_t1, 
                         Il_t1sol, 
                         R_t1, 
                         R_t1sol,
                         t,
                         obs1,
                         obs2,
                         sd_vec, 
                         tau1,
                         tau2
) 
{
  # calculate the log likelihood
  log_lik = log_lik_dens(obs1 = obs1, # global variable
                         obs2 = obs2, # global variable
                         deltal, 
                         gammal,
                         El_t1sol,
                         Il_t1sol,
                         tau1,  # global variable
                         tau2   # global variable
  )
  
  # calculate the prior
  El_log_prior = sum(dnorm(El_t, mean = El_tsol,
                           sd = sd_vec[12], log = TRUE),
                     dnorm(El_t1, mean = El_t1sol, 
                           sd = sd_vec[12], log = TRUE),
                     dnorm(Il_t1, mean = Il_t1sol, 
                           sd = sd_vec[13], log = TRUE),
                     dnorm(R_t1, mean = R_t1sol, 
                           sd = sd_vec[14], log = TRUE)
  )
  
  # calculate the log posterior
  log_postd = log_lik + El_log_prior
  
  return(log_postd)
}
###########################################################################################
El_log_postd_T <- function(El_t,
                           El_tsol,
                           sd_vec
)
{
  El_log_prior = dnorm(El_t, mean = El_tsol, 
                       sd = sd_vec[12], log = TRUE)
  
  # calculate the log posterior
  log_postd = El_log_prior
  
  return(log_postd) 
}
###########################################################################################
###########################################################################################
################ LOG PRIOR, LOG POSTERIOR FUNCTIONS FOR Il_t ##############################
###########################################################################################
###########################################################################################
###########################################################################################
# log prior + log likelihood 
Il_log_postd <- function(deltal,
                         gammal,
                         Il_t, 
                         Il_tsol, 
                         El_t1, 
                         El_t1sol, 
                         Il_t1, 
                         Il_t1sol, 
                         R_t1, 
                         R_t1sol,
                         t,
                         obs1,
                         obs2,
                         sd_vec, 
                         tau1,
                         tau2
) 
{
  # calculate the prior
  Il_log_prior = sum(dnorm(Il_t, mean = Il_tsol, 
                           sd = sd_vec[12], log = TRUE),
                     dnorm(El_t1, mean = El_t1sol, 
                           sd = sd_vec[12], log = TRUE),
                     dnorm(Il_t1, mean = Il_t1sol, 
                           sd = sd_vec[13], log = TRUE),
                     dnorm(R_t1, mean = R_t1sol, 
                           sd = sd_vec[14], log = TRUE)
  )
  
  # calculate the log likelihood
  log_lik = log_lik_dens(obs1 = obs1, # global variable
                         obs2 = obs2, # global variable
                         deltal, 
                         gammal, 
                         El_t1sol,
                         Il_t1sol,
                         tau1,  # global variable
                         tau2   # global variable
  )
  # calculate the log posterior
  log_postd = log_lik + Il_log_prior
  
  return(log_postd)
}

###########################################################################################
Il_log_postd_T <- function(Il_t,
                           Il_tsol,
                           sd_vec
)
{
  Il_log_prior = dnorm(Il_t, mean = Il_tsol, 
                       sd = sd_vec[13], log = TRUE)
  
  # calculate the log posterior
  log_postd = Il_log_prior
  
  return(log_postd) 
}


###########################################################################################
###########################################################################################
################ LOG PRIOR, LOG POSTERIOR FUNCTIONS FOR Rl_t ##############################
###########################################################################################
###########################################################################################
R_log_postd <- function (deltal, 
                         gammal,
                         R_t,
                         R_tsol, 
                         El_t1, 
                         El_t1sol, 
                         Il_t1, 
                         Il_t1sol, 
                         R_t1, 
                         R_t1sol, 
                         t,
                         sd_vec,
                         obs1, 
                         obs2, 
                         tau1,
                         tau2
) 
{
  # calculate the log likelihood
  log_lik = log_lik_dens(obs1, 
                         obs2,
                         deltal, 
                         gammal, 
                         El_t1sol,
                         Il_t1sol,
                         tau1,  # global variable
                         tau2   # global variable
  )
  # calculate the log posterior
  
  # R_log_prior = sum(dnorm(R_t, mean = R_tsol, 
  #                    sd = sd_vec[14], log = TRUE),
  #              dnorm(El_t1, mean = El_t1sol, 
  #                    sd = sd_vec[12], log = TRUE),
  #              dnorm(Il_t1, mean = Il_t1sol, 
  #                    sd = sd_vec[13], log = TRUE),
  #              dnorm(R_t1, mean = R_t1sol, 
  #                    sd = sd_vec[14], log = TRUE)
  # )
  R_log_prior = sum(dnorm(R_t, mean = R_tsol, 
                          sd = sd_vec[14], log = TRUE),
                    dnorm(R_t1, mean = R_t1sol, 
                          sd = sd_vec[14], log = TRUE)
  )
  log_postd = log_lik + R_log_prior
  
  return(log_postd)
}


###########################################################################################
R_log_postd_T <- function(R_t,
                          R_tsol,
                          sd_vec
)
{
  R_log_prior = dnorm(R_t, mean = R_tsol, 
                      sd = sd_vec[14], log = TRUE)
  
  # calculate the log posterior
  log_postd = R_log_prior
  
  return(log_postd) 
}


#######################################################################################
#######################################################################################
#######################################################################################
seirv_for_vec = function(timepoint, 
                         state_values, 
                         parms)
{
  S = state_values[1]
  E = state_values[2]       
  I = state_values[3]       
  R = state_values[4]       
  V = state_values[5]       
  N = S + E + I + R + V
  
  beta_x = parms$beta_x
  gamma_x = parms$gamma_x
  delta_x = parms$delta_x
  v_t = parms$v_t
  
  # beta = beta_x[floor(timepoint+1)]
  # delta = delta_x[floor(timepoint+1)]
  # gamma = gamma_x[floor(timepoint+1)]
  # Ensure we don't use the last NA value
  if (timepoint + 1 > length(beta_x) || is.na(beta_x[timepoint + 1])) {
    beta = beta_x[timepoint]  # Use the last valid value instead
  } else {
    beta = beta_x[timepoint + 1]
  }
  
  if (timepoint + 1 > length(delta_x) || is.na(delta_x[timepoint + 1])) {
    delta = delta_x[timepoint]
  } else {
    delta = delta_x[timepoint + 1]
  }
  
  if (timepoint + 1 > length(gamma_x) || is.na(gamma_x[timepoint + 1])) {
    gamma = gamma_x[timepoint]
  } else {
    gamma = gamma_x[timepoint + 1]
  }
  
  v_t_val = v_t[floor(timepoint+1)]
  
  dS = -beta * S * I/N - v_t_val * parms$v_e * (1-0.0026)
  dE = beta * S * I/N - delta * E
  dI = delta * E - gamma * I
  dR = gamma * I
  dV = v_t_val * parms$v_e * (1-0.0026)
  
  results = c(dS ,dE ,  dI , dR , dV)
  
  return(list(results))
}


# solution of ODE. The state values indices are same for all parameters.
ode_solution_for_vec <- function(state_values,
                                 parms_vec,
                                 time_vec){
  #initial values for ODE
  
  results = lsoda(y = state_values, 
                  times = time_vec, 
                  func = seirv_for_vec, 
                  parms = parms_vec)
  return(results)
}
###########################################################################################
################ LOG PRIOR, LOG POSTERIOR FUNCTIONS FOR tau_1 ##############################
###########################################################################################
tau1_log_postd <- function (tau1_x,
                            shape_p,
                            rate_p,
                            obs1, 
                            deltal, 
                            gammal,
                            El_sol
) 
{
  log_prior = dgamma(tau1_x, shape = shape_p, rate = rate_p, log = TRUE)
  
  mu_obs1 = exp(deltal) * (El_sol)
  
  # Calculate the log likelihood
  log_lik = sum(dnbinom(obs1, size = 1/(tau1_x),
                        mu = mu_obs1, log = TRUE)
  )
  output = log_prior + log_lik
  return(output)
}

###########################################################################################
################ LOG PRIOR, LOG POSTERIOR FUNCTIONS FOR tau_2 ##############################
###########################################################################################
tau2_log_postd <- function (tau2_x,
                            shape_p,
                            rate_p,
                            obs2, 
                            deltal, 
                            gammal,
                            Il_sol
) 
{
  log_prior = dgamma(tau2_x, shape = shape_p, rate = rate_p, log = TRUE)
  
  mu_obs2 = exp(gammal) * (Il_sol)
  
  # Calculate the log likelihood
  log_lik = sum(dnbinom(obs2, size = 1/(tau2_x),
                        mu = mu_obs2, log = TRUE)
  )
  
  output = log_prior + log_lik
  
  return(output)
}
###########################################################################################
###########################################################################################
######################### MAP, mean, burnin for posterior ##################################
###########################################################################################
###########################################################################################
# 
# # This function retrieves the samples for 't', excluding burn-in iterations.
# get_burnin_samples <- function(t, posterior, burnin) {
#   # 'sapply' to loop over each set of samples, excluding the burn-in iterations.
#   # For each set of samples, extract the sample at time 't'.
#   samples <- sapply(posterior[-(1:burnin)], 
#                     function(samples) {
#                       samples[t]
#                     })
#   return(samples)
# }
# 
# ##################################################################################################
# # This function calculates the Maximum A Posteriori (MAP) estimate for a given set of samples.
# calculate_map <- function(samples) {
#   # table of the samples using 'table'.
#   # This groups the samples into bins and counts how many samples fall into each bin.
#   vals <- table(samples)
#   
#   # Identify which bin has the most samples (i.e., find the mode of the samples).
#   # The 'names' function is used to retrieve the names of the bins, which in this case are the sample values.
#   # The 'which.max' function identifies which bin has the most samples.
#   map_estimate <- as.numeric(names(vals)[which.max(vals)]) 
#   return(map_estimate)
# }
# 
# ##################################################################################################
# This function calculates the mean for a given set of samples.
calculate_mean <- function(samples) {
  # Use the 'mean' function to calculate the average of the samples.
  mean_estimate <- mean(samples)
  
  return(mean_estimate)
}


###############################################################################
###############################################################################
###############################################################################
plot_hist <- function(title, 
                      filename_prefix,
                      samp_data) {
  
  # Convert the list of samples into a numeric vector
  samp_vector <- as.numeric(unlist(samp_data))
  
  # Create the histogram with automatic binwidth selection
  p <- ggplot(data.frame(value = samp_vector), aes(x=value)) +
    geom_histogram(fill="steelblue", color="black", alpha=0.7) +  # No binwidth specified
    theme_minimal() +
    theme(text = element_text(size = TITLE_TEXT_SIZE + 10, face="bold"),
          axis.title.x = element_text(size = X_TEXT_SIZE + 10, face="bold"),
          axis.title.y = element_text(size = X_TEXT_SIZE + 10, face="bold"),
          plot.title = element_text(size = TITLE_TEXT_SIZE + 10, face="bold", hjust = 0.5)) +
    labs(title = title, 
         x = "Value", 
         y = "Frequency") +
    scale_x_continuous(limits = c(min(samp_vector), max(samp_vector)))
  
  # Save the plot
  ggsave(paste0(filename_prefix, "_histogram.jpg"), plot = p)
  
  return(p)
}

calc_post_mean_credible <- function(post_samp, timex, lower = 0.025, upper = 0.975) {
  posterior_mean <- sapply(timex, 
                           function(i) {
                             mean(sapply(post_samp, function(sample) sample[i]))
                           }
  )
  
  # Calculate credible intervals (2.5% and 97.5% by default)
  lower_ci <- sapply(timex, 
                     function(i) {
                       quantile(sapply(post_samp, function(sample) sample[i]), probs = lower)
                     }
  )
  
  upper_ci <- sapply(timex, 
                     function(i) {
                       quantile(sapply(post_samp, function(sample) sample[i]), probs = upper)
                     }
  )
  
  return(list(mean = posterior_mean, lower = lower_ci, upper = upper_ci))
}


plot_eda_mean_credible <- function(filename, 
                                   timey, 
                                   y_obs, 
                                   y_mean, 
                                   y_lower, 
                                   y_upper, 
                                   xdate,
                                   name_obs,
                                   name_post,
                                   plot_title) {
  
  # Create a data frame from the input vectors
  df <- data.frame(timey, y_obs, y_mean, y_lower, y_upper, xdate)
  
  # Initialize the plotly figure
  fig <- plot_ly(data = df,
                 width = 1400,
                 height = 700) %>%
    add_trace(x = ~xdate, y = ~y_obs, type = "scatter", mode = "lines", name = name_obs, line = list(color = "red")) %>%
    add_trace(x = ~xdate, y = ~y_mean, type = "scatter", mode = "lines", name = name_post, line = list(color = "blue")) %>%
    # Add the credible interval as a filled region between y_lower and y_upper
    add_ribbons(x = ~xdate, ymin = ~y_lower, ymax = ~y_upper, name = "Credible Interval", fillcolor = 'rgba(100, 100, 255, 0.2)', line = list(color = 'transparent')) 
  
  # Add layout with title and axis configurations
  fig <- layout(fig, 
                title = list(
                  text = plot_title,
                  font = list(family = as.character(x_tf), size = x_ts),  # Customize font and size for the title
                  x = 0.7,  # Center the title
                  y = 0.95,  # Adjust y to move the title down
                  xanchor = 'center',
                  yanchor = 'top'                ),
                xaxis = list(
                  title = paste0("Date starting from  ", start_date, "  to  ", end_date),
                  tickformat = "%b",
                  dtick = "M1",
                  titlefont = list(family = as.character(x_tf), size = x_ts),
                  tickfont = list(family = as.character(x_vf), size = x_vs)),
                
                yaxis = list(
                  title = as.character("Value"),
                  titlefont = list(family = as.character(y_tf), size = y_ts),
                  tickfont = list(family = as.character(y_vf), size = y_vs)),
                
                legend = list(
                  orientation = 'v', 
                  x = 0, 
                  y = 1, 
                  xanchor = 'left', 
                  yanchor = 'top',
                  font = list(size = lg_s, family = as.character(lg_f)))
  )
  
  # Save the plot as an HTML file
  save_path <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(fig, file = save_path, selfcontained = F)
  
  # Convert HTML to PNG
  webshot(save_path, file = filename)
}