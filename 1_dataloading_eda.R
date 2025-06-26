##############################################################################################
############## data loading ##################################################################
##############################################################################################
rm(list = ls())

# install.packages("crayon")
#library loading
{
  # pacman::p_load(deSolve,dplyr,ggplot2,reshape2,e1071,crayon,tidyverse )
  library(deSolve)
  URL_cases <- "https://raw.githubusercontent.com/MoH-Malaysia/covid19-public/main/epidemic/cases_malaysia.csv"
  URL_vax <- "https://raw.githubusercontent.com/MoH-Malaysia/covid19-public/main/vaccination/vax_malaysia.csv"
  cases_data <- read.csv(URL_cases)
  vax_data <- read.csv(URL_vax)
  
  cases_data$date <- as.Date(cases_data$date) # convert to date object
  vax_data$date <- as.Date(vax_data$date) # convert to date object
  merged_data <- merge(vax_data, cases_data, by="date", all=TRUE)
  merged_data[is.na(merged_data)] <- 0
  
  start_date <- as.Date("2021-02-24")
  end_date <- as.Date("2021-10-31")
  end_date_2 <- as.Date("2021-11-15")
    
  merged_data <- merged_data[merged_data$date >= start_date & merged_data$date <= end_date, ]
  cases_new_a <- cases_data$cases_new[cases_data$date >= start_date & cases_data$date <= end_date_2 ]
}
#merged_data <- readxl::read_xlsx("EDA_EXCEL.xlsx")

##############################################################################################
############## Vaccination Efficacy ##################################################################
##############################################################################################
#vax_eff = 1
{
  cases_vax1_total = merged_data$daily_partial
  cases_vax2_total = merged_data$daily_full
  pfizer1 = (merged_data$pfizer1)
  pfizer2 = (merged_data$pfizer2)
  sinovac1 = (merged_data$sinovac1)
  sinovac2 = (merged_data$sinovac2)
  astra1 = (merged_data$astra1)
  astra2 = (merged_data$astra2)
  others1 = merged_data$sinopharm1 + merged_data$cansino
  others2 = merged_data$sinopharm2 + merged_data$cansino3
  
  
  percentage_vax1_brand <- c(
    pfizer1 = 100 * sum(pfizer1)/sum(cases_vax1_total),
    sinovac1 = 100 * sum(sinovac1)/sum(cases_vax1_total),
    astra1 = 100 * sum(astra1)/sum(cases_vax1_total),
    others1 = 100 * sum(others1)/sum(cases_vax1_total)
  )
  
  percentage_vax2_brand <- c(
    pfizer2 = 100 * sum(pfizer2)/sum(cases_vax1_total),
    sinovac2 = 100 * sum(sinovac2)/sum(cases_vax1_total),
    astra2 = 100 * sum(astra2)/sum(cases_vax1_total),
    others2 = 100 * sum(others2)/sum(cases_vax1_total)
  )
  
  
  # Create a data frame
  df <- data.frame(
    brand = rep(c("pfizer", "sinovac", "astra", "others"), 2),
    percentage = c(percentage_vax1_brand, percentage_vax2_brand),
    vax = rep(c("vax1", "vax2"), each = 4)
  )
  
  # # Use ggplot2 to create the plot
  # p <- ggplot(df, aes(x = brand, y = percentage, fill = vax)) +
  #   geom_bar(stat = "identity", position = "dodge") +
  #   geom_text(aes(label = round(percentage, 2)), vjust = -0.3, position = position_dodge(.9)) +
  #   labs(x = "Brand", y = "Percentage", title = "Vaccine brand percentages", fill = "Vaccine") +
  #   theme_minimal()
  # ggsave(paste0("vax_percent", ".png"), plot = p)
}
####################################################################################################################
#from Literature
{
  ve_pfizer1 = 0.80 ;
  ve_pfizer2 = 0.92  ;
  ve_sinovac1 = 0.35 ;
  ve_sinovac2 = 0.74 ;
  ve_astra1 = 0.36 ;
  ve_astra2 = 0.83 ;
  
  ####################################################################################################################
  
  weighted_vax1 = mean((ve_pfizer1 * pfizer1 + ve_sinovac1 * sinovac1 + 
                          ve_astra1 * astra1)/cases_vax1_total)
  
  weighted_vax2 = mean((ve_pfizer2 * pfizer2 + ve_sinovac2 * sinovac2 + 
                          ve_astra2 * astra2)/cases_vax2_total)
  vax_eff = mean(c(weighted_vax1,weighted_vax2))
  print(vax_eff)
}

# GOOOOOOOOOOOOOOOO to line 31 !!!!!!!!!!!!!!!!!!


##############################################################################################
############## Fixed Parameters ##################################################################
##############################################################################################
vax_eff = 0.75
{
  
N = 32600000
delta_period = 5
delta = 1/delta_period
cases_active = merged_data$cases_active
cases_new = merged_data$cases_new
cases_recovered = merged_data$cases_recovered
cases_vax = merged_data$daily_full
}
#cases_vax = merged_data$daily_vaccinated

##############################################################################################
############## New Exposed ##################################################################
##############################################################################################
{
new_exposed <- rep(NA, length(cases_new))
exposed_over_est <- rep(NA, length(cases_new))
exposed_under_est <- rep(NA, length(cases_new))
index <- length(cases_new) 

if (delta_period %% 2 ==0){
  
  #over estimation
  t_0_over = delta_period/2;
  t_n_over = 3 * delta_period / 2;
  period_over <- t_n_over - t_0_over

  for (i in 1:index) {
    ind = i + t_0_over
    exposed_over_est[i] <- delta * sum(cases_new_a[(ind):(ind+period_over)])
  }
  

  #under estimation
  t_0_under = (delta_period/2) + 1
  t_n_under = (3 * delta_period/2) - 1
  under_est_period <- t_n_under - t_0_under
  for (i in 1:index) {
    ind = i + t_0_under 
    exposed_under_est[i] <- delta * sum(cases_new_a[(ind):(ind+under_est_period)])
    new_exposed[i] <- (exposed_over_est[i] + exposed_under_est[i])/2
    
  }
  

  
} else {
  
  #over estimation
  t_0 = (delta_period+1)/2;
  t_n = (3 * delta_period-1) / 2;
  est_period <- t_n - t_0

  for (i in 1:index) {
    ind = i + t_0 - 1
    new_exposed[i] <- delta * sum(cases_new_a[(ind):(ind + est_period)])
  }
  
}

}
##############################################################################################
############## Exposed Compartment ##################################################################
##############################################################################################
{
  
  E_t = c (cases_new[1]/delta, cases_new[2:length(cases_new)]/delta)

  S <- numeric(0)
  S[1] <- N

  for (i in 2:length(merged_data[,1])) {
  S[i] = S[i-1] - new_exposed[i] - (cases_vax[i] * vax_eff)
  }

}
##############################################################################################
############## Beta and Gamma ##################################################################
##############################################################################################
{
  beta_eda <- numeric(0)
  
  n = length(merged_data[,1])
  for (i in 1:n) {
    beta_eda[i] = new_exposed[i] / (cases_active[i] * S[i]/N)
  }
  
  gamma_eda <- numeric(0)
  
  gamma_eda = cases_recovered/cases_active;

}

##############################################################################################
############## Timepoints and data extraction ##################################################################
##############################################################################################
{
  T = length(merged_data$date)
  timepoint <- seq(0, T-1, by = 1)
  indices  <- timepoint + 1
  
  cases_active = merged_data$cases_active[indices]
  cases_new = merged_data$cases_new[indices]
  cases_recovered = merged_data$cases_recovered[indices]
  cases_vax = merged_data$daily_full[indices]
  
  dates <- seq(start_date, end_date, "days")
  xdate <- dates
  options(scipen = 999)  
}


timepoint <- seq(0, length(xdate)-1, by=1)
indices  <- timepoint+1

popsize = 32600000

cases_new = c(0, cases_new)
cases_active = c(0, cases_active)
cases_recovered = c(0, cases_recovered)
cases_vaccinated = c(0,cases_vax)
betal_eda = c(log(beta_eda[1]),log(beta_eda))
deltal_eda = c(log(delta), rep(log(delta), T))
gammal_eda = c(log(gamma_eda[1]), log(gamma_eda))
El_eda = c(E_t[1]-1000 , E_t)
Il_eda = c(cases_active[1]-1000, cases_active)
R_eda = c(0, cumsum(cases_recovered))
S_eda = c(popsize - cases_active[1]-1000 - E_t[1]-1000 , S)

##############################################################################################
##############################################################################################
##############################################################################################

##############################################################################################
############## function loading ##################################################################
##############################################################################################
# 
# seirv_model = function(timepoint, state_values, parms)
#   {
#     S = state_values[1]
#     E = state_values[2]       
#     I = state_values[3]       
#     R = state_values[4]       
#     V = state_values[5]       
#     N = S + E + I + R + V
#     with( 
#       as.list(parms),      
#       {
#         beta = beta_x[floor(timepoint + 1)]
#         gamma = gamma_x[floor(timepoint + 1)]
#         v_t = v_t[floor(timepoint + 1)]
#         
#         dS = -beta * S * I - v_t * v_e
#         dE =  beta * S * I - delta * E
#         dI =  delta * E - gamma * I
#         dR =  gamma * I
#         dV =  v_t * v_e
#         results = c(dS ,dE ,  dI , dR , dV)
#         return(list(results))
#       }
#     )
#   }
#   
#   
# ############## parameters and initial values loading ##################################################################
# ##############################################################################################
# 
# {
#   state_values = c(S = S[1]/N, 
#                    E = E_t[1]/N,
#                    I = cases_active[1]/N, 
#                    R = cases_recovered[1]/N, 
#                    V = 0/N)
#   beta_x = beta_eda
#   gamma_x = gamma_eda
#   v_t = cases_vax/N
# }
# 
# ############## plots ##################################################################
# #######################################################################################
# {
#   parms = c(delta = delta, 
#             beta = beta_x, 
#             gamma = gamma_x, 
#             v_e = vax_eff, 
#             v_t = cases_vax)
#   output = lsoda(state_values, timepoint, seirv_model, parms)
#   SModel = output[,2]
#   EModel = output[,3]
#   IModel = output[,4]
#   RModel = output[,5]
#   VModel = output[,6]
#   
#   
#   cex_value_lab <- 1.5
#   cex_value_axis <- 1.5
#   cex_value_text <- 1.5
#   par(mar = c(5.1, 6.1, 4.1, 2.1))
# }
#   plot(dates, beta_eda, type = "p", col="red", lwd = 3, 
#        ylab = expression(paste("Observed ", beta, " Values")), xlab = "Time period",
#        cex.lab = cex_value_lab, cex.axis = cex_value_axis)
#   text(dates[70], 0.17, 
#        expression(paste("Observed ", beta, " Values")),
#        adj = c(0, 0.5), cex = cex_value_text
#   )
#   title(expression(paste("Observed ", beta, " Values")))
#   
#   par(mar = c(5.1, 6.1, 4.1, 2.1))
#   plot(dates, gamma_eda, type = "p", col="red", lwd = 3, 
#        ylab = expression(paste("Observed ", gamma, " Values")), xlab = "Time period",
#        cex.lab = cex_value_lab, cex.axis = cex_value_axis)
#   text(dates[110], 0.13, 
#        expression(paste("Observed ", gamma, " Values")),
#        adj = c(0, 0.5), cex = cex_value_text
#   )
#   title(expression(paste("Observed ", gamma, " Values")))
#   
#   
#   {  
#   png(file = "EDA_Cases_Active.png")
#   p <- plot(xdate, IModel*N, type = "l", col="blue", lwd = 3, 
#             ylab ="Active Cases", xlab = "Time period", ylim = c(0,400000),
#             cex.lab = cex_value_lab, cex.axis = cex_value_axis)
#   lines(xdate,cases_active , type = "l", col = "red", lw = 3)
#   text(xdate[70], 200000, 
#        paste("Model", "(Blue)","\n" ,"\n", "Observed (Red)"),
#        adj = c(0, 0.5), cex = cex_value_text
#   )
#   title(paste0("Active cases - EDA with vax_eff", vax_eff))
#   dev.off()  # Close the PNG device
#   
#   png(file = "EDA_Cases_New.png")
#   p <- plot(xdate, delta*EModel*N, type = "l", col="blue", lwd = 3, 
#             ylab ="New Cases", xlab = "Time period", ylim = c(0,30000),
#             cex.lab = cex_value_lab, cex.axis = cex_value_axis)
#   lines(xdate,cases_new , type = "l", col = "red", lw = 3)
#   text(xdate[70], 18000, 
#        paste("Model", "(Blue)","\n" ,"\n", "Observed (Red)"),
#        adj = c(0, 0.5), cex = cex_value_text
#   )
#   title(paste0("New Cases EDA with v_e = 0.90"))
#   dev.off()  # Close the PNG device
#   
#   png(file = "EDA_Cum_Cases_Recovered.png")
#   p <-  plot(xdate, RModel*N, type = "l", col="blue", lwd = 3, 
#              ylab ="Cumulative Recovered Cases", xlab = "Time period", ylim = c(0,max(RModel*N)),
#              cex.lab = cex_value_lab, cex.axis = cex_value_axis)
#   lines(xdate,cumsum(cases_recovered) , type = "l", col = "red", lw = 3)
#   text(xdate[70], 2000000, 
#        paste("Model", "(Blue)","\n" ,"\n", "Observed (Red)"),
#        adj = c(0, 0.5), cex = cex_value_text
#   )
#   title(paste0("Cumulative Recovered cases - EDA with v_e = 0.90"))
#   dev.off()  # Close the PNG device
#   
#   png(file = "EDA_Cases_Recovered.png")
#   p <-  plot(xdate[-1], diff(RModel*N), type = "l", col="blue", lwd = 3, 
#              ylab ="Recovered Cases", xlab = "Time period", ylim = c(0,max(diff(RModel*N))),
#              cex.lab = cex_value_lab, cex.axis = cex_value_axis)
#   lines(xdate,(cases_recovered) , type = "l", col = "red", lw = 3)
#   text(xdate[70], 20000, 
#        paste("Model", "(Blue)","\n" ,"\n", "Observed (Red)"),
#        adj = c(0, 0.5), cex = cex_value_text
#   )
#   title(paste0("Recovered cases - EDA, with v_e = 0.90"))
#   dev.off()  # Close the PNG device
# }
