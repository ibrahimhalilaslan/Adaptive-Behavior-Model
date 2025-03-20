#This code estimate aestivation function parameter with MLE

#clear environment 
rm(list = ls())


# Set minimum temperature 
min_tem <- 14;

# Set maximum temperature 
max_tem <- 31;


# Set minimum epsilon 
min_eps <- 0.02;

# Set maximum epsilon 
max_eps <- 0.26;

#set the temperature steps 
temperature_step <- 0.5

#set the epsilon steps 
epsilon_step <- 0.02


# Define temperature and epsilon ranges
temperatures <- seq(min_tem, max_tem, by = temperature_step)
epsilons <- rev(seq(min_eps, max_eps, by = epsilon_step))

# upload gntd data 
prev_percent <- read.csv(file = 'gntd_data_cases_number.csv')

# Filter S. haematobium data directly
sch_haematobium <- subset(prev_percent, parasite_s == "S. haematobium" & !is.na(bio01))

# Convert Bioclimate variables into seasonality
gntd_data <- data.frame(
  MAT = sch_haematobium$bio01,
  Eps = (sch_haematobium$bio10 - sch_haematobium$bio11) * pi / (sch_haematobium$bio01 * 4 * sqrt(2)),
  Prev = sch_haematobium$total_number_pos,
  Cases = sch_haematobium$total_number_exam
)


# Pre-allocate matrix for results
cases_gntd_data <- matrix(NA, nrow = length(epsilons), ncol = length(temperatures))
total_test_gntd_data <- matrix(NA, nrow = length(epsilons), ncol = length(temperatures))

# Use nested loops for temperature and epsilon bins
# This generate GNTD data for each cell grid 
for (i in seq_along(temperatures)) {
  temp_min <- temperatures[i] - temperature_step/2
  temp_max <- temperatures[i] + temperature_step/2
  index_mat <- gntd_data$MAT > temp_min & gntd_data$MAT <= temp_max
  
  for (j in seq_along(epsilons)) {
    eps_min <- epsilons[j] - epsilon_step/2
    eps_max <- epsilons[j] + epsilon_step/2
    index_eps <- gntd_data$Eps > eps_min & gntd_data$Eps <= eps_max
    
    combined_index <- index_mat & index_eps
    if (sum(combined_index) == 0) next  # Skip this iteration if no matching data
    
    cases_gntd_data[j, i] <- sum(gntd_data$Prev[combined_index], na.rm = TRUE)
    total_test_gntd_data[j, i] <- sum(gntd_data$Cases[combined_index], na.rm = TRUE)
  }
}


# run the parameter script
source("par_set_haematobium.R")

# This function calculate the average number of mated pair of worm given MPB 
# if we want to calculate the average number of mated pair of worm we will be using 
phiSh <- function(M) {
  # Define the function to be integrated
  k <- exp(0.5186358*log(M) - 3.253653)
  
  integrand <- function(x) {
    (1 - cos(x)) / ((1 + (M / (M + k)) * cos(x))^(1 + k))
  }
  
  # Perform the integration
  integral_result <- integrate(integrand, lower = 0, upper = 2 * pi)
  
  # Calculate the prevalence
  p <- 1 - (((k / (M + k))^(1 + k))/ (2 * pi)) * integral_result$value
  
  return(p)
}


# convert MPB to prevalence
wormPrevalenceSh <- function(M) {
  k <- exp(0.5186358*log(M) - 3.253653)
  p  <- 1 - (1+M/k)^(-k)    # fraction of humans with at least 1 parasites 
  return(p)
}

# Start to run the model with given aestivation function parameters
objective_function <- function(params) {
  
# number of year to run
year <- 40

step_size <- 0.1
# Set time span for model running 
run_time <- seq(from = 0, to = 365*year, by = step_size)

# The reduction of mortality due to estivation
re_est <- 20
re_est_i <- 2


#set a sample space for temperature 
sample_parameters <- seq(from = 0, to = 365*year+1, by = 1)


# Generate an empty matrix to record outcome of the runs  
out_come_season_mpb <- matrix(nrow = length(epsilons), ncol = length(temperatures))
out_put <- c()

# Generate an empty matrix to record parameter to be used in the model simulation 
preds_mu_i <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))                        
preds_mu <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))                        
preds_nu_s <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))
preds_mu_m <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))
preds_sigma_s <- matrix(nrow = length(temperatures), ncol = length(sample_parameters)) 
preds_nu_c <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))
preds_mu_c <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))
preds_delta_e <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))
preds_beta_s <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))
preds_beta_h <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))
preds_in_aestivation <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))
preds_out_aestivation <- matrix(nrow = length(temperatures), ncol = length(sample_parameters))


  
# We start with setting our aestivation function parameters
# set the temperature of half saturation for snails estivation function at high temperatures 
  aest_temp <- params[1]
  
  # set the rate of snails move to the estivation stages
  aest_rate <- params[2]
  
  # Set the steepness of estivation function for hot temperature 
  steepness <-  params[3]
  
  # Define aestivation function for in aestivation
  in_aest_function <- function(x) {
  y <- aest_rate/(1+exp(-steepness*(x-aest_temp)))
  return(y)
}

# Define aestivation function for out aestivation
out_aest_function <- function(x) {
  y <- aest_rate/(1+exp(steepness*(x-aest_temp)))
  return(y)
}

# start solving ode
for(epsilon in 1:length(epsilons)){
  
  for(j in 1:length(temperatures)){
    
    #set seasonal temperature
    seasonal_temperature <-data.frame(temp = temperatures[j] * (1 + epsilons[epsilon] * sin(2 * pi * sample_parameters/365)))

    # calculate the parameter value for seasonal temperature
    preds_nu_s[j,] <- ifelse(fn_nu_s(seasonal_temperature)$.fitted >= 0, 
                             fn_nu_s(seasonal_temperature)$.fitted,0)
    preds_mu[j,] <- ifelse(seasonal_temperature$temp <= 33, 
                           ifelse(fn_mu_1(seasonal_temperature)$.fitted >= 0,fn_mu_1(seasonal_temperature)$.fitted, 0),
                           ifelse(fn_mu_2(seasonal_temperature)$.fitted >= 0,fn_mu_2(seasonal_temperature)$.fitted, 0))
    preds_mu_m[j,] <- ifelse(fn_mu_m(seasonal_temperature)$.fitted >= 0, 
                             fn_mu_m(seasonal_temperature)$.fitted,0)
    preds_sigma_s[j,] <- ifelse(fn_sigma_s(seasonal_temperature)$.fitted >= 0, 
                                fn_sigma_s(seasonal_temperature)$.fitted,0)
    preds_mu_i[j,] <- ifelse(seasonal_temperature$temp <= 33.5, 
                             ifelse(fn_mu_i_1(seasonal_temperature)$.fitted >= 0,fn_mu_i_1(seasonal_temperature)$.fitted, 0),
                             ifelse(fn_mu_i_2(seasonal_temperature)$.fitted >= 0,fn_mu_i_2(seasonal_temperature)$.fitted, 0))
    preds_nu_c[j,] <- ifelse(fn_nu_c(seasonal_temperature)$.fitted >= 0, 
                             fn_nu_c(seasonal_temperature)$.fitted,0)
    preds_mu_c[j,] <- ifelse(fn_mu_c(seasonal_temperature)$.fitted >= 0, 
                             fn_mu_c(seasonal_temperature)$.fitted,0)
    preds_delta_e[j,] <- ifelse(fn_delta_e(seasonal_temperature)$.fitted >= 0, 
                                fn_delta_e(seasonal_temperature)$.fitted,0)
    preds_beta_s[j,] <- ifelse(fn_beta_s(seasonal_temperature)$.fitted >= 0,
                               fn_beta_s(seasonal_temperature)$.fitted, 0)
    preds_beta_h[j,] <- ifelse(fn_beta_h(seasonal_temperature)$.fitted >= 0,
                               fn_beta_h(seasonal_temperature)$.fitted, 0)
    preds_in_aestivation[j,] <- in_aest_function(seasonal_temperature$temp)
    preds_out_aestivation[j,] <- out_aest_function(seasonal_temperature$temp) 
    
    # Generate linearly interpolate point with temperature dependent parameter function   
    nu_s_afun <- approxfun(x = sample_parameters, y = preds_nu_s[j,])
    mu_afun <- approxfun(x = sample_parameters, y = preds_mu[j,])
    mu_m_afun <- approxfun(x = sample_parameters, y = preds_mu_m[j,])
    sigma_s_afun <- approxfun(x = sample_parameters, y = preds_sigma_s[j,])
    mu_i_afun <- approxfun(x = sample_parameters, y = preds_mu_i[j, ])       
    nu_c_afun <- approxfun(x = sample_parameters, y = preds_nu_c[j,])
    mu_c_afun <- approxfun(x = sample_parameters, y = preds_mu_c[j,])           
    delta_e_afun <- approxfun(x = sample_parameters, y = preds_delta_e[j,])
    beta_s_afun <- approxfun(x = sample_parameters, y = preds_beta_s[j,])
    beta_h_afun <- approxfun(x = sample_parameters, y = preds_beta_h[j,])
    
    in_aestivation_afun <- approxfun(x = sample_parameters, y = preds_in_aestivation[j,])
    out_aestivation_afun <- approxfun(x = sample_parameters, y = preds_out_aestivation[j,])
    
    
    
    #Call the library 
    library(deSolve)
    
    #the diff equation solver 
    thermal_sensitive_model <- function(t, y, parms){
      with(as.list(c(y, parms)),
           {
             nu_s <- nu_s_afun(t)
             mu_m <- mu_m_afun(t)
             mu <- mu_afun(t)
             sigma_s <- sigma_s_afun(t)
             mu_i <- mu_i_afun(t)        
             nu_c <- nu_c_afun(t)             
             mu_c <- mu_c_afun(t)
             delta_e <- delta_e_afun(t)
             beta_s <- beta_s_afun(t)
             beta_h <- beta_h_afun(t)
             
             in_aest <- in_aestivation_afun(t)
             out_aest <- out_aestivation_afun(t)
             
             dS <- (nu_s - (S + E + I) * nu) * (S + r * E) -  lambda * (1-exp(-(beta_s * delta_e * nu_e * h * phiSh(P_m)*P_m)/(mu_m * (S + E + I)))) * S - mu * S - in_aest * S + out_aest * SE
             dSE <- in_aest * S - mu/re_est * SE - out_aest * SE
             dE <- lambda * (1-exp(-(beta_s * delta_e * nu_e *h * phiSh(P_m)*P_m)/(mu_m * (S + E + I)))) * S - (sigma_s + mu_i) * E - in_aest*E + out_aest * EE
             dEE <- in_aest * E - mu_i/re_est_i * EE - out_aest * EE
             dI <- sigma_s * E - mu_i * I - in_aest*I + out_aest * IE
             dIE <- in_aest * I - mu_i/re_est_i * IE - out_aest * IE
             dP <- beta_h * (nu_c/mu_c) * I - sigma_p * P
             dP_m <- sigma_p * P - (mu_h + mu_p) * P_m
             
             return(list(c(S = dS, SE = dSE, E = dE, EE = dEE, I = dI, IE = dIE, P = dP, P_m = dP_m)))
           })
    }
    
    parms0 <- c(nu, lambda)  
    
    # Set the initial conditions 
    y0 <- c(S = (60434 + 33232)/2, SE = 0,  E = 1285, EE = 0, I = 2340, IE = 0, P = 3, P_m = 130)
    
    ## solve the system 
    model_outputs <- ode(y = y0, times = run_time, func = thermal_sensitive_model, parms = parms0)
    
    # record the results 
    out_put[j] <- mean(model_outputs[(length(model_outputs[, 1])-5*366/step_size):length(model_outputs[, 1]),9])
} 
# record the outputs 
out_come_season_mpb[epsilon, ] <- out_put
}  


# Replace entries in out_come_season_mpb with 0 if they are greater than 100 or less than 2
out_come_season_mpb[out_come_season_mpb < 0.4] <- 0

#Convert MPB to the prevalence
prevalence <- wormPrevalenceSh(out_come_season_mpb)

# calculate the Maximumlikelihood value 
likelihood_sum_aestivation <- -2*sum(dbinom(
  x = cases_gntd_data,       # Observed cases
  size = total_test_gntd_data, # Number of trials
  prob = prevalence,        # Probability from model
  log = TRUE
), na.rm = TRUE)  # Handle NA values
return(likelihood_sum_aestivation)

}

 

# Run optimization from different starts
opt_parameter <- optim(
  par = c(23, 0.065, 0.7), 
  fn = objective_function, 
  method = "L-BFGS-B", 
  lower = c(22.5, -Inf, -Inf),  # Lower bound for the first parameter only
  upper = c(30, Inf, Inf)     # Upper bound for the first parameter only
)

