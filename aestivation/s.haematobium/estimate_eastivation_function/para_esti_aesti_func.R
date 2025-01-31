#This code estimate aestivation function parameter

#clear environment 
rm(list = ls())

library(foreach)
library(doParallel)


# We start with setting our eastivation function
# set the temperature of half saturation for snails estivation function at high temperatures 
aest_temp_seq <- seq(from = 22, to = 27, by = .5)  

# set the rate of snails move to the estivation stages
aest_rate_seq <- seq(from = 0.01, to = 0.12, by = 0.05)

# Set the steepness of estivation function for hot temperature 
steepness_seq <-  seq(from = 0.1, to = 1.6, by = 0.5)

# Generate all combinations
combinations <- expand.grid(aest_temp_seq, aest_rate_seq, steepness_seq)

combinations <- combinations[11:20,]

# Set minimum temperature 
min_tem <- 14;

# Set maximum temperature 
max_tem <- 31;


# Set minimum epsilon 
min_eps <- 0.02;

# Set maximum epsilon 
max_eps <- 0.26;


temperature_step <- 0.5
epsilon_step <- 0.02


# Define temperature and epsilon ranges
temperatures <- seq(min_tem, max_tem, by = temperature_step)
epsilons <- seq(min_eps, max_eps, by = epsilon_step)

# upload gntd data 
prev_percent <- read.csv(file = 'gntd_vars_all.csv')

# Filter S. haematobium data directly
sch_haematobium <- subset(prev_percent, parasite_s == "S. haematobium" & percent_pos != 0 & !is.na(bio01))

# Convert Bioclimate variables into seasonality
gntd_data <- data.frame(
  MAT = sch_haematobium$bio01,
  Epsilon = (sch_haematobium$bio10 - sch_haematobium$bio11) * pi / (sch_haematobium$bio01 * 4 * sqrt(2)),
  Prevalence = sch_haematobium$percent_pos
)


# Pre-allocate matrix for results
new_gntd_data <- matrix(NA, nrow = length(epsilons), ncol = length(temperatures))

# Use nested loops for temperature and epsilon bins
for (i in seq_along(temperatures)) {
  temp_min <- temperatures[i]
  temp_max <- temp_min + temperature_step
  index_mat <- gntd_data$MAT > temp_min & gntd_data$MAT <= temp_max
  
  for (j in seq_along(epsilons)) {
    eps_min <- epsilons[j]
    eps_max <- eps_min + epsilon_step
    index_eps <- gntd_data$Epsilon > eps_min & gntd_data$Epsilon <= eps_max
    
    combined_index <- index_mat & index_eps
    new_gntd_data[j, i] <- mean(gntd_data$Prevalence[combined_index], na.rm = TRUE)
  }
}

# Replace NA entries with zero
new_gntd_data[is.na(new_gntd_data)] <- 0

# run the parameter script
source("par_set_haematobium.R")

# This function calculate the average number of mated pair of worm given MPB 
phiSm <- function(M) {
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

# When using estimated clumper parameter for S. heamatobium 
wormPrevalenceSh <- function(M) {
  k <- exp(0.5186358*log(M) - 3.253653)
  p  <- 1 - (1+M/k)^(-k)    # fraction of humans with at least 1 parasites 
  return(p)
}


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


# combine two piece wise function of mu_i
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

square_error <- c()

for (k in 1:nrow(combinations)){
  
  # We start with setting our eastivation function
  # set the temperature of half saturation for snails estivation function at high temperatures 
  aest_temp <- combinations[k, 1]
  
  # set the rate of snails move to the estivation stages
  aest_rate <- combinations[k, 2]
  
  # Set the steepness of estivation function for hot temperature 
  steepness <-  combinations[k, 3]
  
# Define a piecewise function for vector input
in_aest_function <- function(x) {
  y <- aest_rate/(1+exp(-steepness*(x-aest_temp)))
  return(y)
}

# Define a piecewise function for vector input
out_aest_function <- function(x) {
  y <- aest_rate/(1+exp(steepness*(x-aest_temp)))
  return(y)
}

for(epsilon in 1:length(epsilons)){
  
  for(j in 1:length(temperatures)){
    
    #set seasonal temperature
    seasonal_temperature <-data.frame(temp = temperatures[j] * (1 + epsilons[epsilon] * sin(2 * pi * sample_parameters/365)))

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
    
  }
  
  solution_for_temp <- function(index){
    
    # Generate linearly interpolate point with temperature dependent parameter function   
    nu_s_afun <- approxfun(x = sample_parameters, y = preds_nu_s[index,])
    mu_afun <- approxfun(x = sample_parameters, y = preds_mu[index,])
    mu_m_afun <- approxfun(x = sample_parameters, y = preds_mu_m[index,])
    sigma_s_afun <- approxfun(x = sample_parameters, y = preds_sigma_s[index,])
    mu_i_afun <- approxfun(x = sample_parameters, y = preds_mu_i[index, ])       
    nu_c_afun <- approxfun(x = sample_parameters, y = preds_nu_c[index,])
    mu_c_afun <- approxfun(x = sample_parameters, y = preds_mu_c[index,])           
    delta_e_afun <- approxfun(x = sample_parameters, y = preds_delta_e[index,])
    beta_s_afun <- approxfun(x = sample_parameters, y = preds_beta_s[index,])
    beta_h_afun <- approxfun(x = sample_parameters, y = preds_beta_h[index,])
    
    in_aestivation_afun <- approxfun(x = sample_parameters, y = preds_in_aestivation[index,])
    out_aestivation_afun <- approxfun(x = sample_parameters, y = preds_out_aestivation[index,])
    
    
    
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
             
             dS <- (nu_s - (S + E + I) * nu) * (S + r * E) -  lambda * (1-exp(-(beta_s * delta_e * nu_e * h * phiSm(P_m)*P_m)/(mu_m * (S + E + I)))) * S - mu * S - in_aest * S + out_aest * SE
             dSE <- in_aest * S - mu/re_est * SE - out_aest * SE
             dE <- lambda * (1-exp(-(beta_s * delta_e * nu_e *h * phiSm(P_m)*P_m)/(mu_m * (S + E + I)))) * S - (sigma_s + mu_i) * E - in_aest*E + out_aest * EE
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
    out_put <- mean(model_outputs[(length(model_outputs[, 1])-5*366/step_size):length(model_outputs[, 1]),9])
    
    return(out_put)
  }
  
  
  #Setup backend to use many processors
  totalCores = detectCores()
  
  #Leave one core to avoid overload your computer
  cluster <- makeCluster(totalCores[1]-1) 
  registerDoParallel(cluster) 
  
  # Parallel loop using %dopar%
  out_come_for_each_temperature <- foreach(k = 1:length(temperatures), .combine = c) %dopar% {
    
  result <- solution_for_temp(k)
    
    return(result)
    
  }
  
  #Stop cluster
  stopCluster(cluster)
  
  out_come_season_mpb[epsilon, ] <- out_come_for_each_temperature
} 
# Replace entries in out_come_season_mpb with 0 if they are greater than 100 or less than 2
out_come_season_mpb[out_come_season_mpb < .4] <- 0

prevalence_percentage <- wormPrevalenceSh(out_come_season_mpb)*100
square_error[k] <- sum((prevalence_percentage-new_gntd_data)^2)
}


save.image(file = "my_environment_20.RData")

