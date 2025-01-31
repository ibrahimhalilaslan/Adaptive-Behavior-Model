# This  script is analyzing prevalence and MPB with seasonality 
# Mainly we want to see the effect of seasonality on the prevalence and 
# mean burden of parasite and the optimal temperature. 

############### we run the model with multiply many mean temperature and three different seasonality 


rm(list = ls())

# Model parameters are set in this script 

# run the parameter script
source("par_set_haematobium.R")
 
library(foreach)
library(doParallel)


# Load the R environment from the .RData file
load("era5_ts_1990_2023_daily_wide.RData")


# Define the range of rows to read
start_row <- 1
end_row <- 1000

# Read specific rows using 
temperature <- df_daily_wide[start_row:end_row,]




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


# set the temperature of half saturation for snails estivation function at high temperatures 
est_temp <- 23

# set the rate of snails move to the estivation stages
est_rate <- 0.11

# Set the steepness of estivation function for hot temperature 
steepness <-  0.6


# Define a piecewise function for vector input
in_est_function <- function(x) {
  y <- est_rate/(1+exp(-steepness*(x-est_temp)))
  return(y)
}

# Define a piecewise function for vector input
out_est_function <- function(x) {
  y <- est_rate/(1+exp(steepness*(x-est_temp)))
  return(y)
}



# The reduction of mortality due to estivation
re_est <- 20
re_est_i <- 2


#################################################


n <- nrow(temperature)
m <- ncol(temperature)-3
mean_annual_temp <- c()

#set a sample space for temperature 
sample_parameters <- seq(from = 1, to = m, by = 1)

step_size <- 0.1
# Set time span for model running 
run_time <- seq(from = 1, to = m-1, by = step_size)

  #create empty matrix for store the parameters input for the paralel computing 
  preds_mu_i <- matrix(nrow = n, ncol = m)                        
  preds_mu <- matrix(nrow = n, ncol = m)                        
  preds_nu_s <- matrix(nrow = n, ncol = m)
  preds_mu_m <- matrix(nrow = n, ncol = m)
  preds_sigma_s <- matrix(nrow = n, ncol = m) 
  preds_nu_c <- matrix(nrow = n, ncol = m)
  preds_mu_c <- matrix(nrow = n, ncol = m)
  preds_delta_e <- matrix(nrow = n, ncol = m)
  preds_beta_s <- matrix(nrow = n, ncol = m)
  preds_beta_h <- matrix(nrow = n, ncol = m)
  preds_in_estivation <- matrix(nrow = n, ncol = m)
  preds_out_estivation <- matrix(nrow = n, ncol = m)
  
# Start save the parameter value in the matrix for paralel computing   
for(j in 1:n){
    
    #set seasonal temperature
    seasonal_temperature <- data.frame(temp =  as.numeric(temperature[j,4:(m+3)]))
    mean_annual_temp[j] <- mean(as.numeric(temperature[j,4:(m+3)]))
    
    # calculate the parameters values for each corresponding temperature 
    preds_nu_s[j,] <- ifelse(fn_nu_s(seasonal_temperature)$.fitted >= 0, 
                             fn_nu_s(seasonal_temperature)$.fitted,0)
    preds_mu[j,] <- ifelse(seasonal_temperature$temp <= 35.5, 
                           ifelse(fn_mu_1(seasonal_temperature)$.fitted >= 0,fn_mu_1(seasonal_temperature)$.fitted, 0),
                           ifelse(fn_mu_2(seasonal_temperature)$.fitted >= 0,fn_mu_2(seasonal_temperature)$.fitted, 0))
    preds_mu_m[j,] <- ifelse(fn_mu_m(seasonal_temperature)$.fitted >= 0, 
                             fn_mu_m(seasonal_temperature)$.fitted,0)
    preds_sigma_s[j,] <- ifelse(fn_sigma_s(seasonal_temperature)$.fitted >= 0, 
                                fn_sigma_s(seasonal_temperature)$.fitted,0)
    preds_mu_i[j,] <- ifelse(seasonal_temperature$temp <= 34, 
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
    preds_in_estivation[j,] <- in_est_function(seasonal_temperature$temp)
    preds_out_estivation[j,] <- out_est_function(seasonal_temperature$temp) 
    
    
  }
  
  
  
  # This part solve the dynamical system 
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
    
    in_estivation_afun <- approxfun(x = sample_parameters, y = preds_in_estivation[index,])
    out_estivation_afun <- approxfun(x = sample_parameters, y = preds_out_estivation[index,])
    
    
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
             
             in_est <- in_estivation_afun(t)
             out_est <- out_estivation_afun(t)
             
             dS <- (nu_s - (S + E + I) * nu) * (S + r * E) -  lambda * (1-exp(-(beta_s * delta_e * nu_e * h * phiSh(P_m)*P_m)/(mu_m * (S + E + I)))) * S - mu * S - in_est * S + out_est * SE
             dSE <- in_est * S - mu/re_est * SE - out_est * SE
             dE <- lambda * (1-exp(-(beta_s * delta_e * nu_e *h * phiSh(P_m)*P_m)/(mu_m * (S + E + I)))) * S - (sigma_s + mu_i) * E - in_est*E + out_est * EE
             dEE <- in_est * E - mu_i/re_est_i * EE - out_est * EE
             dI <- sigma_s * E - mu_i * I - in_est*I + out_est * IE
             dIE <- in_est * I - mu_i/re_est_i * IE - out_est * IE
             dP <- beta_h * (nu_c/mu_c) * I - sigma_p * P
             dP_m <- sigma_p * P - (mu_h + mu_p) * P_m
             
             return(list(c(S = dS, SE = dSE, E = dE, EE = dEE, I = dI, IE = dIE, P = dP, P_m = dP_m)))
           })
    }
    
    # Specified the parameter value. These might be changed later because they are prevalence dependent. 
    parms0 <- c(nu, lambda)  
    
    # Set the initial conditions 
    y0 <- c(S = (60434 + 33232)/2, SE = 0,  E = 1285, EE = 0, I = 2340, IE = 0, P = 1, P_m = 45)
    
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
  out_come_for_each_temperature <- foreach(k = 1:n, .combine = c) %dopar% {
    
    result <- solution_for_temp(k)
    
    return(result)
    
  }
  
  #Stop cluster
  stopCluster(cluster)
  
  out_come_for_each_temperature[out_come_for_each_temperature < 0.4] <- 0
  
  # When using estimated clumper parameter for S. heamatobium 
  wormPrevalenceSh <- function(M) {
    k <- exp(0.5186358*log(M) - 3.253653)
    p  <- 1 - (1+M/k)^(-k)    # fraction of humans with at least 1 parasites 
    return(p)
  }
  
  
  df_results <- data.frame(Longitude = temperature[, 3],
                     Latitude =temperature[, 2], MAT = mean_annual_temp, MPB = out_come_for_each_temperature, Prevalence = wormPrevalenceSh(out_come_for_each_temperature))

  # Save as CSV with specified column names
  write.csv(df_results,"time_series_1000.csv", row.names = FALSE)
  

  
  
  # Parallel loop using %dopar%
  out_come_for_each_temperature <- foreach(k = 1:n, .combine = c) %dopar% {
    
    result <- solution_for_temp(k)
    
    return(result)
    
  }
  
  #Stop cluster
  stopCluster(cluster)
  
  out_come_for_each_temperature[out_come_for_each_temperature < 0.8] <- 0
  
  wormPrevalenceSm <- function(M) {
    k <- exp(0.61521*log(M) - 4.844146)
    p  <- 1 - (1+M/k)^(-k)    # fraction of humans with at least 1 parasites 
    return(p)
  }
  
  
  df_results <- data.frame(Longitude = temperature[, 3],
                           Latitude = temperature[, 2], MAT = mean_annual_temp, MPB = out_come_for_each_temperature, Prevalence = wormPrevalenceSm(out_come_for_each_temperature))
  
  # Save as CSV with specified column nam/es
  write.csv(df_results,"Sm_time_series_2000.csv", row.names = FALSE)
  
  
  
  
  
  df_1 <- read.csv(file = 'time_series_1000.csv')
  df_2 <- read.csv(file = 'time_series_1500.csv')
  df_3 <- read.csv(file = 'time_series_2000.csv')
  df_4 <- read.csv(file = 'time_series_2500.csv')
  df_5 <- read.csv(file = 'time_series_3000.csv')
  df_6 <- read.csv(file = 'time_series_4000.csv')
  df_7 <- read.csv(file = 'time_series_5000.csv')
  df_8 <- read.csv(file = 'time_series_6000.csv')
  df_9 <- read.csv(file = 'time_series_7000.csv')
  df_10 <- read.csv(file = 'time_series_8000.csv')
  df_11 <- read.csv(file = 'time_series_9000.csv')
  df_12 <- read.csv(file = 'time_series_10000.csv')
  df_13 <- read.csv(file = 'time_series_11000.csv')
  df_14 <- read.csv(file = 'time_series_12000.csv')
  df_15 <- read.csv(file = 'time_series_13000.csv')
  df_16 <- read.csv(file = 'time_series_14000.csv')
  df_17 <- read.csv(file = 'time_series_14630.csv')
  
  combine_df <- bind_rows(df_1, df_2, df_3, df_4, df_5, df_6, df_7, 
                     df_8, df_9, df_10, df_11, df_12, df_13, df_14,
                     df_15, df_16, df_17)
  
  # Save as CSV with specified column names
  write.csv(combine_df,"Sh_simulation_africa_time_series_mpb.csv", row.names = FALSE)
  
  
 
