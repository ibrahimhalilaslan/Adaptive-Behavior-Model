# This  script is analyzing prevalence and MPB with seasonality and adaptive behavior
# Mainly we want to see the effect of seasonality on the prevalence and 
# mean burden of parasite and the optimal temperature. 

############### To create a heat map, we run the model with multiply many mean temperature and three different seasonality 

rm(list = ls())

library(foreach)
library(doParallel)

# run the parameter script
source("par_set_haematobium.R")

# We start with setting the parameter sets for the eastivation function
# set the temperature of half saturation for snails estivation function at high temperatures 
est_temp <- 23

# set the rate of snails move to the estivation stages
est_rate <- 0.065

# Set the steepness of aestivation function for hot temperature 
steepness <-  0.7

# Define aestivation function for in aestivation
in_est_function <- function(x) {
  y <- est_rate/(1+exp(-steepness*(x-est_temp)))
  return(y)
}

# Define aestivation function for out aestivation
out_est_function <- function(x) {
  y <- est_rate/(1+exp(steepness*(x-est_temp)))
  return(y)
}

# The reduction of mortality due to estivation
re_est <- 20
re_est_i <- 2


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

# Set minimum temperature 
min_tem <- 14;

# Set maximum temperature 
max_tem <- 32;

# Set a sequence of temperature  
temperature <- seq(min_tem, max_tem, .2);

# number of year to run
year <- 40

# Set time span for model running 
step_size = 0.1
run_time <- seq(from = 0, to = 365*year, by = step_size)

#set a sample space for temperature 
sample_parameters <- seq(from = 0, to = 365*year+1, by = 1)

# set the min value of seasonality 
min_season <- 0

# set the max value of seasonality 
max_season <- 0.27

#seasonality value
seasonality <- seq(min_season, max_season, 0.01);

# Generate an empty matrix to record outcome of the runs  
out_come_season <- matrix(nrow = length(seasonality), ncol = length(temperature))

# Generate an empty matrix to record all the parameter need in the runs  
preds_mu_i <- matrix(nrow = length(temperature), ncol = length(sample_parameters))                        
preds_mu <- matrix(nrow = length(temperature), ncol = length(sample_parameters))                        
preds_nu_s <- matrix(nrow = length(temperature), ncol = length(sample_parameters))
preds_mu_m <- matrix(nrow = length(temperature), ncol = length(sample_parameters))
preds_sigma_s <- matrix(nrow = length(temperature), ncol = length(sample_parameters)) 
preds_nu_c <- matrix(nrow = length(temperature), ncol = length(sample_parameters))
preds_mu_c <- matrix(nrow = length(temperature), ncol = length(sample_parameters))
preds_delta_e <- matrix(nrow = length(temperature), ncol = length(sample_parameters))
preds_beta_s <- matrix(nrow = length(temperature), ncol = length(sample_parameters))
preds_beta_h <- matrix(nrow = length(temperature), ncol = length(sample_parameters))
preds_in_estivation <- matrix(nrow = length(temperature), ncol = length(sample_parameters))
preds_out_estivation <- matrix(nrow = length(temperature), ncol = length(sample_parameters))


# Start looping the simulaiton for each epsilon and temperature value
for(epsilon in 1:length(seasonality)){
  
  for(j in 1:length(temperature)){
    
    #set seasonal temperature
    seasonal_temperature <-data.frame(temp = temperature[j] * (1 + seasonality[epsilon] * sin(2 * pi * sample_parameters/365)))

    # Set the parameter for seasonal varying temperature 
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
    preds_in_estivation[j,] <- in_est_function(seasonal_temperature$temp)
    preds_out_estivation[j,] <- out_est_function(seasonal_temperature$temp) 
    
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
  out_come_for_each_temperature <- foreach(k = 1:length(temperature), .combine = c) %dopar% {
    
    result <- solution_for_temp(k)
    
    return(result)
    
  }
  
  #Stop cluster
  stopCluster(cluster)
  
  out_come_season[epsilon, ] <- out_come_for_each_temperature
} 



# Repeat the column n times
epsilon_matrix <- matrix(rep(seasonality, length(temperature)), nrow = length(seasonality), byrow = FALSE)
temp_matrix <- matrix(rep(temperature, length(seasonality)), nrow = length(seasonality), byrow = TRUE)



# save the mean parasite burden as data frame for plot 
df <- data.frame(x = as.vector(temp_matrix), y = as.vector(epsilon_matrix), MPB = as.vector(out_come_season))

# upload th library 
library(ggplot2)

# Plot the results and print 
pdf(file="heat_map.pdf", width = 8, height = 5)
par(mfrow=c(2,1))
densit_plot <- ggplot(df, aes(x = x, y = y, colour = MPB,), x ) + geom_point(size = 4, shape = 15)
densit_plot + scale_color_gradient(low="white", high="red") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(x = "Mean annual temperature", y = expression(paste("Seasonality  ", (epsilon))))+ 
  theme(axis.text = element_text(size = 20,color = "black")) + theme(axis.title = element_text(size = 20,color = "black")) 

dev.off()

#Save the environment 
save.image(file = "my_environment.RData")







