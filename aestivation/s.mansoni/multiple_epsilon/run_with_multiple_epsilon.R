# This  script is analyzing the prevalence and MPB in human in a seasonal environment 
# Mainly we want to see the effect of seasonality on the prevalence and 
# mean burden of parasite and the optimal temperature. 

############### we run the model with multiply many mean temperature and three different seasonality 


#We have to runt he model for large epsilon values and we need to try with different seasonality values. 

rm(list = ls())


# upload the necessary packages 
library(foreach)
library(doParallel)

# run the parameter script
source("par_set_mansoni.R")

# We improve the previous model with estivation function which is as follow


# We start with setting the parameters used in the aestivation function
# set the temperature of half saturation for snails estivation function at high temperatures 
est_temp <- 23  

# set the rate of snails move to the aestivation stages
est_rate <- 0.085

# Set the steepness of aestivation function for hot temperature 
steepness <- 0.9


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

# This function calculate the average number of mated pair of worm given MPB 
phiSm <- function(M) {
  # Define the function to be integrated
  k <- exp(0.61521*log(M) - 4.844146)
  
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
temperature <- seq(min_tem, max_tem, 0.1);

# number of year to run
year <- 40


step_size <- 0.1
# Set time span for model running 
run_time <- seq(from = 0, to = 365*year, by = step_size)

#set a sample space for temperature 
sample_parameters <- seq(from = 0, to = 365*year+1, by = 1)

#seasonality three value
sesonality <- c(0, .1, .25)

# create an empty matrix to record the outcome of simulations 
out_come_season <- matrix(nrow = length(sesonality), ncol = length(temperature))



# create empty matrix for the parameters 
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


# start run the model for each seasonality value 

for(epsilon in 1:length(sesonality)){
  
  # run the model for each temperature and corresponding seasonal value 
  
  for(j in 1:length(temperature)){
    
    #set seasonal temperature
    seasonal_temperature <-data.frame(temp = temperature[j] * (1 + sesonality[epsilon] * sin(2 * pi * sample_parameters/365)))
    
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
  
  
  # run the model for each temperature (all parameter values )
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
             
             dS <- (nu_s - (S + E + I) * nu) * (S + r * E) -  lambda * (1-exp(-(beta_s * delta_e * nu_e * h * phiSm(P_m)*P_m)/(mu_m * (S + E + I)))) * S - mu * S - in_est * S + out_est * SE
             dSE <- in_est * S - mu/re_est * SE - out_est * SE
             dE <- lambda * (1-exp(-(beta_s * delta_e * nu_e *h * phiSm(P_m)*P_m)/(mu_m * (S + E + I)))) * S - (sigma_s + mu_i) * E - in_est*E + out_est * EE
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
  out_come_for_each_temperature <- foreach(k = 1:length(temperature), .combine = c) %dopar% {
    
    result <- solution_for_temp(k)
    
    return(result)
    
  }
  
  #Stop cluster
  stopCluster(cluster)
  
  out_come_season[epsilon, ] <- out_come_for_each_temperature
} 



# Plot mean parasite burden for three different seasonality
pdf(file = "mpb_with_estivation.pdf", width = 8, height = 8)
par(mar=c(6, 6, 6, 1), xpd = TRUE)
plot(temperature, out_come_season[1,], col = "blue",  pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     main = "", ylim = c(0, max(out_come_season[1,])+20))  
lines(temperature, out_come_season[2,], col = "orange", pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2)
lines(temperature, out_come_season[3,], col = "red", pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2)

mtext(text = expression(paste("Mean annual temperature (",degree,"C)")), side = 1, line = 4, cex = 2)
mtext(text = "Mean parasite burden", side = 2, line = 4, cex = 2)
mtext(text = expression(paste(italic(S.), " ", italic(mansoni))* " with adaptive behaviors"), side = 3, line = 4, cex = 2)

dev.off()


#Define function to convert the MPB into prevalence 
wormPrevalenceSm <- function(M) {
  k <- exp(0.61521*log(M) - 4.844146)
  p  <- 1 - (1+M/k)^(-k)    # fraction of humans with at least 1 parasites 
  return(p)
}

#convert mean parasite burden to prevalence 
out_come_season[1, ] <- wormPrevalenceSm(out_come_season[1,])
out_come_season[2, ] <- wormPrevalenceSm(out_come_season[2,])
out_come_season[3, ] <- wormPrevalenceSm(out_come_season[3,])



# Plot mean parasite burden for three different seasonality
pdf(file = "prevalence_with_estivation.pdf", width = 8, height = 8)
par(mar=c(6, 6, 6, 1), xpd = TRUE)
plot(temperature, out_come_season[1,], col = "black",  pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     main = "", ylim = c(0, 1))  
lines(temperature, out_come_season[2,], col = "blue", pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2)
lines(temperature, out_come_season[3,], col = "red", pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2)

text(x = 29, y = 0.9, labels = '(23.1, 0.69)', col = "black", cex = 1.5) 
arrows(x0 = temperature[which.max(out_come_season[1,])], y0 = max(out_come_season[1,]), 
       x1 = 27, y1 = 0.9, col = "black", length = 0.1, lty = 3, lwd = 3)


text(x = 29, y = 0.85, labels = '(22.8, 0.63)', col = "blue", cex = 1.5) 
arrows(x0 = temperature[which.max(out_come_season[2,])], y0 = max(out_come_season[2,]), 
       x1 = 27, y1 = 0.85, col = "blue", length = 0.1, lty = 3, lwd = 3)

text(x = 29, y = 0.8, labels = '(20.9, 0.40)', col = "red", cex = 1.5) 
arrows(x0 = temperature[which.max(out_come_season[3,])], y0 = max(out_come_season[3,]), 
       x1 = 27, y1 = 0.8, col = "red", length = 0.1, lty = 3, lwd = 3)

mtext(text = expression(paste("Mean annual temperature (",degree,"C)")), side = 1, line = 4, cex = 2)
mtext(text = expression(paste(italic(Schistosome), " Prevalence")), side = 2, line = 4, cex = 2)
mtext(text = expression(paste(italic(S.), " ", italic(mansoni))* " with adaptive behaviors"), side = 3, line = 4, cex = 2)

# Adding arrow passing over three points for outputs_1
points_to_connect <- c(which.max(out_come_season[1,]), which.max(out_come_season[2,]), which.max(out_come_season[3,]))
arrows(x0 = temperature[points_to_connect[1]], y0 = out_come_season[1,points_to_connect[1]], 
       x1 = temperature[points_to_connect[2]], y1 = out_come_season[2,points_to_connect[2]], 
       col = "black", length = 0.1, lty = 1, lwd = 3)
arrows(x0 = temperature[points_to_connect[2]], y0 = out_come_season[2,points_to_connect[2]], 
       x1 = temperature[points_to_connect[3]], y1 = out_come_season[3,points_to_connect[3]], 
       col = "black", length = 0.1, lty = 1, lwd = 3)

dev.off()






