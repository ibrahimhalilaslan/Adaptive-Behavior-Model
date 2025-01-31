# This  script simulates the mean parasite burden (MPB) with seasonal temperature 
# We want to see the effect of seasonality on MPB and the optimal temperature. 

############### we run the model with multiply many mean temperature and three different seasonality 
rm(list = ls())

library(foreach)
library(doParallel)

# run the parameter script
source("par_set_haematobium.R")

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
min_tem <- 12;

# Set maximum temperature 
max_tem <- 32;

# Set a sequence of temperature  
temperature <- seq(min_tem, max_tem, 0.1);

# Set a matrix to record the outputs
out_puts <- c()

#create a function to run the model with seasonality 
seaonality_effect <- function(epsilon){
   
  # For loop to run the model for each temperature
  for (j in 1:length(temperature)){
  
    # Set the number of year to run the model  
    year <- 40
    
    #Set time span for model running
    step_size = 0.05;
    run_time <- seq(from = 0, to = 365*year, by = step_size)
    
    #set a sample space of parameters for each day 
    sample_parameters <- seq(from = 0, to = 365*year+1, by = 0.1)
    
    
    #set seasonal temperature
    seasonal_temperature <-data.frame(temp = temperature[j] * (1 + epsilon * sin(2 * pi * sample_parameters/365)))

    
    # Generate linearly interpolate point with temperature dependent parameter function   
    nu_s_afun <- approxfun(x = sample_parameters, y = ifelse(fn_nu_s(seasonal_temperature)$.fitted >= 0, 
                                                             fn_nu_s(seasonal_temperature)$.fitted,0))
    mu_m_afun <- approxfun(x = sample_parameters, y = ifelse(fn_mu_m(seasonal_temperature)$.fitted >= 0, 
                                                             fn_mu_m(seasonal_temperature)$.fitted,0))
    mu_afun <- approxfun(x = sample_parameters, y = ifelse(seasonal_temperature$temp <= 33, 
                                                           ifelse(fn_mu_1(seasonal_temperature)$.fitted >= 0,fn_mu_1(seasonal_temperature)$.fitted, 0),
                                                           ifelse(fn_mu_2(seasonal_temperature)$.fitted >= 0,fn_mu_2(seasonal_temperature)$.fitted, 0)))
    sigma_s_afun <- approxfun(x = sample_parameters, y = ifelse(fn_sigma_s(seasonal_temperature)$.fitted >= 0, 
                                                                fn_sigma_s(seasonal_temperature)$.fitted,0))
    mu_i_afun <- approxfun(x = sample_parameters, y = ifelse(seasonal_temperature$temp <= 33.5, 
                                                             ifelse(fn_mu_i_1(seasonal_temperature)$.fitted >= 0,fn_mu_i_1(seasonal_temperature)$.fitted, 0),
                                                             ifelse(fn_mu_i_2(seasonal_temperature)$.fitted >= 0,fn_mu_i_2(seasonal_temperature)$.fitted, 0)))       
    nu_c_afun <- approxfun(x = sample_parameters, y = ifelse(fn_nu_c(seasonal_temperature)$.fitted >= 0, 
                                                             fn_nu_c(seasonal_temperature)$.fitted,0))
    mu_c_afun <- approxfun(x = sample_parameters, y = ifelse(fn_mu_c(seasonal_temperature)$.fitted >= 0, 
                                                             fn_mu_c(seasonal_temperature)$.fitted,0))           
    delta_e_afun <- approxfun(x = sample_parameters, y = ifelse(fn_delta_e(seasonal_temperature)$.fitted >= 0, 
                                                                fn_delta_e(seasonal_temperature)$.fitted,0))
    beta_s_afun <- approxfun(x = sample_parameters, y = ifelse(fn_beta_s(seasonal_temperature)$.fitted >= 0,
                                                               fn_beta_s(seasonal_temperature)$.fitted, 0))
    beta_h_afun <- approxfun(x = sample_parameters, y = ifelse(fn_beta_h(seasonal_temperature)$.fitted >= 0,
                                                               fn_beta_h(seasonal_temperature)$.fitted, 0))
    
    
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
             
             
             dS <- (nu_s - (S + E + I) * nu) * (S + r * E) -  lambda * (1-exp(-(beta_s * delta_e * nu_e * h *  phiSh(P_m)*P_m)/(mu_m * (S + E + I)))) * S - mu * S 
             dE <- lambda * (1-exp(-(beta_s * delta_e * nu_e *h *  phiSh(P_m)*P_m)/(mu_m * (S + E + I)))) * S - (sigma_s + mu_i) * E 
             dI <- sigma_s * E - mu_i * I 
             dP <- beta_h * (nu_c/mu_c) * I - sigma_p * P
             dP_m <- sigma_p * P - (mu_h + mu_p) * P_m
             
             return(list(c(S = dS, E = dE, I = dI, P = dP, P_m = dP_m)))
           })
    }
    
    # Specified the parameter value.
    parms0 <- c(nu, lambda)  
    
    # Set the initial conditions 
    y0 <- c(S = (60434 + 33232)/2,  E = 1285, I = 2340, P = 1, P_m = 45)
    
    # solve the system 
    model_outputs <- ode(y = y0, times = run_time, func = thermal_sensitive_model, parms = parms0)
    
    # record the results 
    out_puts[j] <- mean(model_outputs[(length(model_outputs[, 1])-5*366/step_size):length(model_outputs[, 1]), 6])
    
    }
  return(out_puts)
}

#Run the model with different seasonality 
outputs_1_mpb <- seaonality_effect(0)
outputs_2_mpb <- seaonality_effect(0.1)
outputs_3_mpb <- seaonality_effect(0.25)


# Plot mean parasite burden for three different seasonality
pdf(file = "seasonal_mpb.pdf", width = 8, height = 8)
par(mar=c(6, 6, 6, 1), xpd = TRUE)
plot(temperature, outputs_1_mpb, col = "blue",  pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     main = "", ylim = c(0, max(outputs_1_mpb)))  
lines(temperature, outputs_2_mpb, col = "orange", pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2)
lines(temperature, outputs_3_mpb, col = "red", pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2)

#legend(x = 10, y = 111,  legend =  "No seas. Opt. Temp. =", bty = "n", cex = 1.5, text.col = 1)
#legend(x = 19, y = 111, legend =  round(temperature[which.max(outputs_1[5,])],4), cex = 1.5, bty = "n", text.col = 1)
#legend(x = 10, y = 107,  legend =  "Mild. seas. Opt. Temp. =", bty = "n", cex = 1.5, text.col = 2)
#legend(x = 19, y = 107, legend =  round(temperature[which.max(outputs_2[5,])],4), cex = 1.5, bty = "n", text.col = 2)
#legend(x = 10, y = 103,  legend =  "Sev. seas. Opt. Temp. =", bty = "n", cex = 1.5, text.col = 3)
#legend(x = 19, y = 103, legend =  round(temperature[which.max(outputs_3[5,])],4), cex = 1.5, bty = "n", text.col = 3)


#mtext(text = expression(paste("Mean annual temperature (",degree,"C)")), side = 1, line = 4, cex = 2)
#mtext(text = "Mean parasite burden", side = 2, line = 4, cex = 2)
mtext(text = expression(paste(italic(S.), " ", italic(haematobium))* " without adaptive behaviors"), side = 3, line = 4, cex = 2)

dev.off()

# When using estimated clumper parameter for S. heamatobium 
wormPrevalenceSh <- function(M) {
  k <- exp(0.5186358*log(M) - 3.253653)
  p  <- 1 - (1+M/k)^(-k)    # fraction of humans with at least 1 parasites 
  return(p)
}



outputs_1 <- wormPrevalenceSh(outputs_1_mpb)
outputs_2 <- wormPrevalenceSh(outputs_2_mpb)
outputs_3 <- wormPrevalenceSh(outputs_3_mpb)


# Plot mean parasite burden for three different seasonality
pdf(file = "seasonal_prevalence.pdf", width = 8, height = 8)
par(mar=c(6, 6, 6, 1), xpd = TRUE)
plot(temperature, outputs_1, col = "black",  pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     main = "", ylim = c(0, 1))  
lines(temperature, outputs_2, col = "blue", pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2)
lines(temperature, outputs_3, col = "red", pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2)

#text(x = 29, y = 0.95, labels = expression((T[opt] * "," ~ Prev[max])), col = "black", cex = 1.5)
text(x = 29, y = 0.9, labels = '(25.4, 0.76)', col = "black", cex = 1.5) 
arrows(x0 = temperature[which.max(outputs_1)], y0 = max(outputs_1), 
       x1 = 27, y1 = 0.9, col = "black", length = 0.1, lty = 3, lwd = 3)

text(x = 29, y = 0.85, labels = '(24.9, 0.73)', col = "blue", cex = 1.5) 
arrows(x0 = temperature[which.max(outputs_2)], y0 = max(outputs_2), 
       x1 = 27, y1 = 0.85, col = "blue", length = 0.1, lty = 3, lwd = 3)

text(x = 29, y = 0.8, labels = '(22.8, 0.62)', col = "red", cex = 1.5) 
arrows(x0 = temperature[which.max(outputs_3)], y0 = max(outputs_3), 
       x1 = 27, y1 = 0.8, col = "red", length = 0.1, lty = 3, lwd = 3)

mtext(text = expression(paste(italic(S.), " ", italic(haematobium))* " without adaptive behaviors"), side = 3, line = 4, cex = 2)


# Adding arrow passing over three points for outputs_1
points_to_connect <- c(which.max(outputs_1), which.max(outputs_2), which.max(outputs_3))
arrows(x0 = temperature[points_to_connect[1]], y0 = outputs_1[points_to_connect[1]], 
       x1 = temperature[points_to_connect[2]], y1 = outputs_2[points_to_connect[2]], 
       col = "black", length = 0.1, lty = 1, lwd = 3)
arrows(x0 = temperature[points_to_connect[2]], y0 = outputs_2[points_to_connect[2]], 
       x1 = temperature[points_to_connect[3]], y1 = outputs_3[points_to_connect[3]], 
       col = "black", length = 0.1, lty = 1, lwd = 3)


dev.off()
