#This code estimates the eastivation function with saveral field studies 
# Regarding abundance of snails with temperature. 

# Clear environmet 
rm(list = ls())

library(dplyr)


# set the temperature of half saturation for snails estivation function at high temperatures 
est_temp <- 27


# set the rate of snails move to the estivation stages
est_rate <- 0.01

# Set the steepness of estivation function for hot temperature 
steepness <-  0.4


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
### we set the model to run 


plot(15:35, in_est_function(15:35))

#Call the library 
library(deSolve)


#the diff equation solver 
thermal_sensitive_model <- function(t, y, parms){
  with(as.list(c(y, parms)),
       {
         
         in_est <- in_estivation_afun(t)
         out_est <- out_estivation_afun(t)
         
         dS <- -in_est * S + out_est * SE
         dSE <- in_est * S - out_est * SE
         
         
         return(list(c(S = dS, SE = dSE)))
       })
}



######################## Kenya data ###################
# Read the CSV file for Kenya temperature data 
kenya_data <- read.csv(file = 'kenya_temp_data.csv', nrow = 120)

#Process the Kenya data to group by
kenya_data$system.time_start <- rep(1:12,  10)

# Group the data by 'year' and calculate the mean of temperature 
kenya_grouped_temperature <- kenya_data %>%
  group_by(system.time_start) %>%
  summarise_all(mean, na.rm = TRUE)

#Save the abundance of snails from Kariuki at al. 2004
kenya_snail_aboundance <- c(105,  14,  96,  74,  99, 109, 137,  78,  76,   8,   9, 28,  19,  10,  20, 118, 180, 160, 115,  96, 57, 134,
                            145, 131,   9,  12, 13,  79, 174, 249,   9, 135, 167, 128, 153, 78,   9,   8,  10,  56,   9)

## Kenya month for snails abundance 
kenya_months <- c(3:12, 1:12, 1:12, 1:7)

# Create data frame to group abundance of snails 
kenya_data <- data.frame(months = kenya_months, aboundance = kenya_snail_aboundance)

# Group the data by 'year' and calculate the mean of abundance 
kenya_grouped_aboundance <-kenya_data %>%
  group_by(months) %>%
  summarise_all(mean, na.rm = TRUE)

# We form one year data for snails shows the abundance of snails with varying temperature 
# Now we need to find best epsilon value for the fit of one year temperature 


#set a sample space for temperature 
sample_parameters <- seq(from = 0, to = length(kenya_grouped_aboundance$aboundance)*30+1, by = 1)

temperature <- 26;

epsilon <- 0.07;

m = 70

#set seasonal temperature
seasonal_temperature <-data.frame(temp = temperature * (1 + epsilon * sin(2 * pi * (m+sample_parameters)/360)))


# We plot the the temperature variation over time to see if the epsilon is a good choice 
plot(sample_parameters, seasonal_temperature$temp, type = "l", ylim = c(23, 30))
points(seq(from = 1, to = 360, by = 30), kenya_grouped_temperature$mean)


plot(seasonal_temperature$temp, in_est_function(seasonal_temperature$temp))

#After we find the epsilon and mean annual temperature we run the model with initial conditins 

# Introduce eastivation function
in_estivation_afun <- approxfun(x = sample_parameters, y = in_est_function(seasonal_temperature$temp))
out_estivation_afun <- approxfun(x = sample_parameters, y = out_est_function(seasonal_temperature$temp))



# Set time span for model running 
run_time <- seq(from = 0, to = length(kenya_grouped_aboundance$aboundance)*30, by = 1)
  
## solve the system 
model_outputs <- ode(y = c(S = kenya_grouped_aboundance$aboundance[1], SE = max(kenya_grouped_aboundance$aboundance)- kenya_grouped_aboundance$aboundance[1]), times = run_time, func = thermal_sensitive_model, parms = 0)
  
# Plotting data and simulation
pdf(file = "kenya_bulinus_one_years_abundance.pdf", width = 8, height = 8)
par(mar=c(6, 6, 5, 2), xpd = TRUE)
plot(model_outputs[,1], model_outputs[,2], col = 1,  pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2, main = "", ylim = c(0, 200))  
points(seq(from = 15, to = length(kenya_grouped_aboundance$aboundance)*30-15, by = 30), kenya_grouped_aboundance$aboundance, col = 2, pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, cex.axis=2)

# Adding legends
legend("topleft", legend = c("Simulation", "Data"), col = c(1, 2), pch = 21, lwd = 3)
mtext(text = "Day(s)", side = 1, line = 4, cex = 2)
mtext(text = "Abundance", side = 2, line = 4, cex = 2)
title(main =  expression(paste(italic(Bulinus)," spp. in Msambweni, Kenya for average of three years")), line = 2.5, cex.main = 1.5)

dev.off()

#set a sample space for temperature 
sample_parameters <- seq(from = 0, to = length(kenya_snail_aboundance)*30+1, by = 1)

temperature <- 26;

epsilon <- 0.07;

m = 70

#set seasonal temperature
seasonal_temperature <-data.frame(temp = temperature * (1 + epsilon * sin(2 * pi * (m+sample_parameters)/360)))


# Introduce eastivation function
in_estivation_afun <- approxfun(x = sample_parameters, y = in_est_function(seasonal_temperature$temp))
out_estivation_afun <- approxfun(x = sample_parameters, y = out_est_function(seasonal_temperature$temp))


# Set time span for model running 
run_time <- seq(from = 0, to = length(kenya_snail_aboundance)*30, by = 1)

## solve the system 
model_outputs <- ode(y = c(S = kenya_snail_aboundance[1], SE = max(kenya_snail_aboundance)- kenya_snail_aboundance[1]), times = run_time, func = thermal_sensitive_model, parms = 0)

# Plotting data and simulation
pdf(file = "kenya_bulinus_three_years_abundance.pdf", width = 8, height = 8)
par(mar=c(6, 6, 5, 2), xpd = TRUE)
plot(model_outputs[,1], model_outputs[,2], col = 1,  pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2, main = "", ylim = c(0, 200))  
points(seq(from = 15, to = length(kenya_snail_aboundance)*30-15, by = 30), kenya_snail_aboundance, col = 2, pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, cex.axis=2)

# Adding legends
legend("topleft", legend = c("Simulation", "Data"), col = c(1, 2), pch = 21, lwd = 3)
mtext(text = "Day(s)", side = 1, line = 4, cex = 2)
mtext(text = "Abundance", side = 2, line = 4, cex = 2)
title(main =  expression(paste(italic(Bulinus)," spp. in Msambweni, Kenya for three years")), line = 2.5, cex.main = 1.5)

dev.off()


################################## Tanzania data #######################

# Read the CSV file
tanzania_data <- read.csv(file = 'tanzania_temp_data.csv', nrow = 120)

# Generate class data to group temperature and take mean 
tanzania_data$system.time_start <- rep(1:12,  10)


# Group the data by 'year' and calculate the mean for each variable
tanzania_grouped_data <-tanzania_data %>%
  group_by(system.time_start) %>%
  summarise_all(mean, na.rm = TRUE)

#  These data are from Starkloff at al. 2023

tanzania_sanil_aboundance_non_desica <- c(19, 12, 13, 17, 18, 29, 29, 24, 28, 27, 31, 28)
tanzania_sanil_aboundance_desica <- c(3,  1,  0,  1,  5, 24, 46, 36, 45, 41, 43, 27)

#set a sample space for temperature 
sample_parameters <- seq(from = 0, to = length(kenya_grouped_aboundance$aboundance)*30+1, by = 1)

temperature <- 23.2;

epsilon <- 0.015;

m = 250

#set seasonal temperature
seasonal_temperature <-data.frame(temp = temperature * (1 + epsilon * sin(2 * pi * (m+sample_parameters)/360)))


# We plot the the temperature variation over time to see if the epsilon is a good choice 
plot(sample_parameters, seasonal_temperature$temp, type = "l", ylim = c(22, 25))
points(seq(from = 1, to = 360, by = 30), tanzania_grouped_data$mean)



#After we find the epsilon and mean annual temperature we run the model with initial conditins 

# Introduce eastivation function
in_estivation_afun <- approxfun(x = sample_parameters, y = in_est_function(seasonal_temperature$temp))
out_estivation_afun <- approxfun(x = sample_parameters, y = out_est_function(seasonal_temperature$temp))


####### run with non_desica area

# Set time span for model running 
run_time <- seq(from = 0, to = length(tanzania_sanil_aboundance_non_desica)*30, by = 1)

## solve the system 
model_outputs <- ode(y = c(S = tanzania_sanil_aboundance_non_desica[1], SE = max(tanzania_sanil_aboundance_non_desica)-tanzania_sanil_aboundance_non_desica[1]), times = run_time, func = thermal_sensitive_model, parms = 0)

# Plotting data and simulation
pdf(file = "tanzania_bulinus_abundance.pdf", width = 8, height = 8)
par(mar=c(6, 6, 5, 2), xpd = TRUE)
plot(model_outputs[,1], model_outputs[,2], col = 1,  pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2, main = "", ylim = c(0, 50))  
points(seq(from = 15, to = length(tanzania_sanil_aboundance_non_desica)*30-15, by = 30), tanzania_sanil_aboundance_non_desica, col = 2, pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, cex.axis=2)

# Adding legends
legend("topleft", legend = c("Simulation", "Data"), col = c(1, 2), pch = 21, lwd = 3)
mtext(text = "Day(s)", side = 1, line = 4, cex = 2)
mtext(text = "Abundance", side = 2, line = 4, cex = 2)
title(main =  expression(paste(italic(Bulinus)," spp. in Tanzanian districts of the Lake Victoria")), line = 2.5, cex.main = 1.5)

dev.off()

#run with  desication area

## solve the system 
model_outputs <- ode(y = c(S = tanzania_sanil_aboundance_desica[1], SE = max(tanzania_sanil_aboundance_desica)-tanzania_sanil_aboundance_desica[1]), times = run_time, func = thermal_sensitive_model, parms = 0)


plot(model_outputs[,2], type = "l", ylim = c(0, 100))
points(seq(from = 15, to = length(tanzania_sanil_aboundance_desica)*30-15, by = 30), tanzania_sanil_aboundance_desica)


################################ R.F. sturrock et al. ####################
# We start with seting our eastivation function

# Collection in May 1992, IN Richard Toll lies on the southern bank of the senegal River


time_vary_temp <-  c(25.478, 26.051, 26.433, 27.389, 27.962, 29.108, 29.49, 30.064,  30.828,
                29.873, 29.108, 27.197, 27.197, 25.096, 23.185, 21.274, 19.745, 17.834,
                17.261, 19.554, 19.172, 20.51, 22.229, 22.229, 22.803, 24.331, 27.771,
                30.637, 31.592, 31.975, 32.166, 31.019, 30.446, 29.682, 28.344, 27.197,
                25.86,  24.904, 24.14, 23.376, 21.274, 19.172,  18.408, 17.261, 15.924,
                18.217, 19.554, 20.51, 22.038, 22.803, 23.758, 24.904, 26.433, 27.197,
                27.962, 28.726, 28.344, 29.108, 29.49, 29.49, 29.108, 27.389, 25.096,
                23.567, 22.038, 20.318, 19.363, 18.599, 17.834, 16.306, 15.541, 16.879,
                18.025, 19.172, 20.318, 21.274, 22.611)

time_temp <- c(11.86,   25.54,   40.13,   54.72,   67.49,   89.38,  111.27,  127.69,  151.40,  169.64,  180.59,  197.00,  213.42,  224.36,  228.01,  235.31,
               240.78,  255.37,  264.50,  279.09,  299.15,  326.51,  344.76,  357.52,  379.41,  390.36,  404.95,  410.42,  421.37,  446.91,  472.44,  490.68,
               501.63,  516.22,  530.81,  536.29,  543.58,  549.06,  558.18,  569.12,  587.36,  609.25,  616.55,  629.32,  640.26,  651.21, 662.15,  674.92,
               684.04,  705.93,  727.82,  753.36,  766.12,  778.89,  789.84,  815.37,  829.97,  844.56,  855.50,  879.22,  890.16,  904.76,  923.00,  935.77,
               946.71,  952.18,  964.95,  974.07,  981.37,  995.96, 1008.73, 1019.67, 1037.92, 1052.51, 1063.45, 1078.05, 1086.25)



# make one year seasonal temperature variation 

temperature <- 24.5;

epsilon <- 0.3;

m = -20

sample_parameters <- seq(from = 0, to = 1086, by = 1)


#set seasonal temperature
seasonal_temperature <-data.frame(temp = temperature * (1 + epsilon * sin(2 * pi * (m+sample_parameters)/360)))


# We plot the the temperature variation over time to see if the epsilon is a good choice 
plot(sample_parameters, seasonal_temperature$temp, type = "l", ylim = c(13, 35))
points(time_temp, time_vary_temp)


abundance_snail_year_1 <- c(0.2769230,0.3661538,0.3938461,0.3384615,0.2646153,0.3907692,0.3323076,0.2676923,
                            0.2830769,0.2307692,0.3015384,0.3661538,0.36,0.4061538,0.3076923,0.4338461,0.4953846,
                            0.3384615,0.3538461,0.4061538,0.4984615,0.56,0.5323076,0.6030769,0.6430769)
abundance_snail_year_1 <- (10^abundance_snail_year_1-1)*85

abundance_snail_year_2 <-c(0.5544973,0.4613756,0.5714285,0.5333333,0.5714285,0.5502645,0.4656084,0.4402116,0.3767195,
                           0.2793650,0.1566137,0.2835978,0.2920634,0.3894179,0.3174603,0.4021164,0.435978,0.3894179,
                          0.3301587,0.2751322,0.3894179,0.4190476,0.4613756,0.3386243,0.4190476)
abundance_snail_year_2 <- (10^abundance_snail_year_2-1)*85

abundance_snail_year_3 <-c(0.36258,0.39136,0.46043,0.40863,0.54100,0.25323,0.34532,0.29928,0.18992,0.29928,0.23021,0.28776,
                           0.40287,0.26474,0.17266,0.18992,0.21294,0.25323,0.22446,0.35683,0.37985,0.56978,0.48345,0.56978,0.72517)
abundance_snail_year_3 <- (10^abundance_snail_year_3-1)*85


average_of_three_years <- (abundance_snail_year_1+abundance_snail_year_2+abundance_snail_year_3)/3

sample_parameters <- seq(from = 0, to = length(average_of_three_years)*14, by = 1)


#set seasonal temperature
seasonal_temperature <-data.frame(temp = temperature * (1 + epsilon * sin(2 * pi * (m+sample_parameters)/360)))


# Introduce eastivation function
in_estivation_afun <- approxfun(x = sample_parameters, y = in_est_function(seasonal_temperature$temp))
out_estivation_afun <- approxfun(x = sample_parameters, y = out_est_function(seasonal_temperature$temp))

# Set time span for model running 
run_time <- seq(from = 0, to = length(sample_parameters)-1, by = 1)

## solve the system 
model_outputs <- ode(y = c(S = average_of_three_years[1], SE = max(average_of_three_years)-average_of_three_years[1]), times = run_time, func = thermal_sensitive_model, parms = 0)


# Plotting data and simulation
pdf(file = "richard_toll_bulinus_one_years_abundance.pdf", width = 8, height = 8)
par(mar=c(6, 6, 5, 2), xpd = TRUE)
plot(model_outputs[,1], model_outputs[,2], col = 1,  pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2, main = "", ylim = c(min(average_of_three_years)-70, max(average_of_three_years)+150))  
points(seq(from = 14, to = length(average_of_three_years)*14, by = 14), average_of_three_years, col = 2, pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, cex.axis=2)

# Adding legends
legend("topleft", legend = c("Simulation", "Data"), col = c(1, 2), pch = 21, lwd = 3)
mtext(text = "Day(s)", side = 1, line = 4, cex = 2)
mtext(text = "Abundance", side = 2, line = 4, cex = 2)
title(main =  expression(paste(italic(Bulinus)," spp. in Richard Toll")), line = 2.5, cex.main = 1.5)

dev.off()

###### Run for three years 


sample_parameters <- seq(from = 0, to = length(average_of_three_years)*14*3, by = 1)


#set seasonal temperature
seasonal_temperature <-data.frame(temp = temperature * (1 + epsilon * sin(2 * pi * (m+sample_parameters)/360)))


# Introduce eastivation function
in_estivation_afun <- approxfun(x = sample_parameters, y = in_est_function(seasonal_temperature$temp))
out_estivation_afun <- approxfun(x = sample_parameters, y = out_est_function(seasonal_temperature$temp))

# Set time span for model running 
run_time <- seq(from = 0, to = length(sample_parameters)-1, by = 1)

## solve the system 
model_outputs <- ode(y = c(S = abundance_snail_year_1[1], SE = max(abundance_snail_year_1)-abundance_snail_year_1[1]), times = run_time, func = thermal_sensitive_model, parms = 0)

# Plotting data and simulation
pdf(file = "richard_toll_bulinus_three_years_abundance.pdf", width = 8, height = 8)
par(mar=c(6, 6, 5, 2), xpd = TRUE)
plot(model_outputs[,1], model_outputs[,2], col = 1,  pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2, main = "", ylim = c(min(average_of_three_years)-70, max(average_of_three_years)+150))  
points(seq(from = 14, to = length(average_of_three_years)*14*3, by = 14), c(abundance_snail_year_1, abundance_snail_year_2, abundance_snail_year_3), col = 2, pch = 21, lwd = 3, las = 1, ylab = "",  xlab = "", cex.lab = 2, cex.axis=2)

# Adding legends
legend("topleft", legend = c("Simulation", "Data"), col = c(1, 2), pch = 21, lwd = 3)
mtext(text = "Day(s)", side = 1, line = 4, cex = 2)
mtext(text = "Abundance", side = 2, line = 4, cex = 2)
title(main =  expression(paste(italic(Bulinus)," spp. in Richard Toll")), line = 2.5, cex.main = 1.5)

dev.off()



mid_temp <- 15
epsilon <- 0.1
year <- 3
sample_parameters <- seq(from = 0, to = 360*year, by = 1)


#set seasonal temperature
temp = mid_temp * (1 + epsilon * sin(2 * pi * (sample_parameters)/360))

plot(sample_parameters, temp)
plot(sample_parameters, in_estivation_afun(temp))










