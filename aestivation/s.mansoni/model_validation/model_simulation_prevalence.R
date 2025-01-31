
# upload gntd data 
prev_percent <- read.csv(file = 'gntd_vars_all.csv')


# Extract S.mansoni data set from gntd data set
sch_mansoni_ind <- which(prev_percent$parasite_s == "S. mansoni")

sch_mansoni <- prev_percent[sch_mansoni_ind, ]
sch_mansoni <- sch_mansoni[-which(sch_mansoni$percent_pos == 0),]
sch_mansoni <- sch_mansoni[-c(which(is.na(sch_mansoni$bio01))),]


# Convert Bioclimate variables into seasonality 
gntd_data <- data.frame(
  MAT = sch_mansoni$bio01,
  Epsilon = (sch_mansoni$bio10-sch_mansoni$bio11)*pi/(sch_mansoni$bio01*4*sqrt(2)),
  Prevalence = sch_mansoni$percent_pos
)


load(file ="non_aestivation_heat_map.RData")

gntd_constant_temp_non_aestivation <- c()

for (i in 1:length(gntd_data$MAT)){
  index_mat <- which(gntd_data$MAT[i] <= temp_matrix[1,] & gntd_data$MAT[i]+.2 > temp_matrix[1,]) 
  gntd_constant_temp_non_aestivation[i] <- ifelse(length(index_mat) == 0, out_come_season[1, 1], out_come_season[1, index_mat])
  
}

gntd_non_constant_temp_non_aestivation <- c()

for (i in 1:length(gntd_data$MAT)){
  index_mat <- which(gntd_data$MAT[i] <= temp_matrix[1,] & gntd_data$MAT[i]+.2 > temp_matrix[1,]) 
  index_eps <- which(gntd_data$Epsilon[i] <= epsilon_matrix[,1] & gntd_data$Epsilon[i]+.01 > epsilon_matrix[,1]) 
  gntd_non_constant_temp_non_aestivation[i] <- ifelse(length(index_mat) == 0, out_come_season[1, 1], out_come_season[index_eps, index_mat])
  
}

load(file ="aestivation_heat_map.RData")

gntd_constant_temp_aestivation <- c()

for (i in 1:length(gntd_data$MAT)){
  index_mat <- which(gntd_data$MAT[i] <= temp_matrix[1,] & gntd_data$MAT[i]+.2 > temp_matrix[1,]) 
  gntd_constant_temp_aestivation[i] <- ifelse(length(index_mat) == 0, out_come_season[1,1], out_come_season[1, index_mat])
  
}


gntd_non_constant_temp_aestivation <- c()

for (i in 1:length(gntd_data$MAT)){
  index_mat <- which(gntd_data$MAT[i] <= temp_matrix[1,] & gntd_data$MAT[i]+.2 > temp_matrix[1,]) 
  index_eps <- which(gntd_data$Epsilon[i] <= epsilon_matrix[,1] & gntd_data$Epsilon[i]+.01 > epsilon_matrix[,1]) 
  gntd_non_constant_temp_aestivation[i] <- ifelse(length(index_mat) == 0, out_come_season[1,1], out_come_season[index_eps, index_mat])
  
}



wormPrevalenceSm <- function(M) {
  k <- exp(0.61521*log(M) - 4.844146)
  p  <- 1 - (1+M/k)^(-k)    # fraction of humans with at least 1 parasites 
  return(p)
}


gntd_constant_temp_non_aestivation[is.na(gntd_constant_temp_non_aestivation) | gntd_constant_temp_non_aestivation < .76] <- 0
gntd_non_constant_temp_non_aestivation[is.na(gntd_non_constant_temp_non_aestivation) | gntd_non_constant_temp_non_aestivation < .76] <- 0
gntd_constant_temp_aestivation[is.na(gntd_constant_temp_aestivation) | gntd_constant_temp_aestivation < .74] <- 0
gntd_non_constant_temp_aestivation[is.na(gntd_non_constant_temp_aestivation) | gntd_non_constant_temp_aestivation < .74] <- 0


gntd_constant_temp_non_aestivation <- wormPrevalenceSm(gntd_constant_temp_non_aestivation)
gntd_non_constant_temp_non_aestivation <- wormPrevalenceSm(gntd_non_constant_temp_non_aestivation)
gntd_constant_temp_aestivation <- wormPrevalenceSm(gntd_constant_temp_aestivation)
gntd_non_constant_temp_aestivation <- wormPrevalenceSm(gntd_non_constant_temp_aestivation)




df_frame <- data.frame(
  Lat = sch_mansoni$latitude, 
  Long = sch_mansoni$longitude,
  MAT = gntd_data$MAT,
  Epsilon = gntd_data$Epsilon,
  Prevalence = gntd_data$Prevalence,
  Constant_Temp_Non_Aestivation = gntd_constant_temp_non_aestivation*100,
  Non_Constant_Temp_Non_Aestivation = gntd_non_constant_temp_non_aestivation*100,
  Constant_Temp_Aestivation = gntd_constant_temp_aestivation*100,
  Non_Constant_Temp_Aestivation = gntd_non_constant_temp_aestivation*100)

################# PLOT OBSERVE VERSA PREDICTED ###################

# Scatter plot with customized color scale for a continuous variable
ggplot(data = df_frame) +
  geom_point(aes(x = Prevalence, y = Non_Constant_Temp_Non_Aestivation, color = MAT), size = 0.5, shape = 16) +
  scale_color_gradient(low = "yellow", high = "red") +  # Customize the color gradient
  labs(title = "Seasonality with No Adaptive behavior", 
       x = "Observed", 
       y = "Estimated") +
  theme_minimal()

# Scatter plot with customized color scale for a continuous variable
ggplot(data = df_frame) +
  geom_point(aes(x = Prevalence, y = Non_Constant_Temp_Aestivation, color = MAT), size = 0.5, shape = 16) +
  scale_color_gradient(low = "yellow", high = "red") +  # Customize the color gradient
  labs(title = "Seasonality with Adaptive behavior", 
       x = "Observed", 
       y = "Estimated") +
  theme_minimal()


# Scatter plot with customized color scale for a continuous variable
ggplot(data = df_frame) +
  geom_point(aes(x = Prevalence, y = Non_Constant_Temp_Non_Aestivation, color = Epsilon), size = 0.5, shape = 16) +
  scale_color_gradient(low = "yellow", high = "red") +  # Customize the color gradient
  labs(title = "Seasonality with No Adaptive behavior", 
       x = "Observed", 
       y = "Estimated") +
  theme_minimal()




# Scatter plot with customized color scale for a continuous variable
ggplot(data = df_frame) +
  geom_point(aes(x = Prevalence, y = Non_Constant_Temp_Aestivation, color = Epsilon), size = 0.5, shape = 16) +
  scale_color_gradient(low = "yellow", high = "red") +  # Customize the color gradient
  labs(title = "Seasonality with Adaptive behavior", 
       x = "Observed", 
       y = "Estimated") +
  theme_minimal()

################# PLOT MAT VERSA OBSERVE VERSA PREDICTED ###################

index_eps_1 <- which(df_frame$Epsilon < 0.1)
index_eps_2 <- which(df_frame$Epsilon > 0.1 & df_frame$Epsilon <= 0.15)
index_eps_3 <- which(df_frame$Epsilon > 0.15 & df_frame$Epsilon <= 0.2)
index_eps_4 <- which(df_frame$Epsilon > 0.2 & df_frame$Epsilon <= 0.25)
index_eps_5 <- which(df_frame$Epsilon > 0.25 & df_frame$Epsilon <= 0.3)




# Scatter plot with a line for observed data
ggplot(data = df_frame) +
  geom_point(aes(x = MAT, y = Prevalence), color = "black", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_eps_1, ], aes(x = MAT, y = Non_Constant_Temp_Non_Aestivation), size = 0.5, shape = 16, color = "red") +
  geom_point(data = df_frame[index_eps_3, ], aes(x = MAT, y = Non_Constant_Temp_Non_Aestivation), size = 0.5, shape = 16, color = "blue") +
  geom_point(data = df_frame[index_eps_4, ], aes(x = MAT, y = Non_Constant_Temp_Non_Aestivation), size = 0.5, shape = 16, color = "yellow") +
  labs(title = "Seasonality with No Adaptive Behavior", 
       x = "MAT", 
       y = "Estimated and Observed") +
  theme_minimal()


# Scatter plot with a line for observed data
ggplot(data = df_frame) +
  geom_point(aes(x = MAT, y = Prevalence), color = "black", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_eps_1, ], aes(x = MAT, y = Non_Constant_Temp_Aestivation), size = 0.5, shape = 16, color = "red") +
  geom_point(data = df_frame[index_eps_3, ], aes(x = MAT, y = Non_Constant_Temp_Aestivation), size = 0.5, shape = 16, color = "blue") +
  geom_point(data = df_frame[index_eps_4, ], aes(x = MAT, y = Non_Constant_Temp_Aestivation), size = 0.5, shape = 16, color = "yellow") +
  labs(title = "Seasonality with Adaptive Behavior", 
       x = "MAT", 
       y = "Estimated and Observed") +
  theme_minimal()


################# PLOT OBSERVE VERSA PREDICTED UNDER DIFFERENT EPSILON AND MAT ###################
##### with MAT layer 
index_mat_1 <- which(df_frame$MAT < 15)
index_mat_2 <- which(df_frame$MAT > 15 & df_frame$MAT <= 20)
index_mat_3 <- which(df_frame$MAT > 20 & df_frame$MAT <= 25)
index_mat_4 <- which(df_frame$MAT > 25 & df_frame$MAT <= 30)
index_mat_5 <- which(df_frame$MAT > 30 & df_frame$MAT <= 35)


ggplot() +
  geom_point(data = df_frame[index_mat_2, ], aes(x = Prevalence, y = Non_Constant_Temp_Non_Aestivation), color = "red", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_mat_3, ], aes(x = Prevalence, y = Non_Constant_Temp_Non_Aestivation), color = "blue", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_mat_4, ], aes(x = Prevalence, y = Non_Constant_Temp_Non_Aestivation), color = "yellow", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_mat_5, ], aes(x = Prevalence, y = Non_Constant_Temp_Non_Aestivation), color = "orange", size = 0.5, shape = 16) +
  labs(
    title = "Seasonality with No Adaptive Behavior", 
    x = "Observed", 
    y = "Estimated"
  ) +
  theme_minimal()



ggplot() +
  geom_point(data = df_frame[index_mat_2, ], aes(x = Prevalence, y = Non_Constant_Temp_Aestivation), color = "red", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_mat_3, ], aes(x = Prevalence, y = Non_Constant_Temp_Aestivation), color = "blue", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_mat_4, ], aes(x = Prevalence, y = Non_Constant_Temp_Aestivation), color = "yellow", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_mat_5, ], aes(x = Prevalence, y = Non_Constant_Temp_Aestivation), color = "orange", size = 0.5, shape = 16) +
  labs(
    title = "Seasonality with Adaptive Behavior", 
    x = "Observed", 
    y = "Estimated"
  ) +
  theme_minimal()

####### with Epsilon layer 

ggplot() +
  geom_point(data = df_frame[index_eps_2, ], aes(x = Prevalence, y = Non_Constant_Temp_Non_Aestivation), color = "red", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_eps_3, ], aes(x = Prevalence, y = Non_Constant_Temp_Non_Aestivation), color = "blue", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_eps_4, ], aes(x = Prevalence, y = Non_Constant_Temp_Non_Aestivation), color = "yellow", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_eps_5, ], aes(x = Prevalence, y = Non_Constant_Temp_Non_Aestivation), color = "orange", size = 0.5, shape = 16) +
  labs(
    title = "Seasonality with No Adaptive Behavior", 
    x = "Observed", 
    y = "Estimated"
  ) +
  theme_minimal()



ggplot() +
  geom_point(data = df_frame[index_eps_2, ], aes(x = Prevalence, y = Non_Constant_Temp_Aestivation), color = "red", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_eps_3, ], aes(x = Prevalence, y = Non_Constant_Temp_Aestivation), color = "blue", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_eps_4, ], aes(x = Prevalence, y = Non_Constant_Temp_Aestivation), color = "yellow", size = 0.5, shape = 16) +
  geom_point(data = df_frame[index_eps_5, ], aes(x = Prevalence, y = Non_Constant_Temp_Aestivation), color = "orange", size = 0.5, shape = 16) +
  labs(
    title = "Seasonality with Adaptive Behavior", 
    x = "Observed", 
    y = "Estimated"
  ) +
  theme_minimal()

############################ Write the data in csv file ###############

# Write the data frame to a CSV file
#write.csv(df_frame, file = "s_mansoni_simulation_data.csv", row.names = FALSE)

# Write the data frame to a CSV file
write.csv(df_frame, file = "s_mansoni_simulation_data_without_zeros.csv", row.names = FALSE)






