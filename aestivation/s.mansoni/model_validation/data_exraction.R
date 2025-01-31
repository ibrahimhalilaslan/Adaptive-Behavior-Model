
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


gntd_constant_temp_non_aestivation[is.na(gntd_constant_temp_non_aestivation) | gntd_constant_temp_non_aestivation < .8] <- 0
gntd_non_constant_temp_non_aestivation[is.na(gntd_non_constant_temp_non_aestivation) | gntd_non_constant_temp_non_aestivation < .8] <- 0
gntd_constant_temp_aestivation[is.na(gntd_constant_temp_aestivation) | gntd_constant_temp_aestivation < .8] <- 0
gntd_non_constant_temp_aestivation[is.na(gntd_non_constant_temp_aestivation) | gntd_non_constant_temp_aestivation < .8] <- 0


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


# Load ggplot2
library(ggplot2)

# Scatter plot with points and lines
ggplot(df_frame, aes(x = MAT, y = Prevalence)) +
  geom_point(color = "red", size = 1) +  # Use a reasonable size for visibility
  geom_point(aes(y = Constant_Temp_Non_Aestivation), color = "black") +  # Use geom_line to add a line
  labs(title = "Scatter Plot with Line", x = "MAT", y = "Prevalence / Constant Temp Non-Aestivation") +
  theme_minimal()

# Scatter plot with points and lines
ggplot(df_frame, aes(x = MAT, y = Prevalence)) +
  geom_point(color = "red", size = .5) +  # Use a reasonable size for visibility
  geom_point(aes(y = Non_Constant_Temp_Non_Aestivation), color = "black", size = .5) +  # Use geom_line to add a line
  labs(title = "Scatter Plot with Line", x = "MAT", y = "Prevalence / Constant Temp Non-Aestivation") +
  theme_minimal()


# Scatter plot with points and lines
ggplot(df_frame, aes(x = MAT, y = Prevalence)) +
  geom_point(color = "red", size = .5) +  # Use a reasonable size for visibility
  geom_point(aes(y = Constant_Temp_Aestivation), color = "black", size = .5) +  # Use geom_line to add a line
  labs(title = "Scatter Plot with Line", x = "MAT", y = "Prevalence / Constant Temp Non-Aestivation") +
  theme_minimal()

# Scatter plot with points and lines
ggplot(df_frame, aes(x = MAT, y = Prevalence)) +
  geom_point(color = "red", size = .5) +  # Use a reasonable size for visibility
  geom_point(aes(y = Non_Constant_Temp_Aestivation), color = "black", size = .5) +  # Use geom_line to add a line
  labs(title = "Scatter Plot with Line", x = "MAT", y = "Prevalence / Constant Temp Non-Aestivation") +
  theme_minimal()

# Write the data frame to a CSV file
#write.csv(df_frame, file = "s_mansoni_simulation_data.csv", row.names = FALSE)

# Write the data frame to a CSV file
write.csv(df_frame, file = "s_mansoni_simulation_data_without_zeros.csv", row.names = FALSE)

#Constant_Temp_Non_Aestivation
#Non_Constant_Temp_Non_Aestivation
#Constant_Temp_Aestivation
#Non_Constant_Temp_Aestivation


# Create the ggplot
ggplot() +
  # Plot df data with color based on Prev.
  geom_point(data = df_frame, aes(x = MAT, y = Epsilon, colour = Non_Constant_Temp_Aestivation), size = 1, shape = 15, alpha = 1) + 
  scale_colour_gradient2(aes(colour = Non_Constant_Temp_Aestivation), low = "white", mid = "yellow", high = "red", na.value = "#AED6F1", name = "Simulation Prev.") +
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x = expression(paste("Mean annual temperature (", degree, "C)")), 
       y = expression(paste("Seasonality  ", (epsilon))), 
       title = "") + 
  theme(axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 15))  # Increase legend title size





