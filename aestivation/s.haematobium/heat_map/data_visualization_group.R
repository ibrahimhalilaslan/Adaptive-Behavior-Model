## This script plot heatmap with africa temperature and gntd data set 


# Replace entries in out_come_season with 0 if they are greater than 100 or less than 2
out_come_season[out_come_season < .4] <- 0


# When using estimated clumper parameter for S. heamatobium 
wormPrevalenceSh <- function(M) {
  k <- exp(0.5186358*log(M) - 3.253653)
  p  <- 1 - (1+M/k)^(-k)    # fraction of humans with at least 1 parasites 
  return(p)
}


#convert mean parasite burden to prevalence 
for (i in 1:length(out_come_season[,1])){
  for (j in 1:length(out_come_season[1,])){
    out_come_season[i,j] <- wormPrevalenceSh(out_come_season[i,j])
  }
}


# save the prevalence as data frame for plot 
df <- data.frame(x = as.vector(temp_matrix), y = as.vector(epsilon_matrix), Prevalence = as.vector(out_come_season))

# Set minimum temperature 
min_tem <- 14;

# Set maximum temperature 
max_tem <- 31;

# Set minimum epsilon 
min_eps <- 0.02;

# Set maximum epsilon 
max_eps <- 0.25;

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
  Prevalence = sch_haematobium$percent_pos/100
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
    #new_gntd_data[j, i] <- median(gntd_data$Prevalence[combined_index], na.rm = TRUE)
    new_gntd_data[j, i] <- quantile(gntd_data$Prevalence[combined_index], 
                                   probs = 0.5, 
                                   na.rm = TRUE)
      }
}

# Replace NA entries with zero
new_gntd_data[is.na(new_gntd_data)] <- 0

# Convert the matrix to a data frame of three columns
new_gntd_data_fram <- data.frame(
  MAT = rep(temperatures, times = nrow(new_gntd_data)), # Second column: repeat row indices
  Epsilon = rep(epsilons, each = ncol(new_gntd_data)), # Third column: repeat column indices
  intensity = as.vector(t(new_gntd_data))                # First column: entries of the matrix
  )

new_gntd_data_fram <- new_gntd_data_fram[-which(new_gntd_data_fram$intensity == 0),]


# upload whole africa bioclimate variables   
##africa_temperature <- read.csv(file = 'Africa_BIOCLIM.csv')


# Convert bioclimate variables into seasonality 
#df_africa_temperature <- data.frame(
##  x = africa_temperature$BIO1,
 # y = (africa_temperature$BIO10-africa_temperature$BIO11)*pi/(africa_temperature$BIO1*4*sqrt(2))
#)

# Remove the temperature values out of our interest 
#df_africa_temperature_filtered <- df_africa_temperature[!(df_africa_temperature$x > 35 | df_africa_temperature$x <  12), ]
#df_africa_temperature_filtered <- df_africa_temperature_filtered[!(df_africa_temperature_filtered$y > 0.32 | df_africa_temperature_filtered$y < 0), ]


library(ggplot2)
library(dplyr)

# Your data frame (df) and gntd_data should be defined here

# Create a new variable for coloring
df <- df %>%
  mutate(Prev. = ifelse(Prevalence == 0, NA, Prevalence))


# Print the plot
pdf(file = "bio21_heat_map_haematobium_mean.pdf", width = 8, height = 5)
par(mar = c(6, 6, 6, 1), xpd = TRUE)

# Create the ggplot
ggplot() +
  # Plot df data with color based on Prev.
  geom_point(data = df, aes(x = x, y = y, colour = Prev.), size = 1.25, shape = 15, alpha = 1) + 
  scale_colour_gradient2(aes(colour = Prev.), low = "white", mid = "yellow", high = "red", midpoint = -0.25+ median(df$Prev., na.rm = TRUE), na.value = "#AED6F1", name = "Simulation Prev.") +
  
  # Plot gntd_data with color based on intensity
  geom_point(data = new_gntd_data_fram, aes(x = MAT, y = Epsilon, fill = intensity), size = 2, alpha = .5, shape = 16, color = "black") +
  scale_fill_gradient2(aes(fill = intensity), low = "white", high = "black", midpoint = .5, name = "GNTD Data Prev.") +
  
  # Switch the order of the legends without changing their appearance
  guides(colour = guide_colourbar(order = 1), fill = guide_colourbar(order = 2)) + 
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x = expression(paste("Mean annual temperature (", degree, "C)")), 
       y = "", #expression(paste("Seasonality  ", (epsilon))), 
       title = expression(paste(italic(S.), " ", italic(haematobium))* " with adaptive behaviors"))+
  theme(axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 15))
dev.off()



# Print the plot
pdf(file = "bio21_heat_map_haematobium.pdf", width = 8, height = 5)
par(mar = c(6, 6, 6, 1), xpd = TRUE)

# Create the ggplot
ggplot() +
  # Plot df data with color based on Prev.
  geom_point(data = df, aes(x = x, y = y, colour = Prev.), size = 1.25, shape = 15, alpha = 1) + 
  scale_colour_gradient2(aes(colour = Prev.), low = "white", mid = "yellow", high = "red", midpoint = -0.25+ median(df$Prev., na.rm = TRUE), na.value = "#AED6F1", name = "Simulation Prev.") +
  
  
  # Plot new_gntd_data_fram with intensity using fill aesthetic
  geom_point(data = gntd_data, aes(x = MAT, y = Epsilon, fill = Prevalence), 
             size = .75, shape = 21, stroke = 0, color = "black") + # Shape 21 supports fill
  scale_fill_gradient(low = "white", high = "black", name = "GNTD Data Prev.") +
  
  
  # Switch the order of the legends without changing their appearance
  guides(colour = guide_colourbar(order = 1), fill = guide_colourbar(order = 2)) + 
  
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x = expression(paste("Mean annual temperature (", degree, "C)")), 
       y = "", #expression(paste("Seasonality  ", (epsilon))), 
       title = expression(paste(italic(S.), " ", italic(haematobium))* " with adaptive behaviors"))+
  theme(axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 15))
dev.off()


# Print the plot
pdf(file = "bio21_heat_map_haematobium_trans.pdf", width = 8, height = 5)
par(mar = c(6, 6, 6, 1), xpd = TRUE)

# Create the ggplot
ggplot() +
  # Plot df data with color based on Prev.
  geom_point(data = df, aes(x = x, y = y, colour = Prev.), size = 1.25, shape = 15, alpha = 1) + 
  scale_colour_gradient2(aes(colour = Prev.), low = "white", mid = "yellow", high = "red", midpoint = -0.25+ median(df$Prev., na.rm = TRUE), na.value = "#AED6F1", name = "Simulation Prev.") +
  
  
  # Plot new_gntd_data_fram with intensity using fill aesthetic
  geom_point(data = gntd_data, aes(x = MAT, y = Epsilon, fill = Prevalence), alpha = 0.5,
          size = .75, shape = 21, stroke = 0, color = "black") + # Shape 21 supports fill
  scale_fill_gradient(low = "white", high = "black", name = "GNTD Data Prev.") +
  
  
  # Switch the order of the legends without changing their appearance
  guides(colour = guide_colourbar(order = 1), fill = guide_colourbar(order = 2)) + 
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x = expression(paste("Mean annual temperature (", degree, "C)")), 
       y = "", #expression(paste("Seasonality  ", (epsilon))), 
       title = expression(paste(italic(S.), " ", italic(haematobium))* " with adaptive behaviors"))+
  theme(axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 15))
dev.off()

