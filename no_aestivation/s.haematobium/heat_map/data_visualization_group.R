## This script plots a heatmap with Africa temperature and GNTD dataset 

# Replace small values with zero
out_come_season[out_come_season < .3] <- 0

# define function to convert MPB into prevalence 
wormPrevalenceSh <- function(M) {
  k <- exp(0.5186358 * log(M) - 3.253653)
  p <- 1 - (1 + M / k)^(-k)
  return(p)
}

# Apply prevalence function to outcome matrix
for (i in 1:nrow(out_come_season)) {
  for (j in 1:ncol(out_come_season)) {
    out_come_season[i, j] <- wormPrevalenceSh(out_come_season[i, j])
  }
}

# Convert to data frame for plotting
df <- data.frame(x = as.vector(temp_matrix), 
                 y = as.vector(epsilon_matrix), 
                 Prevalence = as.vector(out_come_season))

# Define temperature and seasonality ranges
min_tem <- 14; max_tem <- 31; min_eps <- 0.02; max_eps <- 0.25
temperature_step <- 0.5; epsilon_step <- 0.02
temperatures <- seq(min_tem, max_tem, by = temperature_step)
epsilons <- seq(min_eps, max_eps, by = epsilon_step)

# Load GNTD data
prev_percent <- read.csv('gntd_vars_all.csv')
sch_haematobium <- subset(prev_percent, parasite_s == "S. haematobium" & 
                            percent_pos != 0 & !is.na(bio01))

# Convert bioclimate variables to seasonality
gntd_data <- data.frame(
  MAT = sch_haematobium$bio01,
  Epsilon = (sch_haematobium$bio10 - sch_haematobium$bio11) * pi / 
    (sch_haematobium$bio01 * 4 * sqrt(2)),
  Prevalence = sch_haematobium$percent_pos / 100
)

# Pre-allocate results matrix
new_gntd_data <- matrix(NA, nrow = length(epsilons), ncol = length(temperatures))

# Bin and compute median prevalence
# we assign each gntd data into grid cell
for (i in seq_along(temperatures)) {
  temp_min <- temperatures[i]
  temp_max <- temp_min + temperature_step
  index_mat <- gntd_data$MAT > temp_min & gntd_data$MAT <= temp_max
  
  for (j in seq_along(epsilons)) {
    eps_min <- epsilons[j]
    eps_max <- eps_min + epsilon_step
    index_eps <- gntd_data$Epsilon > eps_min & gntd_data$Epsilon <= eps_max
    combined_index <- index_mat & index_eps
    new_gntd_data[j, i] <- quantile(gntd_data$Prevalence[combined_index], 0.5, na.rm = TRUE)
  }
}

# Replace NA with zero
new_gntd_data[is.na(new_gntd_data)] <- 0

# Convert matrix to data frame
new_gntd_data_fram <- data.frame(
  MAT = rep(temperatures, times = nrow(new_gntd_data)),
  Epsilon = rep(epsilons, each = ncol(new_gntd_data)),
  intensity = as.vector(t(new_gntd_data))
)
new_gntd_data_fram <- new_gntd_data_fram[new_gntd_data_fram$intensity != 0, ]

# Plotting
library(ggplot2); library(dplyr)

# Add NA for zero prevalence
df <- df %>% mutate(Prev. = ifelse(Prevalence == 0, NA, Prevalence))

# Create heatmap
pdf("bio21_heat_map_haematobium.pdf", width = 8, height = 5)
ggplot() +
  geom_point(data = df, aes(x = x, y = y, colour = Prev.), size = 1.25, shape = 15) + 
  scale_colour_gradient2(low = "white", mid = "yellow", high = "red", 
                         midpoint = median(df$Prev., na.rm = TRUE), 
                         na.value = "#AED6F1", name = "Simulation Prev.") +
  geom_point(data = new_gntd_data_fram, aes(x = MAT, y = Epsilon, fill = intensity), 
             size = 2, alpha = .5, shape = 16, color = "black") +
  scale_fill_gradient2(low = "white", high = "black", midpoint = .5, name = "GNTD Data Prev.") +
  guides(colour = guide_colourbar(order = 1), fill = guide_colourbar(order = 2)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 20, color = "black"), 
        plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
        legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  labs(x = "", y = "", 
       title = expression(paste(italic(S.), " ", italic(haematobium))))
dev.off()