## This code is to plot generate the data frame for constant temperature 
#clear environment 
rm(list = ls())


# Load required libraries
library(maps)
library(ggplot2)
library(sf)
library(tidyr)
library(stars)
library(sfheaders)
library(scales)

#load(file ="december_3.RData")
load(file ="august_2.RData")

prevalence_africa <- read.csv(file = 'Sh_simulation_africa_time_series_mpb.csv')

index <- which(prevalence_africa$MPB < 0.4)

prevalence_africa$MPB[index] <- 0


prevalence_africa$constant_temp_simulation <- numeric(length(prevalence_africa$Prevalence))


#for (i in 1:length(temperature)){
#  index_mat <- which(prevalence_africa$MAT <= temperature[i] + 0.1 & prevalence_africa$MAT > temperature[i]) 
#  prevalence_africa$constant_temp_simulation[index_mat] <- out_come_season[1, i]
#}

for (i in 1:length(temperature)){
  index_mat <- which(prevalence_africa$MAT <= temperature[i] + 0.1 & prevalence_africa$MAT > temperature[i]) 
  prevalence_africa$constant_temp_simulation[index_mat] <- outputs_1_mpb[i]
}


# When using estimated clumper parameter for S. heamatobium 
wormPrevalenceSh <- function(M) {
  k <- exp(0.5186358*log(M) - 3.253653)
  p  <- 1 - (1+M/k)^(-k)    # fraction of humans with at least 1 parasites 
  return(p)
}

prevalence_africa$constant_temp_simulation <- wormPrevalenceSh(prevalence_africa$constant_temp_simulation)
prevalence_africa$Prevalence <- wormPrevalenceSh(prevalence_africa$MPB)

# Get world map data
world_map <- map_data("world")

# Filter data for Africa
africa_map <- subset(world_map, region %in% c("Algeria", "Angola", "Benin", "Botswana", "Burkina Faso", "Western Sahara",
                                              "Burundi", "Cameroon", "Cape Verde", "Central African Republic",
                                              "Chad", "Comoros", "Democratic Republic of the Congo", "Republic of Congo",
                                              "Djibouti", "Egypt", "Equatorial Guinea", "Eritrea", "Eswatini", 
                                              "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
                                              "Ivory Coast", "Kenya", "Lesotho", "Liberia", "Libya", "Madagascar", 
                                              "Malawi", "Mali", "Mauritania", "Mauritius", "Morocco", "Mozambique", 
                                              "Namibia", "Niger", "Nigeria", "Rwanda", "Sao Tome and Principe", 
                                              "Senegal", "Seychelles", "Sierra Leone", "Somalia", "South Africa", 
                                              "South Sudan", "Sudan", "Tanzania", "Togo", "Tunisia", "Uganda", 
                                              "Zambia", "Zimbabwe"))


africa_poly <- st_union(sfheaders::sf_polygon(africa_map, x = "long", y = "lat", polygon_id = "group"))
africa_poly <- st_set_crs(africa_poly, 4326)


# upload prevalence data 
df_africa_mask <- read.csv(file = 'africa_mask_v1.csv')
df_africa_mask <- replace_na(complete(df_africa_mask, lon, lat), list(grid_code = 0))

df_africa_constant <- data.frame(Longitude = prevalence_africa$longitude, Latitude = prevalence_africa$latitude, 
                                constant_temp_simulation = prevalence_africa$constant_temp_simulation)

df_africa_seasonal <- data.frame(Longitude = prevalence_africa$longitude, Latitude = prevalence_africa$latitude, 
                                 seasonal_temperature_simulation = prevalence_africa$Prevalence)


rast_constant <- st_set_crs(st_as_stars(df_africa_constant, coords = c("Longitude", "Latitude")), 4326)
mask_rast <- st_set_crs(st_as_stars(df_africa_mask, coords = c("lat", "lon")), 4326)
mask_rast <- st_warp(mask_rast, rast_constant)
rast_constant$constant_temp_simulation[mask_rast$grid_code == 1] <- NA
rast_constant <- rast_constant[africa_poly]


rast_seasonal <- st_set_crs(st_as_stars(df_africa_seasonal, coords = c("Longitude", "Latitude")), 4326)
mask_rast <- st_set_crs(st_as_stars(df_africa_mask, coords = c("lat", "lon")), 4326)
mask_rast <- st_warp(mask_rast, rast_seasonal)
rast_seasonal$seasonal_temperature_simulation[mask_rast$grid_code == 1] <- NA
rast_seasonal <- rast_seasonal[africa_poly]


rast <- rast_seasonal
rast$prevalence <- rast_constant$constant_temp_simulation - rast_seasonal$seasonal_temperature_simulation



# Load the R environment from the .RData file

# Read the CSV file
my_data <- read.csv("my_vector.csv")

daily_temps <- as.numeric(my_data)

# Create the inset plot
inset_plot_vary_temp <- ggplot(data.frame(1:length(daily_temps), daily_temps), aes(x = 1:length(daily_temps), y = daily_temps)) +
  geom_line(color = "blue") +
  theme_classic() +
  labs(title = "", x = "Days", y = "Temp.") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
    axis.title = element_text(size = 15), 
    axis.text = element_text(size = 10),
    plot.margin = margin(0, 0, 0, 0)) # Remove margins for tight placement

# Create the inset plot
inset_plot_cons_temp <- ggplot(data.frame(1:length(daily_temps), daily_temps), aes(x = 1:length(daily_temps), y = mean(daily_temps))) +
  geom_line(color = "blue") +
  theme_classic() +
  labs(title = " ", x = "Days", y = "Temp.") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
    axis.title = element_text(size = 15), 
    axis.text = element_text(size = 10),
    plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1))# Remove margins for tight placement


# Create the main plot and store it in an object
main_plot <- ggplot() + 
  geom_sf(data = africa_poly, fill = "#AED6F1", color = "#AED6F1", linewidth = 10) + 
  geom_sf(data = africa_poly, fill = "grey40") + 
  geom_stars(aes(fill = prevalence), data = rast, alpha = 1) + 
  scale_fill_gradient2(low = "yellow", mid = "white", high = "red", midpoint = 0, na.value = NA 
                       ) +
  
  geom_map(data = africa_map, map = africa_map,
           aes(map_id = region),
           fill = NA, color = "black", size = 0.5) +
  labs(title = expression(paste(italic(S.), " ", italic(haematobium))), 
       fill = "Prevalence Diff.") +  # Add title and color legend title
  theme_void() +
  coord_sf(ylim = c(-35, max(africa_map$lat))) +
  guides(fill = guide_colorbar(barwidth = unit(10, "lines"), title.position = "top", 
                               title.hjust = 0.5, ticks.colour = "black", frame.colour = "black")) +
  theme(  # Adjust plot margins
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = c(0.225, 0.4), 
    legend.direction = "horizontal",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)) #+  # Center plot title
  
  

# Plot the map with intensity of color based on value
pdf(file = "s_heamatobium_africa_prevalence_difference.pdf", width = 8, height = 8)
print(main_plot)
dev.off()