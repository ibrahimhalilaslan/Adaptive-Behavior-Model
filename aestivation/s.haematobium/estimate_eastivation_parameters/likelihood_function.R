#clear environment 
rm(list = ls())

#Load environment for the model without aestivation 
load(file ="non_aestivation_output.RData")

# Set minimum temperature 
min_tem <- 14;

# Set maximum temperature 
max_tem <- 31;


# Set minimum epsilon 
min_eps <- 0.02;

# Set maximum epsilon 
max_eps <- 0.26;

# set temperature steps 
temperature_step <- 0.5

# set epsilon steps 
epsilon_step <- 0.02


# Define temperature and epsilon ranges
temperatures <- seq(min_tem, max_tem, by = temperature_step)
epsilons <- rev(seq(min_eps, max_eps, by = epsilon_step))

# upload gntd data 
prev_percent <- read.csv(file = 'gntd_data_cases_number.csv')

# Filter S. haematobium data directly
sch_haematobium <- subset(prev_percent, parasite_s == "S. haematobium" & !is.na(bio01))

# calculate seasonality value for gntd data set 
gntd_data <- data.frame(
  MAT = sch_haematobium$bio01,
  Eps = (sch_haematobium$bio10 - sch_haematobium$bio11) * pi / (sch_haematobium$bio01 * 4 * sqrt(2)),
  Prev = sch_haematobium$total_number_pos,
  Cases = sch_haematobium$total_number_exam
)


# Pre-allocate matrix for results
cases_gntd_data <- matrix(NA, nrow = length(epsilons), ncol = length(temperatures))
total_test_gntd_data <- matrix(NA, nrow = length(epsilons), ncol = length(temperatures))

# Use nested loops for temperature and epsilon bins
# This generate GNTD data for each cell grid 
for (i in seq_along(temperatures)) {
  temp_min <- temperatures[i] - temperature_step/2
  temp_max <- temperatures[i] + temperature_step/2
  index_mat <- gntd_data$MAT > temp_min & gntd_data$MAT <= temp_max
  
  for (j in seq_along(epsilons)) {
    eps_min <- epsilons[j] - epsilon_step/2
    eps_max <- epsilons[j] + epsilon_step/2
    index_eps <- gntd_data$Eps > eps_min & gntd_data$Eps <= eps_max
    
    combined_index <- index_mat & index_eps
    if (sum(combined_index) == 0) next  # Skip this iteration if no matching data
    
    cases_gntd_data[j, i] <- sum(gntd_data$Prev[combined_index], na.rm = TRUE)
    total_test_gntd_data[j, i] <- sum(gntd_data$Cases[combined_index], na.rm = TRUE)
  }
}

# set lower bound value for log function to not heat inf
out_come_season_prev <- pmax(pmin(out_come_season_prev, 1 - 1e-15), 1e-15)

# Compute the log-likelihood
log_probs <- dbinom(
  x = cases_gntd_data,       
  size = total_test_gntd_data, 
  prob = out_come_season_prev,        
  log = TRUE
)

# Sum only finite values, ignoring NA
likelihood_sum_seasonal <- sum(log_probs, na.rm = TRUE)

#Likelihood value for eastivation model
likelihood_value_model_seasonal <- -2*likelihood_sum_seasonal

#Load simulation value for the model with constant temperature
load(file ="constant_temp.RData")

# Generate a matrix for likelihood to be calcualted
outputs_prev_mat <- matrix(rep(outputs_prev, times = nrow(cases_gntd_data)), 
                                            nrow = nrow(cases_gntd_data), 
                                            byrow = TRUE)

# set lower bound value for log function to not heat inf
outputs_prev_mat <- pmax(pmin(outputs_prev_mat, 1 - 1e-15), 1e-15)

# Compute the log-likelihood
log_probs <- dbinom(
  x = cases_gntd_data,       
  size = total_test_gntd_data, 
  prob = outputs_prev_mat,        
  log = TRUE
)


# Sum only finite values, ignoring NA
likelihood_sum_constant <- sum(log_probs, na.rm = TRUE)

#Likelihood value for eastivation model
likelihood_value_model_constant <- -2*likelihood_sum_constant


# Load the model simulation with aestivation function
#prevalence_percentage <- readRDS("prevalence_percentage.rds")
prevalence_percentage <- readRDS("prevalence_list.rds")



# Initialize a vector to store likelihood values
likelihood_sum_aestivation <- numeric(length(prevalence_percentage))

# Loop over k matrices
for (k in seq_along(prevalence_percentage)) {
  # Ensure the probability matrix has valid values (0 < p < 1)
  prob_matrix <- prevalence_percentage[[k]]
  prob_matrix[prob_matrix <= 0] <- 1e-15  # Avoid log(0) issues
  prob_matrix[prob_matrix >= 1] <- 1 - 1e-15
  
  # Compute the log-likelihood
  likelihood_sum_aestivation[k] <- sum(dbinom(
    x = cases_gntd_data,       # Observed cases
    size = total_test_gntd_data, # Number of trials
    prob = prob_matrix,        # Probability from model
    log = TRUE
  ), na.rm = TRUE)  # Handle NA values
}

# Recall aestivation parameters
# We start with setting our eastivation function
aest_temp_seq <- seq(from = 22.5, to = 27, by = .5)  

# set the rate of snails move to the estivation stages
aest_rate_seq <- seq(from = 0.005, to = 0.1, by = 0.02)

# Set the steepness of estivation function for hot temperature 
steepness_seq <-  seq(from = 0.1, to = 1.5, by = 0.2)

# Generate all combinations
combinations_1 <- expand.grid(aest_temp_seq, aest_rate_seq, steepness_seq)

# We start with setting our eastivation function
aest_temp_seq <- seq(from = 20, to = 22, by = .5)  

# set the rate of snails move to the estivation stages
aest_rate_seq <- seq(from = 0.005, to = 0.1, by = 0.02)

# Set the steepness of estivation function for hot temperature 
steepness_seq <-  seq(from = 0.1, to = 1.5, by = 0.2)

# Generate all combinations
combinations_2 <- expand.grid(aest_temp_seq, aest_rate_seq, steepness_seq)

# combine the parameter sets 
combinations <- rbind(combinations_1, combinations_2)

#Likelihood value for eastivation model
likelihood_value_model_aestivation <- -2*likelihood_sum_aestivation[which.max(likelihood_sum_aestivation)]+2*3
#The best set of the parameters 
combinations[which.max(likelihood_sum_aestivation),]

#remove the temperature smaller than 22.5
index_remove <- which(combinations[,1] <= 22.5)
likelihood_sum_aestivation_remove <- likelihood_sum_aestivation[-index_remove]
likelihood_value_model_aestivation_remove <- -2*likelihood_sum_aestivation_remove[which.max(likelihood_sum_aestivation_remove)]+2*3
combinations_remove <- combinations[-index_remove,]
combinations_remove[which.max(likelihood_sum_aestivation_remove),]
