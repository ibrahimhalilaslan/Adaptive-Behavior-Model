#clear environment 
rm(list = ls())

load("my_environment_10.RData")
square_error_1 <- square_error

load("my_environment_20.RData")
square_error_2 <- square_error

load("my_environment_30.RData")
square_error_3 <- square_error

load("my_environment_40.RData")
square_error_4 <- square_error

load("my_environment_50.RData")
square_error_5 <- square_error

load("my_environment_60.RData")
square_error_6 <- square_error


load("my_environment_70.RData")
square_error_7 <- square_error


load("my_environment_80.RData")
square_error_8 <- square_error


load("my_environment_90.RData")
square_error_9 <- square_error

load("my_environment_100.RData")
square_error_10 <- square_error

load("my_environment_110.RData")
square_error_11 <- square_error

load("my_environment_120.RData")
square_error_12 <- square_error

load("my_environment_125.RData")
square_error_13 <- square_error

load("my_environment_132.RData")
square_error_14 <- square_error

error <- c(square_error_1, square_error_2, square_error_3, square_error_4, square_error_5,
    square_error_6, square_error_7, square_error_8, square_error_9, square_error_10,
    square_error_11, square_error_12, square_error_13, square_error_14)

# We start with setting our eastivation function
# set the temperature of half saturation for snails estivation function at high temperatures 
aest_temp_seq <- seq(from = 22, to = 27, by = .5)

# set the rate of snails move to the estivation stages
aest_rate_seq <- seq(from = 0.01, to = 0.12, by = 0.05)

# Set the steepness of estivation function for hot temperature 
steepness_seq <-  seq(from = 0.1, to = 1.6, by = 0.5)

# Generate all combinations
combinations <- expand.grid(aest_temp_seq, aest_rate_seq, steepness_seq)
new_data_frame <-  data.frame(combinations,error)

index_temp <- which(combinations$Var1 == 22 | combinations$Var1 == 22.5)

new_data_frame <- new_data_frame[-index_temp,]

new_data_frame[which.min(new_data_frame$error),]



