library(rms)
library(MASS)
library(pROC)
library(ggplot2)

# Read data for both species
df_m <- read.csv('s_mansoni_simulation_data.csv')
df_h <- read.csv('s_haematobium_simulation_data.csv')

# Remove any missing data
df_m <- df_m[complete.cases(df_m),]
df_h <- df_h[complete.cases(df_h),]



#### Calculate RMSE for each simulation

# S. Mansoni
sqrt(mean((df_m$Constant_Temp_Non_Aestivation - df_m$Prevalence)^2)) # Constant temp

#SS_r <- sum((df_m$Constant_Temp_Non_Aestivation - df_m$Prevalence)^2) # Constant temp
#SS_r_m <- mean((df_m$Constant_Temp_Non_Aestivation - df_m$Prevalence)^2) # Constant temp
#SS_t <- sum((df_m$Prevalence - mean(df_m$Prevalence))^2) # Constant temp
#1-SS_r/SS_t

sqrt(mean((df_m$Non_Constant_Temp_Non_Aestivation - df_m$Prevalence)^2)) # Seasonal, no estivation
#SS_r <- sum((df_m$Non_Constant_Temp_Non_Aestivation - df_m$Prevalence)^2) # Seasonal, no estivation
#SS_r_m <- mean((df_m$Non_Constant_Temp_Non_Aestivation - df_m$Prevalence)^2) # Seasonal, no estivation
#1-SS_r/SS_t

sqrt(mean((df_m$Non_Constant_Temp_Aestivation - df_m$Prevalence)^2)) # Seasonal, estivation
#SS_r <- sum((df_m$Non_Constant_Temp_Aestivation - df_m$Prevalence)^2) # Seasonal, estivation
#SS_r_m <- mean((df_m$Non_Constant_Temp_Aestivation - df_m$Prevalence)^2) # Seasonal, estivation
#1-SS_r/SS_t


# S. Haematobium
sqrt(mean((df_h$Constant_Temp_Non_Aestivation - df_h$Prevalence)^2)) # Constant temp
#SS_r_m <- mean((df_h$Constant_Temp_Non_Aestivation - df_h$Prevalence)^2) # Constant temp
#SS_t <- sum((df_h$Prevalence - mean(df_h$Prevalence))^2) # Constant temp
#1-SS_r/SS_t

sqrt(mean((df_h$Non_Constant_Temp_Non_Aestivation - df_h$Prevalence)^2)) # Seasonal, no estivation
#SS_r <- sum((df_h$Non_Constant_Temp_Non_Aestivation - df_h$Prevalence)^2) # Seasonal, no estivation
#SS_r_m <- mean((df_h$Non_Constant_Temp_Non_Aestivation - df_h$Prevalence)^2) # Seasonal, no estivation
#1-SS_r/SS_t

sqrt(mean((df_h$Non_Constant_Temp_Aestivation - df_h$Prevalence)^2)) # Seasonal, estivation
#SS_r <- sum((df_h$Non_Constant_Temp_Aestivation - df_h$Prevalence)^2) # Seasonal, estivation
#SS_r_m <- mean((df_h$Non_Constant_Temp_Aestivation - df_h$Prevalence)^2) # Seasonal, estivation
#1-SS_r/SS_t


#### Are the squared error rates different with seasonality?

# S. Mansoni 
t.test((df_m$Non_Constant_Temp_Non_Aestivation - df_m$Prevalence)^2, 
       (df_m$Constant_Temp_Non_Aestivation - df_m$Prevalence)^2,
       paired = TRUE)

#Maurice suggestion 
ll_const <- as.numeric(logLik(lm(df_m$Prevalence ~ offset(df_m$Constant_Temp_Non_Aestivation))))
ll_seas <- as.numeric(logLik(lm(df_m$Prevalence ~ offset(df_m$Non_Constant_Temp_Non_Aestivation))))
chisq <- -2 * (ll_const - ll_seas)
pchisq(chisq, df = 1, lower.tail = FALSE)

ll_const <- as.numeric(logLik(lm(Prevalence ~ Constant_Temp_Non_Aestivation, data = df_m)))
ll_seas <- as.numeric(logLik(lm(Prevalence ~ Non_Constant_Temp_Non_Aestivation, data = df_m)))
chisq <- -2 * (ll_const - ll_seas)
pchisq(chisq, df = 1, lower.tail = FALSE)


# S. Haematobium
t.test((df_h$Non_Constant_Temp_Non_Aestivation - df_h$Prevalence)^2, 
       (df_h$Constant_Temp_Non_Aestivation - df_h$Prevalence)^2,
       paired = TRUE)
#Maurice suggestion 

ll_const <- as.numeric(logLik(lm(df_h$Prevalence ~ offset(df_h$Constant_Temp_Non_Aestivation))))
ll_seas <- as.numeric(logLik(lm(df_h$Prevalence ~ offset(df_h$Non_Constant_Temp_Non_Aestivation))))
chisq <- -2 * (ll_const - ll_seas)
pchisq(chisq, df = 1, lower.tail = FALSE)



ll_const <- as.numeric(logLik(lm(Prevalence ~ Constant_Temp_Non_Aestivation, data = df_h)))
ll_seas <- as.numeric(logLik(lm(Prevalence ~ Non_Constant_Temp_Non_Aestivation, data = df_h)))
chisq <- -2 * (ll_const - ll_seas)
pchisq(chisq, df = 1, lower.tail = FALSE)


#### Are the error rates different with dormancy (for seasonal model)?

# S. Mansoni 
t.test((df_m$Non_Constant_Temp_Non_Aestivation - df_m$Prevalence)^2, 
       (df_m$Non_Constant_Temp_Aestivation - df_m$Prevalence)^2,
       paired = TRUE)

#Maurice suggestion 

ll_const <- as.numeric(logLik(lm(df_m$Prevalence ~ offset(df_m$Non_Constant_Temp_Non_Aestivation))))
ll_seas <- as.numeric(logLik(lm(df_m$Prevalence ~ offset(df_m$Non_Constant_Temp_Aestivation))))
chisq <- -2 * (ll_const - ll_seas)
pchisq(chisq, df = 3, lower.tail = FALSE)


ll_const <- as.numeric(logLik(lm(Prevalence ~ Non_Constant_Temp_Non_Aestivation, data = df_m)))
ll_seas <- as.numeric(logLik(lm(Prevalence ~ Non_Constant_Temp_Aestivation, data = df_m)))
chisq <- -2 * (ll_const - ll_seas)
pchisq(chisq, df = 3, lower.tail = FALSE)


# S. Haematobium
t.test((df_h$Non_Constant_Temp_Non_Aestivation - df_h$Prevalence)^2, 
       (df_h$Non_Constant_Temp_Aestivation - df_h$Prevalence)^2,
       paired = TRUE)

ll_const <- as.numeric(logLik(lm(df_h$Prevalence ~ offset(df_h$Non_Constant_Temp_Non_Aestivation))))
ll_seas <- as.numeric(logLik(lm(df_h$Prevalence ~ offset(df_h$Non_Constant_Temp_Aestivation))))
chisq <- -2 * (ll_const - ll_seas)
pchisq(chisq, df = 3, lower.tail = FALSE)


ll_const <- as.numeric(logLik(lm(Prevalence ~ Non_Constant_Temp_Non_Aestivation, data = df_h)))
ll_seas <- as.numeric(logLik(lm(Prevalence ~ Non_Constant_Temp_Aestivation, data = df_h)))
chisq <- -2 * (ll_const - ll_seas)
pchisq(chisq, df = 3, lower.tail = FALSE)





#### How many model predictions are false negative? (FN/(FN + TN))

## S. mansoni

# False negative rate for constant temp, no aestivation 
constant_nonA_FN <- sum((df_m$Prevalence > 0) & (df_m$Constant_Temp_Non_Aestivation == 0))/
  (sum((df_m$Prevalence == 0) & (df_m$Constant_Temp_Non_Aestivation == 0)) + sum((df_m$Prevalence > 0) & (df_m$Constant_Temp_Non_Aestivation == 0)))

# False negative rate for seasonal, no aestivation 
non_constant_nonA_FN <- sum((df_m$Prevalence > 0) & (df_m$Non_Constant_Temp_Non_Aestivation == 0))/
  (sum((df_m$Prevalence == 0) & (df_m$Non_Constant_Temp_Non_Aestivation == 0)) + sum((df_m$Prevalence > 0) & (df_m$Non_Constant_Temp_Non_Aestivation == 0)))

# False negative rate for seasonal, aestivation 
non_constant_A_FN <- sum((df_m$Prevalence > 0) & (df_m$Non_Constant_Temp_Aestivation == 0))/
  (sum((df_m$Prevalence == 0) & (df_m$Non_Constant_Temp_Aestivation == 0)) + sum((df_m$Prevalence > 0) & (df_m$Non_Constant_Temp_Aestivation == 0)))

## S. haematobium

# False negative rate for constant temp, no aestivation 
constant_nonA_FN_h <- sum((df_h$Prevalence > 0) & (df_h$Constant_Temp_Non_Aestivation == 0))/
  (sum((df_h$Prevalence == 0) & (df_h$Constant_Temp_Non_Aestivation == 0)) + sum((df_h$Prevalence > 0) & (df_h$Constant_Temp_Non_Aestivation == 0)))

# False negative rate for seasonal, no aestivation 
non_constant_nonA_FN_h <- sum((df_h$Prevalence > 0) & (df_h$Non_Constant_Temp_Non_Aestivation == 0))/
  (sum((df_h$Prevalence == 0) & (df_h$Non_Constant_Temp_Non_Aestivation == 0)) + sum((df_h$Prevalence > 0) & (df_h$Non_Constant_Temp_Non_Aestivation == 0)))

# False negative rate for seasonal, aestivation 
non_constant_A_FN_h <- sum((df_h$Prevalence > 0) & (df_h$Non_Constant_Temp_Aestivation == 0))/
  (sum((df_h$Prevalence == 0) & (df_h$Non_Constant_Temp_Aestivation == 0)) + sum((df_h$Prevalence > 0) & (df_h$Non_Constant_Temp_Aestivation == 0)))


#### Not used currently: Logistic model for mansoni

df_m$Prevalence_binary <- df_m$Prevalence > 0

mean(df_m$Prevalence_binary)

constant_glm <- glm(as.factor(Prevalence_binary) ~ Constant_Temp_Non_Aestivation, data = df_m, family = 'binomial')
summary(constant_glm)
lrm(Prevalence_binary ~ Constant_Temp_Non_Aestivation, data = df_m)

pROC_obj <- roc(as.factor(df_m$Prevalence_binary), df_m$Constant_Temp_Non_Aestivation, plot=TRUE)

seasonal_glm <- glm(Prevalence_binary ~ Non_Constant_Temp_Non_Aestivation, data = df_m, family = 'binomial')
summary(seasonal_glm)
lrm(Prevalence_binary ~ Non_Constant_Temp_Non_Aestivation, data = df_m)

pROC_obj <- roc(as.factor(df_m$Prevalence_binary), df_m$Non_Constant_Temp_Non_Aestivation, plot=TRUE)

seasonal_esti_glm <- glm(Prevalence_binary ~ Non_Constant_Temp_Aestivation, data = df_m, family = 'binomial')
summary(seasonal_esti_glm)
lrm(Prevalence_binary ~ Non_Constant_Temp_Aestivation, data = df_m)

pROC_obj <- roc(as.factor(df_m$Prevalence_binary), df_m$Non_Constant_Temp_Aestivation, plot=TRUE)


#### Not used currently: Logistic model for heamatobium 

df_h$Prevalence_binary <- df_h$Prevalence > 0

mean(df_h$Prevalence_binary)

constant_glm <- glm(as.factor(Prevalence_binary) ~ Constant_Temp_Non_Aestivation, data = df_h, family = 'binomial')
summary(constant_glm)
lrm(Prevalence_binary ~ Constant_Temp_Non_Aestivation, data = df_h)

pROC_obj <- roc(as.factor(df_h$Prevalence_binary), df_h$Constant_Temp_Non_Aestivation, plot=TRUE)

seasonal_glm <- glm(Prevalence_binary ~ Non_Constant_Temp_Non_Aestivation, data = df_h, family = 'binomial')
summary(seasonal_glm)
lrm(Prevalence_binary ~ Non_Constant_Temp_Non_Aestivation, data = df_h)

pROC_obj <- roc(as.factor(df_h$Prevalence_binary), df_h$Non_Constant_Temp_Non_Aestivation, plot=TRUE)

seasonal_esti_glm <- glm(Prevalence_binary ~ Non_Constant_Temp_Aestivation, data = df_h, family = 'binomial')
summary(seasonal_esti_glm)
lrm(Prevalence_binary ~ Non_Constant_Temp_Aestivation, data = df_h)

pROC_obj <- roc(as.factor(df_h$Prevalence_binary), df_h$Non_Constant_Temp_Aestivation, plot=TRUE)



