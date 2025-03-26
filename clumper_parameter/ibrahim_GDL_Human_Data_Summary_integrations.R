### Human Parasitology Data ###

library(tidyverse)
library(lme4)
library(dplyr)

# need file:
human_data <- read_csv("human_data_2018.csv")

human_data$sex <- as.factor(human_data$sex)
human_data$Net.time <- ifelse(human_data$Net.time==1, 'After','Before')
human_data$Intervention.time <- ifelse(human_data$Intervention.time==1, 'After','Before')

PreTreatedKidsID2017_MK_MN <- read_csv("2017_ID_pretreatedKids_MK_NM.csv")
human_data = left_join(human_data, PreTreatedKidsID2017_MK_MN, by=c("ID","year"))


# check if you want if we correctly ID the kids in MK and MN who were pre-treated 
# tmp2 <- subset(tempmat, PreTreatmentDate == "Dec-16")

# make tables 
# overall -- all kids
# by village
# by cluster
# by treatment  -- 3 groups 


##### summary stats for ALL kids #####
prev_data_all <- human_data %>%       ### arithemtic mean
  group_by(year) %>% 
  summarise_at(.vars = c("Sh","Sm"), funs(mean(., na.rm=T), n = sum(!is.na(.)), se = sd(., na.rm=T)/sqrt(sum(!is.na(.))))) %>% 
  mutate(Sh_CI = qnorm(0.975)*Sh_se,
         Sm_CI = qnorm(0.975)*Sm_se)

worm_data_all <- human_data %>%       ### geometric mean
  group_by(year) %>%
  summarize_at(.vars = c("ShW","SmW"), funs(gmean = exp(mean(log(.+1),na.rm=T)),
                                            n = sum(!is.na(.)),
                                            se = sd(exp(log(.+1)), na.rm = T)/sqrt(sum(!is.na(.))))) %>% 
  mutate(ShW_CI = qnorm(0.975)*ShW_se,
         SmW_CI = qnorm(0.975)*SmW_se)

dz_data_all <- left_join(prev_data_all, worm_data_all, by = "year")

other_vars_all <- human_data %>%      ### other important site and intervention variables
  group_by(year) %>% 
  summarize(pop = sum(Pop),
            n = sum(!is.na(sex)),
            nfemale = sum(sex=='F', na.rm = T),
            nTx = sum(!is.na(Tx_early)),
            Tx_early = sum(Tx_early==1, na.rm = T)) %>% 
  mutate(perc_male = 1-(nfemale/n),
         percTx_early = Tx_early/nTx)

human_summary_all <- left_join(dz_data_all, other_vars_all, by = "year") %>% 
  dplyr::select(-ShW_n, -SmW_n, -n, -nfemale, -nTx, -Tx_early) %>% 
  dplyr::select(year, pop, perc_male, percTx_early, Sh_mean, Sh_se, Sh_CI, Sh_n,
                ShW_gmean, ShW_se, ShW_CI, Sm_mean, Sm_se, Sm_CI, Sm_n, SmW_gmean, SmW_se, SmW_CI)


# uncheck it, if you want to write this file
# write.csv(human_summary, file = "human_dz_summary_all.csv")

##### group by village/school #####
prev_data <- human_data %>%       ### arithemtic mean
  group_by(School, year) %>% 
  summarise_at(.vars = c("Sh","Sm"), funs(mean(., na.rm=T), n = sum(!is.na(.)), se = sd(., na.rm=T)/sqrt(sum(!is.na(.))))) %>% 
  mutate(Sh_CI = qnorm(0.975)*Sh_se,
         Sm_CI = qnorm(0.975)*Sm_se)

library(MASS)

worm_data <- human_data %>%       ### geometric mean
  group_by(School, year) %>%
  summarize_at(.vars = c("ShW","SmW"), funs(Gmean = exp(mean(log(.+1),na.rm=T)),
            n = sum(!is.na(.)),
            se = sd(exp(log(.+1)), na.rm = T)/sqrt(sum(!is.na(.))))) %>% 
  mutate(ShW_CI = qnorm(0.975)*ShW_se,
         SmW_CI = qnorm(0.975)*SmW_se)

dz_data <- left_join(prev_data, worm_data, by = c("School" = "School", "year" = "year"))

#################################################################################
# BEGINNING of Giulio's latest integration on Feb 20 2023

# compute clumping parameters

#create a dataframe to store the results
cp_df <- worm_data[, c("School", "year")]
cp_df$Sh_Amean <- NA
cp_df$Sm_Amean <- NA
cp_df$Sh_cp    <- NA
cp_df$Sm_cp    <- NA

LL <- function(pars, data) {
  k <- pars[1]
  m <- pars[2]
  return(-sum(dnbinom(data,size=k,mu=m,log=TRUE)))
}


 for(i in 1:nrow(cp_df)){
    ndf <- filter(human_data, School == cp_df[i,]$School, year ==  cp_df[i,]$year)

    
    opt_sh <-optim(data = round(na.exclude(ndf$ShW)), par = c(1,1), LL, method="L-BFGS-B", lower=c(.000001,.000001))
    opt_sm <-optim(data = round(na.exclude(ndf$SmW)), par = c(1,1), LL, method="L-BFGS-B", lower=c(.000001,.000001))
    
    cp_df[i, 'Sh_Amean'] <- opt_sh$par[2]
    cp_df[i, 'Sm_Amean'] <- opt_sm$par[2]
    cp_df[i, 'Sh_cp'] <- opt_sh$par[1]  
    cp_df[i, 'Sm_cp'] <- opt_sm$par[1] 
    
    
  }

plot(cp_df$Sh_Amean, cp_df$Sh_cp, log = 'xy')

plot(cp_df$Sm_Amean, cp_df$Sm_cp, log = 'xy')

hist(cp_df$Sh_Amean)
hist(cp_df$Sm_Amean)
hist(cp_df$Sh_cp)
hist(cp_df$Sm_cp)

c(mean(cp_df$Sh_Amean), 
mean(cp_df$Sh_Amean[cp_df$year == 2016]), 
mean(cp_df$Sh_Amean[cp_df$year == 2017]), 
mean(cp_df$Sh_Amean[cp_df$year == 2018]))

c(median(cp_df$Sh_Amean), 
  median(cp_df$Sh_Amean[cp_df$year == 2016]), 
  median(cp_df$Sh_Amean[cp_df$year == 2017]), 
  median(cp_df$Sh_Amean[cp_df$year == 2018]))


c(mean(cp_df$Sm_Amean), 
mean(cp_df$Sm_Amean[cp_df$year == 2016]), 
mean(cp_df$Sm_Amean[cp_df$year == 2017]), 
mean(cp_df$Sm_Amean[cp_df$year == 2018]))

c(median(cp_df$Sm_Amean), 
  median(cp_df$Sm_Amean[cp_df$year == 2016]), 
  median(cp_df$Sm_Amean[cp_df$year == 2017]), 
  median(cp_df$Sm_Amean[cp_df$year == 2018]))


ggplot(cp_df, aes(x = Sh_Amean)) + 
  geom_histogram(binwidth=15, color="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sh_Amean)),color="blue", linetype="dashed", size=1) +
  xlab("mean parasite burden") + ylab("count") +
  ggtitle("S. haematobium")

ggplot(na.omit(cp_df, col = "Sm_cp"), aes(x = Sm_Amean)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sm_Amean)),color="blue", linetype="dashed", size=1) +
  xlab("mean parasite burden") + ylab("count") +
  ggtitle("S. mansoni") +
  scale_x_continuous(trans='log10')

  ggplot(cp_df, aes(x = Sh_cp)) + 
    geom_histogram(color="black", fill="white") +
    geom_vline(aes(xintercept=mean(Sh_cp)),color="blue", linetype="dashed", size=1) +
    xlab("clumping paramater") + ylab("count") +
    ggtitle("S. haematobium")
  
  ggplot(na.omit(cp_df, col = "Sm_cp"), aes(x = Sm_cp)) + 
    geom_histogram(color="black", fill="white") +
    geom_vline(aes(xintercept=mean(Sm_cp)),color="blue", linetype="dashed", size=1) +
    xlab("clumping parameter") + ylab("count") +
    ggtitle("S. mansoni")
  
  
  summary(cp_df$Sh_Amean)
  summary(cp_df$Sm_Amean)
  
# linear regression for the two parasites independently, not nested in years
# Here I just selected 2016 as the first year, as no experiment was conducted 
# back then. Note that 2017 is a reinfection year, so it is also interesting
  
regSh <- lm(Sh_cp ~ Sh_Amean, data = cp_df, year == 2016)
summary(regSh)

regSm <- lm(Sm_cp ~ Sm_Amean, data = cp_df, year == 2016)
summary(regSm)

# log-log regression

regSh <- lm(log(Sh_cp) ~ log(Sh_Amean), data = cp_df, year == 2016) #filter(cp_df, year == 2018))
summary(regSh)

regSm <- lm(log(Sm_cp) ~ log(Sm_Amean), data = cp_df, year == 2016)
summary(regSm)



ggplot(data = cp_df,aes(x = log(Sh_cp), y = log(Sh_Amean), color = year, group = 1)) +
  geom_point() +
  geom_smooth(method='lm') +
  xlab("log(parasites)") + ylab("log(clumping parameter 'k')") +
  ggtitle("S. haematobium")

ggplot(data = cp_df,aes(x = log(Sm_cp), y = log(Sm_Amean), color = year, group = 1)) +
  geom_point() +
  geom_smooth(method='lm') +
  xlab("log(parasites)") + ylab("log(clumping parameter 'k')") +
  ggtitle("S. mansoni")




# note that the results are almost identical
# lets us lmer4
library(Matrix)
library(lme4)

cp_df$year <- as.factor(cp_df$year)

regSh <- lmer(log(Sh_cp) ~ log(Sh_Amean) + (1|year), data = cp_df )
summary(regSh)

regSm <- lmer(log(Sm_cp) ~ log(Sm_Amean) + (1|year), data = cp_df )
summary(regSm)

################################

PairPrevShVar <- function(M) {
  k <- M^(summary(regSh)$coefficients[1])*exp(summary(regSh)$coefficients[2])
  p  <-  1 - 2* (1+M/(2*k))^(-k) + (1+M/k)^(-k)    # fraction of humans with at least 2 parasites 
  return(p)
}

PairPrevSmVar <- function(M) {
  k <- M^(summary(regSm)$coefficients[1])*exp(summary(regSm)$coefficients[2])
  p  <-  1 - 2* (1+M/(2*k))^(-k) + (1+M/k)^(-k)    # fraction of humans with at least 2 parasites 
  return(p)
}


# 2 or more worms
1 - 2* (1+cp_df$Sh_Amean[1]/(2 * cp_df$Sh_cp[1]))^(-cp_df$Sh_cp[1]) + (1 + cp_df$Sh_Amean[1]/cp_df$Sh_cp[1])^(-cp_df$Sh_cp[1])

#1 or ore worms
1 - (1/(1 + cp_df$Sh_Amean/cp_df$Sh_cp)^(cp_df$Sh_cp))*(1 + cp_df$Sh_Amean/(1 + cp_df$Sh_Amean/cp_df$Sh_cp)) # fraction of humans with at least 1 

# Estimate prevalence with the values of NBD values 

PairPrevNBD <- function(M, k) {
  p  <-  1 - 2* (1+M/(2*k))^(-k) + (1+M/k)^(-k)    # fraction of humans with at least 2 parasites 
  return(p)
}

WormPrevNBD <- function(M, k) {
  p  <-  1 - (1+M/k)^(-k)    # fraction of humans with at least 2 parasites 
  return(p)
}


cp_df$Sh_Amean
cp_df$Sm_Amean
cp_df$Sh_cp
cp_df$Sm_cp


options(scipen = 9)


#EstPrevSh <- PairPrevShVar(cp_df$Sh_Amean)
#EstPrevSm <- PairPrevSmVar(cp_df$Sm_Amean)


WormEstPrevSh <- WormPrevNBD(cp_df$Sh_Amean,cp_df$Sh_cp)
WormEstPrevSm <- WormPrevNBD(cp_df$Sm_Amean, cp_df$Sm_cp)


EstPrevNBDSh <- PairPrevNBD(cp_df$Sh_Amean,cp_df$Sh_cp)
EstPrevNBDSm <- PairPrevNBD(cp_df$Sm_Amean, cp_df$Sm_cp)

 

ObsPrevSh <- c()
ObsPrevSm <- c()

for(i in 1:nrow(cp_df)){
  ndf <- filter(human_data, School == cp_df[i,]$School, year ==  cp_df[i,]$year)
  ObsPrevSh[i] <- (length(na.exclude(ndf$ShW)) - length(which(na.exclude(ndf$ShW) == 0)))/length(na.exclude(ndf$ShW))
  ObsPrevSm[i] <- (length(na.exclude(ndf$SmW)) - length(which(na.exclude(ndf$SmW) == 0)))/length(na.exclude(ndf$SmW))
  }

ObsPrevSh/EstPrevNBDSh
ObsPrevSm/EstPrevNBDSm

ObsPrevSh/WormEstPrevSh
ObsPrevSm/WormEstPrevSm


##############################################################################
# END of Giulio latest integration on Feb 20 2023
##############################################################################

########################

other_vars <- human_data %>%      ### other important site and intervention variables
  group_by(School, year) %>% 
  summarize(Village = first(Village),
            pop = first(Pop),
            n = sum(!is.na(sex)),
            nfemale = sum(sex=='F', na.rm = T),
            lake = first(as.factor(lake)),
            cluster = first(as.factor(cluster)),
            net = first(as.factor(net)),
            prawn = first(as.factor(prawn)),
            veg_removal = first(as.factor(veg_removal)),
            net.time = first(Net.time),
            intervention.time = first(Intervention.time),
            nTx = sum(!is.na(Tx_early)),
            Tx_early = sum(Tx_early==1, na.rm = T)) %>% 
  mutate(perc_male = 1-(nfemale/n),
         percTx_early = Tx_early/nTx)
  
human_summary_village <- left_join(dz_data, other_vars, by = c("School" = "School", "year" = "year")) %>% 
  dplyr::select(-ShW_n, -SmW_n,-n,-nfemale,-nTx, -Tx_early) %>% 
  dplyr::select(Village, School, year, cluster, lake, pop, perc_male, percTx_early, net, net.time, prawn, veg_removal, intervention.time, Sh_mean, Sh_se, Sh_CI, Sh_n,
         ShW_gmean, ShW_se, ShW_CI,Sm_mean, Sm_se, Sm_CI, Sm_n,SmW_gmean, SmW_se, SmW_CI)

#write.csv(human_summary_village, file = "~/Documents/GitHub/Schisto/Net Effects/human_dz_summary.csv")

# now aggregate to cluster level 
prev_data_cluster <- human_data %>%       ### arithemtic mean
  group_by(cluster, year) %>% 
  summarise_at(.vars = c("Sh","Sm"), funs(mean(., na.rm=T), n = sum(!is.na(.)), se = sd(., na.rm=T)/sqrt(sum(!is.na(.))))) %>% 
  mutate(Sh_CI = qnorm(0.975)*Sh_se,
         Sm_CI = qnorm(0.975)*Sm_se)

worm_data_cluster <- human_data %>%       ### geometric mean
  group_by(cluster, year) %>%
  summarize_at(.vars = c("ShW","SmW"), funs(gmean = exp(mean(log(.+1),na.rm=T)),
                                            n = sum(!is.na(.)),
                                            se = sd(exp(log(.+1)), na.rm = T)/sqrt(sum(!is.na(.))))) %>% 
  mutate(ShW_CI = qnorm(0.975)*ShW_se,
         SmW_CI = qnorm(0.975)*SmW_se)
dz_data_cluster <- left_join(prev_data_cluster, worm_data_cluster, by = c("cluster" = "cluster", "year" = "year"))

other_vars_cluster <- human_data %>%      ### other important site and intervention variables
  group_by(cluster, year) %>% 
  summarize(pop = sum(Pop),
            n = sum(!is.na(sex)),
            nfemale = sum(sex=='F', na.rm = T),
            lake = first(as.factor(lake)),
            net = first(as.factor(net)),
            prawn = first(as.factor(prawn)),
            veg_removal = first(as.factor(veg_removal)),
            nTx = sum(!is.na(Tx_early)),
            Tx_early = sum(Tx_early==1, na.rm = T),
            net.time = first(Net.time),
            intervention.time = first(Intervention.time)) %>% 
  mutate(perc_male = 1-(nfemale/n),
         percTx_early = Tx_early/nTx)

human_summary_cluster <- left_join(dz_data_cluster, other_vars_cluster, by = c("cluster" = "cluster", "year" = "year")) %>% 
  dplyr::select(-ShW_n, -SmW_n,-n,-nfemale,-nTx, -Tx_early) %>% 
  dplyr::select(year, lake, pop, perc_male, percTx_early, net, net.time, prawn, veg_removal, intervention.time, Sh_mean, Sh_se, Sh_CI, Sh_n,
         ShW_gmean, ShW_se, ShW_CI,Sm_mean, Sm_se, Sm_CI, Sm_n,SmW_gmean, SmW_se, SmW_CI)

write.csv(human_summary, file = "G:/.shortcut-targets-by-id/1e2xUwinp8NuntTAkfNCUIW_48Onf6MCd/Belmont Project/Thermal_sensitive_models/r_ibrahim_scripts/mda_model/clumper_parameter")

# now aggregate by intervention  
prev_data_intervention <- human_data %>%       ### arithemtic mean
  group_by(Intervention.type, year) %>% 
  summarise_at(.vars = c("Sh","Sm"), funs(mean(., na.rm=T), n = sum(!is.na(.)), se = sd(., na.rm=T)/sqrt(sum(!is.na(.))))) %>% 
  mutate(Sh_CI = qnorm(0.975)*Sh_se,
         Sm_CI = qnorm(0.975)*Sm_se)

worm_data_intervention <- human_data %>%       ### geometric mean
  group_by(Intervention.type, year) %>%
  summarize_at(.vars = c("ShW","SmW"), funs(gmean = exp(mean(log(.+1),na.rm=T)),
                                            n = sum(!is.na(.)),
                                            se = sd(exp(log(.+1)), na.rm = T)/sqrt(sum(!is.na(.))))) %>% 
  mutate(ShW_CI = qnorm(0.975)*ShW_se,
         SmW_CI = qnorm(0.975)*SmW_se)

dz_data_intervention <- left_join(prev_data_intervention, worm_data_intervention, by = c("Intervention.type" = "Intervention.type", "year" = "year"))

other_vars_intervention <- human_data %>%      ### other important site and intervention variables
  group_by(Intervention.type, year) %>% 
  summarize(pop = sum(Pop),
            n = sum(!is.na(sex)),
            nfemale = sum(sex=='F', na.rm = T),
            lake = first(as.factor(lake)),
            net = first(as.factor(net)),
            prawn = first(as.factor(prawn)),
            veg_removal = first(as.factor(veg_removal)),
            net.time = first(Net.time),
            intervention.time = first(Intervention.time),
            nTx = sum(!is.na(Tx_early)),
            Tx_early = sum(Tx_early==1, na.rm = T)) %>% 
  mutate(perc_male = 1-(nfemale/n),
         percTx_early = Tx_early/nTx)

human_summary_intervention <- left_join(dz_data_intervention, other_vars_intervention, by = c("Intervention.type" = "Intervention.type", "year" = "year")) %>% 
  dplyr::select(-ShW_n, -SmW_n,-n,-nfemale,-nTx, -Tx_early) %>% 
  dplyr::select(year, lake, pop, perc_male, percTx_early, net, net.time, prawn, veg_removal, intervention.time, Sh_mean, Sh_se, Sh_CI, Sh_n,
         ShW_gmean, ShW_se, ShW_CI,Sm_mean, Sm_se, Sm_CI, Sm_n,SmW_gmean, SmW_se, SmW_CI)

write.csv(human_summary, file = "~/Documents/GitHub/Schisto/Net Effects/human_dz_summary_intervention.csv")

