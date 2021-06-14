# Cross-validation: Full model All Data

rm(list=ls()) # reproducibility
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
memory.limit(1e13)
library("nimble")



speciesCode <- "METTOX"
species <- "Metopium toxiferum"
path_to_save <- paste0("E:/", speciesCode, "/") #This should be a location able to store large files

method <- "fullmod_alldata"
load(paste0(path_to_data, "unit_cov2.RData"))
load(paste0(path_to_data, species, "/surv_cov.RData"))
load(paste0(path_to_data, species, "/collectors_of_species.RData"))


for (i in 1:nrow(collectors_of_species)) {
  temp <- subset(surv_cov, surv_cov$collector == collectors_of_species$collector[i])
  collectors_of_species$surveys[i] <- nrow(temp)
  cat("Collector: ", i, " of ", nrow(collectors_of_species), "\n")
}

#The training data set needs to include at least one survey from each:
# * collector (762 unique collectors)
# * year      (124 unique years)
# * month     (12 unique months)

years <- unique(surv_cov$year)
months <- unique(surv_cov$month)
collectors <- unique(surv_cov$collector)

#Divide the 67 sites into five folds using code from Mike Meredith's site
# (http://www.mikemeredith.net/blog/2019/MSOM_CrossVal.htm)
nSites <- 67
nFolds <- 5
( foldSize <- nSites %/% nFolds + 1 )
tmp <- matrix(0, nFolds, foldSize)
for(i in 1:foldSize)
  tmp[, i] <- sample.int(nFolds)
foldID <- as.vector(tmp)[1:nSites]
table(foldID)
#####

#1) Divide the dataset into five groupings (randomly) GROUP 1:5
#2) For each group (e.g. Group 1):
#a) Initially set: Heldout_1 = Group_1,
#                    Training_1 = Group_2 + Group_3 + Group_4 + Group_5
#  b) Check whether all of the collectors that appear in Heldout_1 are in Training_1.
#    *If NOT: Randomly select one record from collector and move out of Heldout_1 and into Training_1
#  c) Check whether all of the years that appear in Heldout_1 are in Training_1
#    *If not: Randomly select one record from that year and move out of Heldout_1 and into Training_1
#  d) Check whether all of the months that appear in Heldout_1 are in Training_1
#    *If not: Randomly select one record from that month and move out of Heldout_1 and into Training_1


heldout1 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==1))
fold1 <- heldout1
training1 <- subset(surv_cov, !as.numeric(surv_cov$county) %in% which(foldID==1))

heldout2 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==2))
fold2 <- heldout2
training2 <- subset(surv_cov, !as.numeric(surv_cov$county) %in% which(foldID==2))

heldout3 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==3))
fold3 <- heldout3
training3 <- subset(surv_cov, !as.numeric(surv_cov$county) %in% which(foldID==3))

heldout4 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==4))
fold4 <- heldout4
training4 <- subset(surv_cov, !as.numeric(surv_cov$county) %in% which(foldID==4))

heldout5 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==5))
fold5 <- heldout5
training5 <- subset(surv_cov, !as.numeric(surv_cov$county) %in% which(foldID==5))




reassign <- function(fold, training) {

  #Check to make sure all collectors in the heldout data are accounted for in the training data
  missing_collectors <- subset(collectors, collectors %in% unique(fold$collector) &
                                 !(collectors %in% unique(training$collector)))
  if(length(missing_collectors)>0) {
    for(i in 1:length(missing_collectors)) {
      temp <- subset(fold, fold$collector == missing_collectors[i])
      rand <- sample(1:nrow(temp), 1)
      training <- rbind(training, temp[rand,])
      fold <- subset(fold, !(fold$index %in% training$index))
    }
  }

  if(sum(!(fold$year %in% training$year)) >0) { #(How many held out years are not in the training set?)
    missing_years <- subset(years, years %in% unique(fold$year) &
                              !(years %in% unique(training$year)))
    for(i in 1:length(missing_years)) {
      temp <- subset(fold, fold$year == missing_years[i])
      rand <- sample(1:nrow(temp), 1)
      training <- rbind(training, temp[rand,])
      fold <- subset(fold, !(fold$index %in% training$index))
    }
  }
  #Check for months
  if(sum(!(fold$month %in% training$month)) >0) { #(How many held out months are not in the training set?)
    missing_months <- subset(months, months %in% unique(fold$month) &
                               !(months %in% unique(training$month)))
    for(i in 1:length(missing_months)) {
      temp <- subset(fold, fold$month == missing_months[i])
      rand <- sample(1:nrow(temp), 1)
      training <- rbind(training, temp[rand,])
      fold <- subset(fold, !(fold$index %in% training$index))
    }
  }

  return(list("training"=training, "fold"=fold))
}

run1 <- reassign(fold1, training1)
run2 <- reassign(fold2, training2)
run3 <- reassign(fold3, training3)
run4 <- reassign(fold4, training4)
run5 <- reassign(fold5, training5)




#First fold
nrow(run1$fold)/(nrow(run1$fold) + nrow(run1$training))
nrow(run2$fold)/(nrow(run2$fold) + nrow(run2$training))
nrow(run3$fold)/(nrow(run3$fold) + nrow(run3$training))
nrow(run4$fold)/(nrow(run4$fold) + nrow(run4$training))
nrow(run5$fold)/(nrow(run5$fold) + nrow(run5$training))

save(run1, run2, run3, run4, run5,
     file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_folds.RData"))
load(paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_folds.RData"))
code <- nimbleCode({

  ############################################################
  ####Level 1: Ecological Process

  #Estimate psi for each site
  for (i in 1:nsites) {
    logit(psi[i]) <- alpha0 + beta*X1[i] + gamma*X2[i] + delta*X3[i]
    z[i] ~ dbern(psi[i])
  }

  #Estimate p for each observation
  for (o in 1:nsites.visits) {
    logit(p[o]) <- obs[IDobs[o]] + zeta*Y1[o] + eta*Y2[o] + mos[IDmos[o]] + yr[IDyr[o]]
    p.temp[o]  <- z[Site[o]] * p[o]
    h[o] ~ dbern(p.temp[o])
  }



  #Priors
  alpha0 ~ dunif(-5, 5)
  beta ~ dunif(-5, 5)
  gamma ~ dunif(-5, 5)
  delta ~ dunif(-5, 5)
  zeta ~ dunif(-5, 5)
  eta ~ dunif(-5, 5)

  #Random effect of collector
  for(j in 1:nobs) {
    obs[j] ~ dnorm(0, tau=tauobs) # added named argument "tau" for NIMBLE
  }

  #Random effect of month
  for(k in 1:nmos){
    mos[k] ~ dnorm(0, tau=taumos) # added named argument "tau" for NIMBLE
  }

  #Random effect of year
  for(l in 1:nyr){
    yr[l] ~ dnorm(0, tau=tauyr) # added named argument "tau" for NIMBLE
  }

  #Hyperpriors

  tauobs ~ dgamma(1, 0.001)
  taumos ~ dgamma(1, 0.001)
  tauyr  ~ dgamma(1, 0.001)

  #Derived Quantities for Cross-Validation
  for (m in 1:nsites_heldout) {
    logit(psi_heldout[Site2_heldout[m]]) <- alpha0 + beta*X1[Site2_heldout[m]] + gamma*X2[Site2_heldout[m]] +
      delta*X3[Site2_heldout[m]]
  }

  for (n in 1:nsites.visits_heldout) { #For each of the held-out datapoints
    logit(p_heldout[n]) <- obs[IDobs_heldout[n]] + zeta*Y1_heldout[n] + eta*Y2_heldout[n] +
      mos[IDmos_heldout[n]] + yr[IDyr_heldout[n]]
    ucpd[n]  <- p_heldout[n]*psi_heldout[Site_heldout[n]]
  }
})



assembleData2 <- function(training, heldout){
  training$year <- factor(training$year)
  training$collector <- factor(training$collector)
  training$month <- factor(training$month)
  heldout$year <- factor(heldout$year)
  heldout$collector <- factor(heldout$collector)
  heldout$month <- factor(heldout$month)
  #heldout$county <- factor(heldout$county)
  #training$county <- factor(heldout$county)



  # constants... never estimated


  constants <- list(
    nsites = 67,
    nsites_heldout = length(unique(heldout$county)),
    nsites.visits = nrow(training),
    nsites.visits_heldout = nrow(heldout),
    nobs = length(unique(training$collector)),
    nmos = length(unique(training$month)),
    nyr = length(unique(training$year)),
    IDobs = as.numeric(training$collector), # convert to numeric factor levels
    IDobs_heldout = as.numeric(heldout$collector),
    IDcounty = as.numeric(training$county), # convert to numeric factor levels
    IDmos = as.numeric(training$month), # convert to numeric factor levels
    IDmos_heldout = as.numeric(heldout$month),
    IDyr = as.numeric(training$year), # convert to numeric factor levels
    IDyr_heldout = as.numeric(heldout$year),
    Site = as.numeric(training$county), # convert to numeric factor levels
    Site_heldout = as.numeric(heldout$county),
    Site2_heldout = unique(as.numeric(heldout$county))
  )


  # data (used in front of "~" statements)
  data <- list(
    X1 = unit_cov2$pop2010,
    X2 = unit_cov2$tmin_min_avg,
    X3 = unit_cov2$area,
    Y1 = training$prev_detected,
    Y1_heldout = heldout$prev_detected,
    Y2 = as.numeric(training$citizen_scientist), # convert to numeric factor levels,
    Y2_heldout = as.numeric(heldout$citizen_scientist),
    h = training$h
  )

  return(list("data"=data, "constants"=constants))
}

###NEW




data_fold1 <- assembleData2(run1$training, run1$fold)
data_fold2 <- assembleData2(run2$training, run2$fold)
data_fold3 <- assembleData2(run3$training, run3$fold)
data_fold4 <- assembleData2(run4$training, run4$fold)
data_fold5 <- assembleData2(run5$training, run5$fold)


monitors <- c("ucpd")

initsFx <- function(df) {

  list(

    z = rep(1, df$constants$nsites),
    beta = rnorm(1),
    gamma = rnorm(1),
    delta = rnorm(1),
    zeta = rnorm(1),
    eta = rnorm(1),
    alpha0 = rnorm(1),
    obs = rnorm(df$constants$nobs),
    mos = rnorm(df$constants$nmos),
    yr = rnorm(df$constants$nyr),
    tauobs = rlnorm(1),
    taumos = rlnorm(1),
    tauyr = rlnorm(1)
  )

}

brier <- function(probs, bool) {
  sum(bool * (1 - probs)^2 + (1 - bool) * probs^2, na.rm=TRUE)
}

timeStart <- Sys.time()
#fold1
mcmc_fold1<- nimbleMCMC(
  code = code, constants = data_fold1$constants,
  data = data_fold1$data,
  monitors = c("ucpd"),
  inits = initsFx(data_fold4),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)
save(mcmc_fold1, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold1.RData"))

br_1 <- rep(NA, 99000)
for(i in 1:99000){
  br_1[i] <- brier(mcmc_fold1[i,], run1$fold$h)
}
br_1 <- br_1/length(run1$fold$h)
rm(mcmc_fold1)
gc()

#fold2
mcmc_fold2<- nimbleMCMC(
  code = code, constants = data_fold2$constants,
  data = data_fold2$data,
  monitors = c("ucpd"),
  inits = initsFx(data_fold2),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)
save(mcmc_fold2, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold2.RData"))

br_2 <- rep(NA, 99000)
for(i in 1:99000){
  br_2[i] <- brier(mcmc_fold2[i,], run2$fold$h)
}
br_2 <- br_2/length(run2$fold$h)
rm(mcmc_fold2)
gc()

#fold3
mcmc_fold3<- nimbleMCMC(
  code = code, constants = data_fold3$constants,
  data = data_fold3$data,
  monitors = c("ucpd"),
  inits = initsFx(data_fold3),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)
save(mcmc_fold3, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold3.RData"))

br_3 <- rep(NA, 99000)
for(i in 1:99000){
  br_3[i] <- brier(mcmc_fold3[i,], run3$fold$h)
}
br_3 <- br_3/length(run3$fold$h)
rm(mcmc_fold3)
gc()

#fold4
mcmc_fold4<- nimbleMCMC(
  code = code, constants = data_fold4$constants,
  data = data_fold4$data,
  monitors = c("ucpd"),
  inits = initsFx(data_fold4),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)
save(mcmc_fold4, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold4.RData"))

br_4 <- rep(NA, 99000)
for(i in 1:99000){
  br_4[i] <- brier(mcmc_fold4[i,], run4$fold$h)
}
br_4 <- br_4/length(run4$fold$h)
rm(mcmc_fold4)
gc()

#fold5
mcmc_fold5<- nimbleMCMC(
  code = code, constants = data_fold5$constants,
  data = data_fold5$data,
  monitors = c("ucpd"),
  inits = initsFx(data_fold5),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)
save(mcmc_fold5, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold5.RData"))

br_5 <- rep(NA, 99000)
for(i in 1:99000){
  br_5[i] <- brier(mcmc_fold5[i,], run5$fold$h)
}
br_5 <- br_5/length(run5$fold$h)
rm(mcmc_fold5)
gc()

br_fullmod_alldata <- mean(mean(br_1), mean(br_2), mean(br_3), mean(br_4), mean(br_5))
save(br_fullmod_alldata, br_1, br_2, br_3, br_4, br_5, file=paste0("G:/occupancy_anacards_data/",
                                                                   speciesCode, "/crossvalidation/", method,
                                                                   "_br.RData"))
timeEnd <- Sys.time()-timeStart




# Cross-validation Full Model Filtered
# Cross-validation

rm(list=ls()) # reproducibility
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
memory.limit(1e13)
library("nimble")




speciesCode <- "METTOX"
species <- "Metopium toxiferum"
workDir <-
  path_to_save <- paste0("E:/", speciesCode, "/") #This should be a location able to store large files

method <- "fullmod_filtered"
load(paste0(path_to_data, "unit_cov2.RData"))
load(paste0(path_to_data, species, "/surv_cov3.RData"))


surv_cov <- subset(surv_cov3, surv_cov3$num_genera > 1 | surv_cov3$h==1)

rm(surv_cov3)
load(paste0(path_to_data, species, "/collectors_of_species.RData"))


for (i in 1:nrow(collectors_of_species)) {
  temp <- subset(surv_cov, surv_cov$collector == collectors_of_species$collector[i])
  collectors_of_species$surveys[i] <- nrow(temp)
  cat("Collector: ", i, " of ", nrow(collectors_of_species), "\n")
}

#The training data set needs to include at least one survey from each:
# * collector (762 unique collectors)
# * year      (124 unique years)
# * month     (12 unique months)

years <- unique(surv_cov$year)
months <- unique(surv_cov$month)
collectors <- unique(surv_cov$collector)

#Divide the 67 sites into five folds using code from Mike Meredith's site
# (http://www.mikemeredith.net/blog/2019/MSOM_CrossVal.htm)
nSites <- 67
nFolds <- 5
( foldSize <- nSites %/% nFolds + 1 )
tmp <- matrix(0, nFolds, foldSize)
for(i in 1:foldSize)
  tmp[, i] <- sample.int(nFolds)
foldID <- as.vector(tmp)[1:nSites]
table(foldID)
#####

#1) Divide the dataset into five groupings (randomly) GROUP 1:5
#2) For each group (e.g. Group 1):
#a) Initially set: Heldout_1 = Group_1,
#                    Training_1 = Group_2 + Group_3 + Group_4 + Group_5
#  b) Check whether all of the collectors that appear in Heldout_1 are in Training_1.
#    *If NOT: Randomly select one record from collector and move out of Heldout_1 and into Training_1
#  c) Check whether all of the years that appear in Heldout_1 are in Training_1
#    *If not: Randomly select one record from that year and move out of Heldout_1 and into Training_1
#  d) Check whether all of the months that appear in Heldout_1 are in Training_1
#    *If not: Randomly select one record from that month and move out of Heldout_1 and into Training_1


heldout1 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==1))
fold1 <- heldout1
training1 <- subset(surv_cov, !as.numeric(surv_cov$county) %in% which(foldID==1))

heldout2 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==2))
fold2 <- heldout2
training2 <- subset(surv_cov, !as.numeric(surv_cov$county) %in% which(foldID==2))

heldout3 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==3))
fold3 <- heldout3
training3 <- subset(surv_cov, !as.numeric(surv_cov$county) %in% which(foldID==3))

heldout4 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==4))
fold4 <- heldout4
training4 <- subset(surv_cov, !as.numeric(surv_cov$county) %in% which(foldID==4))

heldout5 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==5))
fold5 <- heldout5
training5 <- subset(surv_cov, !as.numeric(surv_cov$county) %in% which(foldID==5))




reassign <- function(fold, training) {

  #Check to make sure all collectors in the heldout data are accounted for in the training data
  missing_collectors <- subset(collectors, collectors %in% unique(fold$collector) &
                                 !(collectors %in% unique(training$collector)))
  if(length(missing_collectors)>0) {
    for(i in 1:length(missing_collectors)) {
      temp <- subset(fold, fold$collector == missing_collectors[i])
      rand <- sample(1:nrow(temp), 1)
      training <- rbind(training, temp[rand,])
      fold <- subset(fold, !(fold$index %in% training$index))
    }
  }

  if(sum(!(fold$year %in% training$year)) >0) { #(How many held out years are not in the training set?)
    missing_years <- subset(years, years %in% unique(fold$year) &
                              !(years %in% unique(training$year)))
    for(i in 1:length(missing_years)) {
      temp <- subset(fold, fold$year == missing_years[i])
      rand <- sample(1:nrow(temp), 1)
      training <- rbind(training, temp[rand,])
      fold <- subset(fold, !(fold$index %in% training$index))
    }
  }
  #Check for months
  if(sum(!(fold$month %in% training$month)) >0) { #(How many held out months are not in the training set?)
    missing_months <- subset(months, months %in% unique(fold$month) &
                               !(months %in% unique(training$month)))
    for(i in 1:length(missing_months)) {
      temp <- subset(fold, fold$month == missing_months[i])
      rand <- sample(1:nrow(temp), 1)
      training <- rbind(training, temp[rand,])
      fold <- subset(fold, !(fold$index %in% training$index))
    }
  }

  return(list("training"=training, "fold"=fold))
}

run1 <- reassign(fold1, training1)
run2 <- reassign(fold2, training2)
run3 <- reassign(fold3, training3)
run4 <- reassign(fold4, training4)
run5 <- reassign(fold5, training5)




#First fold
nrow(run1$fold)/(nrow(run1$fold) + nrow(run1$training))
nrow(run2$fold)/(nrow(run2$fold) + nrow(run2$training))
nrow(run3$fold)/(nrow(run3$fold) + nrow(run3$training))
nrow(run4$fold)/(nrow(run4$fold) + nrow(run4$training))
nrow(run5$fold)/(nrow(run5$fold) + nrow(run5$training))

save(run1, run2, run3, run4, run5,
     file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_folds.RData"))
load(paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_folds.RData"))
code <- nimbleCode({

  ############################################################
  ####Level 1: Ecological Process

  #Estimate psi for each site
  for (i in 1:nsites) {
    logit(psi[i]) <- alpha0 + beta*X1[i] + gamma*X2[i] + delta*X3[i]
    z[i] ~ dbern(psi[i])
  }

  #Estimate p for each observation
  for (o in 1:nsites.visits) {
    logit(p[o]) <- obs[IDobs[o]] + zeta*Y1[o] + eta*Y2[o] + mos[IDmos[o]] + yr[IDyr[o]]
    p.temp[o]  <- z[Site[o]] * p[o]
    h[o] ~ dbern(p.temp[o])
  }



  #Priors
  alpha0 ~ dunif(-5, 5)
  beta ~ dunif(-5, 5)
  gamma ~ dunif(-5, 5)
  delta ~ dunif(-5, 5)
  zeta ~ dunif(-5, 5)
  eta ~ dunif(-5, 5)

  #Random effect of collector
  for(j in 1:nobs) {
    obs[j] ~ dnorm(0, tau=tauobs) # added named argument "tau" for NIMBLE
  }

  #Random effect of month
  for(k in 1:nmos){
    mos[k] ~ dnorm(0, tau=taumos) # added named argument "tau" for NIMBLE
  }

  #Random effect of year
  for(l in 1:nyr){
    yr[l] ~ dnorm(0, tau=tauyr) # added named argument "tau" for NIMBLE
  }

  #Hyperpriors

  tauobs ~ dgamma(1, 0.001)
  taumos ~ dgamma(1, 0.001)
  tauyr  ~ dgamma(1, 0.001)

  #Derived Quantities for Cross-Validation
  for (m in 1:nsites_heldout) {
    logit(psi_heldout[Site2_heldout[m]]) <- alpha0 + beta*X1[Site2_heldout[m]] + gamma*X2[Site2_heldout[m]] +
      delta*X3[Site2_heldout[m]]
  }

  for (n in 1:nsites.visits_heldout) { #For each of the held-out datapoints
    logit(p_heldout[n]) <- obs[IDobs_heldout[n]] + zeta*Y1_heldout[n] + eta*Y2_heldout[n] +
      mos[IDmos_heldout[n]] + yr[IDyr_heldout[n]]
    ucpd[n]  <- p_heldout[n]*psi_heldout[Site_heldout[n]]
  }
})



assembleData2 <- function(training, heldout){
  training$year <- factor(training$year)
  training$collector <- factor(training$collector)
  training$month <- factor(training$month)
  heldout$year <- factor(heldout$year)
  heldout$collector <- factor(heldout$collector)
  heldout$month <- factor(heldout$month)
  #heldout$county <- factor(heldout$county)
  #training$county <- factor(heldout$county)



  # constants... never estimated


  constants <- list(
    nsites = 67,
    nsites_heldout = length(unique(heldout$county)),
    nsites.visits = nrow(training),
    nsites.visits_heldout = nrow(heldout),
    nobs = length(unique(training$collector)),
    nmos = length(unique(training$month)),
    nyr = length(unique(training$year)),
    IDobs = as.numeric(training$collector), # convert to numeric factor levels
    IDobs_heldout = as.numeric(heldout$collector),
    IDcounty = as.numeric(training$county), # convert to numeric factor levels
    IDmos = as.numeric(training$month), # convert to numeric factor levels
    IDmos_heldout = as.numeric(heldout$month),
    IDyr = as.numeric(training$year), # convert to numeric factor levels
    IDyr_heldout = as.numeric(heldout$year),
    Site = as.numeric(training$county), # convert to numeric factor levels
    Site_heldout = as.numeric(heldout$county),
    Site2_heldout = unique(as.numeric(heldout$county))
  )


  # data (used in front of "~" statements)
  data <- list(
    X1 = unit_cov2$pop2010,
    X2 = unit_cov2$tmin_min_avg,
    X3 = unit_cov2$area,
    Y1 = training$prev_detected,
    Y1_heldout = heldout$prev_detected,
    Y2 = as.numeric(training$citizen_scientist), # convert to numeric factor levels,
    Y2_heldout = as.numeric(heldout$citizen_scientist),
    h = training$h
  )

  return(list("data"=data, "constants"=constants))
}

###NEW




data_fold1 <- assembleData2(run1$training, run1$fold)
data_fold2 <- assembleData2(run2$training, run2$fold)
data_fold3 <- assembleData2(run3$training, run3$fold)
data_fold4 <- assembleData2(run4$training, run4$fold)
data_fold5 <- assembleData2(run5$training, run5$fold)


monitors <- c("ucpd")

initsFx <- function(df) {

  list(

    z = rep(1, df$constants$nsites),
    beta = rnorm(1),
    gamma = rnorm(1),
    delta = rnorm(1),
    zeta = rnorm(1),
    eta = rnorm(1),
    alpha0 = rnorm(1),
    obs = rnorm(df$constants$nobs),
    mos = rnorm(df$constants$nmos),
    yr = rnorm(df$constants$nyr),
    tauobs = rlnorm(1),
    taumos = rlnorm(1),
    tauyr = rlnorm(1)
  )

}

brier <- function(probs, bool) {
  sum(bool * (1 - probs)^2 + (1 - bool) * probs^2, na.rm=TRUE)
}

timeStart <- Sys.time()
#fold1
mcmc_fold1<- nimbleMCMC(
  code = code, constants = data_fold1$constants,
  data = data_fold1$data,
  monitors = c("ucpd"),
  inits = initsFx(data_fold4),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)
save(mcmc_fold1, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold1.RData"))

br_1 <- rep(NA, 99000)
for(i in 1:99000){
  br_1[i] <- brier(mcmc_fold1[i,], run1$fold$h)
}
br_1 <- br_1/length(run1$fold$h)
rm(mcmc_fold1)
gc()

#fold2
mcmc_fold2<- nimbleMCMC(
  code = code, constants = data_fold2$constants,
  data = data_fold2$data,
  monitors = c("ucpd"),
  inits = initsFx(data_fold2),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)
save(mcmc_fold2, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold2.RData"))

br_2 <- rep(NA, 99000)
for(i in 1:99000){
  br_2[i] <- brier(mcmc_fold2[i,], run2$fold$h)
}
br_2 <- br_2/length(run2$fold$h)
rm(mcmc_fold2)
gc()

#fold3
mcmc_fold3<- nimbleMCMC(
  code = code, constants = data_fold3$constants,
  data = data_fold3$data,
  monitors = c("ucpd"),
  inits = initsFx(data_fold3),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)
save(mcmc_fold3, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold3.RData"))

br_3 <- rep(NA, 99000)
for(i in 1:99000){
  br_3[i] <- brier(mcmc_fold3[i,], run3$fold$h)
}
br_3 <- br_3/length(run3$fold$h)
rm(mcmc_fold3)
gc()

#fold4
mcmc_fold4<- nimbleMCMC(
  code = code, constants = data_fold4$constants,
  data = data_fold4$data,
  monitors = c("ucpd"),
  inits = initsFx(data_fold4),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)
save(mcmc_fold4, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold4.RData"))

br_4 <- rep(NA, 99000)
for(i in 1:99000){
  br_4[i] <- brier(mcmc_fold4[i,], run4$fold$h)
}
br_4 <- br_4/length(run4$fold$h)
rm(mcmc_fold4)
gc()

#fold5
mcmc_fold5<- nimbleMCMC(
  code = code, constants = data_fold5$constants,
  data = data_fold5$data,
  monitors = c("ucpd"),
  inits = initsFx(data_fold5),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)
save(mcmc_fold5, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold5.RData"))

br_5 <- rep(NA, 99000)
for(i in 1:99000){
  br_5[i] <- brier(mcmc_fold5[i,], run5$fold$h)
}
br_5 <- br_5/length(run5$fold$h)
rm(mcmc_fold5)
gc()

br_fullmod_filtered <- mean(mean(br_1), mean(br_2), mean(br_3), mean(br_4), mean(br_5))
save(br_fullmod_filtered, br_1, br_2, br_3, br_4, br_5, file=paste0("G:/occupancy_anacards_data/",
                                                                    speciesCode, "/crossvalidation/", method,
                                                                    "_br.RData"))
timeEnd <- Sys.time()-timeStart






