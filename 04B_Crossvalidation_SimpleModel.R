# Cross validation: Simple model All Data
# Cross-validation

rm(list=ls()) # reproducibility
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
memory.limit(1e13)
library("nimble")



speciesCode <- "METTOX"
species <- "Metopium toxiferum"
method<- "simplemod_alldata"
load(paste0(path_to_data, "unit_cov2.RData"))
load(paste0(path_to_data, "florida.RData"))


load(paste0(path_to_data, species, "/surv_cov.RData"))
counties <- unique(florida@data$NAME_2)
counties <- factor(counties)

h_simple <- data.frame(county=counties)

for(i in 1:length(counties)) {
  temp <- subset(surv_cov, surv_cov$county == counties[i])
  h_simple$h[i] <- sum(temp$h)
  h_simple$nsurveys[i] <- nrow(temp)
}



code <- nimbleCode({

  ############################################################
  ####Level 1: Ecological Process

  #Estimate psi for each site
  for (i in 1:nsites) {

    logit(psi[i]) <- alpha0 + beta*X1[i] + gamma*X2[i] + delta*X3[i]
    z[i] ~ dbern(psi[i])
    p.temp[i] <- z[i]*p
    h[i] ~ dbinom(prob=p.temp[i], size = nsurveys[i])
  }

  #Derived Quantities for Cross-Validation
  for (j in 1:nsites_heldout) {
    logit(psi_heldout[Site_heldout[j]]) <- alpha0 + beta*X1[Site_heldout[j]] + gamma*X2[Site_heldout[j]] + delta*X3[Site_heldout[j]]
  }
  for (m in 1:nsites.visits_heldout) { #For each of the held-out datapoints
    ucpd[m]  <- p*psi_heldout[Site_heldout[m]]
  }

  #Priors
  alpha0 ~ dunif(-5, 5)
  beta ~ dunif(-5, 5)
  gamma ~ dunif(-5, 5)
  delta ~ dunif(-5, 5)
  p ~ dbeta(1,1)


})



#####Cross-validation part:

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


fold1 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==1))
eval1 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==1))
training1 <- subset(h_simple, !as.numeric(h_simple$county) %in% which(foldID==1))

fold2 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==2))
eval2 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==2))
training2 <- subset(h_simple, !as.numeric(h_simple$county) %in% which(foldID==2))

fold3 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==3))
eval3 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==4))
training3 <- subset(h_simple, !as.numeric(h_simple$county) %in% which(foldID==3))

fold4 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==4))
eval4 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==4))
training4 <- subset(h_simple, !as.numeric(h_simple$county) %in% which(foldID==4))

fold5 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==5))
eval5 <- subset(surv_cov, as.numeric(surv_cov$couny) %in% which(foldID==5))
training5 <- subset(h_simple, !as.numeric(h_simple$county) %in% which(foldID==5))


assembleData <- function(h_simple, heldout){

  # constants... never estimated


  constants <- list(
    nsites = length(unique(h_simple$county)),
    nsites_heldout = length(unique(heldout$county)),
    nsurveys = h_simple$nsurveys,
    nsites.visits_heldout = nrow(heldout),
    Site_heldout = as.numeric(heldout$county)
  )


  # data (used in front of "~" statements)
  data <- list(
    X1 = unit_cov2$pop2010,
    X2 = unit_cov2$tmin_min_avg,
    X3 = unit_cov2$area,
    h = h_simple$h
  )

  return(list("data"=data, "constants"=constants))
}

data_fold1 <- assembleData(training1, fold1)
data_fold2 <- assembleData(training2, fold2)
data_fold3 <- assembleData(training3, fold3)
data_fold4 <- assembleData(training4, fold4)
data_fold5 <- assembleData(training5, fold5)




initsFx <- function(df) {
  list(
    p = rbeta(1, 1, 1),
    z = rep(1, df$constants$nsites),
    beta = rnorm(1),
    gamma = rnorm(1),
    delta = rnorm(1),
    alpha0 = rnorm(1)
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
  inits = initsFx(data_fold1),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)

save(mcmc_fold1, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold1.RData"))

br_1 <- rep(NA, 99000)
for(i in 1:99000){
  br_1[i] <- brier(mcmc_fold1[i,], fold1$h)
}
br_1 <- br_1/length(fold1$h)
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
  br_2[i] <- brier(mcmc_fold2[i,], fold2$h)
}
br_2 <- br_2/length(fold2$h)
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
  br_3[i] <- brier(mcmc_fold3[i,], fold3$h)
}
br_3 <- br_3/length(fold3$h)
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
  br_4[i] <- brier(mcmc_fold4[i,], fold4$h)
}
br_4 <- br_4/length(fold4$h)
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
  br_5[i] <- brier(mcmc_fold5[i,], fold5$h)
}
br_5 <- br_5/length(fold5$h)
rm(mcmc_fold5)
gc()

br_simplemod_alldata <- mean(mean(br_1), mean(br_2), mean(br_3), mean(br_4), mean(br_5))
save(br_simplemod_alldata, br_1, br_2, br_3, br_4, br_5, file=paste0("G:/occupancy_anacards_data/", speciesCode,
                                                                     "/crossvalidation/", method, "_br.RData"))
timeEnd <- Sys.time() - timeStart









#Cross validation: Simple model – Filtered Data
# Cross-validation

rm(list=ls()) # reproducibility
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
memory.limit(1e13)
library("nimble")



speciesCode <- "METTOX"
species <- "Metopium toxiferum"
load(paste0(path_to_data, "unit_cov2.RData"))
load(paste0(path_to_data, "florida.RData"))

method <- "simplemod_filtered"

load(paste0(path_to_data, species, "/surv_cov3.RData"))


surv_cov <- subset(surv_cov3, surv_cov3$num_genera > 1 | surv_cov3$h==1)

rm(surv_cov3)
load(paste0(path_to_data, species, "/collectors_of_species.RData"))


counties <- unique(florida@data$NAME_2)
counties <- factor(counties)

h_simple <- data.frame(county=counties)

for(i in 1:length(counties)) {
  temp <- subset(surv_cov, surv_cov$county == counties[i])
  h_simple$h[i] <- sum(temp$h)
  h_simple$nsurveys[i] <- nrow(temp)
}



code <- nimbleCode({

  ############################################################
  ####Level 1: Ecological Process

  #Estimate psi for each site
  for (i in 1:nsites) {

    logit(psi[i]) <- alpha0 + beta*X1[i] + gamma*X2[i] + delta*X3[i]
    z[i] ~ dbern(psi[i])
    p.temp[i] <- z[i]*p
    h[i] ~ dbinom(prob=p.temp[i], size = nsurveys[i])
  }

  #Derived Quantities for Cross-Validation
  for (j in 1:nsites_heldout) {
    logit(psi_heldout[Site_heldout[j]]) <- alpha0 + beta*X1[Site_heldout[j]] + gamma*X2[Site_heldout[j]] + delta*X3[Site_heldout[j]]
  }
  for (m in 1:nsites.visits_heldout) { #For each of the held-out datapoints
    ucpd[m]  <- p*psi_heldout[Site_heldout[m]]
  }

  #Priors
  alpha0 ~ dunif(-5, 5)
  beta ~ dunif(-5, 5)
  gamma ~ dunif(-5, 5)
  delta ~ dunif(-5, 5)
  p ~ dbeta(1,1)


})



#####Cross-validation part:

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


fold1 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==1))
eval1 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==1))
training1 <- subset(h_simple, !as.numeric(h_simple$county) %in% which(foldID==1))

fold2 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==2))
eval2 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==2))
training2 <- subset(h_simple, !as.numeric(h_simple$county) %in% which(foldID==2))

fold3 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==3))
eval3 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==4))
training3 <- subset(h_simple, !as.numeric(h_simple$county) %in% which(foldID==3))

fold4 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==4))
eval4 <- subset(surv_cov, as.numeric(surv_cov$county) %in% which(foldID==4))
training4 <- subset(h_simple, !as.numeric(h_simple$county) %in% which(foldID==4))

fold5 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==5))
eval5 <- subset(surv_cov, as.numeric(surv_cov$couny) %in% which(foldID==5))
training5 <- subset(h_simple, !as.numeric(h_simple$county) %in% which(foldID==5))


assembleData <- function(h_simple, heldout){

  # constants... never estimated


  constants <- list(
    nsites = length(unique(h_simple$county)),
    nsites_heldout = length(unique(heldout$county)),
    nsurveys = h_simple$nsurveys,
    nsites.visits_heldout = nrow(heldout),
    Site_heldout = as.numeric(heldout$county)
  )


  # data (used in front of "~" statements)
  data <- list(
    X1 = unit_cov2$pop2010,
    X2 = unit_cov2$tmin_min_avg,
    X3 = unit_cov2$area,
    h = h_simple$h
  )

  return(list("data"=data, "constants"=constants))
}

data_fold1 <- assembleData(training1, fold1)
data_fold2 <- assembleData(training2, fold2)
data_fold3 <- assembleData(training3, fold3)
data_fold4 <- assembleData(training4, fold4)
data_fold5 <- assembleData(training5, fold5)




initsFx <- function(df) {
  list(
    p = rbeta(1, 1, 1),
    z = rep(1, df$constants$nsites),
    beta = rnorm(1),
    gamma = rnorm(1),
    delta = rnorm(1),
    alpha0 = rnorm(1)
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
  inits = initsFx(data_fold1),
  nchains = 1, niter = 100000, nburnin = 1000, thin = 1,
  summary = F, WAIC = F, samplesAsCodaMCMC = TRUE
)

save(mcmc_fold1, file=paste0("G:/occupancy_anacards_data/", speciesCode, "/crossvalidation/", method, "_mcmc_fold1.RData"))

br_1 <- rep(NA, 99000)
for(i in 1:99000){
  br_1[i] <- brier(mcmc_fold1[i,], fold1$h)
}
br_1 <- br_1/length(fold1$h)
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
  br_2[i] <- brier(mcmc_fold2[i,], fold2$h)
}
br_2 <- br_2/length(fold2$h)
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
  br_3[i] <- brier(mcmc_fold3[i,], fold3$h)
}
br_3 <- br_3/length(fold3$h)
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
  br_4[i] <- brier(mcmc_fold4[i,], fold4$h)
}
br_4 <- br_4/length(fold4$h)
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
  br_5[i] <- brier(mcmc_fold5[i,], fold5$h)
}
br_5 <- br_5/length(fold5$h)
rm(mcmc_fold5)
gc()

br_simplemod_filtered <- mean(mean(br_1), mean(br_2), mean(br_3), mean(br_4), mean(br_5))
save(br_simplemod_filtered, br_1, br_2, br_3, br_4, br_5, file=paste0("G:/occupancy_anacards_data/", speciesCode,
                                                                      "/crossvalidation/", method, "_br.RData"))
timeEnd <- Sys.time() - timeStart































