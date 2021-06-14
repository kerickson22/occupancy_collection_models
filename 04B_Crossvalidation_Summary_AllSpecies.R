#Cross-validation Results - All Species
# NEW CROSS-VALIDATION FORMULAS

rm(list=ls())
memory.limit(1e13)


brier <- function(probs, bool) {
  sum(bool * (1-probs)^2 + (1-bool)*probs^2, na.rm=T)
}

brier_scaled <- function(probs, bool, scale1, scale2) {
  sum((bool*(1-probs)^2)/scale1 + ((1-bool)*probs^2)/scale2)
}

log_lik <- function(ucpd, data) {
  sum(data*log(ucpd) + (1-data)*log(1-ucpd), na.rm=T)
}

log_lik_scaled <- function(ucpd, data, scale1, scale2) {
  sum((data*log(ucpd))/scale1 + ((1-data)*log(1-ucpd))/scale2, na.rm=T)
}




pos <- c(1, 6, 11, 16)
speciesList <- c("MANIND", "METTOX", "RHUCOP", "SCHTER", "TOXPUB", "TOXRAD", "TOXVER")
methods <- c("fullmod_alldata", "fullmod_filtered", "simplemod_alldata", "simplemod_filtered")

#Determined that something was wrong with the simple models. So now I am loading the previously
# completed fullmodel results, and then appending the simple model results to them.
path_to_data <- "E:/"
for(j in 1:1) {

  speciesCode <- speciesList[j]
  # results <- data.frame(species=rep(speciesCode, 20),
  #                     method=rep(methods, each=5),
  #                     fold = rep(1:5, 4),
  #                     br_scale1 = rep(NA, 20),
  #                     br_scale2 = rep(NA, 20),
  #                     log_lik_scale1 = rep(NA, 20),
  #                     log_lik_scale2 = rep(NA, 20))
  load(paste0(path_to_data2, speciesCode, "/crossvalidation/",
              speciesCode, "_crossvalidation_results.RData"))
  if(j>4){path_to_data <- path_to_data2}
  if(j<5) {path_to_data <- "E:/"}
  for(m in 1:2) {
    method <- methods[m]

    #Run1
    load(paste0(path_to_data, speciesCode, "/crossvalidation/", method, "_mcmc_fold1.RData"))
    load(paste0(path_to_data, speciesCode, "/crossvalidation/", method, "_folds.RData"))


    br_1 <- rep(NA, 99000)
    for(i in 1:99000){
      br_1[i] <- brier(mcmc_fold1[i,], run1$fold$h)
    }
    results[pos[m],4 ] <-mean(br_1)

    br_scale1 <- rep(NA, 99000)
    for(i in 1:99000){
      br_scale1 <- brier_scaled(mcmc_fold1[i,], run1$fold$h,
                                sum(run1$fold$h), (length(run1$fold$h)-sum(run1$fold$h)))
    }
    results[pos[m], 5] <- mean(br_scale1)

    loglik_1 <- rep(NA, 99000)
    for(i in 1:99000){
      loglik_1 <- log_lik(mcmc_fold1[i,], run1$fold$h)
    }
    results[pos[m], 6] <- mean(loglik_1)

    loglik_scaled_1 <- rep(NA, 99000)
    for(i in 1:99000){
      loglik_scaled_1 <- log_lik_scaled(mcmc_fold1[i,], run1$fold$h,
                                        sum(run1$fold$h), (length(run1$fold$h)-sum(run1$fold$h)))
    }
    results[pos[m], 7] <- mean(loglik_scaled_1)
    rm(mcmc_fold1)
    cat(1)
    #Run2
    load(paste0(path_to_data, speciesCode, "/crossvalidation/", method, "_mcmc_fold2.RData"))

    br_1 <- rep(NA, 99000)
    for(i in 1:99000){
      br_1[i] <- brier(mcmc_fold2[i,], run2$fold$h)
    }
    results[pos[m]+1,4 ] <-mean(br_1)

    br_scale1 <- rep(NA, 99000)
    for(i in 1:99000){
      br_scale1 <- brier_scaled(mcmc_fold2[i,], run2$fold$h,
                                sum(run2$fold$h), (length(run2$fold$h)-sum(run2$fold$h)))
    }
    results[pos[m]+1, 5] <- mean(br_scale1)

    loglik_1 <- rep(NA, 99000)
    for(i in 1:99000){
      loglik_1 <- log_lik(mcmc_fold2[i,], run2$fold$h)
    }
    results[pos[m]+1, 6] <- mean(loglik_1)

    loglik_scaled_1 <- rep(NA, 99000)
    for(i in 1:99000){
      loglik_scaled_1 <- log_lik_scaled(mcmc_fold2[i,], run2$fold$h,
                                        sum(run2$fold$h), (length(run2$fold$h)-sum(run2$fold$h)))
    }
    results[pos[m]+1, 7] <- mean(loglik_scaled_1)
    rm(mcmc_fold2)
    cat(2)
    #Run3
    load(paste0(path_to_data, speciesCode, "/crossvalidation/", method, "_mcmc_fold3.RData"))

    br_1 <- rep(NA, 99000)
    for(i in 1:99000){
      br_1[i] <- brier(mcmc_fold3[i,], run3$fold$h)
    }
    results[pos[m] +2,4 ] <-mean(br_1)

    br_scale1 <- rep(NA, 99000)
    for(i in 1:99000){
      br_scale1 <- brier_scaled(mcmc_fold3[i,], run3$fold$h,
                                sum(run3$fold$h), (length(run3$fold$h)-sum(run3$fold$h)))
    }
    results[pos[m] +2, 5] <- mean(br_scale1)

    loglik_1 <- rep(NA, 99000)
    for(i in 1:99000){
      loglik_1 <- log_lik(mcmc_fold3[i,], run3$fold$h)
    }
    results[pos[m] +2, 6] <- mean(loglik_1)

    loglik_scaled_1 <- rep(NA, 99000)
    for(i in 1:99000){
      loglik_scaled_1 <- log_lik_scaled(mcmc_fold3[i,], run3$fold$h,
                                        sum(run3$fold$h), (length(run3$fold$h)-sum(run3$fold$h)))
    }
    results[pos[m]+2, 7] <- mean(loglik_scaled_1)
    rm(mcmc_fold3)
    cat(3)
    #Run4
    load(paste0(path_to_data, speciesCode, "/crossvalidation/", method, "_mcmc_fold4.RData"))

    br_1 <- rep(NA, 99000)
    for(i in 1:99000){
      br_1[i] <- brier(mcmc_fold4[i,], run4$fold$h)
    }
    results[pos[m]+3,4 ] <-mean(br_1)

    br_scale1 <- rep(NA, 99000)
    for(i in 1:99000){
      br_scale1 <- brier_scaled(mcmc_fold4[i,], run4$fold$h,
                                sum(run4$fold$h), (length(run4$fold$h)-sum(run4$fold$h)))
    }
    results[pos[m]+3, 5] <- mean(br_scale1)

    loglik_1 <- rep(NA, 99000)
    for(i in 1:99000){
      loglik_1 <- log_lik(mcmc_fold4[i,], run4$fold$h)
    }
    results[pos[m]+3, 6] <- mean(loglik_1)

    loglik_scaled_1 <- rep(NA, 99000)
    for(i in 1:99000){
      loglik_scaled_1 <- log_lik_scaled(mcmc_fold4[i,], run4$fold$h,
                                        sum(run4$fold$h), (length(run4$fold$h)-sum(run4$fold$h)))
    }
    results[pos[m]+3, 7] <- mean(loglik_scaled_1)
    rm(mcmc_fold4)
    cat(4)
    #Run5
    load(paste0(path_to_data, speciesCode, "/crossvalidation/", method, "_mcmc_fold5.RData"))

    br_1 <- rep(NA, 99000)
    for(i in 1:99000){
      br_1[i] <- brier(mcmc_fold5[i,], run5$fold$h)
    }
    results[pos[m]+4,4 ] <-mean(br_1)

    br_scale1 <- rep(NA, 99000)
    for(i in 1:99000){
      br_scale1 <- brier_scaled(mcmc_fold5[i,], run5$fold$h,
                                sum(run5$fold$h), (length(run5$fold$h)-sum(run5$fold$h)))
    }
    results[pos[m]+4, 5] <- mean(br_scale1)

    loglik_1 <- rep(NA, 99000)
    for(i in 1:99000){
      loglik_1 <- log_lik(mcmc_fold5[i,], run5$fold$h)
    }
    results[pos[m]+4, 6] <- mean(loglik_1)

    loglik_scaled_1 <- rep(NA, 99000)
    for(i in 1:99000){
      loglik_scaled_1 <- log_lik_scaled(mcmc_fold5[i,], run5$fold$h,
                                        sum(run5$fold$h), (length(run5$fold$h)-sum(run5$fold$h)))
    }
    results[pos[m]+4, 7] <- mean(loglik_scaled_1)
    rm(mcmc_fold5)
    cat(5)
    cat("Finished with: ", speciesCode, method, "\n")
  }
  save(results, file=paste0(path_to_data2, speciesCode, "/crossvalidation/",
                            speciesCode, "_crossvalidation_results.RData"))
}

##
path_to_data <- "E:/"



for(j in 3:7){
  load(paste0(path_to_data2, speciesList[j], "/crossvalidation/",
              speciesList[j], "_crossvalidation_results.RData"))


  results_summary <- data.frame(species=rep(speciesList[j], 4),
                                method=methods,
                                br = rep(NA,4),
                                br_scale = rep(NA,4),
                                log_lik = rep(NA,4),
                                log_lik_scale = rep(NA,4))
  for (k in 3:6){
    results_summary[1, k] <- mean(results[1:5, k+1], na.rm=T)
    results_summary[2,k]  <- mean(results[6:10, k+1], na.rm=T)
    results_summary[3,k]  <- mean(results[11:15, k+1],na.rm=T)
    results_summary[4,k]  <- mean(results[16:20, k+1], na.rm=T)
  }

  results_summary$log_lik <- -1*results_summary$log_lik
  results_summary$log_lik_scale <- -1*results_summary$log_lik_scale

  save(results_summary, file=paste0(path_to_data2, speciesList[j], "/crossvalidation/",
                                    speciesList[j],
                                    "_crossvalidation_results_summary.RData"))
}




