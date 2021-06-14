# Assemble Simple Model Folds for Crossvalidation
#Assemble simple model folds



load(paste0(path_to_data, "unit_cov2.RData"))
load(paste0(path_to_data, "florida.RData"))



speciesCodes <- c("MANIND", "METTOX", "RHUCOP", "SCHTER", "TOXPUB", "TOXRAD", "TOXVER")
speciesList <- c("Mangifera indica", "Metopium toxiferum", "Rhus copallina", "Schinus terebinthifolia",
                 "Toxicodendron pubescens", "Toxicodendron radicans",
                 "Toxicodendron vernix")
methods <- c("simplemod_alldata", "simplemod_filtered")

counties <- unique(florida@data$NAME_2)
counties <- factor(counties)

for(i in 1:7){


  #simplemod_alldata
  load(paste0(path_to_data, speciesList[i], "/surv_cov.RData"))
  h_simple <- data.frame(county=counties)

  for(j in 1:length(counties)) {
    temp <- subset(surv_cov, surv_cov$county == counties[j])
    h_simple$h[j] <- sum(temp$h)
    h_simple$nsurveys[j] <- nrow(temp)
  }

  nSites <- 67
  nFolds <- 5
  ( foldSize <- nSites %/% nFolds + 1 )
  tmp <- matrix(0, nFolds, foldSize)
  for(k in 1:foldSize)
    tmp[, k] <- sample.int(nFolds)
  foldID <- as.vector(tmp)[1:nSites]
  table(foldID)



  fold1 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==1))
  fold2 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==2))
  fold3 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==3))
  fold4 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==4))
  fold5 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==5))

  reassign <- function(fold) {
    return(list("fold"=fold))
  }

  run1 <- reassign(fold1)
  run2 <- reassign(fold2)
  run3 <- reassign(fold3)
  run4 <- reassign(fold4)
  run5 <- reassign(fold5)
  save(run1, run2, run3, run4, run5,
       file=paste0(path_to_data, speciesCodes[i],
                   "/crossvalidation/", "simplemod_alldata_folds.RData"))

  #simplemod_filtered

  load(paste0(path_to_data, speciesList[i], "/surv_cov3.RData"))
  surv_cov <- subset(surv_cov3, surv_cov3$num_genera > 1 | surv_cov3$h==1)
  rm(surv_cov3)

  h_simple <- data.frame(county=counties)

  for(j in 1:length(counties)) {
    temp <- subset(surv_cov, surv_cov$county == counties[j])
    h_simple$h[j] <- sum(temp$h)
    h_simple$nsurveys[j] <- nrow(temp)
  }

  nSites <- 67
  nFolds <- 5
  ( foldSize <- nSites %/% nFolds + 1 )
  tmp <- matrix(0, nFolds, foldSize)
  for(k in 1:foldSize)
    tmp[, k] <- sample.int(nFolds)
  foldID <- as.vector(tmp)[1:nSites]
  table(foldID)



  fold1 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==1))
  fold2 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==2))
  fold3 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==3))
  fold4 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==4))
  fold5 <- subset(h_simple, as.numeric(h_simple$county) %in% which(foldID==5))

  reassign <- function(fold) {
    return(list("fold"=fold))
  }

  run1 <- reassign(fold1)
  run2 <- reassign(fold2)
  run3 <- reassign(fold3)
  run4 <- reassign(fold4)
  run5 <- reassign(fold5)
  save(run1, run2, run3, run4, run5,
       file=paste0(path_to_data, speciesCodes[i],
                   "/crossvalidation/", "simplemod_filtered_folds.RData"))
}

