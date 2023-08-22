# ***********************************************************************************
# Data
# ***********************************************************************************

rm(list = ls())
setwd("~/Desktop/Sea_Scallops_MWSH/script/")
library(INLA)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(TMB)
library(scales)
options(stringsAsFactors = F)

# Read data
data.covar <- read.csv("../data/20190201-Jessica-temperature_on_tow/tows.2011to2018wTemperature.csv")
data.bio <- read.csv('../data/20181114-Jessica-Data/WeightHeightAge_Strata22_23_24_clean.csv') 

data.all <- inner_join(data.bio, data.covar, by = c("ID", "CRUISE", "TOW_NO", "TOW_DATE", "YEAR", "STRATA_ID")) %>%
    filter(YEAR >= 2012, 
           !is.na(WET_MEAT_WGT),
           !is.na(HEIGHT),
           HEIGHT > 40) %>%
    transmute(weight = WET_MEAT_WGT, 
              height = HEIGHT,
              year = as.factor(YEAR),
              lon=MIDPOINT_LONG, lat=MIDPOINT_LAT, 
              ID_TOW = as.factor(paste(YEAR,TOW_NO,sep = '_')), 
              depth = -DEPTH_M,
              temperature = BOTTOM_TEMP,
              area = case_when(STRATA_ID==22~"SMB",STRATA_ID!=22&lon<=-66.4~"Outside",STRATA_ID!=22&lon>-66.4~"Inside")) 



# # ***********************************************************************************
# # Cross validation: GLMM with year effect
# # ***********************************************************************************
# 
# rm(list = setdiff(ls(),c("data.all")))
# 
# i.model <- "GLMM"
# 
# all_set <- data.all %>% 
#     mutate(height = as.numeric(scale(log(height))))
# 
# 
# # CV settings
# set.seed(111)
# n_fold <- 10
# 
# list_tow <- unique(all_set$ID_TOW)
# folds_tow <- caret::createFolds(list_tow, k = n_fold, list = T, returnTrain = F)
# folds <- lapply(folds_tow, function(x) which(all_set$ID_TOW %in% list_tow[x]))
# 
# obs <- fitted <- resid <- list()
# test_pred <- test_error <- vector("numeric", nrow(all_set))
# 
# for(i_fold in 1:n_fold){
#     # same partitioning of dataset
#     train_set <- all_set[setdiff(1:nrow(all_set), folds[[i_fold]]),]
#     test_set <- all_set[folds[[i_fold]],]
#     
#     # fit model with year as a fixed effect
#     fit.GLMM <- lme4::glmer(
#         data = train_set,
#         formula = weight~0+year+height+(height|ID_TOW),
#         family=Gamma(link = "log"),
#         na.action = na.omit)
#     
#     test_pred[folds[[i_fold]]] <- predict(fit.GLMM, newdata = test_set, re.form=NULL, type="response", allow.new.levels = T)
#     test_error[folds[[i_fold]]]  <- test_set$weight - test_pred[folds[[i_fold]]]
# }
# 
# res <- data.frame("test_pred"=test_pred, "test_resid"=test_error)
# save(file = paste0("TOW_cv/res.log_log.Gamma.", i.model, ".RData"), res)
# 
# 
# 
# # ***********************************************************************************
# # Cross validation: GLMM with year effect, depth and temperature
# # ***********************************************************************************
# 
# rm(list = setdiff(ls(),c("data.all")))
# 
# i.model <- "GLMM-DT"
# 
# all_set <- data.all %>% 
#     mutate(height = as.numeric(scale(log(height))),
#            depth = as.numeric(scale(log(depth))),
#            temperature = as.numeric(scale(log(temperature))))
# 
# 
# # CV settings
# set.seed(111)
# n_fold <- 10
# 
# list_tow <- unique(all_set$ID_TOW)
# folds_tow <- caret::createFolds(list_tow, k = n_fold, list = T, returnTrain = F)
# folds <- lapply(folds_tow, function(x) which(all_set$ID_TOW %in% list_tow[x]))
# 
# obs <- fitted <- resid <- list()
# test_pred <- test_error <- vector("numeric", nrow(all_set))
# 
# for(i_fold in 1:n_fold){
#     # same partitioning of dataset
#     train_set <- all_set[setdiff(1:nrow(all_set), folds[[i_fold]]),]
#     test_set <- all_set[folds[[i_fold]],]
#     
#     # fit model with year as a fixed effect
#     fit.GLMM <- lme4::glmer(
#         data = train_set,
#         formula = weight~0+year+depth+temperature+height+(height|ID_TOW),
#         family=Gamma(link = "log"),
#         na.action = na.omit)
#     
#     test_pred[folds[[i_fold]]] <- predict(fit.GLMM, newdata = test_set, re.form=NULL, type="response", allow.new.levels = T)
#     test_error[folds[[i_fold]]]  <- test_set$weight - test_pred[folds[[i_fold]]]
# }
# 
# res <- data.frame("test_pred"=test_pred, "test_resid"=test_error)
# save(file = paste0("TOW_cv/res.log_log.Gamma.", i.model, ".RData"), res)
# 


# # ***********************************************************************************
# # Cross validation: GLMM with year effect, depth and temperature area effect
# # ***********************************************************************************
# 
# rm(list = setdiff(ls(),c("data.all")))
# 
# i.model <- "GLMM-DTA"
# 
# all_set <- data.all %>%
#     mutate(height = as.numeric(scale(log(height))),
#            depth = as.numeric(scale(log(depth))),
#            temperature = as.numeric(scale(log(temperature))))
# 
# 
# # CV settings
# set.seed(111)
# n_fold <- 10
# 
# list_tow <- unique(all_set$ID_TOW)
# folds_tow <- caret::createFolds(list_tow, k = n_fold, list = T, returnTrain = F)
# folds <- lapply(folds_tow, function(x) which(all_set$ID_TOW %in% list_tow[x]))
# 
# obs <- fitted <- resid <- list()
# test_pred <- test_error <- vector("numeric", nrow(all_set))
# 
# for(i_fold in 1:n_fold){
#     # same partitioning of dataset
#     train_set <- all_set[setdiff(1:nrow(all_set), folds[[i_fold]]),]
#     test_set <- all_set[folds[[i_fold]],]
# 
#     # fit model with year as a fixed effect
#     fit.GLMM <- lme4::glmer(
#         data = train_set,
#         formula = weight~0+year+depth+temperature+height+area+(height|ID_TOW),
#         family=gaussian(link = "log"),
#         na.action = na.omit)
# 
#     test_pred[folds[[i_fold]]] <- predict(fit.GLMM, newdata = test_set, re.form=NULL, type="response", allow.new.levels = T)
#     test_error[folds[[i_fold]]]  <- test_set$weight - test_pred[folds[[i_fold]]]
# }
# 
# res <- data.frame("test_pred"=test_pred, "test_resid"=test_error)
# save(file = paste0("appendix/TOW_cv/res.log_log.Gauss.", i.model, ".RData"), res)


# # ***********************************************************************************
# # Cross validation: spatial model with only static sp effect
# # ***********************************************************************************
# 
# rm(list = setdiff(ls(),c("data.all")))
# 
# i.model <- "SPM"
# 
# all_set <- data.all %>% 
#     mutate(height = as.numeric(scale(log(height)))) %>%
#     mutate(X=lon, Y=lat) %>%
#     `attr<-`("projection", "LL") %>%
#     `attr<-`("zone", "20") %>%
#     PBSmapping::convUL() 
# 
# # create mesh: use UTM
# mesh = inla.mesh.create(all_set[,c("X","Y")], refine = F, extend = F)
# 
# # Create matrices for approximation
# spde <- inla.spde2.matern(mesh, alpha = 2)
# 
# 
# # TMB
# version <- "STM"
# compile(paste0(version,".cpp"))
# dyn.load(dynlib(version))
# 
# # CV settings
# set.seed(111)
# n_fold <- 10
# list_tow <- unique(all_set$ID_TOW)
# folds_tow <- caret::createFolds(list_tow, k = n_fold, list = T, returnTrain = F)
# folds <- lapply(folds_tow, function(x) which(all_set$ID_TOW %in% list_tow[x]))
# 
# # CV
# obs <- fitted <- resid <- list()
# test_pred <- test_error <- vector("numeric", nrow(all_set))
# res_list <- list()
# 
# for(i_fold in 1:n_fold){
#     train_set <- all_set[setdiff(1:nrow(all_set), folds[[i_fold]]),]
#     test_set <- all_set[folds[[i_fold]],]
#     
#     data = list(
#         weight = train_set$weight, 
#         varmat = as.matrix(cbind(1, train_set$height)),
#         ind_loc = train_set$ID_TOW,
#         ind_year = train_set$year,
#         M0 = spde$param.inla$M0,
#         M1 = spde$param.inla$M1,
#         M2 = spde$param.inla$M2)
#     
#     nyear = nlevels(data$ind_year)
#     nloc = nrow(data$M0)
#     nvar = ncol(data$varmat)
#     nobs = nrow(data$varmat)
#     
#     parameters = list(
#         beta = rep(0, nvar),
#         beta_s = matrix(0, nloc, nvar),
#         beta_t = matrix(0, nyear, nvar),
#         beta_st = array(0, dim = c(nloc, nyear, nvar)),
#         log_kappa_s = rep(log(0.1), nvar),
#         log_tau_s = rep(log(0.1), nvar),
#         log_kappa_st = matrix(log(0.1), 1, nvar),
#         log_tau_st = matrix(log(0.1), 1, nvar),
#         log_phi = log(0.1))
#     
#     maps <- list(
#         beta = factor(rep(NA, nvar)),
#         beta_st = factor(array(NA, dim = c(nloc, nyear, nvar))),
#         log_kappa_st = factor(matrix(NA, 1, nvar)),
#         log_tau_st = factor(matrix(NA, 1, nvar))
#     )
#     
#     obj = MakeADFun(data=data, 
#                     parameters=parameters,
#                     map=maps,
#                     random=c("beta_s"),
#                     DLL=version,
#                     silent = F)
#     
#     opt <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=1000, iter.max=1000))
#     
#     rep <- obj$report()
#     
#     obs[[i_fold]] <- train_set$weight
#     fitted[[i_fold]] <- rep$mu
#     resid[[i_fold]] <- obs[[i_fold]] - fitted[[i_fold]]
#     
#     tmp_pred <- tmp_resid <- rep(NA, nrow(test_set))
#     for(i in 1:nrow(test_set)){
#         tmp_pred[i] <- exp((rep$beta + rep$beta_t[as.integer(test_set[i,]$year),] + rep$beta_s[as.integer(test_set[i,]$ID_TOW),] + rep$beta_st[as.integer(test_set[i,]$ID_TOW),as.integer(test_set[i,]$year),]) %*% c(1, test_set[i,]$height))
#         tmp_resid[i] <- test_set$weight[i] - tmp_pred[i]
#     }
#     
#     test_pred[folds[[i_fold]]] <- tmp_pred
#     test_error[folds[[i_fold]]] <- tmp_resid
#     
#     res_list[[i_fold]] <- list(obj = obj, opt = opt, rep = rep)
# }
# 
# res <- data.frame("test_pred"=test_pred, "test_resid"=test_error)
# save(file = paste0("appendix/TOW_cv/res.log_log.Gauss.", i.model, ".RData"), res)
# save(file = paste0("appendix/TOW_cv/res.log_log.Gauss.", i.model, ".res_list.RData"), res_list)



# # ***********************************************************************************
# # Cross validation: spatial temporal model
# # ***********************************************************************************
# 
# rm(list = setdiff(ls(),c("data.all")))
# 
# i.model <- "STM"
# 
# all_set <- data.all %>% 
#     mutate(height = as.numeric(scale(log(height)))) %>%
#     mutate(X=lon, Y=lat) %>%
#     `attr<-`("projection", "LL") %>%
#     `attr<-`("zone", "20") %>%
#     PBSmapping::convUL() 
# 
# # create mesh: use UTM
# mesh = inla.mesh.create(all_set[,c("X","Y")], refine = F, extend = F)
# 
# # Create matrices for approximation
# spde <- inla.spde2.matern(mesh, alpha = 2)
# 
# 
# # TMB
# version <- "STM"
# compile(paste0(version,".cpp"))
# dyn.load(dynlib(version))
# 
# # CV settings
# set.seed(111)
# n_fold <- 10
# list_tow <- unique(all_set$ID_TOW)
# folds_tow <- caret::createFolds(list_tow, k = n_fold, list = T, returnTrain = F)
# folds <- lapply(folds_tow, function(x) which(all_set$ID_TOW %in% list_tow[x]))
# 
# # CV
# obs <- fitted <- resid <- list()
# test_pred <- test_error <- vector("numeric", nrow(all_set))
# res_list <- list()
# 
# for(i_fold in 1:n_fold){
#     train_set <- all_set[setdiff(1:nrow(all_set), folds[[i_fold]]),]
#     test_set <- all_set[folds[[i_fold]],]
#     
#     data = list(
#         weight = train_set$weight, 
#         varmat = as.matrix(cbind(1, train_set$height)),
#         ind_loc = train_set$ID_TOW,
#         ind_year = train_set$year,
#         M0 = spde$param.inla$M0,
#         M1 = spde$param.inla$M1,
#         M2 = spde$param.inla$M2)
#     
#     nyear = nlevels(data$ind_year)
#     nloc = nrow(data$M0)
#     nvar = ncol(data$varmat)
#     nobs = nrow(data$varmat)
#     
#     parameters = list(
#         beta = rep(0, nvar),
#         beta_s = matrix(0, nloc, nvar),
#         beta_t = matrix(0, nyear, nvar),
#         beta_st = array(0, dim = c(nloc, nyear, nvar)),
#         log_kappa_s = rep(log(0.1), nvar),
#         log_tau_s = rep(log(0.1), nvar),
#         log_kappa_st = matrix(log(0.1), 1, nvar),
#         log_tau_st = matrix(log(0.1), 1, nvar),
#         log_phi = log(0.1))
#     
#     maps <- list(
#         beta = factor(rep(NA, nvar))
#     )
#     
#     obj = MakeADFun(data=data, 
#                     parameters=parameters,
#                     map=maps,
#                     random=c("beta_s","beta_st"),
#                     DLL=version,
#                     silent = F)
#     
#     opt <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=1000, iter.max=1000))
#     
#     rep <- obj$report()
#     
#     obs[[i_fold]] <- train_set$weight
#     fitted[[i_fold]] <- rep$mu
#     resid[[i_fold]] <- obs[[i_fold]] - fitted[[i_fold]]
#     
#     tmp_pred <- tmp_resid <- rep(NA, nrow(test_set))
#     for(i in 1:nrow(test_set)){
#         tmp_pred[i] <- exp((rep$beta + rep$beta_t[as.integer(test_set[i,]$year),] + rep$beta_s[as.integer(test_set[i,]$ID_TOW),] + rep$beta_st[as.integer(test_set[i,]$ID_TOW),as.integer(test_set[i,]$year),]) %*% c(1, test_set[i,]$height))
#         tmp_resid[i] <- test_set$weight[i] - tmp_pred[i]
#     }
#     
#     test_pred[folds[[i_fold]]] <- tmp_pred
#     test_error[folds[[i_fold]]] <- tmp_resid
#     
#     res_list[[i_fold]] <- list(obj = obj, opt = opt, rep = rep)
# }
# 
# res <- data.frame("test_pred"=test_pred, "test_resid"=test_error)
# save(file = paste0("TOW_cv/res.log_log.Gauss.", i.model, ".RData"), res)
# save(file = paste0("TOW_cv/res.log_log.Gauss.", i.model, ".res_list.RData"), res_list)



# ***********************************************************************************
# Cross validation: static spatial model + depth temperature on intercept
# ***********************************************************************************

rm(list = setdiff(ls(),c("data.all")))

i.model <- "SPM-DT"

all_set <- data.all %>% 
    mutate(height = as.numeric(scale(log(height))),
           depth = as.numeric(scale(log(depth))),
           temperature = as.numeric(scale(log(temperature)))) %>%
    mutate(X=lon, Y=lat) %>%
    `attr<-`("projection", "LL") %>%
    `attr<-`("zone", "20") %>%
    PBSmapping::convUL() 

# create mesh: use UTM
mesh = inla.mesh.create(all_set[,c("X","Y")], refine = F, extend = F)

# Create matrices for approximation
spde <- inla.spde2.matern(mesh, alpha = 2)


# TMB
version <- "STM_DT"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))

# CV settings
set.seed(111)
n_fold <- 10
list_tow <- unique(all_set$ID_TOW)
folds_tow <- caret::createFolds(list_tow, k = n_fold, list = T, returnTrain = F)
folds <- lapply(folds_tow, function(x) which(all_set$ID_TOW %in% list_tow[x]))

# CV
obs <- fitted <- resid <- list()
test_pred <- test_error <- vector("numeric", nrow(all_set))
res_list <- list()

for(i_fold in 1:n_fold){
    train_set <- all_set[setdiff(1:nrow(all_set), folds[[i_fold]]),]
    test_set <- all_set[folds[[i_fold]],]
    
    data = list(
        weight = train_set$weight, 
        varmat = as.matrix(cbind(1, train_set$height)),
        depth = train_set$depth,
        temperature = train_set$temperature,
        ind_loc = train_set$ID_TOW,
        ind_year = train_set$year,
        M0 = spde$param.inla$M0,
        M1 = spde$param.inla$M1,
        M2 = spde$param.inla$M2)
    
    nyear = nlevels(data$ind_year)
    nloc = nrow(data$M0)
    nvar = ncol(data$varmat)
    nobs = nrow(data$varmat)
    
    parameters = list(
        beta = rep(0, nvar),
        beta_s = matrix(0, nloc, nvar),
        beta_t = matrix(0, nyear, nvar),
        beta_st = array(0, dim = c(nloc, nyear, nvar)),
        beta_depth = matrix(0, nyear, nvar),
        beta_temperature = matrix(0, nyear, nvar),
        log_kappa_s = rep(log(0.1), nvar),
        log_tau_s = rep(log(0.1), nvar),
        log_kappa_st = matrix(log(0.1), 1, nvar),
        log_tau_st = matrix(log(0.1), 1, nvar),
        log_phi = log(0.1))
    
    # one parameter for each of depth and temperature
    maps <- list(
        beta = factor(rep(NA, nvar)),
        beta_st = factor(array(NA, dim = c(nloc, nyear, nvar))),
        log_kappa_st = factor(matrix(NA, 1, nvar)),
        log_tau_st = factor(matrix(NA, 1, nvar)),
        beta_depth = factor(cbind(rep(1,nyear),rep(NA,nyear))),
        beta_temperature = factor(cbind(rep(1,nyear),rep(NA,nyear)))
    )
    
    obj = MakeADFun(data=data, 
                    parameters=parameters,
                    map=maps,
                    random=c("beta_s"),
                    DLL=version,
                    silent = F)
    
    opt <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=1000, iter.max=1000))
    
    rep <- obj$report()
    
    obs[[i_fold]] <- train_set$weight
    fitted[[i_fold]] <- rep$mu
    resid[[i_fold]] <- obs[[i_fold]] - fitted[[i_fold]]
    
    tmp_pred <- tmp_resid <- rep(NA, nrow(test_set))
    for(i in 1:nrow(test_set)){
        tmp_pred[i] <- exp((rep$beta + 
                                rep$beta_t[as.integer(test_set[i,]$year),] + 
                                rep$beta_s[as.integer(test_set[i,]$ID_TOW),] + 
                                rep$beta_st[as.integer(test_set[i,]$ID_TOW),as.integer(test_set[i,]$year),] + 
                                rep$beta_depth[as.integer(test_set[i,]$year),]*test_set[i,]$depth + 
                                rep$beta_temperature[as.integer(test_set[i,]$year),]*test_set[i,]$temperature
        ) %*% c(1, test_set[i,]$height))
        tmp_resid[i] <- test_set$weight[i] - tmp_pred[i]
    }
    
    test_pred[folds[[i_fold]]] <- tmp_pred
    test_error[folds[[i_fold]]] <- tmp_resid
    
    res_list[[i_fold]] <- list(obj = obj, opt = opt, rep = rep)
}

res <- data.frame("test_pred"=test_pred, "test_resid"=test_error)
save(file = paste0("appendix/TOW_cv/res.log_log.Gauss.", i.model, ".RData"), res)
save(file = paste0("appendix/TOW_cv/res.log_log.Gauss.", i.model, ".res_list.RData"), res_list)


# # ***********************************************************************************
# # Cross validation: spatial temporal model + depth temperature on intercept
# # ***********************************************************************************
# 
# rm(list = setdiff(ls(),c("data.all")))
# 
# i.model <- "STM-DT"
# 
# all_set <- data.all %>% 
#     mutate(height = as.numeric(scale(log(height))),
#            depth = as.numeric(scale(log(depth))),
#            temperature = as.numeric(scale(log(temperature)))) %>%
#     mutate(X=lon, Y=lat) %>%
#     `attr<-`("projection", "LL") %>%
#     `attr<-`("zone", "20") %>%
#     PBSmapping::convUL() 
# 
# # create mesh: use UTM
# mesh = inla.mesh.create(all_set[,c("X","Y")], refine = F, extend = F)
# 
# # Create matrices for approximation
# spde <- inla.spde2.matern(mesh, alpha = 2)
# 
# 
# # TMB
# version <- "STM_DT"
# compile(paste0(version,".cpp"))
# dyn.load(dynlib(version))
# 
# # CV settings
# set.seed(111)
# n_fold <- 10
# list_tow <- unique(all_set$ID_TOW)
# folds_tow <- caret::createFolds(list_tow, k = n_fold, list = T, returnTrain = F)
# folds <- lapply(folds_tow, function(x) which(all_set$ID_TOW %in% list_tow[x]))
# 
# # CV
# obs <- fitted <- resid <- list()
# test_pred <- test_error <- vector("numeric", nrow(all_set))
# res_list <- list()
# 
# for(i_fold in 1:n_fold){
#     train_set <- all_set[setdiff(1:nrow(all_set), folds[[i_fold]]),]
#     test_set <- all_set[folds[[i_fold]],]
#     
#     data = list(
#         weight = train_set$weight, 
#         varmat = as.matrix(cbind(1, train_set$height)),
#         depth = train_set$depth,
#         temperature = train_set$temperature,
#         ind_loc = train_set$ID_TOW,
#         ind_year = train_set$year,
#         M0 = spde$param.inla$M0,
#         M1 = spde$param.inla$M1,
#         M2 = spde$param.inla$M2)
#     
#     nyear = nlevels(data$ind_year)
#     nloc = nrow(data$M0)
#     nvar = ncol(data$varmat)
#     nobs = nrow(data$varmat)
#     
#     parameters = list(
#         beta = rep(0, nvar),
#         beta_s = matrix(0, nloc, nvar),
#         beta_t = matrix(0, nyear, nvar),
#         beta_st = array(0, dim = c(nloc, nyear, nvar)),
#         beta_depth = matrix(0, nyear, nvar),
#         beta_temperature = matrix(0, nyear, nvar),
#         log_kappa_s = rep(log(0.1), nvar),
#         log_tau_s = rep(log(0.1), nvar),
#         log_kappa_st = matrix(log(0.1), 1, nvar),
#         log_tau_st = matrix(log(0.1), 1, nvar),
#         log_phi = log(0.1))
#     
#     # one parameter for each of depth and temperature
#     maps <- list(
#         beta = factor(rep(NA, nvar)),
#         beta_depth = factor(cbind(rep(1,nyear),rep(NA,nyear))),
#         beta_temperature = factor(cbind(rep(1,nyear),rep(NA,nyear)))
#     )
#     
#     obj = MakeADFun(data=data, 
#                     parameters=parameters,
#                     map=maps,
#                     random=c("beta_s","beta_st"),
#                     DLL=version,
#                     silent = F)
#     
#     opt <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=1000, iter.max=1000))
#     
#     rep <- obj$report()
#     
#     obs[[i_fold]] <- train_set$weight
#     fitted[[i_fold]] <- rep$mu
#     resid[[i_fold]] <- obs[[i_fold]] - fitted[[i_fold]]
#     
#     tmp_pred <- tmp_resid <- rep(NA, nrow(test_set))
#     for(i in 1:nrow(test_set)){
#         tmp_pred[i] <- exp((rep$beta + 
#                                 rep$beta_t[as.integer(test_set[i,]$year),] + 
#                                 rep$beta_s[as.integer(test_set[i,]$ID_TOW),] + 
#                                 rep$beta_st[as.integer(test_set[i,]$ID_TOW),as.integer(test_set[i,]$year),] + 
#                                 rep$beta_depth[as.integer(test_set[i,]$year),]*test_set[i,]$depth + 
#                                 rep$beta_temperature[as.integer(test_set[i,]$year),]*test_set[i,]$temperature
#         ) %*% c(1, test_set[i,]$height))
#         tmp_resid[i] <- test_set$weight[i] - tmp_pred[i]
#     }
#     
#     test_pred[folds[[i_fold]]] <- tmp_pred
#     test_error[folds[[i_fold]]] <- tmp_resid
#     
#     res_list[[i_fold]] <- list(obj = obj, opt = opt, rep = rep)
# }
# 
# res <- data.frame("test_pred"=test_pred, "test_resid"=test_error)
# save(file = paste0("TOW_cv/res.log_log.Gauss.", i.model, ".RData"), res)
# save(file = paste0("TOW_cv/res.log_log.Gauss.", i.model, ".res_list.RData"), res_list)





# ***********************************************************************************
# Result analysis: selected models
# ***********************************************************************************

rm(list = setdiff(ls(),c("data.all")))


# gather results
x <- list()
load(paste0("appendix/TOW_cv/res.log_log.Gauss.GLMM.RData")); x[["GLMM"]] <- res; rm("res")
load(paste0("appendix/TOW_cv/res.log_log.Gauss.GLMM-DT.RData")); x[["GLMM-DT"]] <- res; rm("res")
load(paste0("appendix/TOW_cv/res.log_log.Gauss.STM.RData")); x[["STM"]] <- res; rm("res")
load(paste0("appendix/TOW_cv/res.log_log.Gauss.STM-DT.RData")); x[["STM-DT"]] <- res; rm("res")
load(paste0("appendix/TOW_cv/res.log_log.Gauss.GLMM-DTA.RData")); x[["GLMM-DTA"]] <- res; rm("res")
load(paste0("appendix/TOW_cv/res.log_log.Gauss.SPM.RData")); x[["SPM"]] <- res; rm("res")
load(paste0("appendix/TOW_cv/res.log_log.Gauss.SPM-DT.RData")); x[["SPM-DT"]] <- res; rm("res")

res <- data.frame(
    "lon"=data.all$lon,
    "lat"=data.all$lat,
    "Year"=data.all$year,
    # "GLMM"=x$`GLMM`$test_resid,
    # "GLMM-DT"=x$`GLMM-DT`$test_resid,
    # "STM"=x$`STM`$test_resid,
    # "STM-DT"=x$`STM-DT`$test_resid,
    "GLMM-DTA"=x$`GLMM-DTA`$test_resid,
    "SPM"=x$`SPM`$test_resid,
    "SPM-DT"=x$`SPM-DT`$test_resid
) %>%
    rename(
        # "GLMM-DT"="GLMM.DT",
        # "STM-DT"="STM.DT",
        "GLMM-DTA"="GLMM.DTA",
        "SPM-DT"="SPM.DT") %>%
    gather(model,resid,-lon,-lat,-Year) %>%
    mutate(model = factor(model, ordered = T))


res.sp <- res %>%
    group_by(Year, lon, lat, model) %>%
    summarise(m.resid = mean(resid), sd.resid = sd(resid), m.abs.resid = mean(abs(resid))) %>% 
    ungroup()


# cv test set resid and sd
bind_rows(
    res %>%
        group_by(Year, model) %>% 
        summarise(indiv.resid.mean = paste0(format(round(mean((resid)),4),nsmall=4, scientific=F)," (", round(sd(resid),4), ")")) %>%
        spread(model, indiv.resid.mean),
    res %>%
        group_by(model) %>% 
        summarise(indiv.resid.mean = paste0(format(round(mean((resid)),4),nsmall=4, scientific=F)," (", round(sd(resid),4), ")")) %>%
        spread(model, indiv.resid.mean) %>% mutate(Year = "2012-2018")
) %>%
    xtable::xtable() %>%
    print(include.rownames=F)
    

# box plot (reviewer suggested)
ggplot(data=res %>%
           group_by(Year, model) %>% 
           summarise(m = mean(resid), s = sd(resid)) %>%
           ungroup() %>%
           mutate(year = as.numeric(as.character(Year)) + as.numeric(as.factor(model))/10 - .25)) +
    geom_pointrange(aes(x = year + as.numeric(factor(model))/10 - .25, y = m, ymin = m-s, ymax = m+s, color = model, shape = model)) +
    geom_line(aes(x = year + as.numeric(factor(model))/10 - .25, y = m, color = model)) +
    geom_hline(yintercept = 0, color = "red", alpha = 0.5) +
    xlab("Year") + ylab("Cross-Validated Residuals") + labs(color = "Model", shape = "Model") +
    ylim(-10,10) +
    scale_colour_manual(values = c("purple", "orange", "blue", "green")) +
    scale_x_continuous(breaks = 2012:2018) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave(filename = "appendix/resid_pointrange_log_log.Gauss.jpg", width = 10, height = 6)



# plot setting
base_map <- ggplot() +
    borders("world", colour="gray50", fill = "gray90", xlim = c(-67,-65), ylim = c(43,45)) +
    coord_map(xlim = c(-67,-65.75), ylim = c(43.5, 44.75)) +
    theme_bw() +
    scale_color_continuous(low = "white", high = "red") +
    scale_size_continuous(guide = FALSE) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(), 
          legend.position = "bottom")

p.resid.tow <- base_map + 
    geom_point(data = res.sp, aes(x=lon, y=lat, color = m.resid),shape = 20) +
    scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits=c(-6,6), oob=squish,
                          guide_legend(title = "Spatial Residual")) +
    facet_grid(model~Year) 
ggsave(filename = "TOW_cv/resid_sp_map_log_log.Gauss.jpg", width = 10, height = 10)






# ***********************************************************************************
# prediction for 100mm
# ***********************************************************************************

rm(list = setdiff(ls(),c("data.all")))

# new data set for prediction
h.sc.a <- attr(scale(log(data.all$height)), "scaled:center")
h.sc.b <- attr(scale(log(data.all$height)), "scaled:scale")

# check scaling
any(!as.numeric(scale(log(data.all$height))) == (log(data.all$height) - h.sc.a)/h.sc.b)

pred.height <- (log(100) - h.sc.a)/h.sc.b

newdata <- data.all %>% 
    mutate(height = as.numeric(scale(log(height))),
           depth = as.numeric(scale(log(depth))),
           temperature = as.numeric(scale(log(temperature)))) %>%
    group_by(ID_TOW) %>% 
    slice(1) %>%
    mutate(height = pred.height) %>%
    select(-weight) %>%
    mutate(GLMM = as.numeric(NA), 
           `GLMM-DT` = as.numeric(NA),
           STM = as.numeric(NA), 
           `STM-DT` = as.numeric(NA))


# ------------------------
# GLMM

all_set <- data.all %>% mutate(height = as.numeric(scale(log(height))))

fit.GLMM <- lme4::glmer(
    data = all_set,
    formula = weight~0+year+height+(height|ID_TOW),
    # family=Gamma(link = "log"),
    family = gaussian(link = "log"),
    na.action = na.omit)
newdata$GLMM <- predict(fit.GLMM, newdata, re.form=NULL, type="response", allow.new.levels = T)

rm("all_set")


# ------------------------
# GLMM with env effect

all_set <-  data.all %>% mutate(height = as.numeric(scale(log(height))),
                                depth = as.numeric(scale(log(depth))),
                                temperature = as.numeric(scale(log(temperature)))) 

fit.GLMM.DT <- lme4::glmer(
    data = all_set,
    formula = weight~0+year+depth+temperature+height+(height|ID_TOW),
    # family=Gamma(link = "log"),
    family = gaussian(link = "log"),
    na.action = na.omit)

newdata$`GLMM-DT` <- predict(fit.GLMM.DT, newdata, re.form=NULL, type="response", allow.new.levels = T)

rm("all_set")


# ------------------------
# spatiotemporal model

all_set <- data.all %>% 
    mutate(height = as.numeric(scale(log(height)))) %>%
    mutate(X=lon, Y=lat) %>%
    `attr<-`("projection", "LL") %>%
    `attr<-`("zone", "20") %>%
    PBSmapping::convUL() 

# create mesh: use UTM
mesh = inla.mesh.create(all_set[,c("X","Y")], refine = F, extend = F)

# Create matrices for approximation
spde <- inla.spde2.matern(mesh, alpha = 2)

# TMB
version <- "STM"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))

data = list(
    weight = all_set$weight, 
    varmat = as.matrix(cbind(1, all_set$height)),
    ind_loc = all_set$ID_TOW,
    ind_year = all_set$year,
    M0 = spde$param.inla$M0,
    M1 = spde$param.inla$M1,
    M2 = spde$param.inla$M2)

nyear = nlevels(data$ind_year)
nloc = nrow(data$M0)
nvar = ncol(data$varmat)
nobs = nrow(data$varmat)

parameters = list(
    beta = rep(0, nvar),
    beta_s = matrix(0, nloc, nvar),
    beta_t = matrix(0, nyear, nvar),
    beta_st = array(0, dim = c(nloc, nyear, nvar)),
    log_kappa_s = rep(log(0.1), nvar),
    log_tau_s = rep(log(0.1), nvar),
    log_kappa_st = matrix(log(0.1), 1, nvar),
    log_tau_st = matrix(log(0.1), 1, nvar),
    log_phi = log(0.1))

maps <- list(
    beta = factor(rep(NA, nvar))
)

obj = MakeADFun(data=data, 
                parameters=parameters,
                map=maps,
                random=c("beta_s","beta_st"),
                DLL=version,
                silent = T)

opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- obj$report()

for(i in 1:nrow(newdata)){
    newdata$STM[i] <- exp((rep$beta + rep$beta_t[as.integer(newdata[i,]$year),] + rep$beta_s[as.integer(newdata[i,]$ID_TOW),] + rep$beta_st[as.integer(newdata[i,]$ID_TOW),as.integer(newdata[i,]$year),]) %*% c(1, newdata[i,]$height))
}


rm("all_set","mesh","spde","data","parameters","maps","obj","opt","rep")


# ------------------------
# spatiotemporal model with depth and temperature

all_set <- data.all %>% 
    mutate(height = as.numeric(scale(log(height))),
           depth = as.numeric(scale(log(depth))),
           temperature = as.numeric(scale(log(temperature)))) %>%
    mutate(X=lon, Y=lat) %>%
    `attr<-`("projection", "LL") %>%
    `attr<-`("zone", "20") %>%
    PBSmapping::convUL() 

# create mesh: use UTM
mesh = inla.mesh.create(all_set[,c("X","Y")], refine = F, extend = F)

# Create matrices for approximation
spde <- inla.spde2.matern(mesh, alpha = 2)


# TMB
version <- "STM_DT"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))

data = list(
    weight = all_set$weight, 
    varmat = as.matrix(cbind(1, all_set$height)),
    depth = all_set$depth,
    temperature = all_set$temperature,
    ind_loc = all_set$ID_TOW,
    ind_year = all_set$year,
    M0 = spde$param.inla$M0,
    M1 = spde$param.inla$M1,
    M2 = spde$param.inla$M2)

nyear = nlevels(data$ind_year)
nloc = nrow(data$M0)
nvar = ncol(data$varmat)
nobs = nrow(data$varmat)

parameters = list(
    beta = rep(0, nvar),
    beta_s = matrix(0, nloc, nvar),
    beta_t = matrix(0, nyear, nvar),
    beta_st = array(0, dim = c(nloc, nyear, nvar)),
    beta_depth = matrix(0, nyear, nvar),
    beta_temperature = matrix(0, nyear, nvar),
    log_kappa_s = rep(log(0.1), nvar),
    log_tau_s = rep(log(0.1), nvar),
    log_kappa_st = matrix(log(0.1), 1, nvar),
    log_tau_st = matrix(log(0.1), 1, nvar),
    log_phi = log(0.1))

# one parameter for each of depth and temperature
maps <- list(
    beta = factor(rep(NA, nvar)),
    beta_depth = factor(cbind(rep(1,nyear),rep(NA,nyear))),
    beta_temperature = factor(cbind(rep(1,nyear),rep(NA,nyear)))
)

obj = MakeADFun(data=data, 
                parameters=parameters,
                map=maps,
                random=c("beta_s","beta_st"),
                DLL=version,
                silent = F)

opt <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=1000, iter.max=1000))

rep <- obj$report()

for(i in 1:nrow(newdata)){
    newdata$`STM-DT`[i]  <- exp((rep$beta + 
                            rep$beta_t[as.integer(newdata[i,]$year),] + 
                            rep$beta_s[as.integer(newdata[i,]$ID_TOW),] + 
                            rep$beta_st[as.integer(newdata[i,]$ID_TOW),as.integer(newdata[i,]$year),] + 
                            rep$beta_depth[as.integer(newdata[i,]$year),]*newdata[i,]$depth + 
                            rep$beta_temperature[as.integer(newdata[i,]$year),]*newdata[i,]$temperature
    ) %*% c(1, newdata[i,]$height))
}
    


save(file = "model_predict/predict_100mm_log_log.Gauss.rda",newdata)


# ------------------------
# generate figures 

load("appendix/model_predict/predict_100mm_log_log.Gauss.rda")

# plots

# plot setting
base_map <- ggplot() +
    borders("world", colour="gray50", fill = "gray90", xlim = c(-67,-65), ylim = c(43,45)) +
    coord_map(xlim = c(-67,-65.75), ylim = c(43.5, 44.75)) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(), 
          legend.position = "bottom")



pred.100 <- newdata %>% 
    select(year, lon, lat, ID_TOW, `GLMM`, `GLMM-DT`, STM, `STM-DT`) %>%
    gather("model", "pred.mw", -year, -lon, -lat, -ID_TOW) %>%
    rename(Year = year) %>%
    mutate(Model = factor(model, ordered = T, levels = c("GLMM", "GLMM-DT", "STM", "STM-DT")))


# # plots
# p.pred.100 <- base_map + 
#     geom_point(data = pred.100, aes(x=lon, y=lat, color=pred.mw, alpha = pred.mw)) +
#     scale_color_continuous(low = "white", high = "red", guide_legend(title = "Predicted Meat Weight/g")) +
#     facet_grid(Model~Year) +
#     guides(alpha = FALSE)
# ggsave(filename = "model_predict/pred_100_map_log.jpg", width = 10, height = 8)

# terrain.colors(10)
# rainbow(10)

base_map + 
    geom_point(data = pred.100, aes(x=lon, y=lat, color = pred.mw), shape=20) +
    scale_color_gradientn(colours = rainbow(5), limits=c(5,21),oob=scales::squish) +
    facet_grid(Model~Year) +
    labs(colour = "Predicted Meat Weight")
ggsave(filename = "model_predict/pred_100_map_log_log.Gauss.jpg", width = 10, height = 10)


p.pred.100_long <- ggplot(pred.100) +
    geom_point(aes(x=lon,pred.mw, color = Model)) +
    facet_wrap(~Year, nrow = 2) +
    theme_bw() +
    xlab("Longitude")+
    ylab("Predicted Meat Weight/g")+
    theme(legend.position = c(0.9, 0.2))
ggsave(filename = "model_predict/pred_100_long_log_log.Gauss.jpg", p.pred.100_long, width = 10, height = 8)


# box plot
ggplot(data=pred.100 %>%
           group_by(Year, model) %>% 
           summarise(m = mean(pred.mw), s = sd(pred.mw)) %>%
           ungroup() %>%
           mutate(year = as.numeric(as.character(Year)) + as.numeric(as.factor(model))/10 - .25)) +
    geom_pointrange(aes(x = year + as.numeric(factor(model))/10 - .25, y = m, ymin = m-s, ymax = m+s, color = model, shape = model)) +
    geom_line(aes(x = year + as.numeric(factor(model))/10 - .25, y = m, color = model)) +
    xlab("Year") + ylab("Predicted Meat Weight") + labs(color = "Model", shape = "Model") +
    scale_colour_manual(values = c("purple", "orange", "blue", "green")) +
    scale_x_continuous(breaks = 2012:2018) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          panel.grid.major.y = element_line(linetype = 2),
          legend.position = "bottom") + 
    scale_y_continuous(breaks = c(8:17))
ggsave(filename = "appendix/pred_100_pointrange_log_log.Gauss.jpg", width = 10, height = 6)







# ***********************************************************************************
# extract model parameters
# ***********************************************************************************

rm(list = setdiff(ls(),c("data.all")))


# ------------------------
# GLMM

all_set <- data.all %>% mutate(height = as.numeric(scale(log(height))))

fit.GLMM <- lme4::glmer(
    data = all_set,
    formula = weight~0+year+height+(height|ID_TOW),
    # family=Gamma(link = "log"),
    family = gaussian(link = "log"),
    na.action = na.omit)

save(file = "appendix/model_predict/fit.glmm.rda", fit.GLMM)

# ------------------------
# GLMM with year effect

all_set <-  data.all %>% mutate(height = as.numeric(scale(log(height))),
                                depth = as.numeric(scale(log(depth))),
                                temperature = as.numeric(scale(log(temperature)))) 

fit.GLMM.DT <- lme4::glmer(
    data = all_set,
    formula = weight~0+year+depth+temperature+height+(height|ID_TOW),
    # family=Gamma(link = "log"),
    family = gaussian(link = "log"),
    na.action = na.omit)

save(file = "appendix/model_predict/fit.glmm.dt.rda", fit.GLMM.DT)


# ------------------------
# spatiotemporal model

all_set <- data.all %>% 
    mutate(height = as.numeric(scale(log(height)))) %>%
    mutate(X=lon, Y=lat) %>%
    `attr<-`("projection", "LL") %>%
    `attr<-`("zone", "20") %>%
    PBSmapping::convUL() 

# create mesh: use UTM
mesh = inla.mesh.create(all_set[,c("X","Y")], refine = F, extend = F)

# Create matrices for approximation
spde <- inla.spde2.matern(mesh, alpha = 2)

# TMB
version <- "STM"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))

data = list(
    weight = all_set$weight, 
    varmat = as.matrix(cbind(1, all_set$height)),
    ind_loc = all_set$ID_TOW,
    ind_year = all_set$year,
    M0 = spde$param.inla$M0,
    M1 = spde$param.inla$M1,
    M2 = spde$param.inla$M2)

nyear = nlevels(data$ind_year)
nloc = nrow(data$M0)
nvar = ncol(data$varmat)
nobs = nrow(data$varmat)

parameters = list(
    beta = rep(0, nvar),
    beta_s = matrix(0, nloc, nvar),
    beta_t = matrix(2.5, nyear, nvar),
    beta_st = array(0, dim = c(nloc, nyear, nvar)),
    log_kappa_s = rep(-2, nvar),
    log_tau_s = rep(2, nvar),
    log_kappa_st = matrix(0, 1, nvar),
    log_tau_st = matrix(1, 1, nvar),
    log_phi = 1)

maps <- list(
    beta = factor(rep(NA, nvar))
)

obj = MakeADFun(data=data, 
                parameters=parameters,
                map=maps,
                random=c("beta_s","beta_st"),
                DLL=version,
                silent = T)

opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- obj$report()
sd.rep <- sdreport(obj)
summary(sd.rep, "fixed")

fit.STM <- list(obj=obj,opt=opt,rep=rep,sd.rep=sd.rep)

save(file = "appendix/model_predict/fit.stm.rda", fit.STM)

rm("all_set","mesh","spde","data","parameters","maps","obj","opt","rep")


# ------------------------
# spatiotemporal model with depth and temperature

all_set <- data.all %>% 
    mutate(height = as.numeric(scale(log(height))),
           depth = as.numeric(scale(log(depth))),
           temperature = as.numeric(scale(log(temperature)))) %>%
    mutate(X=lon, Y=lat) %>%
    `attr<-`("projection", "LL") %>%
    `attr<-`("zone", "20") %>%
    PBSmapping::convUL() 

# create mesh: use UTM
mesh = inla.mesh.create(all_set[,c("X","Y")], refine = F, extend = F)

# Create matrices for approximation
spde <- inla.spde2.matern(mesh, alpha = 2)


# TMB
version <- "STM_DT"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))

data = list(
    weight = all_set$weight, 
    varmat = as.matrix(cbind(1, all_set$height)),
    depth = all_set$depth,
    temperature = all_set$temperature,
    ind_loc = all_set$ID_TOW,
    ind_year = all_set$year,
    M0 = spde$param.inla$M0,
    M1 = spde$param.inla$M1,
    M2 = spde$param.inla$M2)

nyear = nlevels(data$ind_year)
nloc = nrow(data$M0)
nvar = ncol(data$varmat)
nobs = nrow(data$varmat)

parameters = list(
    beta = rep(0, nvar),
    beta_s = matrix(0, nloc, nvar),
    beta_t = matrix(0, nyear, nvar),
    beta_st = array(0, dim = c(nloc, nyear, nvar)),
    beta_depth = matrix(0, nyear, nvar),
    beta_temperature = matrix(0, nyear, nvar),
    log_kappa_s = rep(log(0.1), nvar),
    log_tau_s = rep(log(0.1), nvar),
    log_kappa_st = matrix(log(0.1), 1, nvar),
    log_tau_st = matrix(log(0.1), 1, nvar),
    log_phi = log(0.1))

# one parameter for each of depth and temperature
maps <- list(
    beta = factor(rep(NA, nvar)),
    beta_depth = factor(cbind(rep(1,nyear),rep(NA,nyear))),
    beta_temperature = factor(cbind(rep(1,nyear),rep(NA,nyear)))
)

obj = MakeADFun(data=data, 
                parameters=parameters,
                map=maps,
                random=c("beta_s","beta_st"),
                DLL=version,
                silent = F)

opt <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=1000, iter.max=1000))
rep <- obj$report()
sd.rep <- sdreport(obj)

fit.STM.DT <- list(obj=obj,opt=opt,rep=rep,sd.rep=sd.rep)

save(file = "appendix/model_predict/fit.stm.dt.rda", fit.STM.DT)

rm("all_set","mesh","spde","data","parameters","maps","obj","opt","rep")





# ------------------------
# GLMM with area effect


all_set <-  data.all %>% 
    mutate(height = as.numeric(scale(log(height))),
           depth = as.numeric(scale(log(depth))),
           temperature = as.numeric(scale(log(temperature))))

fit.GLMM.DTA <- lme4::glmer(
    data = all_set,
    formula = weight~0+year+depth+temperature+height+(height|ID_TOW)+area,
    family = gaussian(link = "log"),
    na.action = na.omit)

summary(fit.GLMM.DTA)

save(file = "appendix/model_predict/fit.glmm.dta.rda",fit.GLMM.DTA)

newdata$`GLMM-DT` <- predict(fit.GLMM.DT, newdata, re.form=NULL, type="response", allow.new.levels = T)

rm("all_set")


# plot subs-areas
ggplot() +
    borders("world", colour="gray50", fill = "gray90", xlim = c(-67,-65), ylim = c(43,45)) +
    coord_map(xlim = c(-67,-65.75), ylim = c(43.5, 44.75)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    geom_point(data=data.all,aes(lon,lat,color=area,shape=area)) + 
    scale_color_discrete() +
    xlab("Longitude") + ylab("Latitude") + 
    labs(color = "Sub-Area", shape = "Sub-Area")
ggsave(filename = "appendix/area_effect.jpg", width = 6, height = 8)



# ------------------------
# AIC and BIC
# GLMM
# AIC      BIC   logLik deviance df.resid 
# 79962.9  80053.4 -39969.4  79938.9    13999 

# GLMM DT
# AIC      BIC   logLik deviance df.resid 
# 79669.9  79775.5 -39820.9  79641.9    13997 

# GLMM-DTA
# AIC      BIC   logLik deviance df.resid 
# 78729.7  78850.4 -39348.8  78697.7    13995 

# STM
# nll
# 38805.18

# STM DT
# nll
# 38629.44


# BIC
# -2 * LL + log(N) * k


# ------------------------
# tabulate parameters
load("appendix/model_predict/fit.glmm.rda")
load("appendix/model_predict/fit.glmm.dt.rda")
load("appendix/model_predict/fit.stm.rda")
load("appendix/model_predict/fit.stm.dt.rda")

# sjPlot::tab_model(fit.GLMM,fit.GLMM.DT)



summary(fit.STM.DT$sd.rep, "random")
summary(fit.STM.DT$sd.rep, "fixed", p.value = TRUE)

# year effect, with CI
data.frame(
    Year = levels(all_set$year),
    GLMM = paste0(round(fit.GLMM@beta[1:nlevels(all_set$year)],2), " (", round(summary(fit.GLMM)$coefficients[1:nlevels(all_set$year),"Std. Error"],2), ")"), 
    GLMM.DT = paste0(round(fit.GLMM.DT@beta[1:nlevels(all_set$year)],2), " (", round(summary(fit.GLMM.DT)$coefficients[1:nlevels(all_set$year),"Std. Error"],2), ")"),
    STM = paste0(round(fit.STM$rep$beta_t[,1],2), " (", round(summary(fit.STM$sd.rep, "fixed")[1:nlevels(all_set$year),"Std. Error"],2), ")"),
    STM.DT = paste0(round(fit.STM.DT$rep$beta_t[,1],2), " (", round(summary(fit.STM.DT$sd.rep, "fixed")[1:nlevels(all_set$year),"Std. Error"],2), ")")
) %>% 
    xtable::xtable() %>%
    print(include.rownames=F)


# environmental var
# depth and temperature: p-value

# GLMM.DT
# Estimate  Std. Error    t value     Pr(>|z|)
# depth       -0.07352694 0.025820073  -2.847666 4.404113e-03
# temperature  0.16360482 0.033625931   4.865436 1.142046e-06

# STM.DT
# Estimate  Std. Error     z value    Pr(>|z^2|)
# beta_depth       -0.05321325 0.023571998  -2.2574774  2.397826e-02
# beta_temperature  0.19035525 0.030741805   6.1920649  5.938108e-10

# random effect: parameter estimates

# response variance
# # GLMM: 
# phi^2 = 11.98
# 
# # GLMM DT
# phi^2 = 11.97

# # STM: 
# phi^2 = 12.73685
# 
# log_kappa_s   0.63039735 0.29661458   2.1253081  3.356092e-02
# log_kappa_s  -4.45118749 6.12010084  -0.7273062  4.670384e-01
# log_tau_s    -0.49053584 0.41513158  -1.1816394  2.373488e-01
# log_tau_s    11.69933418 8.89049942   1.3159367  1.881953e-01
# log_kappa_st -1.51293085 0.14409538 -10.4995094  8.683019e-26
# log_kappa_st -0.08503343 0.27146206  -0.3132424  7.540965e-01
# log_tau_st    1.64581645 0.16313664  10.0885766  6.206290e-24
# log_tau_st    0.67635259 0.35312402   1.9153400  5.544915e-02

# 
# # STM DT
# phi^2 = 12.82726
# 
# log_kappa_s      -3.08971725 1.042711476  -2.9631565  3.045017e-03
# log_kappa_s       0.92283225 0.238588591   3.8678809  1.097853e-04
# log_tau_s         5.23460730 1.519556579   3.4448255  5.714281e-04
# log_tau_s        -0.32473696 0.377033983  -0.8612936  3.890764e-01
# log_kappa_st     -0.55702939 0.119524356  -4.6603840  3.156201e-06
# log_kappa_st     -1.57786907 0.242185427  -6.5151281  7.262767e-11
# log_tau_st        0.96390926 0.135950437   7.0901520  1.339648e-12
# log_tau_st        2.65752890 0.294106192   9.0359502  1.625837e-19


# plot random intercept and slope







