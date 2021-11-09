### Preliminary model comparisons using simple GLM(M) models for length-weight relationships, used as starting place for the development of Spatial Length-weight models.
#  Analysis to justify approaches taken in YYin paper 
#  J.Sameoto, June 2021 
#  Questions of interest to answer: 
#  * Are depth and temperature informative covariates? 
#  * What is appropriate random effects structure, if any?
#  
#
#
############################################## Step 1... load functions and data + tidy things up for model ############################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(bbmle)
library(predictmeans)
options(stringsAsFactors = F)

# Need to set up a variance inflation factor function, these 2 functions from
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

# Along with this helper function...
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

# Read data from Github
dats <- c("https://raw.githubusercontent.com/Dave-Keith/SPA3_Spatial_LWRs/master/Data/WeightHeightAge_Strata22_23_24_clean.csv",
         "https://raw.githubusercontent.com/Dave-Keith/SPA3_Spatial_LWRs/master/Data/tows.2011to2018wTemperature.csv")
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
csvs <- NULL
for(dat in dats) 
{
  download.file(dat,destfile = basename(dat))
  csvs[[dat]] <- read.csv(paste0(getwd(),"/",basename(dat)))
  file.remove(paste0(getwd(),"/",basename(dat)))
} # end for(fun in funs)

data.bio <- csvs[[1]]
data.covar <- csvs[[2]]

# Merge the datasets, remove NA's, data before 2012, and any data with SH < 40 mm, then give everything nicer names.
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
              temperature = BOTTOM_TEMP) 

# Standardize some of the variables for the model.
all_set <- data.all %>% 
    mutate(height = as.numeric(scale(log(height))),
           depth = as.numeric(scale(log(depth))),
           temperature = as.numeric(scale(log(temperature))))

# Plot the standardized depth vs temp 
ggplot(data = all_set, aes(x=depth, y=temperature)) + geom_point()




##########################################   Step 2 - Start data exploration and model analysis. ##############################################################   

####
## Test if depth AND temperature could both be use in modelling
#  see how strong collinearity is; use threshold of VIF>=3 to drop covariates as per Zuur et al. 2010 
#Variance inflation factors are < 3 so we should be able to retain both of these for the models.
corvif(all_set[,c("depth","temperature")])


####
## Test for best random effect structure - following Zuur et al 2009 p 120+ 
#
all_set.weightlogctr <- data.all %>% 
  mutate(height = as.numeric(scale(log(height))),
         depth = as.numeric(scale(log(depth))),
         temperature = as.numeric(scale(log(temperature))),
         weight.log.ctr = as.numeric(scale(log(weight))))

model.gls <- nlme::gls(weight.log.ctr~0+year*height,
    na.action = na.omit, method = "REML", 
    data = all_set.weightlogctr)
summary(model.gls)

LME.1 <- nlme::lme(weight.log.ctr~ 0+year*height, random = ~ 1|ID_TOW,
                   data = all_set.weightlogctr,
  na.action = na.omit, method = "REML")
summary(LME.1)

LME.2 <- nlme::lme(weight.log.ctr~ 0+year*height, random = ~ height|ID_TOW,
                   data = all_set.weightlogctr,
                   na.action = na.omit, method = "REML")
summary(LME.2)

# LME.2 looks best of these
AIC(model.gls,LME.1,LME.2)
anova(model.gls,LME.1,LME.2) 

# Diagnostics plots 
 residplot(model.gls)
 residplot(LME.1)
 residplot(LME.2)
 
# Conclude random intercept and slope is best random effect structure for GLMM model 
 

########
#  Next Test Random effect model structure when model has all covariates, max interaction, note moving to lmer package from here out
# Run the same test with lme4::lmer 

#random intercept 
GLMM.DT.1.lmer <- lme4::lmer(
  data = all_set.weightlogctr,
  REML = TRUE,
  formula = weight.log.ctr~ 0+year*height*depth*temperature+(1|ID_TOW),
  na.action = na.omit)

summary(GLMM.DT.1.lmer)

# These residuals appear acceptable if not brilliant
plot(GLMM.DT.1.lmer)
plot(GLMM.DT.1.lmer,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1)

# But Not particularly lovely residuals here...
lattice::qqmath(GLMM.DT.1.lmer)
plot(GLMM.DT.1.lmer, rstudent(.) ~ hatvalues(.))

#random slope and intercept 
GLMM.DT.2.lmer <- lme4::lmer(
  data = all_set.weightlogctr,
  REML = TRUE,
  formula = weight.log.ctr~ 0+year*height*depth*temperature+(height|ID_TOW),
  na.action = na.omit)

summary(GLMM.DT.2.lmer)

# This look decent but not brilliant
plot(GLMM.DT.2.lmer)
plot(GLMM.DT.2.lmer,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1)

# But Again pretty ugly
lattice::qqmath(GLMM.DT.2.lmer)
plot(GLMM.DT.2.lmer, rstudent(.) ~ hatvalues(.))

#compare models 
# When use all fixed effects in model, model with random slope and random intercept is still the best random effects structure 
# Note that below in the fixed effect comparisons we only look at additive models
AIC(GLMM.DT.1.lmer, GLMM.DT.2.lmer)


############ 
## Using best random effect structure we can compare a number of additive models that include year, temperature, and depth
# Using gaussian family ie. using lmer  because of convergence issues with Gamma & glmer (this was tested outside this script) 
# Build from simple models out 

#1. Simplest Model, height only 
GLMM.DT.lmer.1 <- lme4::lmer(
  data = all_set.weightlogctr,
  REML = FALSE,
  formula = weight.log.ctr~ 0+height+ (height|ID_TOW),
  na.action = na.omit)

summary(GLMM.DT.lmer.1)

#2. height and year, does relationship change each year  
GLMM.DT.lmer.2 <- lme4::lmer(
  data = all_set.weightlogctr,
  REML = FALSE,
  formula = weight.log.ctr~ 0+height+year + (height|ID_TOW),
  na.action = na.omit)

summary(GLMM.DT.lmer.2)

#3. height and temperature, does temperature affect relationship 
GLMM.DT.lmer.3 <- lme4::lmer(
  data = all_set.weightlogctr,
  REML = FALSE,
  formula = weight.log.ctr~ 0+height+temperature+(height|ID_TOW),
  na.action = na.omit)

summary(GLMM.DT.lmer.3)


#4. height and depth 
GLMM.DT.lmer.4 <- lme4::lmer(
  data = all_set.weightlogctr,
  REML = FALSE,
  formula = weight.log.ctr~ 0+height+depth+(height|ID_TOW),
  na.action = na.omit)

summary(GLMM.DT.lmer.4)


#5. height, year, temperature 
GLMM.DT.lmer.5 <- lme4::lmer(
  data = all_set.weightlogctr,
  REML = FALSE,
  formula = weight.log.ctr~ 0+height+year+temperature+(height|ID_TOW),
  na.action = na.omit)

summary(GLMM.DT.lmer.5)

#6. height, year, depth 
GLMM.DT.lmer.6 <- lme4::lmer(
  data = all_set.weightlogctr,
  REML = FALSE,
  formula = weight.log.ctr~ 0+height+year+depth+(height|ID_TOW),
  na.action = na.omit)

summary(GLMM.DT.lmer.6)

#7. height + depth + temperature ; Maybe the temperature is a surrogate for year and we don't need both
GLMM.DT.lmer.7 <- lme4::lmer(
  data = all_set.weightlogctr,
  REML = FALSE,
  formula = weight.log.ctr~ 0+height+depth+temperature+(height|ID_TOW),
  na.action = na.omit)

summary(GLMM.DT.lmer.7)

#8. height + year + depth + temperature 
GLMM.DT.lmer.8 <- lme4::lmer(
  data = all_set.weightlogctr,
  REML = FALSE,
  formula = weight.log.ctr~ 0+height+year+depth+temperature+(height|ID_TOW),
  na.action = na.omit)

summary(GLMM.DT.lmer.8)

# Here we can compare these models together.... and looks like model 8 is our best model
AICctab(GLMM.DT.lmer.1,GLMM.DT.lmer.2,GLMM.DT.lmer.3,GLMM.DT.lmer.4,
        GLMM.DT.lmer.5,GLMM.DT.lmer.6,GLMM.DT.lmer.7,
        GLMM.DT.lmer.8,logLik=T,weights=T,base=T)


## Model is 8 the best model of the bunch, has a year, depth and temperature effect.  Based on the initial models (2-4) depth is the most informative piece of information.  
# Notice that models 4 and 7 are essentially equivalent, that says temperature is not providing us much new information compared to the depth only model.
# You see something similar with models 5 and 6 again the only difference is the switch between depth and temperature there.  

