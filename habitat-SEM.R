library(tidyverse)
library(lavaan)
library(piecewiseSEM)
library(mgcv)
.pardefault <- par()

## load the data ----
dat = read.csv("./data/divers_prod_hab.csv") %>% 
  dplyr::rename(prod = 'production_sum_mean', habH = 'shannon_habitat',
                bioH_richness = 'RICHNESS', bioH_shannon = 'SHANNON_B',
                bioH_simpson = 'SIMPSON_B', bioH_evenness = 'EVENNESS_B', is = 'median_nolog', sediment_size = 'seg_weighted_size') %>% 
  dplyr::mutate(across(c(prod, habH,bioH_richness, bioH_shannon, bioH_simpson, bioH_evenness, sediment_size, is), ~as.numeric(scale(.x)), .names = "{.col}_scaled")) %>% 
  dplyr::select(prod, prod_scaled, habH, habH_scaled, bioH_richness, bioH_richness_scaled, bioH_shannon, bioH_shannon_scaled, bioH_simpson, bioH_simpson_scaled, bioH_evenness, bioH_evenness_scaled,sediment_size, sediment_size_scaled, is, is_scaled)

# confirmatory factor analysis ----
## richness
model_cfa_richness <- '
is_scaled ~ 1 + bioH_richness_scaled + prod_scaled
prod_scaled ~ 1 + bioH_richness_scaled + habH_scaled +  sediment_size_scaled
bioH_richness_scaled ~ 1 +  habH_scaled + sediment_size_scaled
'
## shannon
model_cfa_shannon <- '
is_scaled ~ 1 + bioH_shannon_scaled + prod_scaled
prod_scaled ~ 1 + bioH_shannon_scaled + habH_scaled + sediment_size_scaled
bioH_shannon_scaled ~ 1 + habH_scaled + sediment_size_scaled
'
## simpson
model_cfa_simpson <- '
is_scaled ~ 1+bioH_simpson_scaled + prod_scaled
prod_scaled ~ 1+bioH_simpson_scaled + habH_scaled
bioH_simpson_scaled ~ 1+habH_scaled +sediment_size_scaled
'
## evenness
model_cfa_evenness <- '
is_scaled ~ 1+bioH_evenness_scaled + prod_scaled
prod_scaled ~ 1+bioH_evenness_scaled + habH_scaled
bioH_evenness_scaled ~ 1+habH_scaled+ sediment_size_scaled
'
# fit cfa
fit_cfa_richness =cfa(model_cfa_richness, data = dat)
fit_cfa_shannon =cfa(model_cfa_shannon, data = dat)
fit_cfa_simpson = cfa(model_cfa_simpson, data = dat)
fit_cfa_evenness = cfa(model_cfa_evenness, data = dat)

# get summarise of cfa fits
summary(fit_cfa_richness, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_cfa_shannon, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_cfa_simpson, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_cfa_evenness, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

# SEMs ----
# fit with full path analysis/sems
# read in the path analysis model syntax
source("habitat-model-specs.R")
## raw path analysis/SEMs ----
## the raw, unscaled models all have issues fitting due to differences in scales in production values vs other variables.
fit_sem_shannon <- sem(model_shannon, data = dat, auto.var = TRUE)
fit_sem_simpson <- sem(model_simpson, data = dat, auto.var = TRUE)
fit_sem_evenness <- sem(model_evenness, data = dat, auto.var = TRUE)

## scaled SEMs ----
### these fit fine
fit_sem_richness_scaled <- sem(model_richness_scaled, data = dat, auto.var = TRUE)
fit_sem_shannon_scaled <- sem(model_shannon_scaled, data = dat, auto.var = TRUE)
fit_sem_simpson_scaled <- sem(model_simpson_scaled, data = dat, auto.var = TRUE)
fit_sem_evenness_scaled <- sem(model_evenness_scaled, data = dat, auto.var = TRUE)
## do some weird stuff to estimate total effects
# I am still working on perfecting this
fit_sem_simpson_eff = sem(model_simpson_effects, data = dat, auto.var = TRUE)

# get the summaries for SEM fits
## no standard error estimates because model variances are on very different scales
summary(fit_sem_shannon, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_sem_simpson, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_sem_evenness, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

##
summary(fit_sem_richness_scaled, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_sem_shannon_scaled, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_sem_simpson_scaled, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_sem_evenness_scaled, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
#effects summary
summary(fit_sem_simpson_eff, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

# quick look at AIC for comm diversity measures
# really not that interesting just exploratory
# next, fit different model structures 
AIC(fit_sem_shannon_scaled, fit_sem_simpson_scaled, fit_sem_evenness_scaled)

# Get estimates of confidence intervals
parameterestimates(fit_sem_simpson_scaled)


# piecewiseSEM ----
# using piecewiseSEM to do CFA allowing for different data distributions!

# Here it will make sense to look at the data and think about the data generation processes and error structures
plot(dat$prod_scaled, dat$bioH_simpson_scaled)
plot(dat$prod_scaled, dat$habH_scaled)# non-linear issues
plot(dat$habH_scaled, dat$bioH_shannon_scaled)
# plot(lm(bioH_shannon_scaled ~ habH_scaled, data = dat))

# fitting piecewise 
model_psem_shannon <- psem(
  glm(is ~ prod + bioH_shannon, family = gaussian(link = 'log'), data = dat),
  glm(prod ~ habH + bioH_shannon, family = gaussian(link = 'log'), data = dat),
  lm(bioH_shannon ~ habH, data = dat),
  data = dat)

summary(model_psem_shannon)

model_psem_shannon_scaled <- psem(
  glm(is_scaled ~ prod_scaled + bioH_shannon_scaled, family = gaussian(link = "identity"), data = dat),
  lm(prod_scaled ~ habH_scaled + bioH_shannon_scaled + sediment_size_scaled, data = dat),
  lm(bioH_shannon_scaled ~ habH_scaled+ sediment_size_scaled, data = dat),
  data = dat)

summary(model_psem_shannon_scaled)

model_psem_simpson_scaled <- psem(
  glm(is_scaled ~ prod_scaled + bioH_simpson_scaled, family = gaussian(link = "identity"), data = dat),
  glm(prod_scaled ~ bioH_simpson_scaled + habH_scaled+ sediment_size_scaled, family = gaussian(link = "identity"), data = dat),
  lm(bioH_simpson_scaled ~ habH_scaled + sediment_size_scaled, data = dat),
  data = dat)

summary(model_psem_simpson_scaled)
dSep(model_psem_simpson_scaled)
# all plotting of these are meh
plot(model_psem_simpson_scaled)


## playing with other model structures e.g., GAMs
model_psem_simpson_scaled_gam <- psem(
  gam(is_scaled ~ s(prod_scaled) + s(bioH_simpson_scaled), family = gaussian(link = "identity"), data = dat),
  gam(prod_scaled ~ s(bioH_simpson_scaled) + s(habH_scaled), family = gaussian(link = "identity"), data = dat),
  gam(bioH_simpson_scaled ~ s(habH_scaled) + s(sediment_size_scaled), family = gaussian(link = "identity"),data = dat),
  data = dat)

summary(model_psem_simpson_scaled_gam)

model_psem_evenness_scaled_gam <- psem(
  glm(is_scaled ~ prod_scaled + bioH_evenness_scaled, family = gaussian(link = "identity"), data = dat),
  glm(prod_scaled ~ habH_scaled + bioH_evenness_scaled, family = gaussian(link = "identity"), data = dat),
  lm(bioH_evenness_scaled ~ habH_scaled, data = dat),
  data = dat)

model_psem_richness_scaled <- psem(
  glm(is_scaled ~ prod_scaled + bioH_richness_scaled, family = gaussian(link = "identity"), data = dat),
  glm(prod_scaled ~ habH_scaled + bioH_richness_scaled, family = gaussian(link = "identity"), data = dat),
  lm(bioH_richness_scaled ~ habH_scaled, data = dat),
  data = dat)
dSep(model_psem_simpson_scaled)
plot(model_psem_simpson_scaled)


## attempts to visualize...they all pretty much suck

library(semPlot)
semPaths(
  object = fit_sem_simpson,
  what = "path",
  whatLabels = "par",
  rotation = 2
)






