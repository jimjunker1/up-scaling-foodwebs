library(tidyverse)
library(lavaan)
library(piecewiseSEM)
library(mgcv)
library(MASS)
library(smwrStats)
remotes::install_github("USGS-R/smwrStats")
.pardefault <- par()

## load the data ----
dat = read.csv("./data/divers_prod_hab.csv") %>% 
  dplyr::mutate(production_sum_mean = (production_sum_mean/1000),
                seg_weighted_size = (seg_weighted_size/10),
                log_prod = log10(production_sum_mean)) %>% 
  dplyr::mutate(prod_2 = production_sum_mean^2,
                prod_3 = production_sum_mean^3,
                habH_2 = shannon_habitat^2,
                habH_3 = shannon_habitat^3,
                is10 = median_nolog*10) %>% 
  dplyr::rename(prod = 'production_sum_mean', habH = 'shannon_habitat',
                bioH_richness = 'RICHNESS', bioH_shannon = 'SHANNON_B',
                bioH_simpson = 'SIMPSON_B', bioH_evenness = 'EVENNESS_B', is = 'median_nolog', sediment_size = 'seg_weighted_size') %>% 
  # scale the variables because production is on a vastly different scale than rest
  dplyr::mutate(across(c(prod, prod_2, prod_3, habH, habH_2, habH_3, bioH_richness, bioH_shannon, bioH_simpson, bioH_evenness, sediment_size, is), ~as.numeric(scale(.x)), .names = "{.col}_scaled")) %>%
  dplyr::mutate(d0 = bioH_richness,
                d1 = exp(bioH_shannon),
                d2 = 1/bioH_simpson) %>% 
  dplyr::select(log_prod, prod,prod_2,prod_3, prod_scaled, prod_2_scaled, prod_3_scaled, habH, habH_2, habH_2_scaled, habH_3, habH_scaled, habH_2_scaled, bioH_richness, bioH_richness_scaled, bioH_shannon, bioH_shannon_scaled, bioH_simpson, bioH_simpson_scaled, bioH_evenness, bioH_evenness_scaled,sediment_size, sediment_size_scaled, is, is10, is_scaled, d0, d1, d2)

dat_mod = dat %>%
  dplyr::select(is,prod,bioH_shannon, bioH_simpson, bioH_richness, habH, sediment_size) %>% 
  dplyr::mutate(across(everything(), ~.x + 0.001))

# SEMs ----
# fit with full path analysis/sems
mvnormtest::mshapiro.test(t(dat[,c('is_scaled','log_prod','bioH_shannon','habH','sediment_size')]))
smwrStats::optimBoxCox(dat_mod[,c('is','prod','bioH_shannon','bioH_simpson','bioH_richness','habH','sediment_size')])
dat_mod = dat_mod %>% 
  dplyr::mutate(is = ((is^0.9)-1)/0.9,
                prod = ((prod^-1.1)-1)/-1.1,
                bioH_shannon = log(bioH_shannon),
                bioH_simpson = ((bioH_simpson^0.9)-1)/0.9,
                bioH_richness = ((bioH_richness^0.8)-1)/0.8,
                habH = sqrt(habH),
                sediment_size = ((sediment_size^0.3)-1)/0.3)

mvnormtest::mshapiro.test(t(dat_mod[,c('is','prod','bioH_shannon','bioH_simpson','bioH_richness','habH','sediment_size')]))
model_sem_shannon <-'
# the full model of interaction strengths
# is10 is IS * 10 to get variance on similar scales influenced by production and diversity

is_scaled ~ log_prod + bioH_shannon

# here prod itself has dependence on habitat Diversity
# and biodiversity

log_prod ~ habH + bioH_shannon + sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity and sediment size

bioH_shannon ~ habH + sediment_size

sediment_size ~ habH
'
fit_sem_bioH_shannon <- sem(model_sem_shannon, data = dat, auto.var = TRUE)
summary(fit_sem_bioH_shannon,fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
standardizedsolution(fit_sem_bioH_shannon)

fit_psem_bioH_shannon = psem(
  lm(is10 ~ log_prod + bioH_shannon, dat),
  lm(log_prod ~ habH + bioH_shannon + sediment_size, dat),
  lm(bioH_shannon ~ habH + sediment_size, dat),
  lm(sediment_size ~ habH, dat)
)
summary(fit_psem_bioH_shannon)
AIC(fit_psem_bioH_shannon, AIC.type = "dsep")

model_latentBio_simple <-'
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is_scaled ~ log_prod + bioH

# here prod itself has dependence on habitat Diversity
# and biodiversity

log_prod ~ habH + bioH + sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity and sediment size

bioH ~ habH + sediment_size

sediment_size ~ habH

## identify the latent variables

bioH =~ bioH_shannon_scaled + bioH_richness_scaled + bioH_simpson_scaled

## covariances
bioH_shannon_scaled ~~ bioH_richness_scaled
bioH_richness_scaled ~~ bioH_simpson_scaled
'
fit_sem_bioH_simple <- sem(model_latentBio_simple, data = dat, auto.var = TRUE)
summary(fit_sem_bioH_simple, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
standardizedsolution(fit_sem_bioH_simple)
### models with latent variable for biodiversity (richness, shannon, evenness)
# production composite f
dat_mod = dat[-which(dat$prod > 5),]
is_prod_mod = lm(is_scaled ~ prod_scaled + prod_2_scaled, data = dat)
qqnorm(resid(is_prod_mod));shapiro.test(resid(is_prod_mod))
summary(is_prod_mod)$coefficients
# prod beta: -0.092290891
# prod_2 beta: 0.031447172
# prod_3 beta: -0.003274193
# prod_scaled beta: -0.10343031
# prod_2_scaled beta: 0.17780430
# prod_3_scaled beta: -0.09021509
# is_scale~prod_scaled beta: -6.79282355
# is_scale~prod_2_scaled beta: 12.87363157
# is_scale~prod_3_scaled beta: -7.30716041
prod_composite = -1.0343031 * dat$prod_scaled + 1.7780430 * dat$prod_2_scaled + -0.9021509 * dat$prod_3_scaled

## habitat composite
prod_habH_mod = lm(prod_composite ~ habH_scaled + habH_2_scaled, data = data.frame(dat,prod_composite))
plot(prod_habH_mod)
shapiro.test(resid(prod_habH_mod))
summary(prod_habH_mod)$coefficients
# habH beta: -0.12291289
# habH_2 beta: 0.08859962
# habH_scaled beta: -0.032259388
# habH_2_scaled beta: 0.018491807
# is_scale~habH_scaled beta: -1.79139440
#is_scale~habH_2_Scaled beta: 1.02686784
habH_composite = -0.12684885 * dat$habH + 0.09386793 * dat$habH_2

# more structure 
model_latentBio <-'
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is_scaled ~ prod_composite + bioH

# production has a non-linear response so we can create a composite variable

prod_composite <~ -0.12232501 * prod_scaled + 0.23182805 * prod_2_scaled + -0.13158717 * prod_3_scaled

# prod_composite <~ -6.79282355 * prod_scaled + 12.87363157 * prod_2_scaled + -7.30716041 * prod_3_scaled
# habH has a non-linear response with production so another composite is necessary

habH_composite <~ -0.032259388 * habH_scaled + 0.018491807 * habH_2_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

# prod_composite ~ habH_scaled + bioH + sediment_size
prod_composite ~ habH_composite + bioH + sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity and sediment size

# bioH ~ habH_scaled + sediment_size
bioH ~ habH_composite + sediment_size

# sediment_size ~ habH_scaled
sediment_size ~ habH_composite

## identify the latent variables

bioH =~ bioH_shannon_scaled + bioH_richness_scaled + bioH_evenness_scaled

## covariances
bioH_shannon_scaled ~~ bioH_richness_scaled
bioH_shannon_scaled ~~ bioH_evenness_scaled
bioH_richness_scaled ~~ bioH_evenness_scaled
'
fit_sem_bioH <- sem(model_latentBio, data = dat, auto.var = TRUE, fixed.x = FALSE)#, estimator = 'MLM', se = 'robust')
varTable(fit_sem_bioH)
modificationindices(fit_sem_bioH)
summary(fit_sem_bioH,fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
standardizedsolution(fit_sem_bioH)
# quick look at AIC for comm diversity measures
# really not that interesting just exploratory
# next, fit different model structures 
AIC(fit_sem_bioH_scaled, fit_sem_bioH_scaled_2)

# Get estimates of confidence intervals
parameterestimates(fit_sem_bioH_scaled, standardized = TRUE)

# piecewiseSEM ----
# using piecewiseSEM to do SEM allowing for different data distributions!

# Here it will make sense to look at the data and think about the data generation processes and error structures
## is

is_prod_mod = lm(is_scaled ~ prod_scaled + prod_2_scaled + prod_3_scaled, data = dat)
summary(is_prod_mod)
prod_composite = -5.744e+00 * dat$prod_scaled + 9.874e+00 * dat$prod_2_scaled + -5.010e+00 * dat$prod_3_scaled
dat$prod_composite <- prod_composite


## prod
plot(dat$prod_composite, y = dat$is_scaled)
plot(dat$prod_composite, x = (dat$prod))
plot(dat$prod_composite, x = dat$sediment_size)# non-linear issues
plot(dat$prod_composite, x = dat$habH)
# check residuals from individual glms
prod_lm = lm(prod_composite ~ habH + bioH_shannon + sediment_size, data = dat)
shapiro.test(resid(prod_lm));hist(resid(prod_lm));summary(prod_lm);qqnorm(resid(prod_lm))

## bioH
plot(dat$bioH_shannon, x = dat$habH)
plot(dat$bioH_shannon, x = dat$sediment_size)# non-linear issues
# check residuals from individual glms
bioH_lm = lm(bioH_shannon ~ habH + sediment_size, data = dat)
shapiro.test(resid(bioH_lm));hist(resid(bioH_lm));summary(bioH_lm);qqnorm(resid(bioH_glm))

# fitting piecewise models
model_psem_shannon <- psem(
  lm(is_scaled~ prod_composite + bioH_shannon_scaled,  data = dat),
  lm(prod_composite ~ habH_scaled + bioH_shannon_scaled + sediment_size_scaled, data = dat),
  lm(bioH_shannon_scaled ~ habH_scaled + sediment_size_scaled, data = dat),
  lm(sediment_size_scaled ~ habH_scaled, data = dat)
);
basisSet(model_psem_shannon)
dSep(model_psem_shannon)
summary(model_psem_shannon)

LLchisq(model_psem_shannon)

# fitting simpson
is_prod_mod = lm(is_scaled ~ prod_scaled + prod_2_scaled + prod_3_scaled, data = dat)
summary(is_prod_mod)
prod_composite = -5.744e+00 * dat$prod_scaled + 9.874e+00 * dat$prod_2_scaled + -5.010e+00 * dat$prod_3_scaled
dat$prod_composite <- prod_composite


## prod
plot(dat$prod_composite, y = dat$is_scaled)
plot(dat$prod_composite, x = (dat$prod))
plot(dat$prod_composite, x = dat$sediment_size)# non-linear issues
plot(dat$prod_composite, x = dat$habH)
# check residuals from individual glms
prod_lm = lm(prod_composite ~ habH + bioH_simpson + sediment_size, data = dat)
shapiro.test(resid(prod_lm));hist(resid(prod_lm));summary(prod_lm);qqnorm(resid(prod_lm))

## bioH
plot(dat$bioH_simpson, x = dat$habH)
plot(dat$bioH_simpson, x = dat$sediment_size)# non-linear issues
# check residuals from individual glms
bioH_lm = lm(bioH_simpson ~ habH + sediment_size, data = dat)
shapiro.test(resid(bioH_lm));hist(resid(bioH_lm));summary(bioH_lm);qqnorm(resid(bioH_glm))

# fitting piecewise models
model_psem_simpson <- psem(
  lm(is_scaled~ prod_composite + bioH_simpson,  data = dat),
  lm(prod_composite ~ habH + bioH_simpson + sediment_size, data = dat),
  lm(bioH_simpson ~ habH + sediment_size, data = dat),
  lm(sediment_size ~ habH, data = dat)
);
basisSet(model_psem_simpson)
dSep(model_psem_simpson)
summary(model_psem_simpson)
# Spare(d) code ----
## confirmatory factor analysis ----
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
prod_scaled ~ 1+bioH_simpson_scaled + habH_scaled + sediment_size_scaled
bioH_simpson_scaled ~ 1+habH_scaled +sediment_size_scaled
'
## evenness
model_cfa_evenness <- '
is_scaled ~ 1+bioH_evenness_scaled + prod_scaled
prod_scaled ~ 1+bioH_evenness_scaled + habH_scaled + sediment_size_scaled
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





## attempts to visualize...they all pretty much suck

library(semPlot)
semPaths(
  object = fit_sem_simpson,
  what = "path",
  whatLabels = "par",
  rotation = 2
)



# read in the path analysis model syntax
source("habitat-model-specs.R")
## scaled SEMs ----
### these converge fine, but still show "poor" model fits CFI < 0.9
fit_sem_richness <- sem(model_richness, data = dat, auto.var = TRUE)
varTable(fit_sem_richness)
fit_sem_shannon <- sem(model_shannon, data = dat, auto.var = TRUE)
fit_sem_simpson <- sem(model_simpson, data = dat, auto.var = TRUE)
fit_sem_evenness <- sem(model_evenness, data = dat, auto.var = TRUE)
varTable(fit_sem_evenness)
# get the summaries for SEM fits
summary(fit_sem_richness, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_sem_shannon, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_sem_simpson, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
summary(fit_sem_evenness_scaled, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

## modification indices
modificationindices(fit_sem_shannon)
# lhs           op    rhs       mi      epc     sepc.lv sepc.all sepc.nox
# is_scaled     ~ sediment_size  2.629  0.257   0.257    0.132    0.132
# is_scaled     ~          habH 74.103 -3.054  -3.054   -0.777   -3.070
# prod          ~     is_scaled 62.835  1.037   1.037    0.925    0.925
# bioH_shannon  ~     is_scaled 10.811 -0.171  -0.171   -0.500   -0.500
# sediment_size ~     is_scaled 11.062  0.821   0.821    1.596    1.596
# habH          ~     is_scaled 74.103 -1.258  -1.258   -4.945   -4.945

## just focus on Shannon index and add in options to make linear with CFA



