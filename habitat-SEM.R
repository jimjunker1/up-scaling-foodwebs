library(knitr)
library(tidyverse)
library(lavaan)
library(MASS)
library(tidySEM)
library(kableExtra)
library(flextable)
## load the data ----
dat = read.csv("./data/divers_prod_hab.csv") %>% 
  dplyr::mutate(production_sum_mean = (production_sum_mean/1000),
                seg_weighted_size = (seg_weighted_size/10),
                log_prod = log10(production_sum_mean),
                is10 = median_nolog*10) %>% 
  dplyr::rename(prod = 'production_sum_mean', habH = 'shannon_habitat',
                bioH_richness = 'RICHNESS', bioH_shannon = 'SHANNON_B',
                bioH_simpson = 'SIMPSON_B', bioH_evenness = 'EVENNESS_B', is = 'median_nolog', sediment_size = 'seg_weighted_size') %>% 
  # scale the variables because production is on a vastly different scale than rest
  dplyr::mutate(across(c(prod, habH, bioH_richness, bioH_shannon, bioH_simpson, bioH_evenness, sediment_size, is), ~as.numeric(scale(.x)), .names = "{.col}_scaled")) %>%
  dplyr::mutate(d0 = bioH_richness,
                d1 = exp(bioH_shannon),
                d2 = 1/bioH_simpson) %>% 
  dplyr::select(log_prod, prod, prod_scaled, habH,habH_scaled, bioH_richness, bioH_richness_scaled, bioH_shannon, bioH_shannon_scaled, bioH_simpson, bioH_simpson_scaled, bioH_evenness, bioH_evenness_scaled,sediment_size, sediment_size_scaled, is, is10, is_scaled, d0, d1, d2)


# SEMs ----
# fit with full path analysis/sems
## shannon
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
fit_sem_bioH_shannon_summ = summary(fit_sem_bioH_shannon,fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
fit_sem_bioH_shannon_stand = standardizedsolution(fit_sem_bioH_shannon)
## richness
model_sem_richness <-'
# the full model of interaction strengths
# is10 is IS * 10 to get variance on similar scales influenced by production and diversity

is_scaled ~ log_prod + bioH_richness

# here prod itself has dependence on habitat Diversity
# and biodiversity

log_prod ~ habH + bioH_richness + sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity and sediment size

bioH_richness ~ habH + sediment_size

sediment_size ~ habH
'
fit_sem_bioH_richness <- sem(model_sem_richness, data = dat, auto.var = TRUE)
fit_sem_bioH_richness_summ = summary(fit_sem_bioH_richness,fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
fit_sem_bioH_richness_stand = standardizedsolution(fit_sem_bioH_richness)
## simpson
model_sem_simpson <-'
# the full model of interaction strengths
# is10 is IS * 10 to get variance on similar scales influenced by production and diversity

is_scaled ~ log_prod + bioH_simpson

# here prod itself has dependence on habitat Diversity
# and biodiversity

log_prod ~ habH + bioH_simpson + sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity and sediment size

bioH_simpson ~ habH + sediment_size

sediment_size ~ habH
'
fit_sem_bioH_simpson <- sem(model_sem_simpson, data = dat, auto.var = TRUE)
fit_sem_bioH_simpson_summ = summary(fit_sem_bioH_simpson,fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
fit_sem_bioH_simpson_stand = standardizedsolution(fit_sem_bioH_simpson)
## latent variable for biodiversity
model_latentBio_simple <-'
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is_scaled ~ g*log_prod + h*bioH

# here prod itself has dependence on habitat Diversity
# and biodiversity

log_prod ~ b * habH + f* bioH + d * sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity and sediment size

bioH ~ c*habH + e*sediment_size

sediment_size ~ a*habH

## identify the latent variables

bioH =~ bioH_shannon_scaled + bioH_richness_scaled + bioH_simpson_scaled

## covariances
bioH_shannon_scaled ~~ bioH_richness_scaled
# bioH_shannon_scaled ~~ bioH_simpson_scaled
bioH_richness_scaled ~~ bioH_simpson_scaled

# calculating the effects of habitat diversity through production
habH_prodTE := (b*g)+(c*f*g)+(a*d*g)+(a*e*f*g)
# calculating the effects of habitat diversity through biodiversity
habH_bioTE := (c*h)+(a*e*h)
# calculating the total effect of habitat diversity mediated through both production and biodiversity
habH_TE := (b*g)+(c*f*g)+(a*d*g)+(a*e*f*g)+(c*h)+(a*e*h)

'
fit_sem_bioH_simple <- sem(model_latentBio_simple, data = dat, auto.var = TRUE)

fit_sem_summ = summary(fit_sem_bioH_simple, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

fit_sem_bioH_stand = standardizedsolution(fit_sem_bioH_simple)

