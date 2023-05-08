# scaled path models with single biodiversity measures
## scaled models ----
model_shannon <-'
# set out the model structure
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is_scaled ~ 1+prod + bioH_shannon

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ 1+habH + bioH_shannon+ sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity and sediment size

bioH_shannon ~ 1+habH + sediment_size

sediment_size ~ 1+habH
'

model_simpson <-'
# set out the model structure
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is_scaled ~ 1+prod + bioH_simpson

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ 1+habH + bioH_simpson+ sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity

bioH_simpson ~ 1+habH+ sediment_size

sediment_size ~ 1+habH
'

model_evenness <-'
# set out the model structure
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is_scaled ~ 1+prod + bioH_evenness_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ 1+habH + bioH_evenness_scaled+ sediment_size

# OTOH, biodiversity is_scaled only influenced by habitat Diversity

bioH_evenness_scaled ~ 1+habH+ sediment_size

sediment_size ~ 1+habH
'

model_richness <-'
# set out the model structure
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is_scaled ~ 1+prod + bioH_richness_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ 1+habH + bioH_richness_scaled+ sediment_size

# OTOH, biodiversity is_scaled only influenced by habitat Diversity

bioH_richness_scaled ~ 1+habH+ sediment_size

sediment_size ~ 1+habH
'

## SEM with latent variables for biodiversity

# exercis_scalede to determine the structure of model
model_latentBio_min <- '
# set out the model structure
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is_scaled ~ prod + bioH

# here prod itself has dependence on habitat Diversity + sediment size

prod ~ habH + sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity and sediment size

bioH ~ habH + sediment_size

sediment_size ~ habH

## identify the latent variables

bioH =~ bioH_shannon + bioH_richness_scaled + bioH_evenness_scaled
'
model_latentBio <-'
# set out the model structure
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is ~ log_prod + bioH

# here prod itself has dependence on habitat Diversity
# and biodiversity

log_prod ~ habH + bioH + sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity and sediment size

bioH ~ habH + sediment_size

sediment_size ~ habH

## identify the latent variables

bioH =~ bioH_shannon + bioH_richness_scaled + bioH_evenness_scaled

## covariances
bioH_shannon ~~ bioH_richness_scaled
bioH_shannon ~~ bioH_evenness_scaled
bioH_richness_scaled ~~ bioH_evenness_scaled

'

model_latentBio_2 <-'
# set out the model structure
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is_scaled ~ 1+prod + bioH

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ 1+habH + bioH + sediment_size

# OTOH, biodiversity is_scaled influenced by habitat Diversity and sediment size

bioH ~ habH + sediment_size

sediment_size ~ habH

## identify the latent variables

bioH =~ bioH_shannon + bioH_richness_scaled + bioH_evenness_scaled
'

## effects models ----
model_simpson_effects <-'
# set out the model structure
# the full model of interaction strengths
# is_scaled is_scaled influenced by production and diversity

is_scaled ~ b * prod + d * bioH_simpson

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ a * habH + e * bioH_simpson+ sediment_size

# OTOH, biodiversity is_scaled only influenced by habitat Diversity

bioH_simpson ~ c * habH + f * sediment_size

# calculating the effects of habitat diversity through biodiversity
habH_TIE := a*b*c*d*e
bioH_TE := d*b*e
'


