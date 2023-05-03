# scaled path models with single biodiversity measures
## scaled models ----
model_shannon_scaled <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is_scaled ~ 1+prod_scaled + bioH_shannon_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod_scaled ~ 1+habH_scaled + bioH_shannon_scaled+ sediment_size_scaled

# OTOH, biodiversity is influenced by habitat Diversity and sediment size

bioH_shannon_scaled ~ 1+habH_scaled + sediment_size_scaled

sediment_size_scaled ~ 1+habH_scaled
'

model_simpson_scaled <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is_scaled ~ 1+prod_scaled + bioH_simpson_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod_scaled ~ 1+habH_scaled + bioH_simpson_scaled+ sediment_size_scaled

# OTOH, biodiversity is influenced by habitat Diversity

bioH_simpson_scaled ~ 1+habH_scaled+ sediment_size_scaled

sediment_size_scaled ~ 1+habH_scaled
'

model_evenness_scaled <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is_scaled ~ 1+prod_scaled + bioH_evenness_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod_scaled ~ 1+habH_scaled + bioH_evenness_scaled+ sediment_size_scaled

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_evenness_scaled ~ 1+habH_scaled+ sediment_size_scaled

sediment_size_scaled ~ 1+habH_scaled
'

model_richness_scaled <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is_scaled ~ 1+prod_scaled + bioH_richness_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod_scaled ~ 1+habH_scaled + bioH_richness_scaled+ sediment_size_scaled

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_richness_scaled ~ 1+habH_scaled+ sediment_size_scaled

sediment_size_scaled ~ 1+habH_scaled
'

## SEM with latent variables for biodiversity
model_latentBio_scaled <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is_scaled ~ 1+prod_scaled + bioH_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod_scaled ~ 1+habH_scaled + bioH_scaled+ sediment_size_scaled

# OTOH, biodiversity is influenced by habitat Diversity and sediment size

bioH_scaled ~ habH_scaled + sediment_size_scaled

sediment_size_scaled ~ habH_scaled

## identify the latent variables

bioH_scaled =~ bioH_shannon_scaled + bioH_richness_scaled + bioH_evenness_scaled
'

model_latentBio_scaled_2 <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is_scaled ~ prod_scaled + bioH_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod_scaled ~ habH_scaled + sediment_size_scaled

# OTOH, biodiversity is influenced by habitat Diversity and sediment size

bioH_scaled ~ habH_scaled + + prod_scaled + sediment_size_scaled

sediment_size_scaled ~ habH_scaled

## identify the latent variables

bioH_scaled =~ bioH_shannon_scaled + bioH_richness_scaled + bioH_evenness_scaled
'

 
## effects models ----
model_simpson_effects <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is ~ b * prod_scaled + d * bioH_simpson_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ a * habH_scaled + e * bioH_simpson_scaled+ sediment_size_scaled

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_simpson ~ c * habH_scaled + f * sediment_size_scaled

# calculating the effects of habitat diversity through biodiversity
habH_TIE := a*b*c*d*e
bioH_TE := d*b*e
'


