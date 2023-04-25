# raw models ----
model_shannon <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is ~ 1 + prod + bioH_shannon

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ 1 + habH + bioH_shannon

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_shannon ~ 1 + habH

# # estimating the variance of habitat diversity
# 
# habH ~~ habH
# 
# # estimating the residual variances
# 
# bioH ~~ bioH
# prod ~~ prod

# estimating the covariance for residuals

#bioH ~~ prod
'

model_simpson <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is ~ 1 + prod + bioH_simpson

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ 1 + habH + bioH_simpson

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_simpson ~ 1 + habH
'

model_evenness <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is ~ 1 + prod + bioH_evenness

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ 1 + habH + bioH_evenness

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_evenness ~ 1 + habH
'

# scaled models ----
model_shannon_scaled <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is_scaled ~ 1 + prod_scaled + bioH_shannon_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod_scaled ~ 1 + habH_scaled + bioH_shannon_scaled+ sediment_size_scaled

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_shannon_scaled ~ 1 + habH_scaled + sediment_size_scaled

# # estimating the variance of habitat diversity
# 
# habH ~~ habH
# 
# # estimating the residual variances
# 
# bioH ~~ bioH
# prod ~~ prod

# estimating the covariance for residuals

#bioH ~~ prod
'

model_simpson_scaled <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is_scaled ~ 1 + prod_scaled + bioH_simpson_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod_scaled ~ 1 + habH_scaled + bioH_simpson_scaled+ sediment_size_scaled

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_simpson_scaled ~ 1 + habH_scaled+ sediment_size_scaled
'

model_evenness_scaled <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is_scaled ~ 1 + prod_scaled + bioH_evenness_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod_scaled ~ 1 + habH_scaled + bioH_evenness_scaled+ sediment_size_scaled

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_evenness_scaled ~ 1 + habH_scaled+ sediment_size_scaled
'

model_richness_scaled <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is_scaled ~ 1 + prod_scaled + bioH_richness_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod_scaled ~ 1 + habH_scaled + bioH_richness_scaled+ sediment_size_scaled

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_richness_scaled ~ 1 + habH_scaled+ sediment_size_scaled
'
## effects models ----
model_simpson_effects <-'
# set out the model structure
# the full model of interaction strengths
# IS is influenced by production and diversity

is ~ 1 + b * prod_scaled + d * bioH_simpson_scaled

# here prod itself has dependence on habitat Diversity
# and biodiversity

prod ~ 1 + a * habH_scaled + e * bioH_simpson_scaled+ sediment_size_scaled

# OTOH, biodiversity is only influenced by habitat Diversity

bioH_simpson ~ 1 + c * habH_scaled + f * sediment_size_scaled

# calculating the effects of habitat diversity through biodiversity
habH_TIE := a*b*c*d*e
bioH_TE := d*b*e
'


