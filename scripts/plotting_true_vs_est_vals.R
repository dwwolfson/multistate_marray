# November 11, 2025
# The point of this script is to get the code together to run model scenarios 1 and 2, which were written 
# by Todd shortly after the arctic goose conference.

# Model 1 is very similar to the model I presented on (a Burnham joint live-encounter model), 
# where only 2 age classes can be recognizedat time of banding. An update to this version is that
# the pi parameter (proportion of marked AHY that are SY at time of banding) has annual variation instead of being fixed over the whole 20 yrs.

# Model 2 is the same as Model 1 except with a stage-based Leftkovich formulation instead of age classes.

# The JAGS code for model 1 is also in a txt file at "~/scripts/burnham_marray.Rmd"
# The JAGS code for model 2 is also in a txt file at "~/scripts/burnham_lefkovitch_marray.Rmd"


# Libraries used for both scenarios
library(IPMbook)
library(jagsUI)
library(here)


