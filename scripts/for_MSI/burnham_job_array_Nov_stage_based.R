# November 13, 2025
# The point of this script is to get the code together to run model scenarios 1 and 2, which were written 
# by Todd shortly after the arctic goose conference.

# Model 1 is very similar to the model I presented on (a Burnham joint live-encounter model), 
# where only 2 age classes can be recognizedat time of banding. An update to this version is that
# the pi parameter (proportion of marked AHY that are SY at time of banding) has annual variation instead of being fixed over the whole 20 yrs.

# Model 2 is the same as Model 1 except with a stage-based Leftkovich formulation instead of age classes.


# This script will be for model 2

# Libraries used for both scenarios
packages<-c("here", "jagsUI", "readr", "glue", "tibble", "IPMbook")

# load packages
invisible(lapply(packages, library, character.only = TRUE))
#source(here("MSI/simulations_head/job_array_test/scripts/kery_schaub_helper_functions.R"))  #make sure that this is in the here:here directory location
# I'll switch to just using the IPMbook package instead of referencing a separate file with the Kery and Schaub functions written out
options(scipen = 999)


# cat(file = "burnham_lefkovitch_marray.txt", "
# model {
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 1. SPECIFY PRIORS
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # pi is proportion of AHY that are subadults at time of banding
# pi.mu ~ dnorm(0,0.33)
# pi.sd ~ dunif(0.01, 2); pi.tau <- pow(pi.sd, -2)
# # g is proportion of surviving subadults that transition to become adult breeders
# g.mu ~ dnorm(0,0.33)
# g.sd ~ dunif(0.01, 2); g.tau <- pow(g.sd, -2)
# for (t in 1:(n_yrs)){
#   logit.pi[t] ~ dnorm(pi.mu, pi.tau)
#   logit(pi[t]) <- logit.pi[t]
#   logit.g[t] ~ dnorm(g.mu, g.tau)
#   logit(g[t]) <- logit.g[t]
# }
# 
# # age-specific hyperparameters 
# # include code for temporal covariation across age classes?
# for (a in 1:n_ages){
#   s.mu[a] ~ dnorm(0, 0.33) # survival
#   s.sd[a] ~ dunif(0.01, 2); s.tau[a] <- pow(s.sd[a], -2) # annual variation, convert to precision
#   F.mu[a] ~ dnorm(0, 0.33) # fidelity to breeding site
#   F.sd[a] ~ dunif(0.01, 2); F.tau[a] <- pow(F.sd[a], -2)
#   r.mu[a] ~ dnorm(0, 0.33) # seber recovery probability
#   r.sd[a] ~ dunif(0.01, 2); r.tau[a] <- pow(r.sd[a], -2)
#   p.mu[a] ~ dnorm(0, 0.33)  # recapture probability, same site
#   p.sd[a] ~ dunif(0.01, 2); p.tau[a] <- pow(p.sd[a], -2)
#   for (t in 1:(n_yrs)){
#     logit.s[a,t] ~ dnorm(s.mu[a], s.tau[a])
#     logit(s[a,t]) <- logit.s[a,t] # backtransform to real scale 
#     logit.F[a,t] ~ dnorm(F.mu[a], F.tau[a])
#     logit(F[a,t]) <- logit.F[a,t]  
#     logit.r[a,t] ~ dnorm(r.mu[a], r.tau[a])
#     logit(r[a,t]) <- logit.r[a,t]  
#     logit.p[a,t] ~ dnorm(p.mu[a], p.tau[a])
#     logit(p[a,t]) <- logit.p[a,t]  
#   } # t
# } # a
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 2. Multistate state-transition and observation matrix (a vector suffices)
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # alive resident 1: HY, 2: SY, 3: ASY, 4: AHY
# # recovered dead, all ages: 5 (thanks Thomas!)
# # alive emigrant 6: SY, 7: ASY (unobservable states)
# # dead or unobserved (absorbing)
# 
# for (t in 1:(n_yrs-1)){
#   # departure HY release
#    psi[1,t,1] <- 0
#    psi[1,t,2] <- s[1,t]*F[1,t] # survive to subadult, remain
#    psi[1,t,3] <- 0
#    psi[1,t,4] <- 0
#    psi[1,t,5] <- (1-s[1,t])*r[1,t] # die, recovered
#    psi[1,t,6] <- s[1,t]*(1-F[1,t]) # survive, emigrate
#    psi[1,t,7] <- 0
#    psi[1,t,8] <- (1-s[1,t])*(1-r[1,t]) # die, not recovered
#   # departure subadult resident
#    psi[2,t,1] <- 0
#    psi[2,t,2] <- s[2,t]*F[2,t]*(1-g[t]) # persist as resident subadults
#    psi[2,t,3] <- s[2,t]*F[2,t]*g[t] # graduate to resident adults
#    psi[2,t,4] <- 0
#    psi[2,t,5] <- (1-s[2,t])*r[2,t]
#    psi[2,t,6] <- s[2,t]*(1-F[2,t])*(1-g[t]) # persist as dispersed subadults
#    psi[2,t,7] <- s[2,t]*(1-F[2,t])*g[t] # graduate to dispersed adults
#    psi[2,t,8] <- (1-s[2,t])*(1-r[2,t])
#   # departure adult resident
#    psi[3,t,1] <- 0
#    psi[3,t,2] <- 0
#    psi[3,t,3] <- s[3,t]*F[3,t]
#    psi[3,t,4] <- 0
#    psi[3,t,5] <- (1-s[3,t])*r[3,t]
#    psi[3,t,6] <- 0
#    psi[3,t,7] <- s[3,t]*(1-F[3,t])
#    psi[3,t,8] <- (1-s[3,t])*(1-r[3,t])
#   # departure AHY resident (individuals starting as subadult can persist or graduate)
#    psi[4,t,1] <- 0
#    psi[4,t,2] <- pi[t]*s[2,t]*F[2,t]*(1-g[t]) # persist as subadult
#    psi[4,t,3] <- pi[t]*s[2,t]*F[2,t]*g[t] + (1-pi[t])*s[3,t]*F[3,t] # graduate or start as adult
#    psi[4,t,4] <- 0
#    psi[4,t,5] <- pi[t]*(1-s[2,t])*r[2,t] + (1-pi[t])*(1-s[3,t])*r[3,t]
#    psi[4,t,6] <- pi[t]*s[2,t]*(1-F[2,t])*(1-g[t])
#    psi[4,t,7] <- pi[t]*s[2,t]*(1-F[2,t])*g[t] + (1-pi[t])*s[3,t]*(1-F[3,t])
#    psi[4,t,8] <- pi[t]*(1-s[2,t])*(1-r[2,t]) + (1-pi[t])*(1-s[3,t])*(1-r[3,t])
#   # recently dead (simplified coding, all transition to long dead) 
#    psi[5,t,1:8] <- c(0, 0, 0, 0, 0, 0, 0, 1)
#   # departure subadult emigrant
#    psi[6,t,1] <- 0
#    psi[6,t,2] <- 0
#    psi[6,t,3] <- 0
#    psi[6,t,4] <- 0
#    psi[6,t,5] <- (1-s[2,t])*r[2,t]
#    psi[6,t,6] <- s[2,t]*(1-g[t])
#    psi[6,t,7] <- s[2,t]*g[t]
#    psi[6,t,8] <- (1-s[2,t])*(1-r[2,t])
#   # departure adult emigrant
#    psi[7,t,1] <- 0
#    psi[7,t,2] <- 0
#    psi[7,t,3] <- 0
#    psi[7,t,4] <- 0
#    psi[7,t,5] <- (1-s[3,t])*r[3,t]
#    psi[7,t,6] <- 0
#    psi[7,t,7] <- s[3,t]
#    psi[7,t,8] <- (1-s[3,t])*(1-r[3,t])
#   # already dead (simplified coding) 
#    psi[8,t,1:8] <- c(0, 0, 0, 0, 0, 0, 0, 1)
# 
#  # observation matrix 
#    po[1,t] <- 1 # HY releases, all seen (but condition on capture so ignored)
#    po[2,t] <- p[2,t] # captured as SY where last captured
#    po[3,t] <- p[3,t] # captured as ASY where last captured
#    po[4,t] <- 1 # AHY releases
#    po[5,t] <- 1 # dead recoveries, r in state-transition matrix
#    po[6,t] <- 0 # SY emigrants, unobservable
#    po[7,t] <- 0 # ASY emigrant, unobservable
#    po[8,t] <- 0 # dead unrecovered or long dead
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Miche's magic code, never needs to be modified
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#  # Calculate probability of non-encounter (dq) and reshape the array for the encounter probabilities
#       for (k in 1:ns){
#          dp[k,t,k] <- po[k,t]
#          dq[k,t,k] <- 1-po[k,t]} #k
#       for (k in 1:(ns-1)){
#          for (m in (k+1):ns){
#             dp[k,t,m] <- 0
#             dq[k,t,m] <- 0} #m
#       } #k
#       for (k in 2:ns){
#          for (m in 1:(k-1)){
#             dp[k,t,m] <- 0
#             dq[k,t,m] <- 0} #m
#       } #k
#    } #t
# 
# # Define the multinomial likelihood
# for (t in 1:((n_yrs-1)*ns)){
#    marr[t,1:(n_yrs*ns-(ns-1))] ~ dmulti(pri[t,], rel[t])
# } # t
# 
# # Define the cell probabilities of the multistate m-array
# for (t in 1:(n_yrs-2)){
#    U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
#    for (j in (t+1):(n_yrs-1)){
#    ## fixed based on published erratum of book code
#       U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,j-1,] %*% dq[,j-1,] # fixed
#    } #j
# } #t
# U[(n_yrs-2)*ns+(1:ns), (n_yrs-2)*ns+(1:ns)] <- ones
# 
# # Diagonal
# for (t in 1:(n_yrs-2)){
#    pri[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*% dp[,t,]
#    # Above main diagonal
#    for (j in (t+1):(n_yrs-1)){
#       pri[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% psi[,j,] %*% dp[,j,]
#    } #j
# } #t
# pri[(n_yrs-2)*ns+(1:ns), (n_yrs-2)*ns+(1:ns)] <- psi[,n_yrs-1,] %*% dp[,n_yrs-1,]
# 
# # Below main diagonal
# for (t in 2:(n_yrs-1)){
#    for (j in 1:(t-1)){
#       pri[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
#    } #j
# } #t
# 
# # Last column: probability of non-recapture
# for (t in 1:((n_yrs-1)*ns)){
#    pri[t,(n_yrs*ns-(ns-1))] <- 1-sum(pri[t,1:((n_yrs-1)*ns)])
# } #t
# } # end model
# ")

parms2 <- c("s.mu", "s.sd", "F.mu", "F.sd", "r.mu", "r.sd", "p.mu", "p.sd", 
            "pi.mu", "pi.sd", "g.mu", "g.sd",   
            "pi", "g", "s", "F", "r", "p")



##########################################
# static model settings
###########################################


n_yrs <- 20 # length of banding-recovery
n_stages <- 3
n_states <- 8

# set annual releases
HY_banded <- 1000 # reasonable real-world example
AHY_banded <- 1000 
rel_HY <- rep(HY_banded, n_yrs)
rel_AHY <- rep(AHY_banded, n_yrs)
rel <- cbind(rel_HY, rel_AHY) # modify if SY and ASY releases

### simulation code doesn't seem to work if binding columns of all 0s but that could likely be fixed
n_cohorts <- dim(rel)[2]
rel_cohort <- c(1,4) # modify to c(1:4) if releasing known SY and ASY too

###########################################################################
# Start simulation 
##########################################################################

xx<- Sys.getenv('SLURM_ARRAY_TASK_ID')

#############################
# simulate true parameters
############################

# specify age-specific means and variances
# note the code allows biologically unlikely scenarios like juvs surviving better than adults
r.mu <- qlogis(runif(3, 0.05, 0.3))
r.sd <- runif(3, 0.1, 1)
p.mu <- qlogis(runif(3, 0.05, 0.3))
p.sd <- runif(3, 0.1, 1)
s.mu <- qlogis(runif(3, 0.6, 0.9))
s.sd <- runif(3, 0.1, 1)
F.mu <- qlogis(runif(3, 0.7, 0.95))
F.sd <- runif(3, 0.1, 0.5) # allowing less annual variation for F
pi.mu <- qlogis(runif(1, 0.1, 0.3))
pi.sd <- runif(1, 0.1, 1)
pi <- plogis(rnorm(n_yrs, pi.mu, pi.sd))
g.mu <- qlogis(runif(1, 0.25, 0.95))
g.sd <- runif(1, 0.1, 1)
g <- plogis(rnorm(n_yrs, g.mu, g.sd))

# generate matrices of age-specific vital rates
r <- p <- s <- F <- matrix(NA, nrow = n_stages, ncol = n_yrs)
for (a in 1:n_stages){
  r[a,1:n_yrs] <- plogis(rnorm(n_yrs, r.mu[a], r.sd[a]))  
  p[a,1:n_yrs] <- plogis(rnorm(n_yrs, p.mu[a], p.sd[a]))  
  s[a,1:n_yrs] <- plogis(rnorm(n_yrs, s.mu[a], s.sd[a]))  
  F[a,1:n_yrs] <- plogis(rnorm(n_yrs, F.mu[a], F.sd[a]))  
}

####################################
# Store true parameters
###################################

# path for true parameters
# ~/output/true_params/

# store the same hyperparameter for pi mu and sd for all 3 age classes
param_list<-list(F, p, r, s, pi, g,
                 F.mu, F.sd, r.mu, r.sd, p.mu, p.sd, s.mu, s.sd, 
                 rep(pi.mu,3), rep(pi.sd, 3), rep(g.mu, 3), rep(g.sd, 3)) #tripling to keep dimensions same

names(param_list)<-c("F", "p", "r", "s", "pi",
                     "F.mu", "F.sd", "r.mu", "r.sd", "p.mu", "p.sd", "s.mu", "s.sd", 
                     "pi.mu", "pi.sd", "g.mu", "g.sd")

# save to file
saveRDS(param_list, here(glue::glue("MSI/simulations_head/Nov_2025_lefkovitch/output/true_params/param_set_simulation_run_{xx}.rds")))

########################################
# Simulate dataset from true parameters
########################################

# create array of transition probabilities
psi <- array(NA, dim = c(n_states, n_yrs, n_states))
for (i in 1:n_yrs){
  psi[1,i,1:8] <- c(0, s[1,i]*F[1,i], 0, 0, (1-s[1,i])*r[1,i], s[1,i]*(1-F[1,i]), 0, (1-s[1,i])*(1-r[1,i]))
  psi[2,i,1:8] <- c(0, s[2,i]*F[2,i]*(1-g[i]), s[2,i]*F[2,i]*g[i], 0, (1-s[2,i])*r[2,i], 
                    s[2,i]*(1-F[2,i])*(1-g[i]), s[2,i]*(1-F[2,i])*g[i], (1-s[2,i])*(1-r[2,i]))
  psi[3,i,1:8] <- c(0, 0, s[3,i]*F[3,i], 0, (1-s[3,i])*r[3,i], 0, s[3,i]*(1-F[3,i]), (1-s[3,i])*(1-r[3,i]))
  psi[4,i,1:8] <- c(0, pi[i]*s[2,i]*F[2,i]*(1-g[i]), pi[i]*s[2,i]*F[2,i]*g[i] + (1-pi[i])*s[3,i]*F[3,i], 0, 
                    pi[i]*(1-s[2,i])*r[2,i] + (1-pi[i])*(1-s[3,i])*r[3,i], 
                    pi[i]*s[2,i]*(1-F[2,i])*(1-g[i]), pi[i]*s[2,i]*(1-F[2,i])*g[i] + (1-pi[i])*s[3,i]*(1-F[3,i]), 
                    pi[i]*(1-s[2,i])*(1-r[2,i]) + (1-pi[i])*(1-s[3,i])*(1-r[3,i]))
  psi[5,i,1:8] <- c(0, 0, 0, 0, 0, 0, 0, 1)
  psi[6,i,1:8] <- c(0, 0, 0, 0, (1-s[2,i])*r[2,i], s[2,i]*(1-g[i]), s[2,i]*g[i], (1-s[2,i])*(1-r[2,i]))
  psi[7,i,1:8] <- c(0, 0, 0, 0, (1-s[3,i])*r[3,i], 0, s[3,i], (1-s[3,i])*(1-r[3,i]))
  psi[8,i,1:8] <- c(0, 0, 0, 0, 0, 0, 0, 1)
}
# verify all calculations sum to 1
rowSums(psi[1:8,n_yrs,])

# create empty matrix to hold truth and CH data
truth <- CH <- matrix(0, nrow = sum(rel), ncol = n_yrs)

c <- 1 # set row counter to 1
for (i in 1:(n_yrs-1)){
  for (j in 1:n_cohorts){ # release cohorts
    for (k in 1:rel[i,j]){
      truth[c,i] <- CH[c,i] <- rel_cohort[j] 
      for (m in (i+1):n_yrs){
        truth[c,m] <- sample(1:n_states, size = 1, prob = psi[truth[c,m-1], m-1, 1:n_states])
        CH[c,m] <- ifelse(truth[c,m] == 2,
                          truth[c,m]*rbinom(1, 1, p[2, m]), # chance of being seen if SY
                          ifelse(truth[c,m] == 3,
                                 truth[c,m]*rbinom(1, 1, p[3, m]), # chance of being seen if ASY
                                 truth[c,m]))
      } # m
      c <- c + 1
    } # k
  } # j
} # i
# no information in final year except release status so no Brownie's f in final year
for (i in n_yrs:n_yrs){
  for (j in 1:2){ # release cohorts
    for (k in 1:rel[i,j]){
      CH[c,i] <- rel_cohort[j]
      c <- c + 1
    } # k
  } # j
} # i

# convert 6:8 to 0 (unobservable states)
ch <- CH
ch[ch[,] >= 6] <- 0
# inspect
#head(ch)
#tail(ch)

# convert to multistate CH using M Schaub IPM book code
marr2 <- marray(ch, unobs = 3)


####################################
# Save simulated datasets used
###################################

# path for simulated datasets
# ~/output/sim_datasets/

# save truth capture history to file
write_csv(as.data.frame(truth), here(glue::glue("MSI/simulations_head/Nov_2025_lefkovitch/output/sim_datasets/simulated_truth_ch_simulation_run_{xx}.csv")))

# save m-array to file
write_csv(as.data.frame(marr2), here(glue::glue("MSI/simulations_head/Nov_2025_lefkovitch/output/sim_datasets/simulated_marray_simulation_run_{xx}.csv")))

###############################
# Run JAGS model
###############################

# bundle data for jags
jags.data <- list(n_yrs = n_yrs,
                  n_ages = n_stages,
                  ns = n_states, # number of states
                  marr = marr2, 
                  rel = rowSums(marr2),
                  ones = diag(n_states), 
                  zero = matrix(0, ncol = n_states, nrow = n_states))

# MCMC settings
ni <- 15000; nb <- 5000; nt <- 1; nc <- 4; na <- 1000 # (~140 min, all except penultimate converged and well mixed)
#ni <- 7500; nb <- 2500; nt <- 1; nc <- 4; na <- 1000 # (~70 min, most parameters converged)
#ni <- 1500; nb <- 500; nt <- 1; nc <- 4; na <- 1000 # (exploratory, ~15 min, most converge, hypers need more time)
iterations <- ((ni - nb)/nt)*nc

# track model run time
start_time<-Sys.time()

# Call JAGS from R and check convergence
out2 <- jagsUI(jags.data, inits = NULL, parms2, "burnham_lefkovitch_marray.txt", 
               n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, 
               parallel = TRUE)

end_time<-Sys.time()

# run time in minutes
run_min<-as.numeric(difftime(end_time, start_time, units = "min"))

cat("Simulation run number ", xx, "took ", run_min, " minutes to fit the jags model.", "\n",
    file=here("MSI/simulations_head/Nov_2025_lefkovitch/output/model_results/model_runtimes.txt"), append=T)


####################################
# Save JAGS model
###################################

# path for simulated datasets
# ~/output/sim_datasets/

# don't think it's necessary to save fit jags model full output
#saveRDS(out, here(glue::glue("MSI/simulations_head/Nov_2025/output/model_results/fit_model_simulation_run_{xx}.rds")))

# save model summary
df<-as.data.frame(out2$summary)
df<-rownames_to_column(df, var="param")
write_csv(df, here(glue::glue("MSI/simulations_head/Nov_2025_lefkovitch/output/model_results/model_summary_simulation_run_{xx}.csv")))





























