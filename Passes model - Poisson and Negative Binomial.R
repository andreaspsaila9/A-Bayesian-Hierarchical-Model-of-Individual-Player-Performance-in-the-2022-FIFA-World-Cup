# necessary packages
library(dplyr)
library(rjags)
library(coda)


# Preparing Data
set.seed(88)

# players <- read.csv("C:/Users/psail/OneDrive/Desktop/UOM 4th yr/Thesis/Dataset/Proper Dataset/Players minutes without gks.csv")
# passing <- read.csv("C:/Users/psail/OneDrive/Desktop/UOM 4th yr/Thesis/Dataset/Proper Dataset/player_passingandshootingwithoutgks.csv")
# defense <- read.csv("C:/Users/psail/OneDrive/Desktop/UOM 4th yr/Thesis/Dataset/Proper Dataset/player_defensewithoutgks.csv")

players <- read.csv("C:/Users/psail/OneDrive/Desktop/UOM 4th yr/Thesis/Dataset/Proper Dataset/Players minutes without gks.csv",
                    fileEncoding = "UTF-8")

passing <- read.csv("C:/Users/psail/OneDrive/Desktop/UOM 4th yr/Thesis/Dataset/Proper Dataset/player_passingandshootingwithoutgks.csv",
                    fileEncoding = "UTF-8")

defense <- read.csv("C:/Users/psail/OneDrive/Desktop/UOM 4th yr/Thesis/Dataset/Proper Dataset/player_defensewithoutgks.csv",
                    fileEncoding = "UTF-8")

df <- players %>%
  inner_join(passing,  by = c("player", "team")) %>%
  inner_join(defense, by = c("player", "team")) # merging the 3 excel files into 1


# this code tells us if and who has missing passing data - Ismaila Mohamad
df %>% filter(is.na(passes_completed)) %>% select(player, team, minutes, minutes_90s)

# NOW we remove the player with missing passing data before creating IDs
df_model <- df %>%
  filter(!is.na(passes_completed)) %>%                    # removed player (sample becomes 310 players, need to mention in write up)
  mutate(
    player_id = as.numeric(factor(player)),               # converts player name into numeric ID
    team_id   = as.numeric(factor(team)),                 # converts team name into numeric ID
    tau       = minutes_90s,              # total minutes in 90-min units
    tau_star  = minutes_90s / games,    # avg minutes per match in 90-min units
    minutes_90 = minutes_90s,                             # minutes per 90 for all players since some players played less
    is_DF = ifelse(position.x == "DF", 1, 0),             # dummy for defenders
    is_MF = ifelse(position.x == "MF", 1, 0)              # dummy for midfielders (FW is baseline by elimination)
  )


delta <- as.matrix(df_model[, 7:38]) # teams columns in dataset: from Argentina to Wales

# checks for clean data 
sapply(df_model[, c("passes_completed", "interceptions","shots" ,"minutes_90",
                    "player_id", "team_id", "is_DF", "is_MF")],
       function(x) sum(is.na(x)))

#Checks if passing data is suitable for a poisson and a negative binomial
summary(df_model$passes_completed) 
hist(df_model$passes_completed,
     main = "Histogram of Passes Completed",
     xlab = "Passes Completed",
     ylab = "Frequency")


# jags data list
data_jags_pass <- list(
  y_pass = df_model$passes_completed, # response variable 1 (att)
  player_id = df_model$player_id, # player id 
  team_id = df_model$team_id, # team id
  tau = df_model$tau,
  tau_star = df_model$tau_star,
  delta = delta,
  N = nrow(df_model), # amount of observations
  N_players = max(df_model$player_id), # amount of players
  N_teams = ncol(delta) # amount of teams
)

model_pois_pass <- "
model {

  # Priors on player effects
  for (p in 1:N_players) {
    Delta_star[p] ~ dnorm(m, s)
    #trick to sum-to-zero constraints
    Delta[p] <- Delta_star[p] - mean(Delta_star[])
  }

  # Priors on team effects
  for (k in 1:N_teams) {
  #  lambda_e_star[k] ~ dnorm(0, 0.01) # like team effects in whitaker and baio and blangiardo
  #  lambda_Edive_star[k] ~ dnorm(0, 0.01) # like team effects in whitaker and baio and blangiardo
    
    lambda_e_star[k] ~ dnorm(mu.lambda_e, tau.lambda_e)  # like team effects in whitaker and baio and blangiardo
lambda_Edive_star[k] ~ dnorm(mu.lambda_Edive, tau.lambda_Edive) # like team effects in whitaker and baio and blangiardo

  # sum-to-zero constraints
    lambda_e[k] <- lambda_e_star[k] - mean(lambda_e_star[])
    lambda_Edive[k] <- lambda_Edive_star[k] - mean(lambda_Edive_star[])
  }

# baio and blangiardo impose sum-to-zero constraints to their team effects, so:

  # Likelihood
  for (i in 1:N) {

    y_pass[i] ~ dpois(eta[i] * tau[i])

    opp_sum[i] <- inprod(lambda_Edive[], delta[i,]) # total opponent ability

    log(eta[i]) <- Delta[player_id[i]] +
                   tau[i] * lambda_e[team_id[i]] -
                   tau_star[i] * opp_sum[i]
  }
  # inspired by Baio and Blangiardo code, priors in the random effects
  m~dnorm(0,0.01)
  s~dgamma(1,0.1) # controls spread of player effects
  
# priors on the random effects
mu.lambda_e ~ dnorm(0,0.01)
mu.lambda_Edive ~ dnorm(0,0.01) # mean 0 with precision very small and weak, we do not 
# know anything about team effects apart from the fact that they 
# are around zero, but very weakly-informative
tau.lambda_e ~ dgamma(1,1)
tau.lambda_Edive ~ dgamma(1,1)
# priors are very broad, we let the data decide if teams are similar or different
}
"



# now the jags model is initiated with 3 parallel markov chains which help us 
# check for convergence
jags_model <- jags.model( 
  textConnection(model_pois_pass),
  data = data_jags_pass,
  n.chains = 3,
  n.adapt = 5000
)

update(jags_model, 30000) # account for burn-in

params_pass_pois <- c(
  "Delta", "lambda_e", "lambda_Edive",
   "tau.lambda_e", "tau.lambda_Edive" # "mu.lambda_e", "mu.lambda_Edive", "m", "s"
)


samples_pass_pois <- coda.samples(jags_model, variable.names = params_pass_pois, n.iter = 170000, thin = 1)




# Get summary stats and quantiles
summary_stats <- summary(samples_pass_pois)$statistics
summary_quants <- summary(samples_pass_pois)$quantiles

# Identify row indices corresponding to player ability parameters (Deltas)
Delta_idx <- grep("^Delta\\[", rownames(summary_stats)) # delta indices
Delta_stats <- summary_stats[Delta_idx, ] #  Extract posterior summary statistics for player abilities

player_lookup_pass <- df_model %>%
  distinct(player_id, player) %>%
  arrange(player_id)

#player names instead of delta[1]..., in the appropriate order
player_names_pass <- player_lookup_pass$player 

player_table_poisson <- data.frame(
  player        = player_names, # extracts player names
  ability_mean  = Delta_stats[,"Mean"], # posterior mean for player ability
  ability_sd    = Delta_stats[,"SD"], # posterior standard deviation - measures uncertainty
  ability_lower = summary_quants[Delta_idx, "2.5%"],  # lower bound of the 95% credible interval 
  ability_upper = summary_quants[Delta_idx, "97.5%"]  # upper bound of the 95% credible interval
)

# Data is viewed
View(player_table_poisson)

# CONVERGENCE:

# 1. Gelman-Rubin R-hat overall
gd_uni_pass_pois <- gelman.diag(samples_pass_pois, autoburnin = FALSE, multivariate = FALSE)

# R-hat worst case
max_psrf_pass_pois <- max(gd_uni_pass_pois$psrf[, "Point est."], na.rm = TRUE)
max_psrf_pass_pois


psrf <- as.data.frame(gd_uni_pass_pois$psrf)        
psrf$param <- rownames(psrf)

# 1) Which parameters exceed 1.05 R hat?
bad_point <- psrf %>%
  filter(`Point est.` > 1.05) %>%
  arrange(desc(`Point est.`))

bad_point


# 2. Heidelberger-Welch stationarity 
hw <- heidel.diag(samples_pass_pois)
hw

# 3. ESS
ess_pass_pois <- effectiveSize(samples_pass_pois)
summary(ess_pass_pois)

# 4. ACF
autocorr.plot(samples_pass_pois)

# 5. Geweke for all
g_pass_pois <- geweke.diag(samples_pass_pois)
geweke.plot(samples_pass_pois)

# we check for australia defensive team effect after update:
param_def2_pass_pois <- samples_pass_pois[, "lambda_team_def[2]"]

gelman.diag(param_def2_pass_pois, autoburnin = FALSE, multivariate = FALSE)

heidel.diag(param_def2_pass_pois)

plot(samples_pass_pois[, "lambda_team_def[2]"])

autocorr.plot(samples_pass_pois[, "lambda_team_def[2]"])

# amrabat, hakimi, ronaldo convergence, CHECK DELTA NUMBERS

RonaldoPoisPass<- samples_pass_pois[, "Delta[113]"]
heidel.diag(RonaldoPoisPass)
plot(RonaldoPoisPass, main = expression("Plot for " * Delta[113] * " (Cristiano Ronaldo)"))
gelman.diag(RonaldoPoisPass)
geweke.diag(RonaldoPoisPass)
effectiveSize(RonaldoPoisPass)

AmrabatPoisPass <- samples_pass_pois[, "Delta[103]"]
heidel.diag(AmrabatPoisPass)
plot(AmrabatPoisPass, main = expression("Plot for " * Delta[103] * " (Sofyan Amrabat)"))
gelman.diag(AmrabatPoisPass)
geweke.diag(AmrabatPoisPass)
effectiveSize(AmrabatPoisPass)

HakimiPoisPass <- samples_pass_pois[, "Delta[103]"]
heidel.diag(HakimiPoisPass)
plot(HakimiPoisPass, main = expression("Plot for " * Delta[103] * " (Sofyan Amrabat)"))
gelman.diag(HakimiPoisPass)
geweke.diag(HakimiPoisPass)
effectiveSize(HakimiPoisPass)

# ---------------------------- PREDICTED PASSES OVER WHOLE TOURNAMENT

library(coda)
library(dplyr)

# samples is the list with 3 chains. we combined all these into one matrix
S_pass_poi <- as.matrix(samples_pass_pois) 

# split matrix into parameter blocks: matrix for players abilities, draws for team attacking and defensive effects
Delta_star_draws <- S_pass_poi[, grep("^Delta\\[", colnames(S_pass_poi)), drop = FALSE]                 
team_star_draws  <- S_pass_poi[, grep("^lambda_e\\[", colnames(S_pass_poi)), drop = FALSE]       
oppteam_star_draws <- S_pass_poi[, grep("^lambda_Edive\\[", colnames(S_pass_poi)), drop = FALSE]       

# for each observation i we need player name, team name, and how many minutes theY played
pid  <- df_model$player_id
tid  <- df_model$team_id
tau      <- df_shots$tau
tau_star <- df_shots$tau_star

# this computes all opponent sums at once - opponent defensive sum for each observation in each posterior draw 
opp_sum_draws_pass <- oppteam_star_draws %*% t(delta)   

# build eta 
log_eta_draws <- Delta_star_draws[, pid, drop = FALSE] +
  sweep(team_star_draws[, tid, drop = FALSE], 2, tau, `*`) -
  sweep(opp_sum_draws_pass, 2, tau_star, `*`)

eta_draws_pass <- exp(log_eta_draws)


# extract predicted expected passes per observation using η
mu_draws_pass <- sweep(eta_draws_pass, 2, tau, `*`)  

# extend to tournament totals per player
N_players_pass <- data_jags_pass$N_players

mu_player_draws <- sapply(1:N_players_pass, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_pass[, cols, drop = FALSE])
})


# summarise posterior per player - 2.5% and 97.5% quantiles
player_pred_total_pass_poi <- data.frame(
  player = player_names,
  pred_pass_total_mean  = apply(mu_player_draws, 2, mean),
  pred_pass_total_lower = apply(mu_player_draws, 2, quantile, probs = 0.025),
  pred_pass_total_upper = apply(mu_player_draws, 2, quantile, probs = 0.975)
)

# add observed totals + minutes totals for comparison
# analyze model and how good it is by comapring actual vs predicted
obs_pass_totals <- df_model %>%
  group_by(player_id) %>%
  summarise(
    obs_pass_total = sum(passes_completed),
    minutes_90_total = sum(minutes_90),
    .groups = "drop"
  )



player_pred_total_pass_poi <- player_pred_total_pass_poi %>%
  mutate(player_id = 1:N_players_pass) %>%
  left_join(obs_pass_totals, by = "player_id") %>%
  arrange(desc(pred_pass_total_mean))

errors_pass_poi <- player_pred_total_pass_poi$obs_pass_total - 
  player_pred_total_pass_poi$pred_pass_total_mean

MAE_model_pass_poi  <- mean(abs(errors_pass_poi))
RMSE_model_pass_poi <- sqrt(mean(errors_pass_poi^2))

MAE_model_pass_poi
RMSE_model_pass_poi

View(player_pred_total_pass_poi)


# now the NEGATIVE BINOMIAL jags model is initiated with 3 parallel markov chains 
# which help us check for convergence. neg bin should have better results than poisson
model_nb_pass <- "
model {

  # Priors on player effects
  for (p in 1:N_players) {
    Delta_star[p] ~ dnorm(m, s)
    #trick to sum-to-zero constraints
    Delta[p] <- Delta_star[p] - mean(Delta_star[])
  }

  # Priors on team effects
  for (k in 1:N_teams) {
  #  lambda_e_star[k] ~ dnorm(0, 0.1) # like team effects in whitaker and baio and blangiardo
  #  lambda_Edive_star[k] ~ dnorm(0, 0.1) # like team effects in whitaker and baio and blangiardo
    
    lambda_e_star[k] ~ dnorm(mu.lambda_e, tau.lambda_e)  # like team effects in whitaker and baio and blangiardo
    lambda_Edive_star[k] ~ dnorm(mu.lambda_Edive, tau.lambda_Edive) # like team effects in whitaker and baio and blangiardo

  # sum-to-zero constraints
    lambda_e[k] <- lambda_e_star[k] - mean(lambda_e_star[])
    lambda_Edive[k] <- lambda_Edive_star[k] - mean(lambda_Edive_star[])
  }

# baio and blangiardo impose sum-to-zero constraints to their team effects, so:

  # Likelihood
  for (i in 1:N) {

    opp_sum[i] <- inprod(lambda_Edive[], delta[i,]) # total opponent ability

    log(eta[i]) <- Delta[player_id[i]] +
                   tau[i] * lambda_e[team_id[i]] -
                   tau_star[i] * opp_sum[i]

    mu[i] <- eta[i] * tau[i]

    p[i] <- r / (r + mu[i]) #  since jags needs the first element in dnegbin to be 
    # a probability parameter we introduce p which is a probability parameter
    # which uses mu. explain this
    # The negative binomial model was specified in terms of the mean μ_i = η_i τ_i
    # and dispersion parameter r. For implementation in JAGS, this was 
    # reparameterised as p_i = r / (r + μ_i), since the dnegbin distribution
    # is defined using a probability and size parameter.

    y_pass[i] ~ dnegbin(p[i], r)
  }
  # inspired by Baio and Blangiardo code, priors in the random effects
  m~dnorm(0,1)
  s~dgamma(1,1) # controls spread of player effects
  
# priors on the random effects
mu.lambda_e ~ dnorm(0,1)
mu.lambda_Edive ~ dnorm(0,1) # mean 0 with precision very small and weak, we do not 
# know anything about team effects apart from the fact that they 
# are around zero, but very weakly-informative
tau.lambda_e ~ dgamma(2,2)
tau.lambda_Edive ~ dgamma(1,1)
# priors are very broad, we let the data decide if teams are similar or different

r ~ dgamma(1,1) # since r must be positive and controls overdispersion
}
" 



jags_model_nb <- jags.model( 
  textConnection(model_nb_pass),
  data = data_jags_pass,
  n.chains = 3,
  n.adapt = 2000
)

update(jags_model_nb, 20000) # account for burn-in # 2400

params_pass_nb <- c(
  "Delta", "lambda_e", "lambda_Edive",
  "tau.lambda_e", "tau.lambda_Edive","r"
)


samples_pass_nb <- coda.samples(jags_model_nb, variable.names = params_pass_nb, n.iter = 50000, thin = 1) # 16000

# Get summary stats and quantiles for NB model
summary_nb_stats  <- summary(samples_pass_nb)$statistics
summary_nb_quants <- summary(samples_pass_nb)$quantiles

# extracts Delta parameters (player abilities) from NB summary(stats)
Delta_idx_nb <- grep("^Delta\\[", rownames(sum_stats_nb))
Delta_stats_nb <- sum_stats_nb[Delta_idx_nb, ] # Extract Delta rows

# player names instead of Delta[1]..., in the appropriate order 
player_names <- df_model %>%
  distinct(player_id, player) %>%
  arrange(player_id) %>%
  pull(player)

player_table_nb <- data.frame(
  player        = player_names,
  ability_mean  = Delta_stats_nb[,"Mean"],
  ability_sd    = Delta_stats_nb[,"SD"],
  ability_lower = summary_nb_quants[Delta_idx_nb, "2.5%"],
  ability_upper = summary_nb_quants[Delta_stats_nb, "97.5%"]
)

# NB abilities are presented
View(player_table_nb)

# save to disk
saveRDS(samples_pass_nb, file = "samples_pass_nb.rds")
saveRDS(jags_model_nb, file = "jags_model_nb.rds")

# checking for convergence

# checking for convergence
# 2. Gelman-Rubin R-hat overall
gd_uni_pass_nb <- gelman.diag(samples_pass_nb, autoburnin = FALSE, multivariate = FALSE)

# Worst-case R-hat
max_psrf_pass_nb <- max(gd_uni_pass_nb$psrf[, "Point est."], na.rm = TRUE)
max_psrf_pass_nb

psrf <- as.data.frame(gd_uni_pass_nb$psrf)        
psrf$param <- rownames(psrf)

# 1) Which parameters exceed 1.05 R hat?
bad_point <- psrf %>%
  filter(`Point est.` > 1.05) %>%
  arrange(desc(`Point est.`))

bad_point

# 3. ESS
ess_pass_nb <- effectiveSize(samples_pass_nb)
summary(ess_pass_nb)

# 4. ACF
autocorr.plot(samples_pass_nb)

# 5. Geweke for all
g <- geweke.diag(samples_pass_nb)
geweke.plot(samples_pass_nb)

# amrabat, hakimi, ronaldo convergence, CHECK DELTA NUMBERS

RonaldoNBPass<- samples_pass_nb[, "Delta[113]"]
heidel.diag(RonaldoNBPass)
plot(RonaldoNBPass, main = expression("Plot for " * Delta[113] * " (Cristiano Ronaldo)"))
gelman.diag(RonaldoNBPass)
geweke.diag(RonaldoNBPass)
effectiveSize(RonaldoNBPass)

AmrabatNBPass <- samples_pass_nb[, "Delta[103]"]
heidel.diag(AmrabatNBPass)
plot(AmrabatNBPass, main = expression("Plot for " * Delta[103] * " (Sofyan Amrabat)"))
gelman.diag(AmrabatNBPass)
geweke.diag(AmrabatNBPass)
effectiveSize(AmrabatNBPass)

HakimiNBPass <- samples_pass_nb[, "Delta[103]"]
heidel.diag(HakimiNBPass)
plot(HakimiNBPass, main = expression("Plot for " * Delta[103] * " (Sofyan Amrabat)"))
gelman.diag(HakimiNBPass)
geweke.diag(HakimiNBPass)
effectiveSize(HakimiNBPass)


# PREDICTED PASSES OVER WHOLE TOURNAMENT TO USE RMSE AND MAE (NEG BIN)

library(coda)
library(dplyr)

# samples_pass_nb is the list with 3 chains. combine into one matrix
S_pass_nb <- as.matrix(samples_pass_nb)

# split matrix into parameter blocks
Delta_star_draws <- S_pass_nb[, grep("^Delta\\[", colnames(S_pass_nb)), drop = FALSE]
team_star_draws   <- S_pass_nb[, grep("^lambda_e\\[", colnames(S_pass_nb)), drop = FALSE]
oppteam_star_draws  <- S_pass_nb[, grep("^lambda_Edive\\[", colnames(S_pass_nb)), drop = FALSE]
r_draws_pass     <- S_pass_nb[, "r"]   # optional (only needed if you simulate counts)

# for each observation i we need player id, team id, and minutes
pid  <- df_model$player_id
tid  <- df_model$team_id
tau      <- df_model$tau
tau_star <- df_model$tau_star

# opponent defensive sum for each observation in each posterior draw
opp_sum_draws_nb <- oppteam_draws_pass %*% t(delta)

# build eta 
log_eta_draws <- Delta_star_draws[, pid, drop = FALSE] +
  sweep(team_star_draws[, tid, drop = FALSE], 2, tau, `*`) -
  sweep(opp_sum_draws_nb, 2, tau_star, `*`)

eta_draws_pass <- exp(log_eta_draws)

# expected passes per observation (NB mean is still mu)
mu_draws_nb <- sweep(eta_draws_pass, 2, tau, `*`)

# extend to tournament totals per player (expected totals; excludes NB randomness)
N_players_pass_nb <- data_jags_pass$N_players

mu_player_draws_nb <- sapply(1:N_players_pass_nb, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_nb[, cols, drop = FALSE])
})

# summarise posterior per player
player_pred_total_pass_nb <- data.frame(
  player = player_names,
  pred_pass_total_mean  = apply(mu_player_draws_nb, 2, mean),
  pred_pass_total_lower = apply(mu_player_draws_nb, 2, quantile, probs = 0.025),
  pred_pass_total_upper = apply(mu_player_draws_nb, 2, quantile, probs = 0.975)
)

# add observed totals + minutes totals for comparison
obs_pass_totals_nb <- df_model %>%
  group_by(player_id) %>%
  summarise(
    obs_pass_total = sum(passes_completed),
    minutes_90_total = sum(minutes_90),
    .groups = "drop"
  )

player_pred_total_pass_nb <- player_pred_total_pass_nb %>%
  mutate(player_id = 1:N_players_pass_nb) %>%
  left_join(obs_pass_totals_nb, by = "player_id") %>%
  arrange(desc(pred_pass_total_mean))


# compute overall model MAE and RMSE (single values for model)
errors_pass_nb <- player_pred_total_pass_nb$obs_pass_totals_nb - 
  player_pred_total_pass_nb$pred_pass_total_mean

MAE_model_pass_nb  <- mean(abs(errors_pass_nb))
RMSE_model_pass_nb <- sqrt(mean(errors_pass_nb^2))

MAE_model_pass_nb 
RMSE_model_pass_nb 


View(player_pred_total_pass_nb)

load("my_workspace.RData")
save.image("my_workspace.RData")