# necessary packages
library(dplyr)
library(rjags)
library(coda)


# Preparing Data
set.seed(88)

players <- read.csv("C:/Users/psail/OneDrive/Desktop/UOM 4th yr/Thesis/Dataset/Proper Dataset/Players minutes without gks.csv")
passing <- read.csv("C:/Users/psail/OneDrive/Desktop/UOM 4th yr/Thesis/Dataset/Proper Dataset/player_passingandshootingwithoutgks.csv")
defense <- read.csv("C:/Users/psail/OneDrive/Desktop/UOM 4th yr/Thesis/Dataset/Proper Dataset/player_defensewithoutgks.csv")

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
    minutes_90 = minutes_90s,                             # minutes per 90 for all players since some players played less
    is_DF = ifelse(position.x == "DF", 1, 0),             # dummy for defenders
    is_MF = ifelse(position.x == "MF", 1, 0)              # dummy for midfielders (FW is baseline by elimination)
  )

delta <- as.matrix(df_model[, 7:38]) # teams columns in dataset: from Argentina to Wales

# checks for clean data 
sapply(df_model[, c("passes_completed", "interceptions", "minutes_90",
                    "player_id", "team_id", "is_DF", "is_MF")],
       function(x) sum(is.na(x)))

#Checks if passing data is suitable for a poisson and a negative binomial
summary(df_model$passes_completed) 
hist(df_model$passes_completed,
     main = "Histogram of Passes Completed",
     xlab = "Passes Completed",
     ylab = "Frequency")


# jags data list
data_jags <- list(
  y_pass = df_model$passes_completed, # response variable 1 (att)
  minutes = df_model$minutes_90, # values like 3.4, 5.1 etc
  player_id = df_model$player_id, # player id 
  team_id = df_model$team_id, # team id
  is_DF = df_model$is_DF, # position is defender?
  is_MF = df_model$is_MF, # position is midfielder?
  delta = delta,
  N = nrow(df_model), # amount of observations
  N_players = max(df_model$player_id), # amount of players
  N_teams = ncol(delta) # amount of teams
)

model_pois_att <- "
model {
  
  # Player ability (centered)
  for (p in 1:N_players) {
    Delta_raw[p] ~ dnorm(0, 1)
  }
  mean_Delta <- mean(Delta_raw[])
  for (p in 1:N_players) {
    Delta[p] <- Delta_raw[p] - mean_Delta
  }
  
  # Team attack/defense (centered)
  for (k in 1:N_teams) {
    att_raw[k] ~ dnorm(0, 1)
    def_raw[k] ~ dnorm(0, 4)
  }
  mean_att <- mean(att_raw[])
  mean_def <- mean(def_raw[])
  for (k in 1:N_teams) {
    lambda_team_att[k] <- att_raw[k] - mean_att
    lambda_team_def[k] <- def_raw[k] - mean_def
  }
  
  # Likelihood
  for (i in 1:N) {
    
    # opponent sum
  opponent_sum[i] <- inprod(lambda_team_def[], delta[i,])

    
    # linear predictor for passing rate
    eta[i] <- Delta[player_id[i]] +
      lambda_team_att[team_id[i]] -
      opponent_sum[i]
    
    # exposure offset
    log(mu[i]) <- log(minutes[i]) + eta[i]
    
    y_pass[i] ~ dpois(mu[i])
  }
}
"


# now the jags model is initiated with 3 parallel markov chains which help us 
# check for convergence
jags_model <- jags.model( 
  textConnection(model_pois_att),
  data = data_jags,
  n.chains = 3,
  n.adapt = 6000
)

update(jags_model, 30000) # account for burn-in

params_pass_pois <- c("Delta", "lambda_team_att", "lambda_team_def") # , "mu")

samples_pass_pois <- coda.samples(jags_model, variable.names = params_pass_pois, n.iter = 170000, thin = 1)




# Get summary stats and quantiles
summary_stats <- summary(samples_pass_pois)$statistics
summary_quants <- summary(samples_pass_pois)$quantiles

# Identify row indices corresponding to player ability parameters (Deltas)
Delta_indices <- grep("^Delta", rownames(summary_stats))

#  Extract posterior summary statistics for player abilities
Delta_stats <- summary_stats[Delta_indices, ]

#player names instead of delta[1]..., in the appropriate order
player_names <- levels(factor(df_model$player))

player_table_poisson <- data.frame(
  player        = player_names, # extracts player names
  ability_mean  = Delta_stats[,"Mean"], # posterior mean for player ability
  ability_sd    = Delta_stats[,"SD"], # posterior standard deviation - measures uncertainty
  ability_lower = summary(samples_pass_pois)$quantiles[Delta_indices, "2.5%"], # lower bound of the 95% credible interval
  ability_upper = summary(samples_pass_pois)$quantiles[Delta_indices, "97.5%"] # upper bound of the 95% credible interval
)

# Data is viewed
View(player_table_poisson)
plot(samples_pass_pois) # takes long since it gives traceplot for each player and team
 


#extracts posterior for latent player ability
Delta_stats[, ]

mean(as.matrix(samples_pass_pois)[, "lambda_team_att[3]"])
# the above calculates the attacking contribution to passes - team 3,
# measured on the log rate scale

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

# ---------------------------- PREDICTED PASSES OVER WHOLE TOURNAMENT

library(coda)
library(dplyr)

# samples is the list with 3 chains. we combined all these into one matrix
S <- as.matrix(samples_pass_pois) 

# split matrix into parameter blocks: matrix for players abilities, draws for team attacking and defensive effects
Delta_draws <- S[, grep("^Delta\\[", colnames(S)), drop = FALSE]                 
att_draws   <- S[, grep("^lambda_team_att\\[", colnames(S)), drop = FALSE]       
def_draws   <- S[, grep("^lambda_team_def\\[", colnames(S)), drop = FALSE]       

# for each observation i we need player name, team name, and how many minutes theY played
pid  <- df_model$player_id
tid  <- df_model$team_id
mins <- df_model$minutes_90

# this computes all opponent sums at once - opponent defensive sum for each observation in each posterior draw 
opp_sum_draws_poi <- def_draws %*% t(delta)   

# we build eta for each posterior draw and each observation
eta_draws <- Delta_draws[, pid, drop = FALSE] +
  att_draws[, tid, drop = FALSE] - opp_sum_draws_poi                  

# extract predicted expected passes per observation using η
mu_draws <- sweep(exp(eta_draws), 2, mins, `*`)  

# extend to tournament totals per player
N_players <- data_jags$N_players

mu_player_draws <- sapply(1:N_players, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws[, cols, drop = FALSE])
})

# rebuild player names in the same order as player_id (1..N_players)
player_names <- df_model %>%
  distinct(player_id, player) %>%
  arrange(player_id) %>%
  pull(player)

# summarise posterior per player - 2.5% and 97.5% quantiles
player_pred_total_poi <- data.frame(
  player = player_names,
  pred_pass_total_mean  = apply(mu_player_draws, 2, mean),
  pred_pass_total_lower = apply(mu_player_draws, 2, quantile, probs = 0.025),
  pred_pass_total_upper = apply(mu_player_draws, 2, quantile, probs = 0.975)
)

# add observed totals + minutes totals for comparison
# analyze model and how good it is by comapring actual vs predicted
obs_totals <- df_model %>%
  group_by(player_id) %>%
  summarise(
    obs_pass_total = sum(passes_completed),
    minutes_90_total = sum(minutes_90),
    .groups = "drop"
  )



player_pred_total_poi <- player_pred_total_poi %>%
  mutate(player_id = 1:N_players) %>%
  left_join(obs_totals, by = "player_id") %>%
  arrange(desc(pred_pass_total_mean))

errors_pass_poi <- player_pred_total_poi$obs_pass_total - 
  player_pred_total_poi$pred_pass_total_mean

MAE_model_pass_poi  <- mean(abs(errors_pass_poi))
RMSE_model_pass_poi <- sqrt(mean(errors_pass_poi^2))

MAE_model_pass_poi
RMSE_model_pass_poi

View(player_pred_total_poi)


# now the NEGATIVE BINOMIAL jags model is initiated with 3 parallel markov chains 
# which help us check for convergence. neg bin should have better results than poisson
model_nb_att <- "
model {

  # Player ability (Centred)
  for (p in 1:N_players) {
    Delta_raw[p] ~ dnorm(0, 4)
  }
  mean_Delta <- mean(Delta_raw[])
  for (p in 1:N_players) {
    Delta[p] <- Delta_raw[p] - mean_Delta
  }

  # Team attack/defense (Centred)
  for (k in 1:N_teams) {
    att_raw[k] ~ dnorm(0, 4)
    def_raw[k] ~ dnorm(0, 4)
  }
  mean_att <- mean(att_raw[])
  mean_def <- mean(def_raw[])
  for (k in 1:N_teams) {
    lambda_team_att[k] <- att_raw[k] - mean_att
    lambda_team_def[k] <- def_raw[k] - mean_def
  }

  # Overdispersion parameter (NB size/shape)
  r ~ dgamma(2, 1)

  # Likelihood
  for (i in 1:N) {

    # opponent sum
    opponent_sum[i] <- inprod(lambda_team_def[], delta[i,])

    # linear predictor
    eta[i] <- Delta[player_id[i]] +
              lambda_team_att[team_id[i]] -
              opponent_sum[i]

    # exposure offset
    log(mu[i]) <- log(minutes[i]) + eta[i]

    # Negative binomial with mean mu[i]
    p_nb[i] <- r / (r + mu[i]) # here jags doesnt parameterise NB by mean, so wesolve for p in terms of miu (simple subject of formula)
    y_pass[i] ~ dnegbin(p_nb[i], r) # neg bin (mean, var)
  }
}
"


jags_model_nb <- jags.model( 
  textConnection(model_nb_att),
  data = data_jags,
  n.chains = 3,
  n.adapt = 3000
)

update(jags_model_nb, 15000) # account for burn-in # 2400

params_pass_nb <- c("Delta", "lambda_team_att", "lambda_team_def", "r") # "mu"

samples_pass_nb <- coda.samples(jags_model_nb, variable.names = params_pass_nb, n.iter = 40000, thin = 1) # 16000

# Get summary stats and quantiles for NB model
summary_nb_stats  <- summary(samples_pass_nb)$statistics
summary_nb_quants <- summary(samples_pass_nb)$quantiles

# extracts Delta parameters (player abilities) from NB summary(stats)
Delta_nb_indices <- grep("^Delta", rownames(summary_nb_stats))

# Extract Delta rows
Delta_nb_stats <- summary_nb_stats[Delta_nb_indices, ]

# player names instead of Delta[1]..., in the appropriate order 
player_names <- df_model %>%
  distinct(player_id, player) %>%
  arrange(player_id) %>%
  pull(player)

player_table_nb <- data.frame(
  player        = player_names,
  ability_mean  = Delta_nb_stats[,"Mean"],
  ability_sd    = Delta_nb_stats[,"SD"],
  ability_lower = summary_nb_quants[Delta_nb_indices, "2.5%"],
  ability_upper = summary_nb_quants[Delta_nb_indices, "97.5%"]
)

# NB abilities are presented
View(player_table_nb)

# save to disk
saveRDS(samples_pass_nb, file = "samples_pass_nb.rds")
saveRDS(jags_model_nb, file = "jags_model_nb.rds")

# checking for convergence

# 1. plots for all players - takes long
plot(samples_pass_nb)

# checking for convergence
# 2. Gelman-Rubin R-hat overall
gd_uni_pass_nb <- gelman.diag(samples_pass_nb, autoburnin = FALSE, multivariate = FALSE)

# Worst-case R-hat
max_psrf_pass_nb <- max(gd_uni_pass_nb$psrf[, "Point est."], na.rm = TRUE)
max_psrf_pass_nb


# 3. ESS
ess_pass_nb <- effectiveSize(samples_pass_nb)
summary(ess_pass_nb)

# 4. ACF
autocorr.plot(samples_pass_nb)

# 5. Geweke for all
g <- geweke.diag(samples_pass_nb)
geweke.plot(samples_pass_nb)


# PREDICTED PASSES OVER WHOLE TOURNAMENT TO USE RMSE AND MAE (NEG BIN)

library(coda)
library(dplyr)

# samples_pass_nb is the list with 3 chains. combine into one matrix
S_nb <- as.matrix(samples_pass_nb)

# split matrix into parameter blocks
Delta_draws_nb <- S_nb[, grep("^Delta\\[", colnames(S_nb)), drop = FALSE]
att_draws_nb   <- S_nb[, grep("^lambda_team_att\\[", colnames(S_nb)), drop = FALSE]
def_draws_nb   <- S_nb[, grep("^lambda_team_def\\[", colnames(S_nb)), drop = FALSE]
r_draws_nb     <- S_nb[, "r"]   # optional (only needed if you simulate counts)

# for each observation i we need player id, team id, and minutes
pid  <- df_model$player_id
tid  <- df_model$team_id
mins <- df_model$minutes_90

# opponent defensive sum for each observation in each posterior draw
opp_sum_draws_nb <- def_draws_nb %*% t(delta)

# build eta for each posterior draw and each observation
eta_draws_nb <- Delta_draws_nb[, pid, drop = FALSE] +
  att_draws_nb[, tid, drop = FALSE] - opp_sum_draws_nb

# expected passes per observation (NB mean is still mu)
mu_draws_nb <- sweep(exp(eta_draws_nb), 2, mins, `*`)

# extend to tournament totals per player (expected totals; excludes NB randomness)
N_players <- data_jags$N_players

mu_player_draws_nb <- sapply(1:N_players, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_nb[, cols, drop = FALSE])
})

# summarise posterior per player
player_pred_total_nb <- data.frame(
  player = player_names,
  pred_pass_total_mean  = apply(mu_player_draws_nb, 2, mean),
  pred_pass_total_lower = apply(mu_player_draws_nb, 2, quantile, probs = 0.025),
  pred_pass_total_upper = apply(mu_player_draws_nb, 2, quantile, probs = 0.975)
)

# add observed totals + minutes totals for comparison
obs_totals <- df_model %>%
  group_by(player_id) %>%
  summarise(
    obs_pass_total = sum(passes_completed),
    minutes_90_total = sum(minutes_90),
    .groups = "drop"
  )

player_pred_total_nb <- player_pred_total_nb %>%
  mutate(player_id = 1:N_players) %>%
  left_join(obs_totals, by = "player_id") %>%
  arrange(desc(pred_pass_total_mean))

# compute per-player absolute error (useful diagnostic)
player_pred_total_nb <- player_pred_total_nb %>%
  mutate(abs_error = abs(obs_pass_total - pred_pass_total_mean))

# compute overall model MAE and RMSE (single values for model)
errors_pass_nb <- player_pred_total_nb$obs_pass_total - 
  player_pred_total_nb$pred_pass_total_mean

MAE_model_pass_nb  <- mean(abs(errors_pass_nb))
RMSE_model_pass_nb <- sqrt(mean(errors_pass_nb^2))

MAE_model_pass_nb 
RMSE_model_pass_nb 


View(player_pred_total_nb)

# -------------------------------------------------------------------------

# Model for interceptions only for defenders, since does not make sense for
# mids and forwards
# Defenders only (temporary subset)
df_defmid <- df_model %>%
  filter(is_DF == 1 | is_MF == 1) %>%           # keep DF + MF only
  mutate(defmid_player_id = as.numeric(factor(player)))
# defender player IDs


delta_defmid <- delta[df_model$is_DF == 1 | df_model$is_MF == 1, ] # latent player ability only defenders and mids

# delta needs to have same no of rows as no of defenders
# we check this in the code: 
nrow(df_defmid)
nrow(delta_defmid)   # MUST match

# checks for clean data (after filtering and creating clean IDs)
sapply(df_defmid[, c("interceptions", "minutes_90", "player_id", "team_id", "is_DF","is_MF")],
       function(x) sum(is.na(x)))

#Checks if interception data is suitable for a poisson and a negative binomial
summary(df_defmid$interceptions) 
hist(df_defmid$interceptions,
     main = "Histogram of Interceptions Completed (Defenders and Midfielders Only)",
     xlab = "Interceptions Completed",
     ylab = "Frequency")


# Now N_teams will be a proper integer
N_teams <- ncol(delta_defmid)

# Rebuild the JAGS data list 
data_jags_int_poi <- list(
  y_int     = df_defmid$interceptions,
  minutes   = df_defmid$minutes_90,
  player_id = df_defmid$defmid_player_id,
  team_id   = df_defmid$team_id,
  delta     = delta_defmid,
  N         = nrow(df_defmid),
  N_players = max(df_defmid$defmid_player_id),
  N_teams   = ncol(delta_defmid)
)

model_pois_int <- "
model {

  # Player interception ability (defenders and midfielders) (centred)
  for (p in 1:N_players) {
    Delta_raw[p] ~ dnorm(0, 1)
  }
   mean_Delta <- mean(Delta_raw[]) 
  for (p in 1:N_players) {
    Delta[p] <- Delta_raw[p] - mean_Delta
  }


  # Team defense - interceptions (centred)
for (k in 1:N_teams) {
    team_int_raw[k] ~ dnorm(0, 1)
  }
  mean_team_int <- mean(team_int_raw[])
  for (k in 1:N_teams) {
    lambda_team_int[k] <- team_int_raw[k] - mean_team_int
  }


  # Opponent attacking quality (centered)
  for (k in 1:N_teams) {
    opp_att_raw[k] ~ dnorm(0, 1)
  }
  mean_opp_att <- mean(opp_att_raw[])
  for (k in 1:N_teams) {
    lambda_opp_att[k] <- opp_att_raw[k] - mean_opp_att
  }
  
  # Likelihood
  for (i in 1:N) {
    
    # opponent sum
  opponent_sum[i] <- inprod(lambda_opp_att[], delta[i,])

    # linear predictor for interception rate
    eta[i] <- Delta[player_id[i]] +
              lambda_team_int[team_id[i]] +
              opponent_sum[i]  #(stronger opponent attack quality, more interceptions for defending team, mention in thesis)

    # exposure offset - account for minutes
    log(mu[i]) <- log(minutes[i]) + eta[i]
    
  # observation
  y_int[i] ~ dpois(mu[i])

   
  }

}
"

jags_model_int_poi <- jags.model(
  textConnection(model_pois_int),
  data = data_jags_int_poi,
  n.chains = 3,
  n.adapt = 320 
)

update(jags_model_int_poi, 3800)  
 
params_int_poi <- c("Delta", "lambda_team_int", "lambda_opp_att")

samples_int_poi <- coda.samples(jags_model_int_poi, variable.names = params_int_poi, n.iter = 25000, thin = 1) 
# 6) Summaries + defender ability table
summary_stats_int  <- summary(samples_int_poi)$statistics
summary_quants_int <- summary(samples_int_poi)$quantiles

Delta_idx <- grep("^Delta\\[", rownames(summary_stats_int))
Delta_stats <- summary_stats_int[Delta_idx, ]

defendermidfielder_names <- levels(factor(df_defmid$player)) # extract only defenders

player_table_poisson_int <- data.frame(
  player        = defendermidfielder_names,
  ability_mean  = Delta_stats[, "Mean"],
  ability_sd    = Delta_stats[, "SD"],
  ability_lower = summary_quants_int[Delta_idx, "2.5%"],
  ability_upper = summary_quants_int[Delta_idx, "97.5%"]
)

View(player_table_poisson_int)

# traceplots
plot(samples_int_poi)

#Achraf Hakimi plot
plot(samples_int_poi[, 13],
     main = expression("Plot for " * Delta[13] * " (Achraf Hakimi)"))

#Sofyan Amrabat plot
plot(samples_int_poi[, 213],
     main = expression("Plot for " * Delta[213] * " (Sofyan Amrabat)"))

HakimiPoisInt <- samples_int_poi[, "Delta[13]"]
heidel.diag(HakimiPoisInt)

AmrabatPoisInt <- samples_int_poi[, "Delta[213]"]
heidel.diag(AmrabatPoisInt)

#extracts posterior for latent player ability
Delta_stats[, ]

# checking for convergence

# 1. traceplots for all players - takes long
plot(samples_int_poi)

# checking for convergence
# 2. Gelman-Rubin R-hat overall
gd_uni_int <- gelman.diag(samples_int_poi, autoburnin = FALSE, multivariate = FALSE)

# Worst-case R-hat
max_psrf <- max(gd_uni_int$psrf[, "Point est."], na.rm = TRUE)
max_psrf


# 3. ESS
ess <- effectiveSize(samples_int_poi)
summary(ess)

# 4. ACF
autocorr.plot(samples_int_poi)

# 5. Geweke for all
g <- geweke.diag(samples_int_poi)
geweke.plot(samples_int_poi)

# PRECIDTED INTERCEPTIONS OVER WHOLE TOURNAMENT TO USE RMSE AND MAE (POISS)
library(coda)
library(dplyr)

# combine chains into one matrix 
S_int <- as.matrix(samples_int_poi)

# split into player effects, each column is a defender, each row a posterior draw
Delta_draws_int <- S_int[, grep("^Delta\\[", colnames(S_int)), drop = FALSE]
team_draws_int  <- S_int[, grep("^lambda_team_int\\[", colnames(S_int)), drop = FALSE]
opp_draws_int   <- S_int[, grep("^lambda_opp_att\\[", colnames(S_int)), drop = FALSE]

# player and team ids + minutes for defenders
pid  <- df_defmid$defmid_player_id
tid  <- df_defmid$team_id
mins <- df_defmid$minutes_90

# matrix mult between posterior draw and opponent parameterto obtain opponent contribution
opp_sum_draws_int <- opp_draws_int %*% t(delta_defmid)

# build eta for each posterior draw and observation
# we add attacking opponent strength because the more attacks the more chance
# of interceptions by defenders
eta_draws_int <- Delta_draws_int[, pid, drop = FALSE] +
  team_draws_int[, tid, drop = FALSE] + opp_sum_draws_int

# expected interceptions per observation
# achieved by mutiplying exp(eta) by minutes for exposure
mu_draws_int <- sweep(exp(eta_draws_int), 2, mins, `*`)

# total expected interceptions over tournament per defender/midfielder
N_players_defmid <- max(df_defmid$defmid_player_id)

#posterior distributions for posterior draws for each defender
mu_player_draws_int <- sapply(1:N_players_defmid, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_int[, cols, drop = FALSE])
})

# summarise posterior per defender
player_pred_total_int <- data.frame(
  player = defendermidfielder_names,
  pred_int_total_mean  = apply(mu_player_draws_int, 2, mean),
  pred_int_total_lower = apply(mu_player_draws_int, 2, quantile, probs = 0.025),
  pred_int_total_upper = apply(mu_player_draws_int, 2, quantile, probs = 0.975)
)

# observed totals + minutes totals for comparison
obs_totals_int <- df_defmid %>%
  group_by(defmid_player_id) %>%
  summarise(
    obs_int_total = sum(interceptions),
    minutes_90_total = sum(minutes_90),
    .groups = "drop"
  )

#join predicted and observed data
player_pred_total_int <- player_pred_total_int %>%
  mutate(defmid_player_id = 1:N_players_defmid) %>%
  left_join(obs_totals_int, by = "defmid_player_id") %>%
  arrange(desc(pred_int_total_mean))

# compute overall model MAE and RMSE (single values for model)
errors <- player_pred_total_int$obs_int_total - 
  player_pred_total_int$pred_int_total_mean

MAE_model_int_poi  <- mean(abs(errors))
RMSE_model_int_poi <- sqrt(mean(errors^2))

MAE_model_int_poi
RMSE_model_int_poi

View(player_pred_total_int)

# INTERCEPTIONS - NEGATIVE BINOMIAL

library(dplyr)
library(rjags)
library(coda)

 
# JAGS data
data_jags_int_negbin <- list(
  y_int     = df_defmid$interceptions,
  minutes   = df_defmid$minutes_90,
  player_id = df_defmid$defmid_player_id,
  team_id   = df_defmid$team_id,
  delta     = delta_defmid,
  N         = nrow(df_defmid),
  N_players = max(df_defmid$defmid_player_id),
  N_teams   = ncol(delta_defmid)
)

# Centered NB model 
model_nb_int <- "
model {

  # Player interception ability (defenders and midfielders) (centred)
  for (p in 1:N_players) {
    Delta_raw[p] ~ dnorm(0, 1)
  }
   mean_Delta <- mean(Delta_raw[]) 
  for (p in 1:N_players) {
    Delta[p] <- Delta_raw[p] - mean_Delta
  }

  # Team defense - interceptions (centred)
  for (k in 1:N_teams) {
    team_int_raw[k] ~ dnorm(0, 1)
  }
  mean_team_int <- mean(team_int_raw[])
  for (k in 1:N_teams) {
    lambda_team_int[k] <- team_int_raw[k] - mean_team_int
  }

  # Opponent attacking quality (centered)
  for (k in 1:N_teams) {
    opp_att_raw[k] ~ dnorm(0, 1)
  }
  mean_opp_att <- mean(opp_att_raw[])
  for (k in 1:N_teams) {
    lambda_opp_att[k] <- opp_att_raw[k] - mean_opp_att
  }

  # Overdispersion parameter (Negative Binomial)
  r ~ dgamma(0.01, 0.01)

  # Likelihood
  for (i in 1:N) {
    
    # opponent sum
    opponent_sum[i] <- inprod(lambda_opp_att[], delta[i,])

    # linear predictor for interception rate
    eta[i] <- Delta[player_id[i]] +
              lambda_team_int[team_id[i]] +
              opponent_sum[i]  #(stronger opponent attack quality, more interceptions for defending team, mention in thesis)

    # exposure offset - account for minutes
    log(mu[i]) <- log(minutes[i]) + eta[i]

    # Negative binomial with mean mu[i]
    p_nb[i] <- r / (r + mu[i])
    y_int[i] ~ dnegbin(p_nb[i], r)
  }

}
"


# 4) Fit NB model
jags_model_int_nb <- jags.model(
  textConnection(model_nb_int),
  data = data_jags_int_negbin,
  n.chains = 3,
  n.adapt = 1200 
)

update(jags_model_int_nb, 3000) 

params_int_nb <- c("Delta", "lambda_team_int", "lambda_opp_att", "r")
samples_int_nb <- coda.samples(
  jags_model_int_nb,
  variable.names = params_int_nb,
  n.iter = 22000, 
  thin = 1
)

# 5) Defender ability table (NB)
summary_stats_int_nb  <- summary(samples_int_nb)$statistics
summary_quants_int_nb <- summary(samples_int_nb)$quantiles

Delta_idx_nb   <- grep("^Delta\\[", rownames(summary_stats_int_nb))
Delta_stats_nb <- summary_stats_int_nb[Delta_idx_nb, ]

defendermidfielder_names <- levels(factor(df_defmid$player))

player_table_nb_int <- data.frame(
  player        = defendermidfielder_names,
  ability_mean  = Delta_stats_nb[, "Mean"],
  ability_sd    = Delta_stats_nb[, "SD"],
  ability_lower = summary_quants_int_nb[Delta_idx_nb, "2.5%"],
  ability_upper = summary_quants_int_nb[Delta_idx_nb, "97.5%"]
)

View(player_table_nb_int)

# Convergence checks:
# 1) Gelman-Rubin:
gd_uni_int_nb <- gelman.diag(samples_int_nb, autoburnin = FALSE, multivariate = FALSE)
max_psrf_int_nb <- max(gd_uni_int_nb$psrf[, "Point est."], na.rm = TRUE)
max_psrf_int_nb

#2) Effective Sample Size:
ess_int_nb <- effectiveSize(samples_int_nb)
summary(ess_int_nb)

#3) Autocorrelation plot:
autocorr.plot(samples_int_nb)

#4) Geweke plot:
g_int_nb <- geweke.diag(samples_int_nb)
geweke.plot(samples_int_nb)

#5) Traceplot:
traceplot(samples_int_nb)

# PREDICTED INTERCEPTIONS OVER WHOLE TOURNAMENT (NB)
# expected totals (mu), plus overall RMSE/MAE
# =========================
# NEGATIVE BINOMIAL VERSION
# =========================

# combine chains into one matrix
S_int_nb <- as.matrix(samples_int_nb)

# split into parameter blocks (each row is a posterior draw)
Delta_draws_int_nb <- S_int_nb[, grep("^Delta\\[", colnames(S_int_nb)), drop = FALSE]
team_draws_int_nb  <- S_int_nb[, grep("^lambda_team_int\\[", colnames(S_int_nb)), drop = FALSE]
opp_draws_int_nb   <- S_int_nb[, grep("^lambda_opp_att\\[", colnames(S_int_nb)), drop = FALSE]
r_draws_int_nb     <- S_int_nb[, "r"]  

# player and team ids + minutes for defenders
pid  <- df_defmid$player_index_defmid
tid  <- df_defmid$team_id
mins <- df_defmid$minutes_90

# matrix mult, posterior draws x opponent design -> opponent contribution per obs
opp_sum_draws_int_nb <- opp_draws_int_nb %*% t(delta_defmid)

# build eta for each posterior draw and observation
# add opponent attacking strength: more attacks -> more interception opportunities
eta_draws_int_nb <- Delta_draws_int_nb[, pid, drop = FALSE] +
  team_draws_int_nb[, tid, drop = FALSE] +
  opp_sum_draws_int_nb

# expected interceptions per observation (NB mean is mu = exp(eta)*exposure)
mu_draws_int_nb <- sweep(exp(eta_draws_int_nb), 2, mins, `*`)

# total expected interceptions over tournament per defender
N_players_defmid <- max(df_defmid$player_index_defmid)

# posterior distribution of tournament totals for each defender (sum across their obs)
mu_player_draws_int_nb <- sapply(1:N_players_defmid, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_int_nb[, cols, drop = FALSE])
})

# summarise posterior per defender
player_pred_total_int_nb <- data.frame(
  player = defendermidfielder_names,
  pred_int_total_mean  = apply(mu_player_draws_int_nb, 2, mean),
  pred_int_total_lower = apply(mu_player_draws_int_nb, 2, quantile, probs = 0.025),
  pred_int_total_upper = apply(mu_player_draws_int_nb, 2, quantile, probs = 0.975)
)

# observed totals + minutes totals for comparison
obs_totals_int <- df_defmid %>%
  group_by(player_index_defmid) %>%
  summarise(
    obs_int_total     = sum(interceptions),
    minutes_90_total  = sum(minutes_90),
    .groups = "drop"
  )

# join predicted and observed data
player_pred_total_int_nb <- player_pred_total_int_nb %>%
  mutate(player_index_defmid = 1:N_players_defmid) %>%
  left_join(obs_totals_int, by = "player_index_defmid") %>%
  arrange(desc(pred_int_total_mean))

# compute overall model MAE and RMSE (single values for model)
errors <- player_pred_total_int_nb$obs_int_total - 
  player_pred_total_int_nb$pred_int_total_mean

MAE_model  <- mean(abs(errors))
RMSE_model <- sqrt(mean(errors^2))

MAE_model_int_negbin
RMSE_model_int_negbin

View(player_pred_total_int_nb)



# SHOTS MODEL for forwards and midfielders
library(dplyr)
library(rjags)
library(coda)

set.seed(88)

# 1) Subset to midfielders + forwards, build player IDs 

# MF+FW sample
keep_idx <- which(
  (df_model$is_MF == 1 | df_model$position.x == "FW") &
    !is.na(df_model$shots) &
    df_model$minutes_90 > 0
)

# Subset df_model using the SAME sample
df_shots <- df_model[keep_idx, ] %>%
  mutate(
    shot_player_id = as.numeric(factor(player)),
    is_MF_shots = ifelse(position.x == "MF", 1, 0)
  )
nrow(df_shots)


# Subset delta using the SAME sample
delta_shots <- as.matrix(delta[keep_idx, ])

# check for equality
stopifnot(nrow(df_shots) == nrow(delta_shots))
N_teams <- ncol(delta_shots)
stopifnot(max(df_shots$team_id) <= N_teams)

# Check player number 
player_lookup <- df_shots %>%
  distinct(shot_player_id, player) %>%
  arrange(shot_player_id)

player_lookup

#Amrabat 156, Ronaldo 40


# quick NA check
sapply(df_shots[, c("shots", "minutes_90", "shot_player_id", "team_id", "is_MF_shots")],
       function(x) sum(is.na(x)))

#Checks if shot data is suitable for a poisson and a negative binomial
summary(df_shots$shots) 
hist(df_shots$shots,
     main = "Histogram of Shots Taken (Attackers and Midfielders Only)",
     xlab = "Shots Taken",
     ylab = "Frequency")

# 2) JAGS data
data_jags_shots <- list(
  y_shots   = df_shots$shots,
  minutes   = df_shots$minutes_90,
  player_id = df_shots$shot_player_id,
  team_id   = df_shots$team_id,
  is_MF     = df_shots$is_MF_shots,
  delta     = delta_shots,
  N         = nrow(df_shots),
  N_players = max(df_shots$shot_player_id),
  N_teams   = N_teams
)

# ---------------------------- POISSON SHOTS (CENTERED)
model_pois_shots <- "
model {

  # Player shot volume ability (centered)
  for (p in 1:N_players) {
    Delta_raw[p] ~ dnorm(0, 1)
  }
  mean_Delta <- mean(Delta_raw[])
  for (p in 1:N_players) {
    Delta[p] <- Delta_raw[p] - mean_Delta
  }

  # Team shot-creation (attack) and opponent shot-suppression (defense) (centered)
  for (k in 1:N_teams) {
    att_raw[k] ~ dnorm(0, 1)
    def_raw[k] ~ dnorm(0, 4)
  }
  mean_att <- mean(att_raw[])
  mean_def <- mean(def_raw[])
  for (k in 1:N_teams) {
    lambda_team_att[k] <- att_raw[k] - mean_att
    lambda_team_def[k] <- def_raw[k] - mean_def
  }

  # Position effect: MF vs FW baseline
  beta_MF ~ dnorm(0, 1)

  # Likelihood
  for (i in 1:N) {

    # opponent sum (same style as your passing model)
    opponent_sum[i] <- inprod(lambda_team_def[], delta[i,])

    # linear predictor for shot rate
    eta[i] <- Delta[player_id[i]] +
              lambda_team_att[team_id[i]] -
              opponent_sum[i] +
              beta_MF * is_MF[i]

    # exposure offset - account for minutes
    log(mu[i]) <- log(minutes[i]) + eta[i]

    # observation
    y_shots[i] ~ dpois(mu[i])
  }
}
"

jags_model_shots_poi <- jags.model(
  textConnection(model_pois_shots),
  data = data_jags_shots,
  n.chains = 3,
  n.adapt = 170 # 2000
)

update(jags_model_shots_poi, 2300)  # 30000

params_shots_poi <- c("Delta", "lambda_team_att", "lambda_team_def", "beta_MF")
samples_shots_poi <- coda.samples(
  jags_model_shots_poi,
  variable.names = params_shots_poi,
  n.iter = 15000, # 170000
  thin = 1
) 

# player ability table (Poisson)
sum_stats_poi <- summary(samples_shots_poi)$statistics
sum_quants_poi <- summary(samples_shots_poi)$quantiles

Delta_idx <- grep("^Delta\\[", rownames(sum_stats_poi))
Delta_stats <- sum_stats_poi[Delta_idx, ]

player_names_shots <- levels(factor(df_shots$player))

player_table_shots_poi <- data.frame(
  player        = player_names_shots,
  ability_mean  = Delta_stats[, "Mean"],
  ability_sd    = Delta_stats[, "SD"],
  ability_lower = sum_quants_poi[Delta_idx, "2.5%"],
  ability_upper = sum_quants_poi[Delta_idx, "97.5%"]
)

View(player_table_shots_poi)

# Convergence checks:
# 1) Gelman-Rubin:
gd_uni_shots_poi <- gelman.diag(samples_shots_poi, autoburnin = FALSE, multivariate = FALSE)
max_psrf_shots_poi <- max(gd_uni_shots_poi$psrf[, "Point est."], na.rm = TRUE)
max_psrf_shots_poi

#2) Effective Sample Size:
ess_shots_poi <- effectiveSize(samples_shots_poi)
summary(ess_shots_poi)

#3) Autocorrelation plot:
autocorr.plot(samples_shots_poi)

#4) Geweke plot:
g_shots_poi <- geweke.diag(samples_shots_poi)
geweke.plot(samples_shots_poi)

#5) Traceplot:
plot(samples_shots_poi)

RonaldoPoisShots<- samples_shots_poi[, "Delta[40]"]
heidel.diag(RonaldoPoisShots)

AmrabatPoisShots <- samples_shots_poi[, "Delta[156]"]
heidel.diag(AmrabatPoisShots)

#Amrabat 156, Ronaldo 40

# combine the above code with the observed data for comparison
# =========================

# combine chains into one matrix
S_shots_poi <- as.matrix(samples_shots_poi)

# split into parameter blocks
Delta_draws_shots <- S_shots_poi[, grep("^Delta\\[", colnames(S_shots_poi)), drop = FALSE]
att_draws_shots   <- S_shots_poi[, grep("^lambda_team_att\\[", colnames(S_shots_poi)), drop = FALSE]
def_draws_shots   <- S_shots_poi[, grep("^lambda_team_def\\[", colnames(S_shots_poi)), drop = FALSE]
beta_MF_draws     <- S_shots_poi[, "beta_MF"]

# ids + exposure
pid   <- df_shots$shot_player_id
tid   <- df_shots$team_id
mins  <- df_shots$minutes_90
is_MF <- df_shots$is_MF_shots

# opponent defensive contribution
opp_sum_draws_shots <- def_draws_shots %*% t(delta_shots)

# build eta (beta_MF replicated across observations)
eta_draws_shots <- Delta_draws_shots[, pid, drop = FALSE] +
  att_draws_shots[, tid, drop = FALSE] -
  opp_sum_draws_shots +
  beta_MF_draws %*% t(is_MF)

# expected shots per observation
mu_draws_shots <- sweep(exp(eta_draws_shots), 2, mins, `*`)

# tournament totals per player
N_players_shots <- max(df_shots$shot_player_id)

mu_player_draws_shots <- sapply(1:N_players_shots, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_shots[, cols, drop = FALSE])
})

# summarise posterior per player
player_pred_total_shots_poi <- data.frame(
  player = N_players_shots,
  pred_shots_total_mean  = apply(mu_player_draws_shots, 2, mean),
  pred_shots_total_lower = apply(mu_player_draws_shots, 2, quantile, probs = 0.025),
  pred_shots_total_upper = apply(mu_player_draws_shots, 2, quantile, probs = 0.975)
)

# observed totals + minutes totals for comparison
obs_shots_total <- df_shots %>%
  group_by(shot_player_id) %>%
  summarise(
    obs_shots_total = sum(shots),
    minutes_90_total = sum(minutes_90),
    .groups = "drop"
  )



# join predicted and observed data + overall RMSE/MAE
player_pred_total_shots_poi <- player_pred_total_shots_poi %>%
  mutate(shot_player_id = 1:N_players_shots) %>%
  left_join(obs_shots_total, by = "shot_player_id") %>%
  arrange(desc(pred_shots_total_mean)) %>%


# compute overall model MAE and RMSE (single values for model)
errors_shots_poi <- player_pred_total_int_nb$obs_int_total - 
  player_pred_total_int_nb$pred_int_total_mean

MAE_model_shots_poi  <- mean(abs(errors_shots_poi))
RMSE_model_shots_poi <- sqrt(mean(errors_shots_poi^2))

MAE_model_shots_poi
RMSE_model_shots_poi

View(player_pred_total_shots_poi)

# ---------------------------- NEGATIVE BINOMIAL SHOTS (CENTERED)
model_nb_shots <- "
model {

  # Player shot volume ability (centered)
  for (p in 1:N_players) {
    Delta_raw[p] ~ dnorm(0, 1)
  }
  mean_Delta <- mean(Delta_raw[])
  for (p in 1:N_players) {
    Delta[p] <- Delta_raw[p] - mean_Delta
  }

  # Team shot-creation (attack) and opponent shot-suppression (defense) (centered)
  for (k in 1:N_teams) {
    att_raw[k] ~ dnorm(0, 1)
    def_raw[k] ~ dnorm(0, 4)
  }
  mean_att <- mean(att_raw[])
  mean_def <- mean(def_raw[])
  for (k in 1:N_teams) {
    lambda_team_att[k] <- att_raw[k] - mean_att
    lambda_team_def[k] <- def_raw[k] - mean_def
  }

  # Position effect: MF vs FW baseline
  beta_MF ~ dnorm(0, 1)

  # Overdispersion parameter (Negative Binomial)
  r ~ dgamma(0.01, 0.01)

  # Likelihood
  for (i in 1:N) {

    # opponent sum
    opponent_sum[i] <- inprod(lambda_team_def[], delta[i,])

    # linear predictor for shot rate
    eta[i] <- Delta[player_id[i]] +
              lambda_team_att[team_id[i]] -
              opponent_sum[i] +
              beta_MF * is_MF[i]

    # exposure offset - account for minutes
    log(mu[i]) <- log(minutes[i]) + eta[i]

    # Negative binomial with mean mu[i]
    p_nb[i] <- r / (r + mu[i])
    y_shots[i] ~ dnegbin(p_nb[i], r)
  }
}
"

jags_model_shots_nb <- jags.model(
  textConnection(model_nb_shots),
  data = data_jags_shots,
  n.chains = 3,
  n.adapt = 1000
)

update(jags_model_shots_nb, 2400)

params_shots_nb <- c("Delta", "lambda_team_att", "lambda_team_def", "beta_MF", "r")
samples_shots_nb <- coda.samples(
  jags_model_shots_nb,
  variable.names = params_shots_nb,
  n.iter = 16000,
  thin = 1
)

# player ability table (NB)
sum_stats_nb <- summary(samples_shots_nb)$statistics
sum_quants_nb <- summary(samples_shots_nb)$quantiles

Delta_idx_nb <- grep("^Delta\\[", rownames(sum_stats_nb))
Delta_stats_nb <- sum_stats_nb[Delta_idx_nb, ]

player_names_shots <- levels(factor(df_shots$player))

player_table_shots_nb <- data.frame(
  player        = player_names_shots,
  ability_mean  = Delta_stats_nb[, "Mean"],
  ability_sd    = Delta_stats_nb[, "SD"],
  ability_lower = sum_quants_nb[Delta_idx_nb, "2.5%"],
  ability_upper = sum_quants_nb[Delta_idx_nb, "97.5%"]
)

View(player_table_shots_nb)


# compare woth observed data NEG BIN
# =========================

S_shots_nb <- as.matrix(samples_shots_nb)

Delta_draws_shots_nb <- S_shots_nb[, grep("^Delta\\[", colnames(S_shots_nb)), drop = FALSE]
att_draws_shots_nb   <- S_shots_nb[, grep("^lambda_team_att\\[", colnames(S_shots_nb)), drop = FALSE]
def_draws_shots_nb   <- S_shots_nb[, grep("^lambda_team_def\\[", colnames(S_shots_nb)), drop = FALSE]
beta_MF_draws_nb     <- S_shots_nb[, "beta_MF"]
r_draws_shots_nb     <- S_shots_nb[, "r"]  # not used unless simulating

# opponent defensive contribution
opp_sum_draws_shots_nb <- def_draws_shots_nb %*% t(delta_shots)

# build eta
eta_draws_shots_nb <- Delta_draws_shots_nb[, pid, drop = FALSE] +
  att_draws_shots_nb[, tid, drop = FALSE] -
  opp_sum_draws_shots_nb +
  beta_MF_draws_nb %*% t(is_MF)

# expected shots per observation
mu_draws_shots_nb <- sweep(exp(eta_draws_shots_nb), 2, mins, `*`)

# Number of players for shots
N_players_shots <- max(df_shots$shot_player_id)


# tournament totals per player
mu_player_draws_shots_nb <- sapply(1:N_players_shots, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_shots_nb[, cols, drop = FALSE])
})

player_names_shots <- levels(factor(df_shots$player))


# summarise posterior per player
player_pred_total_shots_nb <- data.frame(
  player = player_names_shots,
  pred_shots_total_mean  = apply(mu_player_draws_shots_nb, 2, mean),
  pred_shots_total_lower = apply(mu_player_draws_shots_nb, 2, quantile, probs = 0.025),
  pred_shots_total_upper = apply(mu_player_draws_shots_nb, 2, quantile, probs = 0.975)
)


# join observed + RMSE/MAE
player_pred_total_shots_nb <- player_pred_total_shots_nb %>%
  mutate(shot_player_id = 1:N_players_shots) %>%
  left_join(obs_shots_total, by = "shot_player_id") %>%
  arrange(desc(pred_shots_total_mean)) 

  
# compute overall model MAE and RMSE (single values for model)
errors_shots_nb <- player_pred_total_int_nb$obs_int_total - 
player_pred_total_int_nb$pred_int_total_mean

MAE_model_shots_nb  <- mean(abs(errors_shots_nb))
RMSE_model_shots_nb <- sqrt(mean(errors_shots_nb^2))

MAE_model_shots_nb
RMSE_model_shots_nb


View(player_pred_total_shots_nb)

# =========================
# CONVERGENCE: SHOTS NB MODEL
# =========================

# 1. Gelman-Rubin R-hat overall
gd_uni <- gelman.diag(samples_shots_nb, autoburnin = FALSE, multivariate = FALSE)

# R-hat worst case
max_psrf <- max(gd_uni$psrf[, "Point est."], na.rm = TRUE)
max_psrf

psrf <- as.data.frame(gd_uni$psrf)        
psrf$param <- rownames(psrf)

# 1) Which parameters exceed 1.05 R-hat?
bad_point <- psrf %>%
  filter(`Point est.` > 1.05) %>%
  arrange(desc(`Point est.`))

bad_point

# 2) (stricter) Which parameters exceed 1.05 on Upper CI?
bad_upper <- psrf %>%
  filter(`Upper C.I.` > 1.05) %>%
  arrange(desc(`Upper C.I.`))

bad_upper

# 2. Heidelberger–Welch stationarity
hw <- heidel.diag(samples_shots_nb)
table(hw[,1])   # stationarity pass/fail
table(hw[,2])   # halfwidth pass/fail

# 3. Effective Sample Size
ess <- effectiveSize(samples_shots_nb)
summary(ess)

# 4. Autocorrelation plots
autocorr.plot(samples_shots_nb)

# 5. Geweke diagnostics
g <- geweke.diag(samples_shots_nb)

# optional plot
geweke.plot(samples_shots_nb)

hw <- heidel.diag(samples_nb)

# collect stationarity + halfwidth results
stationarity <- unlist(lapply(hw, function(x) x[, "stest"]))
halfwidth    <- unlist(lapply(hw, function(x) x[, "htest"]))

# which failed stationarity?
fail_stationarity <- names(stationarity[stationarity != "passed"])

# which failed halfwidth?
fail_halfwidth <- names(halfwidth[halfwidth != "passed"])

fail_stationarity
fail_halfwidth

# what exactly is the ability mean?
# bayes factor between the 2 models, and if it is done only for
# nested models. it is calculated using the marginal likelihood
# may need to use bridge sampling to achieve bayes factor between the 2 models
# interacting abilities of teams (lambda_team_att?)
 
load("my_workspace.RData")
save.image("my_workspace.RData")
