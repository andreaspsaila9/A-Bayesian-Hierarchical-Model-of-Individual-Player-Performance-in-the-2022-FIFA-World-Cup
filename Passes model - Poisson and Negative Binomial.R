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
    minutes_90 = minutes_90s,                             # minutes per 90 for all players since some players played less
    is_DF = ifelse(position.x == "DF", 1, 0),             # dummy for defenders
    is_MF = ifelse(position.x == "MF", 1, 0)              # dummy for midfielders (FW is baseline by elimination)
  )

delta <- as.matrix(df_model[, 7:38]) # teams columns in dataset: from Argentina to Wales

# checks for clean data 
sapply(df_model[, c("passes_completed", "interceptions","shots" "minutes_90",
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
# RE RUN CAUSE OF UTF-8 ENCODING WRONG NAMES



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

load("my_workspace.RData")
save.image("my_workspace.RData")