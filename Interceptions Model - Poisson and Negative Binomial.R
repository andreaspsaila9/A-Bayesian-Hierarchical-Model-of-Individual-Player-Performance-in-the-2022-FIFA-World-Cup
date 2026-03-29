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
sapply(df_model[, c("passes_completed", "interceptions","shots", "minutes_90",
                    "player_id", "team_id", "is_DF", "is_MF")],
       function(x) sum(is.na(x)))
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
  n.adapt = 420 
)

update(jags_model_int_poi, 4200)  

params_int_poi <- c("Delta", "lambda_team_int", "lambda_opp_att")

samples_int_poi <- coda.samples(jags_model_int_poi, variable.names = params_int_poi, n.iter = 30000, thin = 1) 
# 6) Summaries + defender ability table
summary_stats_int  <- summary(samples_int_poi)$statistics
summary_quants_int <- summary(samples_int_poi)$quantiles

Delta_idx <- grep("^Delta\\[", rownames(summary_stats_int))
Delta_stats <- summary_stats_int[Delta_idx, ]

player_lookup_defmid <- df_defmid %>%
  distinct(defmid_player_id, player) %>%
  arrange(defmid_player_id)

player_table_poisson_int <- data.frame(
  player        = player_lookup_defmid$player,
  ability_mean  = Delta_stats[, "Mean"],
  ability_sd    = Delta_stats[, "SD"],
  ability_lower = summary_quants_int[Delta_idx, "2.5%"],
  ability_upper = summary_quants_int[Delta_idx, "97.5%"]
)

View(player_table_poisson_int)

# traceplots
plot(samples_int_poi)



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

HakimiPoisInt <- samples_int_poi[, "Delta[8]"] # "Delta[131]"
heidel.diag(HakimiPoisInt)
plot(HakimiPoisInt, main = expression("Plot for " * Delta[8] * " (Achraf Hakimi)"))
gelman.diag(HakimiPoisInt)
geweke.diag(HakimiPoisInt)
effectiveSize(HakimiPoisInt)

AmrabatPoisInt <- samples_int_poi[, "Delta[212]"]
heidel.diag(AmrabatPoisInt)
plot(AmrabatPoisInt, main = expression("Plot for " * Delta[212] * " (Sofyan Amrabat)"))
gelman.diag(AmrabatPoisInt)
geweke.diag(AmrabatPoisInt)
effectiveSize(AmrabatPoisInt)
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
  player = player_lookup_defmid$player,
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
  n.adapt = 200 
)

update(jags_model_int_nb, 1000) 

params_int_nb <- c("Delta", "lambda_team_int", "lambda_opp_att", "r")
samples_int_nb <- coda.samples(
  jags_model_int_nb,
  variable.names = params_int_nb,
  n.iter = 10000, 
  thin = 1
)

# 5) Defender ability table (NB)
summary_stats_int_nb  <- summary(samples_int_nb)$statistics
summary_quants_int_nb <- summary(samples_int_nb)$quantiles

Delta_idx_nb   <- grep("^Delta\\[", rownames(summary_stats_int_nb))
Delta_stats_nb <- summary_stats_int_nb[Delta_idx_nb, ]

player_lookup_defmid <- df_defmid %>%
  distinct(defmid_player_id, player) %>%
  arrange(defmid_player_id)


player_table_nb_int <- data.frame(
  player        = player_lookup_defmid$player,
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

# Negative Binomial model: Achraf Hakimi
HakimiNBInt <- samples_int_nb[, "Delta[8]"]
heidel.diag(HakimiNBInt)
plot(HakimiNBInt, main = expression("Plot for " * Delta[8] * " (Achraf Hakimi)"))
gelman.diag(HakimiNBInt)
geweke.diag(HakimiNBInt)
effectiveSize(HakimiNBInt)

# Negative Binomial model: Sofyan Amrabat
AmrabatNBInt <- samples_int_nb[, "Delta[212]"]
heidel.diag(AmrabatNBInt)
plot(AmrabatNBInt, main = expression("Plot for " * Delta[212] * " (Sofyan Amrabat)"))
gelman.diag(AmrabatNBInt)
geweke.diag(AmrabatNBInt)
effectiveSize(AmrabatNBInt)

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
pid  <- df_defmid$defmid_player_id
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
N_players_defmid <- max(df_defmid$defmid_player_id)

# posterior distribution of tournament totals for each defender (sum across their obs)
mu_player_draws_int_nb <- sapply(1:N_players_defmid, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_int_nb[, cols, drop = FALSE])
})

# summarise posterior per defender
player_pred_total_int_nb <- data.frame(
  player = player_lookup_defmid$player,
  pred_int_total_mean  = apply(mu_player_draws_int_nb, 2, mean),
  pred_int_total_lower = apply(mu_player_draws_int_nb, 2, quantile, probs = 0.025),
  pred_int_total_upper = apply(mu_player_draws_int_nb, 2, quantile, probs = 0.975)
)

# observed totals + minutes totals for comparison
obs_totals_int <- df_defmid %>%
  group_by(defmid_player_id) %>%
  summarise(
    obs_int_total     = sum(interceptions),
    minutes_90_total  = sum(minutes_90),
    .groups = "drop"
  )

# join predicted and observed data
player_pred_total_int_nb <- player_pred_total_int_nb %>%
  mutate(defmid_player_id = 1:N_players_defmid) %>%
  left_join(obs_totals_int, by = "defmid_player_id") %>%
  arrange(desc(pred_int_total_mean))

# compute overall model MAE and RMSE (single values for model)
errors <- player_pred_total_int_nb$obs_int_total - 
  player_pred_total_int_nb$pred_int_total_mean

MAE_model_int_negbin  <- mean(abs(errors))
RMSE_model_int_negbin <- sqrt(mean(errors^2))

MAE_model_int_negbin
RMSE_model_int_negbin

View(player_pred_total_int_nb)

load("my_workspace.RData")
save.image("my_workspace.RData")