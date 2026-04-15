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
sapply(df_model[, c("passes_completed", "interceptions","shots", "minutes_90",
                    "player_id", "team_id", "is_DF", "is_MF")],
       function(x) sum(is.na(x)))

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

# Subset df_model using the SAME 
# Remember that shots model uses a subset so different player indices
df_shots <- df_model[keep_idx, ] %>%
  mutate(
    shot_player_id = match(player, unique(player)), # = as.numeric(factor(player)),
    is_MF_shots = ifelse(position.x == "MF", 1, 0),
    player_id = as.numeric(factor(player)),
    team_id   = as.numeric(factor(team))
  )
nrow(df_shots)


# Subset delta using the SAME sample
delta_shots <- as.matrix(delta[keep_idx, ])

# check for equality
stopifnot(nrow(df_shots) == nrow(delta_shots))
N_teams <- ncol(delta_shots)
stopifnot(max(df_shots$team_id) <= N_teams)

# Check player number 
player_lookup_attmid <- df_shots %>%
  distinct(shot_player_id, player) %>%
  arrange(shot_player_id)

player_lookup_attmid

#Amrabat 103, Ronaldo 113


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
  player_id = df_shots$shot_player_id,
  team_id   = df_shots$team_id,
  delta     = delta_shots,
  tau = df_shots$tau,
  tau_star = df_shots$tau_star,
  N         = nrow(df_shots),
  N_players = max(df_shots$shot_player_id),
  N_teams   = N_teams
)




model_pois_shots <- "
model {

  # Priors on player effects
  for (p in 1:N_players) {
    Delta_star[p] ~ dnorm(m, s)
    #trick to sum-to-zero constraints
    Delta[p] <- Delta_star[p] - mean(Delta_star[])
  }

  # Priors on team effects
  for (k in 1:N_teams) {
  #  lambda_e_star[k] ~ dnorm(0, 0.0001) # like team effects in whitaker and baio and blangiardo
  #  lambda_Edive_star[k] ~ dnorm(0, 0.0001) # like team effects in whitaker and baio and blangiardo
    
    lambda_e_star[k] ~ dnorm(mu.lambda_e, tau.lambda_e)  # like team effects in whitaker and baio and blangiardo
lambda_Edive_star[k] ~ dnorm(mu.lambda_Edive, tau.lambda_Edive) # like team effects in whitaker and baio and blangiardo

  # sum-to-zero constraints
    lambda_e[k] <- lambda_e_star[k] - mean(lambda_e_star[])
    lambda_Edive[k] <- lambda_Edive_star[k] - mean(lambda_Edive_star[])
  }

# baio and blangiardo impose sum-to-zero constraints to their team effects, so:

  # Likelihood
  for (i in 1:N) {

    y_shots[i] ~ dpois(eta[i] * tau[i])

    opp_sum[i] <- inprod(lambda_Edive[], delta[i,]) # total opponent ability

    log(eta[i]) <- Delta[player_id[i]] +
                   tau[i] * lambda_e[team_id[i]] -
                   tau_star[i] * opp_sum[i]
  }
  # inspired by Baio and Blangiardo code, priors in the random effects
  m~dnorm(0,0.0001)
  s~dgamma(0.1,0.01) # controls spread of player effects
  
# priors on the random effects
mu.lambda_e ~ dnorm(0,0.0001)
mu.lambda_Edive ~ dnorm(0,0.0001) # mean 0 with precision very small and weak, we do not 
# know anything about team effects apart from the fact that they 
# are around zero, but very weakly-informative
tau.lambda_e ~ dgamma(.01,.01)
tau.lambda_Edive ~ dgamma(.01,.01)
# priors are very broad, we let the data decide if teams are similar or different
}
"

# now the jags model is initiated with 3 parallel markov chains which help us 
# check for convergence
jags_model <- jags.model( 
  textConnection(model_pois_shots),
  data = data_jags_shots,
  n.chains = 3,
  n.adapt = 1000
)

update(jags_model, 6000) # account for burn-in

 # params_shots_pois <- c("Delta", "lambda_e", "lambda_Edive") # , "mu")
# params_shots_pois <- c("Delta", "Delta_star", "lambda_e", "lambda_Edive", "m", "s")
params_shots_pois <- c(
  "Delta", "lambda_e", "lambda_Edive",
  "mu.lambda_e", "mu.lambda_Edive", "tau.lambda_e", "tau.lambda_Edive","m", "s"
)

samples_shots_pois <- coda.samples(jags_model, variable.names = params_shots_pois, n.iter = 30000, thin = 1)



# player ability table (Poisson)
sum_stats_poi <- summary(samples_shots_pois)$statistics
sum_quants_poi <- summary(samples_shots_pois)$quantiles

Delta_idx <- grep("^Delta\\[", rownames(sum_stats_poi))
Delta_stats <- sum_stats_poi[Delta_idx, ]

player_names_shots <- player_lookup_attmid$player#levels(factor(df_shots$player))

# player_table_shots_poi <- data.frame(
 # player        = player_names_shots,
#  ability_mean  = Delta_stats[, "Mean"],
 # ability_sd    = Delta_stats[, "SD"],
  #ability_lower = sum_quants_poi[Delta_idx, "2.5%"],
  #ability_upper = sum_quants_poi[Delta_idx, "97.5%"]
# )

samples_shots_pois_mcmc <- as.mcmc(as.matrix(samples_shots_pois)) # combine all chains into one

Delta_hpd <- HPDinterval(samples_shots_pois_mcmc, prob = 0.95)
Delta_hpd <- Delta_hpd[grep("^Delta\\[", rownames(Delta_hpd)), ]

player_table_shots_poi <- data.frame(
  player        = player_names_shots,
  ability_mean  = Delta_stats[, "Mean"],
  ability_sd    = Delta_stats[, "SD"],
  ability_lower = Delta_hpd[, "lower"],
  ability_upper = Delta_hpd[, "upper"]
)

View(player_table_shots_poi)

# Convergence checks:
# 1) Gelman-Rubin:
gd_uni_shots_poi <- gelman.diag(samples_shots_pois, autoburnin = FALSE, multivariate = FALSE)
max_psrf_shots_poi <- max(gd_uni_shots_poi$psrf[, "Point est."], na.rm = TRUE)
max_psrf_shots_poi

psrf_vals <- gd_uni_shots_poi$psrf[, "Point est."]

# all parameters above 1.05
psrf_vals[psrf_vals > 1.05]

#2) Effective Sample Size:
ess_shots_pois <- effectiveSize(samples_shots_pois)
summary(ess_shots_pois)

#3) Autocorrelation plot:
autocorr.plot(samples_shots_pois)

#4) Geweke plot:
g_shots_poi <- geweke.diag(samples_shots_pois)
geweke.plot(samples_shots_pois)

#5) Traceplot:
plot(samples_shots_pois)

#Amrabat 103, Ronaldo 113, mbappe 56

MbappePoisShots <- samples_shots_pois[, "Delta[56]"]

heidel.diag(MbappePoisShots)
plot(MbappePoisShots, main = expression("Plot for " * Delta[56] * " (Kylian Mbappé)"))
gelman.diag(MbappePoisShots)
geweke.diag(MbappePoisShots)
effectiveSize(MbappePoisShots)

AmrabatPoisShots <- samples_shots_pois[, "Delta[103]"]
heidel.diag(AmrabatPoisShots)
plot(AmrabatPoisShots, main = expression("Plot for " * Delta[103] * " (Sofyan Amrabat)"))
gelman.diag(AmrabatPoisShots)
geweke.diag(AmrabatPoisShots)
effectiveSize(AmrabatPoisShots)



# combine the above code with the observed data for comparison
# =========================

# combine chains into one matrix
S_shots_poi <- as.matrix(samples_shots_pois)

# split into parameter blocks

# Delta_draws_shots <- S_shots_poi[, grep("^Delta\\[", colnames(S_shots_poi)), drop = FALSE]
# att_star_draws   <- S_shots_poi[, grep("^lambda_e_star\\[", colnames(S_shots_poi)), drop = FALSE]
# def_star_draws   <- S_shots_poi[, grep("^lambda_Edive_star\\[", colnames(S_shots_poi)), drop = FALSE]

Delta_draws_shots <- S_shots_poi[, grep("^Delta\\[", colnames(S_shots_poi)), drop = FALSE]
team_star_draws   <- S_shots_poi[, grep("^lambda_e\\[", colnames(S_shots_poi)), drop = FALSE]
oppteam_star_draws   <- S_shots_poi[, grep("^lambda_Edive\\[", colnames(S_shots_poi)), drop = FALSE]

# ids + exposure
pid      <- df_shots$shot_player_id
tid      <- df_shots$team_id
tau      <- df_shots$tau
tau_star <- df_shots$tau_star

# opponent defence
opp_sum_draws_shots <- oppteam_star_draws %*% t(delta_shots)

# build eta 

log_eta_draws <- Delta_draws_shots[, pid, drop = FALSE] +
  sweep(team_star_draws[, tid, drop = FALSE], 2, tau, `*`) -
  sweep(opp_sum_draws_shots, 2, tau_star, `*`)

eta_draws_shots <- exp(log_eta_draws)

# expected shots per observation
mu_draws_shots <- sweep(eta_draws_shots, 2, tau, `*`) 

# tournament totals per player
N_players_shots <- data_jags_shots$N_players

mu_player_draws_shots <- sapply(1:N_players_shots, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_shots[, cols, drop = FALSE])
})

# summarise posterior per player
 #player_pred_total_shots_poi <- data.frame(
#  player = player_names_shots,
#  pred_shots_total_mean  = apply(mu_player_draws_shots, 2, mean),
#  pred_shots_total_lower = apply(mu_player_draws_shots, 2, quantile, probs = 0.025),
#  pred_shots_total_upper = apply(mu_player_draws_shots, 2, quantile, probs = 0.975)
# )
# summarise posterior per player using 95% HPD intervals
hpd_pred_shots <- t(sapply(1:ncol(mu_player_draws_shots), function(j) {
  hpd <- HPDinterval(as.mcmc(mu_player_draws_shots[, j]), prob = 0.95)
  c(lower = hpd[1, "lower"], upper = hpd[1, "upper"])
}))

player_pred_total_shots_poi <- data.frame(
  player = player_names_shots,
  pred_shots_total_mean  = apply(mu_player_draws_shots, 2, mean),
  pred_shots_total_lower = hpd_pred_shots[, "lower"],
  pred_shots_total_upper = hpd_pred_shots[, "upper"]
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
  arrange(desc(pred_shots_total_mean)) 

# compute overall model MAE and RMSE (single values for model)
errors_shots_poi <- player_pred_total_shots_poi$obs_shots_total - 
  player_pred_total_shots_poi$pred_shots_total_mean

MAE_model_shots_poi  <- mean(abs(errors_shots_poi))
RMSE_model_shots_poi <- sqrt(mean(errors_shots_poi^2))

MAE_model_shots_poi
RMSE_model_shots_poi

View(player_pred_total_shots_poi)

# ---------------------------- NEGATIVE BINOMIAL SHOTS (CENTERED)
model_nb_shots <- "
model {

  # Priors on player effects
  for (p in 1:N_players) {
    Delta_star[p] ~ dnorm(m, s)
    #trick to sum-to-zero constraints
    Delta[p] <- Delta_star[p] - mean(Delta_star[])
  }

  # Priors on team effects
  for (k in 1:N_teams) {
  #  lambda_e_star[k] ~ dnorm(0, 0.0001) # like team effects in whitaker and baio and blangiardo
  #  lambda_Edive_star[k] ~ dnorm(0, 0.0001) # like team effects in whitaker and baio and blangiardo
    
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

    y_shots[i] ~ dnegbin(p[i], r)
  }
  # inspired by Baio and Blangiardo code, priors in the random effects
  m~dnorm(0,0.0001)
  s~dgamma(0.1,0.01) # controls spread of player effects
  
# priors on the random effects
mu.lambda_e ~ dnorm(0,0.0001)
mu.lambda_Edive ~ dnorm(0,0.0001) # mean 0 with precision very small and weak, we do not 
# know anything about team effects apart from the fact that they 
# are around zero, but very weakly-informative
tau.lambda_e ~ dgamma(.01,.01)
tau.lambda_Edive ~ dgamma(.01,.01)
# priors are very broad, we let the data decide if teams are similar or different

r ~ dgamma(.01,.01) # since r must be positive and controls overdispersion
}
" 

jags_model_shots_nb <- jags.model(
  textConnection(model_nb_shots),
  data = data_jags_shots,
  n.chains = 3,
  n.adapt = 1000
)


update(jags_model_shots_nb, 6000)

params_shots_nb <- c(
  "Delta", "lambda_e", "lambda_Edive",
 "tau.lambda_e", "tau.lambda_Edive","r"
)

samples_shots_nb <- coda.samples(
  jags_model_shots_nb,
  variable.names = params_shots_nb,
  n.iter = 30000, 
  thin = 1
)

# player ability table (NB)
sum_stats_nb <- summary(samples_shots_nb)$statistics
sum_quants_nb <- summary(samples_shots_nb)$quantiles

samples_shots_nb_mcmc <- as.mcmc(as.matrix(samples_shots_nb))

Delta_idx_nb <- grep("^Delta\\[", rownames(sum_stats_nb))
Delta_stats_nb <- sum_stats_nb[Delta_idx_nb, ]

Delta_hpd_shots_nb <- HPDinterval(samples_shots_nb_mcmc, prob = 0.95)
Delta_hpd_shots_nb <- Delta_hpd_shots_nb[grep("^Delta\\[", rownames(Delta_hpd_shots_nb)), ]

player_table_shots_nb <- data.frame(
  player        = player_names_shots,
  ability_mean  = Delta_stats_nb[, "Mean"],
  ability_sd    = Delta_stats_nb[, "SD"],
  ability_lower = Delta_hpd_shots_nb[, "lower"],
  ability_upper = Delta_hpd_shots_nb[, "upper"]
)

View(player_table_shots_nb)


# compare woth observed data NEG BIN
# =========================

S_shots_nb <- as.matrix(samples_shots_nb)

Delta_draws_shots_nb <- S_shots_nb[, grep("^Delta\\[", colnames(S_shots_nb)), drop = FALSE]
team_draws_shots_nb   <- S_shots_nb[, grep("^lambda_e\\[", colnames(S_shots_nb)), drop = FALSE]
oppteam_draws_shots_nb   <- S_shots_nb[, grep("^lambda_Edive\\[", colnames(S_shots_nb)), drop = FALSE]

r_draws_shots_nb     <- S_shots_nb[, "r"]  

# ids + exposure
pid      <- df_shots$shot_player_id
tid      <- df_shots$team_id
tau      <- df_shots$tau
tau_star <- df_shots$tau_star


# opponent defence, opponent team defensive contribution for each draw and observation
opp_sum_draws_shots_nb <- oppteam_draws_shots_nb %*% t(delta_shots)

# build eta 

log_eta_draws_nb <- Delta_draws_shots_nb[, pid, drop = FALSE] +
  sweep(team_draws_shots_nb[, tid, drop = FALSE], 2, tau, `*`) -
  sweep(opp_sum_draws_shots_nb, 2, tau_star, `*`)

eta_draws_shots_nb <- exp(log_eta_draws_nb)

# expected shots per observation  or posterior expected mean at the observation level.
mu_observed_draws_shots <- sweep(eta_draws_shots_nb, 2, tau, `*`) 

# Number of players for shots
# tournament totals per player
N_players_shots <- data_jags_shots$N_players


# tournament totals per player
mu_player_total_draws_shots_nb <- sapply(1:N_players_shots, function(p) {
  cols <- which(pid == p)
  rowSums(mu_observed_draws_shots[, cols, drop = FALSE])
})


# summarise posterior per player using 95% HPD intervals
hpd_pred_shots_nb <- t(sapply(1:ncol(mu_player_total_draws_shots_nb), function(j) {
  hpd <- HPDinterval(as.mcmc(mu_player_total_draws_shots_nb[, j]), prob = 0.95)
  c(lower = hpd[1, "lower"], upper = hpd[1, "upper"])
}))

player_pred_total_shots_nb <- data.frame(
  player = player_names_shots,
  pred_shots_total_mean  = apply(mu_player_total_draws_shots_nb, 2, mean),
  pred_shots_total_lower = hpd_pred_shots_nb[, "lower"],
  pred_shots_total_upper = hpd_pred_shots_nb[, "upper"]
)

# join observed + RMSE/MAE
player_pred_total_shots_nb <- player_pred_total_shots_nb %>%
  mutate(shot_player_id = 1:N_players_shots) %>%
  left_join(obs_shots_total, by = "shot_player_id") %>%
  arrange(desc(pred_shots_total_mean)) 

# compute overall model MAE and RMSE (single values for model)
errors_shots_nb <- player_pred_total_shots_nb$obs_shots_total - 
  player_pred_total_shots_nb$pred_shots_total_mean

MAE_model_shots_nb  <- mean(abs(errors_shots_nb))
RMSE_model_shots_nb <- sqrt(mean(errors_shots_nb^2))

MAE_model_shots_nb
RMSE_model_shots_nb

View(player_pred_total_shots_nb)

# =========================
# CONVERGENCE: SHOTS NB MODEL
# =========================

# 1. Gelman-Rubin R-hat overall
gd_uni_shots_nb <- gelman.diag(samples_shots_nb, autoburnin = FALSE, multivariate = FALSE)

# R-hat worst case
max_psrf_shots_nb <- max(gd_uni_shots_nb$psrf[, "Point est."], na.rm = TRUE)
max_psrf_shots_nb

psrf <- as.data.frame(gd_uni_shots_nb$psrf)        
psrf$param <- rownames(psrf)

# 1) Which parameters exceed 1.05 R-hat?
bad_point <- psrf %>%
  filter(`Point est.` > 1.05) %>%
  arrange(desc(`Point est.`))

bad_point



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

hw <- heidel.diag(samples_shots_nb)


MbappeNBGShots <- samples_shots_nb[, "Delta[56]"]

heidel.diag(MbappeNBGShots)
plot(MbappeNBGShots, main = expression("Plot for " * Delta[56] * " (Kylian Mbappé)"))
gelman.diag(MbappeNBGShots)
geweke.diag(MbappeNBGShots)
effectiveSize(MbappeNBGShots)

AmrabatNegBinShots <- samples_shots_nb[, "Delta[103]"]
heidel.diag(AmrabatNegBinShots)
plot(AmrabatNegBinShots, main = expression("Plot for " * Delta[103] * " (Sofyan Amrabat)"))
gelman.diag(AmrabatNegBinShots)
geweke.diag(AmrabatNegBinShots)
effectiveSize(AmrabatNegBinShots)

# what exactly is the ability mean?
# bayes factor between the 2 models, and if it is done only for
# nested models. it is calculated using the marginal likelihood
# may need to use bridge sampling to achieve bayes factor between the 2 models
# interacting abilities of teams (lambda_team_att?)

load("my_workspace.RData")
save.image("my_workspace.RData")
