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


# ---------------------------- POISSON INTERCEPTIONS (CENTERED)

# 2) JAGS data
data_jags_int <- list(
  y_int     = df_defmid$interceptions,
  player_id = df_defmid$defmid_player_id,
  team_id   = df_defmid$team_id,
  delta     = delta_defmid,
  tau       = df_defmid$minutes_90,
  tau_star = df_defmid$tau_star,
  N         = nrow(df_defmid),
  N_players = max(df_defmid$defmid_player_id),
  N_teams   = ncol(delta_defmid)
)

model_pois_int <- "
model {
  
  # Priors on player effects
  for (p in 1:N_players) {
    Delta_star[p] ~ dnorm(m, s)
    #trick to sum-to-zero constraints
    Delta[p] <- Delta_star[p] - mean(Delta_star[])
  }
  
  # Priors on team effects
  for (k in 1:N_teams) {
    
    lambda_e_star[k] ~ dnorm(mu.lambda_e, tau.lambda_e)  # like team effects in whitaker and baio and blangiardo
    lambda_Edive_star[k] ~ dnorm(mu.lambda_Edive, tau.lambda_Edive) # like team effects in whitaker and baio and blangiardo
    
    # sum-to-zero constraints
    lambda_e[k] <- lambda_e_star[k] - mean(lambda_e_star[])
    lambda_Edive[k] <- lambda_Edive_star[k] - mean(lambda_Edive_star[])
  }
  
  # baio and blangiardo impose sum-to-zero constraints to their team effects, so:
  
  # Likelihood
  for (i in 1:N) {
    
    y_int[i] ~ dpois(eta[i] * tau[i])
    
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

jags_model_int_poi <- jags.model(
  textConnection(model_pois_int),
  data = data_jags_int,
  n.chains = 3,
  n.adapt = 1000 #1000
)

update(jags_model_int_poi, 9000) #9000?

params_int_poi <- c(
  "Delta", "lambda_e", "lambda_Edive",
  "tau.lambda_e", "tau.lambda_Edive" # "mu.lambda_e", "mu.lambda_Edive","m", "s"
)

samples_int_poi <- coda.samples(
  jags_model_int_poi,
  variable.names = params_int_poi,
  n.iter = 40000, # 40000?
  thin = 1
)

# player ability table (Poisson interceptions)
sum_stats_int_poi <- summary(samples_int_poi)$statistics
sum_quants_int_poi <- summary(samples_int_poi)$quantiles

Delta_idx_int <- grep("^Delta\\[", rownames(sum_stats_int_poi))
Delta_stats_int <- sum_stats_int_poi[Delta_idx_int, ]

player_lookup_defmid <- df_defmid %>%
  distinct(defmid_player_id, player) %>%
  arrange(defmid_player_id)

# player_table_poisson_int <- data.frame(
#  player        = player_lookup_defmid$player,
#  ability_mean  = Delta_stats_int[, "Mean"],
# ability_sd    = Delta_stats_int[, "SD"],
#ability_lower = sum_quants_int_poi[Delta_idx_int, "2.5%"],
#ability_upper = sum_quants_int_poi[Delta_idx_int, "97.5%"]
# )


samples_int_poi_mcmc <- as.mcmc(as.matrix(samples_int_poi))

Delta_hpd_poi_int <- HPDinterval(samples_int_poi_mcmc, prob = 0.95)
Delta_hpd_poi_int <- Delta_hpd_poi_int[grep("^Delta\\[", rownames(Delta_hpd_poi_int)), ]

player_table_poisson_int <- data.frame(
  player        = player_lookup_defmid$player,
  ability_mean  = Delta_stats_int[, "Mean"],
  ability_sd    = Delta_stats_int[, "SD"],
  ability_lower = Delta_hpd_poi_int[, "lower"],
  ability_upper = Delta_hpd_poi_int[, "upper"]
)

View(player_table_poisson_int)




#extracts posterior for latent player ability
Delta_stats_int[, ]

# checking for convergence

# 1. traceplots for all players - takes long
plot(samples_int_poi)

# checking for convergence
# 2. Gelman-Rubin R-hat overall
gd_uni_int <- gelman.diag(samples_int_poi, autoburnin = FALSE, multivariate = FALSE)

# Worst-case R-hat
max_psrf <- max(gd_uni_int$psrf[, "Point est."], na.rm = TRUE)
max_psrf

psrf <- as.data.frame(gd_uni_int$psrf)        
psrf$param <- rownames(psrf)

# 1) Which parameters exceed 1.05 R-hat?
bad_point <- psrf %>%
  filter(`Point est.` > 1.05) %>%
  arrange(desc(`Point est.`))

bad_point


# 3. ESS
ess <- effectiveSize(samples_int_poi)
summary(ess)

HakimiPoisInt <- samples_int_poi[, "Delta[8]"] # "Delta[131]"
heidel.diag(HakimiPoisInt)
plot(HakimiPoisInt, main = expression("Plot for " * Delta[8] * " (Achraf Hakimi)"))
gelman.diag(HakimiPoisInt)
geweke.diag(HakimiPoisInt)
effectiveSize(HakimiPoisInt)

AmrabatPoisInt <- samples_int_poi[, "Delta[213]"]
heidel.diag(AmrabatPoisInt)
plot(AmrabatPoisInt, main = expression("Plot for " * Delta[213] * " (Sofyan Amrabat)"))
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
team_draws_int  <- S_int[, grep("^lambda_e\\[", colnames(S_int)), drop = FALSE]
oppteam_draws_int   <- S_int[, grep("^lambda_Edive\\[", colnames(S_int)), drop = FALSE]

# player and team ids + minutes for defenders
pid  <- df_defmid$defmid_player_id
tid  <- df_defmid$team_id
tau      <- df_defmid$tau
tau_star <- df_defmid$tau_star


# matrix mult between posterior draw and opponent parameterto obtain opponent contribution
opp_sum_draws_int <- oppteam_draws_int %*% t(delta_defmid)

# build eta for each posterior draw and observation
# we add attacking opponent strength because the more attacks the more chance
# of interceptions by defenders
log_eta_draws <- Delta_draws_int[, pid, drop = FALSE] +
  sweep(team_draws_int[, tid, drop = FALSE], 2, tau, `*`) -
  sweep(opp_sum_draws_int, 2, tau_star, `*`)

eta_draws_int <- exp(log_eta_draws)

# expected interceptions per observation
mu_draws_int <- sweep(eta_draws_int, 2, tau, `*`) 


# total expected interceptions over tournament per defender/midfielder
N_players_int <- data_jags_int$N_players

#posterior distributions for posterior draws for each defender
mu_player_draws_int <- sapply(1:N_players_int, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_int[, cols, drop = FALSE])
})

# summarise posterior per defender
# player_pred_total_int <- data.frame(
#  player = player_lookup_defmid$player,
#  pred_int_total_mean  = apply(mu_player_draws_int, 2, mean),
# pred_int_total_lower = apply(mu_player_draws_int, 2, quantile, probs = 0.025),
#  pred_int_total_upper = apply(mu_player_draws_int, 2, quantile, probs = 0.975)
# 

# summarise posterior per player using 95% HPD intervals
hpd_pred_int_poi <- t(sapply(1:ncol(mu_player_draws_int), function(j) {
  hpd <- HPDinterval(as.mcmc(mu_player_draws_int[, j]), prob = 0.95)
  c(lower = hpd[1, "lower"], upper = hpd[1, "upper"])
}))

player_pred_total_int_poi <- data.frame(
  player = player_lookup_defmid$player,
  pred_int_total_mean  = apply(mu_player_draws_int, 2, mean),
  pred_int_total_lower = hpd_pred_int_poi[, "lower"],
  pred_int_total_upper = hpd_pred_int_poi[, "upper"]
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
player_pred_total_int_poi <- player_pred_total_int_poi %>%
  mutate(defmid_player_id = 1:N_players_int) %>%
  left_join(obs_totals_int, by = "defmid_player_id") %>%
  arrange(desc(pred_int_total_mean))

# compute overall model MAE and RMSE (single values for model)
errors <- player_pred_total_int_poi$obs_int_total - 
  player_pred_total_int_poi$pred_int_total_mean

MAE_model_int_poi  <- mean(abs(errors))
RMSE_model_int_poi <- sqrt(mean(errors^2))

MAE_model_int_poi
RMSE_model_int_poi

View(player_pred_total_int_poi)

# INTERCEPTIONS - NEGATIVE BINOMIAL

library(dplyr)
library(rjags)
library(coda)



# Centered NB model 
model_nb_int <- "
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
    
    y_int[i] ~ dnegbin(p[i], r)
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


# 4) Fit NB model
jags_model_int_nb <- jags.model(
  textConnection(model_nb_int),
  data = data_jags_int,
  n.chains = 3,
  n.adapt = 1000 
)

update(jags_model_int_nb, 9000) 

params_int_nb <- c(
  "Delta", "lambda_e", "lambda_Edive",
  "tau.lambda_e", "tau.lambda_Edive","r" # "mu.lambda_e", "mu.lambda_Edive","m", "s"
)

samples_int_nb <- coda.samples(
  jags_model_int_nb,
  variable.names = params_int_nb,
  n.iter = 40000, 
  thin = 1
)

# player ability: Neg Bin interceptions
summary_stats_int_nb  <- summary(samples_int_nb)$statistics
summary_quants_int_nb <- summary(samples_int_nb)$quantiles

Delta_idx_nb   <- grep("^Delta\\[", rownames(summary_stats_int_nb))
Delta_stats_nb <- summary_stats_int_nb[Delta_idx_nb, ]

player_lookup_defmid <- df_defmid %>%
  distinct(defmid_player_id, player) %>%
  arrange(defmid_player_id)


#player_table_nb_int <- data.frame(
#  player        = player_lookup_defmid$player,
#  ability_mean  = Delta_stats_nb[, "Mean"],
#  ability_sd    = Delta_stats_nb[, "SD"],
#  ability_lower = summary_quants_int_nb[Delta_idx_nb, "2.5%"],
#  ability_upper = summary_quants_int_nb[Delta_idx_nb, "97.5%"]
# )
samples_int_nb_mcmc <- as.mcmc(as.matrix(samples_int_nb))

Delta_hpd_nb_int <- HPDinterval(samples_int_nb_mcmc, prob = 0.95)
Delta_hpd_nb_int <- Delta_hpd_nb_int[grep("^Delta\\[", rownames(Delta_hpd_nb_int)), ]

player_table_nb_int <- data.frame(
  player        = player_lookup_defmid$player,
  ability_mean  = Delta_stats_nb[, "Mean"],
  ability_sd    = Delta_stats_nb[, "SD"],
  ability_lower = Delta_hpd_nb_int[, "lower"],
  ability_upper = Delta_hpd_nb_int[, "upper"]
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
AmrabatNBInt <- samples_int_nb[, "Delta[213]"]
heidel.diag(AmrabatNBInt)
plot(AmrabatNBInt, main = expression("Plot for " * Delta[213] * " (Sofyan Amrabat)"))
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
team_draws_int_nb  <- S_int_nb[, grep("^lambda_e\\[", colnames(S_int_nb)), drop = FALSE]
oppteam_draws_int_nb   <- S_int_nb[, grep("^lambda_Edive\\[", colnames(S_int_nb)), drop = FALSE]
r_draws_int_nb     <- S_int_nb[, "r"]  

# player and team ids + minutes for defenders
pid  <- df_defmid$defmid_player_id
tid  <- df_defmid$team_id
tau      <- df_defmid$tau
tau_star <- df_defmid$tau_star


# matrix mult, posterior draws x opponent design -> opponent contribution per obs
opp_sum_draws_int_nb <- oppteam_draws_int_nb %*% t(delta_defmid)

# build eta for each posterior draw and observation
# we add attacking opponent strength because the more attacks the more chance
# of interceptions by defenders
log_eta_draws <- Delta_draws_int_nb[, pid, drop = FALSE] +
  sweep(team_draws_int_nb[, tid, drop = FALSE], 2, tau, `*`) -
  sweep(opp_sum_draws_int_nb, 2, tau_star, `*`)

eta_draws_int_nb <- exp(log_eta_draws)


# expected interceptions per observation (NB mean is mu = exp(eta)*exposure)
mu_draws_int_nb <- sweep(eta_draws_int_nb, 2, tau, `*`) 

# total expected interceptions over tournament per defender
N_players_int <- data_jags_int$N_players

# posterior distribution of tournament totals for each defender (sum across their obs)
mu_player_draws_int_nb <- sapply(1:N_players_int, function(p) {
  cols <- which(pid == p)
  rowSums(mu_draws_int_nb[, cols, drop = FALSE])
})

# summarise posterior per defender
# player_pred_total_int_nb <- data.frame(
#  player = player_lookup_defmid$player,
#  pred_int_total_mean  = apply(mu_player_draws_int_nb, 2, mean),
#  pred_int_total_lower = apply(mu_player_draws_int_nb, 2, quantile, probs = 0.025),
#  pred_int_total_upper = apply(mu_player_draws_int_nb, 2, quantile, probs = 0.975)
# )

# summarise posterior per defender using 95% HPD intervals
hpd_pred_int_nb <- t(sapply(1:ncol(mu_player_draws_int_nb), function(j) {
  hpd <- HPDinterval(as.mcmc(mu_player_draws_int_nb[, j]), prob = 0.95)
  c(lower = hpd[1, "lower"], upper = hpd[1, "upper"])
}))

player_pred_total_int_nb <- data.frame(
  player = player_lookup_defmid$player,
  pred_int_total_mean  = apply(mu_player_draws_int_nb, 2, mean),
  pred_int_total_lower = hpd_pred_int_nb[, "lower"],
  pred_int_total_upper = hpd_pred_int_nb[, "upper"]
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
  mutate(defmid_player_id = 1:N_players_int) %>%
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