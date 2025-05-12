
# ==============================================
# Simulating and Modeling multinomial mixed model
# Set up is investors invest in 4 'deal types' D1, D2, D3, and D4
#
# - Simulate investor preferences using multinomial distributions
# - Estimate covariance empirically
# - Estimate covariance structure via Bayesian multinomial regression
# - Posterior predictive checks for model validation
# ==============================================



library(tidyverse)
library(brms)

set.seed(34)

# Parameters for fake data
num_investors <- 200
dealtypes <- c("D1", "D2", "D3", "D4")
num_dealtypes <- length(dealtypes)

# Covariance of relative log odds (relative to D1)
sigma = matrix(c(.5, 0.1, 0.0,
                 0.1, .1, 0.1,
                 0.0, 0.1, .4), nrow=3)

# Generate 'true' investor preferences in logodds space
investor_prefs <- MASS::mvrnorm(
  n = num_investors,
  mu = rep(0, num_dealtypes - 1),  # baseline first category is 0
  Sigma = sigma
)

# add zero baseline for the first category explicitly (Pivot)
investor_prefs <- cbind(0, investor_prefs)

# Calculate probabilities from logodds ratios (Softmax)
# i.e.  -> exp(-33.0)/(1 + exp(.77) + exp(-33.0) + exp(4))
investor_prefs_prob <- exp(investor_prefs)
investor_prefs_prob <- investor_prefs_prob / rowSums(investor_prefs_prob)


# simulate observed deal counts per investor
total_deals <- sample(20:50, num_investors, replace=TRUE)

# Now do draws from multinomial using the probabilities
observed_counts <- t(sapply(1:num_investors, function(i) {
  rmultinom(1, total_deals[i], investor_prefs_prob[i,])
}))

colnames(observed_counts) <- dealtypes

investor_data <- data.frame(
  investor = factor(1:num_investors),
  total_deals = total_deals,
  observed_counts
)

###############################################################
#
# Analyzing the fake data
#
##############################################################

## First, lets try a simple calculation of covariance of relative odds

# The data probably contains zeros, so we will add a pseudo count to each observation
# This is loosely equivalent to a Dirichlet prior with alpha_i = pseudo_count

pseudo_count = 1 # 0.5 is also reasonable. Bigger counts biases toward less preference
inv_log_odds = investor_data |>
           mutate(
             log_odds_D2 = log((D2 + pseudo_count) / (D1 + pseudo_count)),
             log_odds_D3 = log((D3 + pseudo_count) / (D1 + pseudo_count)),
             log_odds_D4 = log((D4 + pseudo_count) / (D1 + pseudo_count)),
             weight = total_deals
           ) |>
            select(
              log_odds_D2, log_odds_D3, log_odds_D4, weight
            )




weighted_cov <- cov.wt(inv_log_odds[,1:3],wt = inv_log_odds$weight)
print(weighted_cov$cov)

# These are a bit too big
# The reason is it is measuring both the true covariance
# AND extra multinomial sampling noise. This gets better with more data.
# Could also increase pseudo count
#
#             log_odds_D2 log_odds_D3 log_odds_D4
#log_odds_D2   0.7424503   0.2522264   0.1560024
#log_odds_D3   0.2522264   0.3121256   0.1820508
#log_odds_D4   0.1560024   0.1820508   0.4966619

################################################
#  BRMS Multinomial fit
################################

## Ok lets do it using logistic regression directly
# Advantages are
# - uncertainty quantification and generative statistical framework
# - Separate covariance variation from binomial variation
# - ability to incorporate other predictors
# - ability to incorporate priors


# Combine counts for brms

investor_data$counts = with(investor_data, cbind(D1,D2,D3,D4))



fit <- brm(
  # 'i' is a dummy to indicate we want correlations between dealtypes
  # Also note that the trials(total_deals) is required even though redundant
  counts | trials(total_deals)  ~  1 + (1  | i | investor),
  data = investor_data,
  family = multinomial(),
  iter = 2000, chains = 4, cores = 4
)



# generate covariance
cov_fit = VarCorr(fit)$investor$cov
mean_cov_matrix <- matrix(
  c(
    cov_fit[1, 1, "muD2_Intercept"],
    cov_fit[1, 1, "muD3_Intercept"],
    cov_fit[1, 1, "muD4_Intercept"],

    cov_fit[2, 1, "muD2_Intercept"],
    cov_fit[2, 1, "muD3_Intercept"],
    cov_fit[2, 1, "muD4_Intercept"],

    cov_fit[3, 1, "muD2_Intercept"],
    cov_fit[3, 1, "muD3_Intercept"],
    cov_fit[3, 1, "muD4_Intercept"]
  ),
  nrow = 3, byrow = TRUE,
  dimnames = list(
    c("muD2", "muD3", "muD4"),
    c("muD2", "muD3", "muD4")
  )
)
print(mean_cov_matrix)

# Reproduces the gen model pretty well
#
#        muD2       muD3       muD4
#muD2 0.56984367 0.14426463 0.03653532
#muD3 0.14426463 0.08785884 0.07581221
#muD4 0.03653532 0.07581221 0.34067843

# Summary includes errors
summary(fit)




################################################
#
#  some posterior checks, under development
########################################
library(bayesplot)
# Generate posterior predictive samples explicitly
pp_samples <- posterior_predict(fit, draws = 500)  # shape: draws x investors x categories

# Convert posterior predictive counts to proportions first (by draw)
predicted_proportions_draws <- sweep(pp_samples, c(1,2), investor_data$total_deals, "/")

# Then average these proportions across draws
predicted_proportions_mean <- apply(predicted_proportions_draws, c(2,3), mean)
predicted_proportions_lower <- apply(predicted_proportions_draws, c(2,3), quantile, 0.025)
predicted_proportions_upper <- apply(predicted_proportions_draws, c(2,3), quantile, 0.975)

# Observed proportions
observed_proportions <- investor_data$counts / investor_data$total_deals

# Plot observed vs predicted proportions
par(mfrow = c(2,2))
dealtypes <- c("D1", "D2", "D3", "D4")
for (i in 1:4) {
  plot(
    observed_proportions[, i],
    predicted_proportions_mean[, i],
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = paste("Observed proportion", dealtypes[i]),
    ylab = paste("Predicted proportion", dealtypes[i]),
    main = paste("Observed vs Predicted:", dealtypes[i]),
    pch = 16, col = rgb(0,0,1,0.5)
  )
  arrows(
    x0 = observed_proportions[, i],
    y0 = predicted_proportions_lower[, i],
    y1 = predicted_proportions_upper[, i],
    angle = 90, length = 0.05, code = 3, col = rgb(0,0,0,0.3)
  )
  abline(0,1,col='red', lwd=2)
}

# check zeros separately
observed_zeros <- colSums(investor_data$counts == 0)

predicted_zeros_per_draw <- apply(pp_samples == 0, c(1,3), sum) # draws x categories



# TODO make this a quad

for (dealtype_idx in 1:4) {
  p <- bayesplot::ppc_stat(
    y = observed_zeros[dealtype_idx],
    yrep = matrix(predicted_zeros_per_draw[, dealtype_idx], ncol = 1),
    stat = "identity",
    binwidth = 1
  ) +
    ggtitle(paste("Zero-count Investors for", dealtypes[dealtype_idx])) +
    theme_minimal()
  print(p)

}




# Choose a few random investors for detailed checks
set.seed(123)
investor_ids <- sample(1:num_investors, 4)

par(mfrow=c(2,2))
for (inv in investor_ids) {
  observed <- investor_data$counts[inv,]
  predicted <- pp_samples[, inv, ]  # draws x categories

  boxplot(
    predicted,
    names = dealtypes,
    main = paste("Investor", inv, "Counts"),
    ylab = "Deal Counts"
  )
  points(1:4, observed, col = "red", pch=19)
  legend("topright", legend = c("Observed", "Predicted"), col = c("red", "black"), pch=c(19, NA), lty=c(NA, 1))
}

######################################
#
# Alternative fit directly using stan
#########################################
if(FALSE){
  library(cmdstanr)
  investor_data$counts = with(investor_data, cbind(D1,D2,D3,D4))
  stan_data <- list(
    N = nrow(investor_data),
    K = length(dealtypes),
    counts = as.matrix(investor_data$counts)
  )
  # compare to stan version of logistic normal
  model_lm = cmdstan_model("logisticNorm.stan")
  fit <- model_lm$sample(data = stan_data, chains =4, cores=4, iter_sampling = 1000, iter_warmup =500)

  draws_df <- fit$draws(format = "draws_df")
  mcmc_trace(draws_df, pars = vars(starts_with("sigma")))
  mcmc_trace(draws_df, pars = vars(starts_with("Omega")))
  # Extract means for Omega and sigma
  Omega_means <- draws_df |>
    select(starts_with("Omega")) |>
    summarise(across(everything(), mean)) |>
    unlist()

  sigma_means <- draws_df |>
    select(starts_with("sigma")) |>
    summarise(across(everything(), mean)) |>
    unlist()

  K_minus_1 <- length(sigma_means)
  Omega_mean_matrix <- matrix(Omega_means, nrow = K_minus_1, byrow = FALSE)
  cov_matrix <- diag(sigma_means) %*% Omega_mean_matrix %*% diag(sigma_means)
  colnames(cov_matrix) <- rownames(cov_matrix) <- paste0("logodds_", 2:(K_minus_1 + 1))
  print(cov_matrix)  # Exactly the same as brms

}


################################
# Alternative- Dirichlet_multinomial
#################################
# Note that this is workable but actualy less flexible.
# All correlations must be negative and symmetric

if(FALSE){ #disable this section for now
library(cmdstanr)

# Prepare data
stan_data <- list(
  N = nrow(investor_data),
  K = length(dealtypes),
  counts = as.matrix(investor_data$counts)
)

# Compile the model
model <- cmdstan_model("dirichletMn.stan")

# Run the model
fit <- model$sample(data = stan_data, chains = 4, iter_sampling = 1000, iter_warmup = 500)

# Summary
fit$summary("alpha")
library(bayesplot)
draws <- fit$draws()
mcmc_trace(draws)

mcmc_dens(draws, regex_pars = "alpha")



#####
# MLE
#####
#install.packages("MGLM") # old ?
library(MGLM)
dm_fit <- MGLMfit(as.matrix(investor_data$counts), dist="DM")
dm_fit
}
