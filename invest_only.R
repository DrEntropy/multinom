library(tidyverse)
library(brms)

set.seed(42)

# Parameters for fake data
num_investors <- 200
dealtypes <- c("D1", "D2", "D3", "D4")
num_dealtypes <- length(dealtypes)

# Generate true investor preferences
investor_prefs <- MASS::mvrnorm(
  n = num_investors,
  mu = rep(0, num_dealtypes - 1),  # baseline first category is 0
  Sigma = matrix(c(1.0, -0.3, 0.2,
                   -0.3, 1.5, -0.4,
                   0.2, -0.4, 2.0), nrow=3)
)

# add zero baseline for the first category explicitly
investor_prefs <- cbind(0, investor_prefs)

investor_prefs_prob <- apply(investor_prefs, 1, brms::inv_logit_scaled) %>% t()
investor_prefs_prob <- investor_prefs_prob / rowSums(investor_prefs_prob)

# simulate observed deal counts per investor
total_deals <- sample(20:50, num_investors, replace=TRUE)

observed_counts <- t(sapply(1:num_investors, function(i) {
  rmultinom(1, total_deals[i], investor_prefs_prob[i,])
}))

colnames(observed_counts) <- dealtypes

investor_data <- data.frame(
  investor = factor(1:num_investors),
  total_deals = total_deals,
  observed_counts
)
investor_data$counts = with(investor_data, cbind(D1,D2,D3,D4))


fit <- brm(
  counts | trials(total_deals)  ~ 0 + (1 | i | investor),
  data = investor_data,
  family = multinomial(),
  iter = 2000, chains = 4, cores = 4
)

summary(fit)
library(bayesplot)
# Generate posterior predictive samples explicitly
pp_samples <- posterior_predict(fit, draws = 500)  # shape: draws x investors x categories

# Convert posterior predictive counts to proportions first (by draw)
predicted_proportions_draws <- sweep(pp_samples, c(1,2), investor_data$total_deals, "/")

# Then average these proportions across draws
predicted_proportions_mean <- apply(predicted_proportions_draws, c(2,3), mean)

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
  abline(0,1,col='red', lwd=2)
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




