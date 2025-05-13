library(dplyr)
library(tidyr)
library(MASS) # For multivariate normal sampling
library(brms)

set.seed(42)

# Parameters
num_groups <- 100
num_categories <- 3 # Categories: 0 (pivot), 1, 2 . Minimally interesting



mean_effects <- c(0, 0) # categories 1,2,3; 0 is baseline
cov_matrix <- matrix(c(
  1.0, 0.5,
  0.5, 2.0
), nrow = 2, byrow = TRUE)

# Random effects for group / subject
group_effects <- mvrnorm(n = num_groups, mu = mean_effects, Sigma = cov_matrix)
group_effects <- cbind(0, group_effects) # Add pivot logits = 0
colnames(group_effects) <- paste0("logit_", 0:(num_categories - 1))

# simulated_df contains the underlying logits
simulated_df <- data.frame(
  group= paste0("group_", 1:num_companies),
  group_effects
)

#  generate observations for each group
data_list <- lapply(1:num_groups, function(i) {
  n_obs <- rpois(lambda = 30, n=1)
  # Generate observations
  obs <- lapply(1:n_obs, function(j) {
    logits <- group_effects[i, ]
    probs <- exp(logits) / sum(exp(logits))
    type <- sample(0:(num_categories - 1), size = 1, prob = probs)

    data.frame(
      group = simulated_df$group[i],

      type = paste0("type_", type)
    )
  })

  bind_rows(obs)
})

df <- bind_rows(data_list)

# Check simulated data
head(df)

# 3. Fit categorical model using brms
fit <- brm(
  formula = type ~   (1 | i | group),
  data = df,
  family = categorical(),
  chains = 4, cores = 4, iter = 2000
)

summary(fit)
print(VarCorr(fit)$group$cov)

