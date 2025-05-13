library(dplyr)
library(tidyr)
library(brms)
library(bayesplot)
library(posterior)

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
group_effects <- MASS::mvrnorm(n = num_groups, mu = mean_effects, Sigma = cov_matrix)
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

#  Fit categorical model using brms
fit <- brm(
  formula = type ~   (1 | i | group),
  data = df,
  family = categorical(),
  chains = 4, cores = 4, iter = 2000
)

summary(fit)
print(VarCorr(fit)$group$cov)

##
# Alternative, convert to counts

agg_df <- df |>
  group_by(group,type) |> summarize(count = n()) |>
  pivot_wider(id_cols = group, names_from = type, values_from = count, values_fill = 0) |>
  mutate( counts = cbind(type_0,type_1,type_2)) |> ungroup() |>
  mutate( total = type_0 + type_1 + type_2)

fit2 <- brm(
  formula = counts | trials(total) ~   (1 | i | group),
  data =agg_df,
  family = multinomial(),
  chains = 4, cores = 4, iter = 2000
)

summary(fit2)
print(VarCorr(fit2)$group$cov)

## Use stan directly
library(cmdstanr)

stan_data <- list(
  N = nrow(agg_df),
  K = 3,
  counts = as.matrix(agg_df$counts)
)
model_lm = cmdstan_model("logisticNorm.stan")
fit <- model_lm$sample(data = stan_data, chains =4, cores=4)

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
