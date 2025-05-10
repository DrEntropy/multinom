library(dplyr)
library(tidyr)
library(MASS) # For multivariate normal sampling
library(brms)

set.seed(42)

# Parameters
num_companies <- 40
num_categories <- 4 # Categories: 0 (pivot), 1, 2, 3
num_prev_categories <- 4

# 1. Generate company random effects (pivot category 0 logit is 0)
mean_effects <- c(0, 0, 0) # categories 1,2,3; 0 is baseline
cov_matrix <- matrix(c(
  1.0, 0.5, 0.3,
  0.5, 1.0, 0.4,
  0.3, 0.4, 1.0
), nrow = 3, byrow = TRUE)

# Random effects for companies
company_effects <- mvrnorm(n = num_companies, mu = mean_effects, Sigma = cov_matrix)
company_effects <- cbind(0, company_effects) # Add pivot logits = 0
colnames(company_effects) <- paste0("logit_", 0:(num_categories - 1))

company_df <- data.frame(
  company = paste0("company_", 1:num_companies),
  company_effects
)

# 2. Generate observations for each company
data_list <- lapply(1:num_companies, function(i) {
  n_obs <- sample(1:10, 1)

  # Random prev_type for each observation
  prev_type <- sample(0:(num_prev_categories - 1), n_obs, replace = TRUE)

  # Fixed effects of prev_type (pivot=0)
  prev_type_effects <- matrix(c(
    0,    0,    0,    # prev_type = 0 (reference)
    0.5, -0.5,  0.2,  # prev_type = 1
    -0.3,  0.6, -0.4,  # prev_type = 2
    0.2, -0.2,  0.5   # prev_type = 3
  ), byrow = TRUE, nrow = 4)
  prev_type_effects <- cbind(0, prev_type_effects)

  # Generate observations
  obs <- lapply(1:n_obs, function(j) {
    logits <- company_effects[i, ] + prev_type_effects[prev_type[j] + 1, ]
    probs <- exp(logits) / sum(exp(logits))
    type <- sample(0:(num_categories - 1), size = 1, prob = probs)

    data.frame(
      company = company_df$company[i],
      prev_type = paste0("type_", prev_type[j]),
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
  formula = type ~ prev_type + (1 | type | company),
  data = df,
  family = categorical(),
  chains = 2, cores = 2, iter = 2000
)

summary(fit)
