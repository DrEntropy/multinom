library(dplyr)
library(tidyr)
library(MASS) # For multivariate normal sampling
library(brms)

set.seed(42)

# Parameters
num_companies <- 40
num_categories <- 3 # Categories: 0 (pivot), 1, 2 . Minimally interesting


# 1. Generate company random effects (pivot category 0 logit is 0)
mean_effects <- c(0, 0) # categories 1,2,3; 0 is baseline
cov_matrix <- matrix(c(
  1.0, 0.5,
  0.5, 2.0
), nrow = 2, byrow = TRUE)

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




  # Generate observations
  obs <- lapply(1:n_obs, function(j) {
    logits <- company_effects[i, ]
    probs <- exp(logits) / sum(exp(logits))
    type <- sample(0:(num_categories - 1), size = 1, prob = probs)

    data.frame(
      company = company_df$company[i],

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
  formula = type ~0 +  (1 | i | company),
  data = df,
  family = categorical(),
  chains = 4, cores = 4, iter = 2000
)

summary(fit)
