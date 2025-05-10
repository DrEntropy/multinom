library(brms)
library(dplyr)

set.seed(123) # for reproducibility

# Parameters
n <- 500                     # number of observations
n_groups <- 10               # number of groups for random effects
groups <- paste0("G", 1:n_groups)
predictor_levels <- c("A", "B")

# Simulate predictor
df <- data.frame(
  predictor = sample(predictor_levels, size = n, replace = TRUE),
  group = sample(groups, size = n, replace = TRUE)
)

# Simulate random effects (multivariate normal)
# 2 non-baseline categories => 2-dimensional random effects
mu_rand <- c(0, 0)
Sigma_rand <- matrix(c(1.0, 0.5,
                       0.5, 1.5), nrow = 2)

rand_effects <- MASS::mvrnorm(n = n_groups, mu = mu_rand, Sigma = Sigma_rand)
rownames(rand_effects) <- groups
colnames(rand_effects) <- c("cat2", "cat3")

# Fixed effect (predictor)
beta <- matrix(c(0.5, -0.5,   # predictor = B vs. A for cat2, cat3
                 -1.0, 1.0),  # intercepts for cat2, cat3
               ncol = 2, byrow = TRUE)
rownames(beta) <- c("predictorB", "Intercept")
colnames(beta) <- c("cat2", "cat3")

# Generate linear predictor (logits)
df <- df %>%
  mutate(
    intercept_cat2 = beta["Intercept", "cat2"] + rand_effects[group, "cat2"],
    intercept_cat3 = beta["Intercept", "cat3"] + rand_effects[group, "cat3"],
    predictor_effect_cat2 = ifelse(predictor == "B", beta["predictorB", "cat2"], 0),
    predictor_effect_cat3 = ifelse(predictor == "B", beta["predictorB", "cat3"], 0),
    logit_cat2 = intercept_cat2 + predictor_effect_cat2,
    logit_cat3 = intercept_cat3 + predictor_effect_cat3
  )

# Baseline logits (cat1) = 0
logits <- with(df, cbind(cat1 = 0, cat2 = logit_cat2, cat3 = logit_cat3))

# Convert logits to probabilities and simulate categorical response
prob <- exp(logits) / rowSums(exp(logits))
df$response <- apply(prob, 1, function(p) sample(c("cat1", "cat2", "cat3"), 1, prob = p))

# Convert to factors explicitly
df$response <- factor(df$response, levels = c("cat1", "cat2", "cat3"))
df$predictor <- factor(df$predictor)
df$group <- factor(df$group)

# Fit the BRMS multinomial model with random intercepts
fit <- brm(
  response ~ predictor + (1 | response | group),
  data = df,
  family = categorical()
)

summary(fit)
