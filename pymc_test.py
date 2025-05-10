import pymc as pm
import arviz as az
import numpy as np

# Generate synthetic data
np.random.seed(123)
data = np.random.normal(loc=0, scale=1, size=100)
import pytensor 
#pytensor.config.cxx=“”
# Simple PyMC model
with pm.Model():
    mu = pm.Normal('mu', mu=0, sigma=1)
    sigma = pm.HalfNormal('sigma', sigma=1)
    obs = pm.Normal('obs', mu=mu, sigma=sigma, observed=data)

    # Sample from the posterior
    trace = pm.sample(500, tune=500, cores = 1)

# Summarize and print results
summary = az.summary(trace)
print(summary)