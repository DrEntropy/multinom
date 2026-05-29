"""Smoke test: PyMC 6 sampling with nutpie (default when installed)."""

import arviz as az
import numpy as np
import pymc as pm

print(f"PyMC {pm.__version__}")
# Pytensor cxx not used anymore by default. New version uses numba and nutpie.
#import pytensor
#print(f"Pytensor config.cxx: {pytensor.config.cxx}")

np.random.seed(123)
data = np.random.normal(loc=0, scale=1, size=100)

with pm.Model():
    mu = pm.Normal("mu", mu=0, sigma=1)
    sigma = pm.HalfNormal("sigma", sigma=1)
    pm.Normal("obs", mu=mu, sigma=sigma, observed=data)
    trace = pm.sample(500, tune=500, cores=1, chains=1)

print(az.summary(trace))
