# Multinomial Regression with PyMC

This repository is a small project designed to explore multinomial regression using [PyMC](https://www.pymc.io/), a probabilistic programming library for Bayesian modeling.

## GOALS:

### Data

Focus on multinomial regression.  Simplest cases are probably better,  it needs to involved groups however.    One categorical predictor but treating them as groups is probably a good idea (see https://doingbayesiandataanalysis.blogspot.com/2026/01/problem-in-multilevel-hierarchical.html for example using fake data)

Data options:
* fake data
* Ice cream flavor, movie genre, IRIS? 
* See also Advanced Regression Course.

### Cross validation / model comparison

**Core problem:** After aggregation, what is the right way to do cross-validation? Does hierarchical structure matter?

**Context:** The answer depends on the use case:

- If predicting for a **new single observation**, we need observation-level CV
- If comparing models with **grouped data** (unique combinations of predictors), leaving out at group level may be better (this is what brms does by default)

For example if the data is customer purchases among N-items,  if you want to predict next purchase, thats individual, if you want to predict a new customers buying patterns that would be at group level.

TO do observational level k-fold:

1.  Disaggregate -> split -> aggregate
    * Take the data, disaggreate it into individual observations.  Do k-split at that level, then aggregate and fit , predicting on the held out 1/K fold (disaggregated?).
    * Straightforward, but perhaps slow with large data sets.
2. Data thinking / splitting per cell.
    * I need to learn about this, but apparenetly it is possible to do the above without creating the intermediate individual level data.  
    * This applies to "Convolution-closed distributions" [Data thinning for convolution-closed distributions](https://arxiv.org/abs/2301.07276)

See also [Stan discourse discussion](https://discourse.mc-stan.org/t/understanding-loo-and-binomial-models/23500) 

See also Advanced Regression Course.


### ZSN
Discuss zerosum normal - in particular that it should also sum to zero on outcome which avoids issues with choosing a pivot/refecence category when there is no suitable choice.  Note further that with random effects, if you have a reference, there are no random effects for that category, breaking symmetry. 


### using varying intercepts to emulate overdispersion

This is sort of extra note.. not sure how it fits in if at all. 

If you have one row per ‘group’ of counts , you can use random intecepts (one per row) to emulate over dispersion.  Except when doing posterior predictions dont use the fit group level intercepts, draw new. (this simulated the idea that the group level interecepts are ‘latent’ and not known in such a case.)

### Misc refs
Other examples: [Nominal Regression in Stan](https://quantscience.rbind.io/posts/2022-07-30-multinomial-regression-in-stan/)


## TODO:

[] review claude-assisted article (multinomial_regression_article.ipynb), or start over. This does capture a lot of my thoughts though
[] Review data thinning paper
[] Find a good dataset!
[] Compare BRMS / span / kfold. 
[] move this as subdirectory in NotesonDatasci

## Overview

Multinomial regression is a statistical technique used to model outcomes with more than two categories. This repository demonstrates how to implement and analyze multinomial regression models using PyMC.

* multinom.ipynb: This looks at the Iris data set and uses a multinomial regression model to predict the species of iris based on the features. 

* random_effects_mn.ipynb: This explores random effects in multinomial regression.  In this case we look at the simplest possible model: Three categories and just one random effect. In the future may add covariats.

* multinomial_regression_article.ipyb  -  Claude created document, possible hot garbage. 


## Requirements

see yaml file. Note that os of May 10 2025, there are issues with pymc installs on windows, see 'working_env.yml' for a workaround.

## References

- [PyMC Documentation](https://www.pymc.io/projects/docs/en/stable/)
- [Multinomial Regression Overview](https://en.wikipedia.org/wiki/Multinomial_logistic_regression)
- [Stan manual](https://mc-stan.org/docs/stan-users-guide/regression.html#multi-logit.section)
- [Bambi example](https://bambinos.github.io/bambi/notebooks/categorical_regression.html)

Also this old pymc discourse post:
- [PYMC Discourse](https://discourse.pymc.io/t/multivariate-multinomial-logistic-regression/5242)

And the old version of Bayesian Analysis in Python covered this (Chapter 4)