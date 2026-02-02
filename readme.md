# Multinomial Regression with PyMC

This repository is a small project designed to explore multinomial regression using [PyMC](https://www.pymc.io/), a probabilistic programming library for Bayesian modeling.

## GOALS:

### Data

Focus on multinomial regression.  Simplest cases are probably better,  it needs to involved groups however.    One categorical predictor but treating them as groups is probably a good idea (see https://doingbayesiandataanalysis.blogspot.com/2026/01/problem-in-multilevel-hierarchical.html for example using fake data)

Data options:
* fake data
* Ice cream flavor, movie genre, IRIS? 
* See also Advanced Regression Course.



### ZSN
Discuss zerosum normal - in particular that it should also sum to zero on outcome which avoids issues with choosing a pivot/refecence category when there is no suitable choice.  Note further that with random effects, if you have a reference, there are no random effects for that category, breaking symmetry. 

### Cross validation / model comparison
Discuss how to do cross validation . The right way to do this depends on the use case: If you actually want to assess prediction quality for the model for a new observation, that observation is going to be a single one. So is there an efficient way to do this? Sure but i dont think this is automated in any of the packages. For model comparison,  if you leave out at the group level (i.e. unique combination of predictors) then it is harder problem and maybe the better way? This is what brms would do, for example. But see thsi: [https://discourse.mc-stan.org/t/understanding-loo-and-binomial-models/23500](https://discourse.mc-stan.org/t/understanding-loo-and-binomial-models/23500)    

In any event to do the proper individual observation level k-fold requires manually doing it i believe.  That is,  'disagregate' the data so that you have one row per observation.  Do the k-fold split and then re-aggregate each of the splits, and fit as usual.  I think this might be possible without the actual intermediary steps (at least for binomial case) with something called binomial thinning / splitting per cell?   [Data thinning for convolution-closed distributions](https://arxiv.org/abs/2301.07276)


But Clarify: This doesnt apply to all cases, often teh observations are the whole , for example how many times a particular customer purchases for N visits (not 'next visit').  This leaving out whole rows is what BRMS does when using the kfold function.  

 

### using varying intercepts to emulate overdispersion

This is sort of extra note.. not sure how it fits in if at all. 

If you have one row per ‘group’ of counts , you can use random intecepts (one per row) to emulate over dispersion.  Except when doing posterior predictions dont use the fit group level intercepts, draw new. (this simulated the idea that the group level interecepts are ‘latent’ and not known in such a case.)

### Misc refs
Other examples: [Nominal Regression in Stan](https://quantscience.rbind.io/posts/2022-07-30-multinomial-regression-in-stan/)


## TODO:

[] review claude coded article, or start over. This does capture a lot of my thoughts though
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