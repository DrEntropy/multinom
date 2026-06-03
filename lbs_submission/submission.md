# LBS Contest Submission — Bayesian Workflow

**Notebook:** [`aggregated_cv.ipynb`](aggregated_cv.ipynb)

**Problem:**

In multinomial choice modeling, individual-level data gets large fast, so I often need to 
aggregate it into per-group count vectors to make a hierarchical fit tractable. For example:
customer 17 chose categories (3,1,0,5) across 9 purchases. (This arose with proprietary
investor/deal-type data; the example is synthetic but structurally identical.) Handed this
data, you might just fit independent random effects per category — modeling the correlations between
products is less obvious, and isn't even offered in Bambi. Do they earn their keep? On
synthetic data a targeted posterior predictive check separates the two models cleanly,
but cross-validating that comparison is where I get stuck — it seems to depend on what I'm
predicting: an existing customer's next purchase, or a new customer's whole basket.
az.loo answers the second, but isn't trustworthy here (bad Pareto-k). After
aggregating for computation, what's the best workflow to compare these kinds of models
reliably, and at scale?
