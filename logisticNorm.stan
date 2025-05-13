data {
  int<lower=1> N;              // number of sources / rows
  int<lower=2> K;              // number of categories
  array[N, K] int counts;      // observed counts per investor/category
}

parameters {
  vector[K - 1] mu;                        // mean log-odds
  cholesky_factor_corr[K - 1] L_Omega;     // correlation structure
  vector<lower=0>[K - 1] sigma;            // scale (std dev) for log-odds
  array[N] vector[K - 1] z;                // standardized latent variables
}



model {
  // Priors
  mu ~ normal(0, 3);
  sigma ~ student_t(3, 0, 2.5);
  L_Omega ~ lkj_corr_cholesky(1);

  // Non centered parameterization
  for (n in 1:N) {
    z[n] ~ normal(0, 1);
  }

  // Likelihood
  for (n in 1:N) {
    vector[K - 1] log_odds = mu + diag_pre_multiply(sigma, L_Omega) * z[n];
    vector[K] log_p = append_row(0, log_odds); // reference (first) category at baseline
    counts[n] ~ multinomial_logit(log_p);
  }
}

//   retrieve correlation matrix explicitly

generated quantities {
  corr_matrix[K - 1] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
