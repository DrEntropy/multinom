data {
  int<lower=1> N;             // number of investors
  int<lower=2> K;             // number of categories
  array[N,K] int counts;           // observed counts
}

parameters {
  vector<lower=0>[K] alpha;   // Dirichlet parameters
}

model {
  alpha ~ gamma(2, 0.1);  // weakly informative prior
  for (n in 1:N)
    counts[n] ~ dirichlet_multinomial(alpha);
}
