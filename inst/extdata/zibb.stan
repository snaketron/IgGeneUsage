functions {
  // skeleton from brms -> binomial_lpmf changed to beta_binomial_lpmf
  real zibb_lpmf(int y, int trials, real a, real b, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         beta_binomial_lpmf(0 | trials, a, b));
    } else {
      return bernoulli_lpmf(0 | zi) +
             beta_binomial_lpmf(y | trials, a, b);
    }
  }

  // TODO: test rng extensively
  // Affects predictions at repertoire (sample) level if zi >> 0
  int zibb_rng(int y, int trials, real a, real b, real zi) {
    if (bernoulli_rng(zi) == 1) {
      return (0);
    } else {
      return (beta_binomial_rng(trials, a, b));
    }
  }
}


data {
  int <lower = 0> N_sample; // number of samples (repertoires)
  int <lower = 0> N_gene; // number of genes
  int Y [N_gene, N_sample]; // number of successes (cells) in samples x gene
  int N [N_sample]; // number of total tries (repertoire size)
  int <lower = -1, upper = 1> X[N_sample]; // condition
}


transformed data {
  // convert int N -> real N fo convenient division
  // in generated quantities block
  real Nreal [N_sample];
  Nreal = N;
}


parameters {
  real alpha_grand;
  real beta_grand;
  real <lower = 0> alpha_sigma;
  real <lower = 0> beta_sigma;
  real <lower = 0> beta_gene_sigma;
  vector [N_gene] alpha_raw;
  vector [N_gene] beta_raw [N_sample];
  vector [N_gene] beta_gene_raw;
  real <lower = 0> phi;
  real<lower = 0> tau;
  real <lower = 0, upper = 1> z;  // zero-inflation probability
}



transformed parameters {
  vector [N_gene] beta [N_sample];
  vector [N_gene] alpha_gene;
  vector [N_gene] beta_gene;
  vector <lower = 0> [N_gene] a [N_sample];
  vector <lower = 0> [N_gene] b [N_sample];

  // non-centered params (at repertoire level)
  alpha_gene = alpha_grand + alpha_sigma * alpha_raw;
  beta_gene = beta_grand + beta_gene_sigma * beta_gene_raw;

  // non-centered params (at gene level)
  for(i in 1:N_sample) {
    beta[i] = beta_gene + beta_sigma * beta_raw[i];
    a[i] = inv_logit(alpha_gene + beta[i]*X[i]) * phi;
    b[i] = phi - a[i];
  }
}


model {
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      // likelihood: TODO speedup
      target += zibb_lpmf(Y[j, i] | N[i], a[i][j], b[i][j], z);
    }
  }


  // priors
  alpha_grand ~ normal(0, 10);
  beta_grand ~ normal(0, 5);

  alpha_sigma ~ cauchy(0, 1);
  beta_sigma ~ cauchy(0, 1);
  beta_gene_sigma ~ cauchy(0, 1);

  alpha_raw ~ normal(0, 1);
  for(i in 1:N_sample) {
    beta_raw[i] ~ normal(0, 1);
  }
  beta_gene_raw ~ normal(0, 1);

  phi ~ exponential(tau); //pareto ||
  tau ~ gamma(3, 0.1);
  z ~ beta(1, 3);
}


generated quantities {
  // PPC: count usage
  int Yhat [N_gene, N_sample];

  // PPC: proportion usage
  real Yhat_individual [N_gene, N_sample];

  // PPC: proportion usage at a gene level
  matrix [2, N_gene] Yhat_gene;

  // LOG-LIK
  matrix [N_gene, N_sample] log_lik;

  //TODO: speedup, run in C++ not big factor on performance
  for(j in 1:N_gene) {
    for(i in 1:N_sample) {
      log_lik[j,i] = zibb_lpmf(Y[j, i] | N[i], a[i][j], b[i][j], z);
      Yhat[j, i] = zibb_rng(Y[j, i], N[i], a[i][j], b[i][j], z);

      if(Nreal[i] == 0.0) {
        Yhat_individual[j, i] = 0;
      }
      else {
        Yhat_individual[j, i] = Yhat[j,i]/Nreal[i];
      }
    }
    Yhat_gene[1, j] = inv_logit(alpha_gene[j]+beta_gene[j]*1.0);
    Yhat_gene[2, j] = inv_logit(alpha_gene[j]+beta_gene[j]*(-1.0));
  }
}
