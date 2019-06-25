// generated with brms 2.9.0
functions {

  /* zero-inflated binomial log-PDF of a single response
   * Args:
   *   y: the response value
   *   trials: number of trials of the binomial part
   *   theta: probability parameter of the binomial part
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_binomial_lpmf(int y, int trials,
                                   real theta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         binomial_lpmf(0 | trials, theta));
    } else {
      return bernoulli_lpmf(0 | zi) +
             binomial_lpmf(y | trials, theta);
    }
  }
  /* zero-inflated binomial log-PDF of a single response
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   trials: number of trials of the binomial part
   *   theta: probability parameter of the binomial part
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_binomial_logit_lpmf(int y, int trials,
                                         real theta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         binomial_lpmf(0 | trials, theta));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             binomial_lpmf(y | trials, theta);
    }
  }
  /* zero-inflated binomial log-PDF of a single response
   * logit parameterization of the binomial part
   * Args:
   *   y: the response value
   *   trials: number of trials of the binomial part
   *   eta: linear predictor for binomial part
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_binomial_blogit_lpmf(int y, int trials,
                                          real eta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         binomial_logit_lpmf(0 | trials, eta));
    } else {
      return bernoulli_lpmf(0 | zi) +
             binomial_logit_lpmf(y | trials, eta);
    }
  }
  /* zero-inflated binomial log-PDF of a single response
   * logit parameterization of the binomial part
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   trials: number of trials of the binomial part
   *   eta: linear predictor for binomial part
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_binomial_blogit_logit_lpmf(int y, int trials,
                                                real eta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         binomial_logit_lpmf(0 | trials, eta));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             binomial_logit_lpmf(y | trials, eta);
    }
  }
  // zero-inflated binomial log-CCDF and log-CDF functions
  real zero_inflated_binomial_lccdf(int y, int trials, real theta, real zi) {
    return bernoulli_lpmf(0 | zi) + binomial_lccdf(y | trials, theta);
  }
  real zero_inflated_binomial_lcdf(int y, int trials, real theta, real zi) {
    return log1m_exp(zero_inflated_binomial_lccdf(y | trials, theta, zi));
  }
}


data {
  int <lower = 0> N_sample; // number of samples
  int <lower = 0> N_gene; // number of genes
  int Y [N_gene, N_sample]; // number of cells in patient x gene
  int N [N_sample]; // number of total cells
  int <lower = -1, upper = 1> X[N_sample]; // design variable = condition
}

parameters {
  real alpha_grand;
  real beta_grand;
  real <lower = 0> alpha_sigma;
  real <lower = 0> beta_sigma;
  real <lower = 0> alpha_gene_sigma;
  real <lower = 0> beta_gene_sigma;
  vector [N_gene] alpha_raw [N_sample];
  vector [N_gene] beta_raw [N_sample];
  vector [N_gene] alpha_gene_raw;
  vector [N_gene] beta_gene_raw;
}



transformed parameters {
  vector [N_gene] alpha [N_sample];
  vector [N_gene] alpha_gene;
  vector [N_gene] beta [N_sample];
  vector [N_gene] beta_gene;

  // non-centered params (top)
  alpha_gene = alpha_grand + alpha_gene_sigma * alpha_gene_raw;
  beta_gene = beta_grand + beta_gene_sigma * beta_gene_raw;

  // non-centered params (low)
  for(i in 1:N_sample) {
    alpha[i] = alpha_gene + alpha_sigma * alpha_raw[i];
    beta[i] = beta_gene + beta_sigma * beta_raw[i];
  }
}


model {
  for(i in 1:N_sample) {
    Y[, i] ~ binomial_logit(N[i], alpha[i] + beta[i]*X[i]);
    // target += zero_inflated_binomial_blogit_lpmf(Y[n] | trials[n], mu[n], zi);
  }

  alpha_grand ~ normal(0, 20);
  beta_grand ~ normal(0, 5);

  alpha_sigma ~ cauchy(0, 3);
  alpha_gene_sigma ~ cauchy(0, 3);
  beta_sigma ~ cauchy(0, 3);
  beta_gene_sigma ~ cauchy(0, 3);

  for(i in 1:N_sample) {
    alpha_raw[i] ~ normal(0, 1);
    beta_raw[i] ~ normal(0, 1);
  }
  alpha_gene_raw ~ normal(0, 1);
  beta_gene_raw ~ normal(0, 1);
}


generated quantities {
  int Yhat [N_gene, N_sample];
  real Yhat_individual [N_gene, N_sample];
  matrix [2, N_gene] Yhat_gene;
  matrix [N_gene, N_sample] log_lik;
  real temp;

  for(j in 1:N_gene) {
    for(i in 1:N_sample) {
      temp = N[i];

      log_lik[j,i] = binomial_logit_lpmf(Y[j,i] | N[i], alpha[i][j] + beta[i][j] * X[i]);
      Yhat[j, i] = binomial_rng(N[i], inv_logit(alpha[i][j]+beta[i][j]*X[i]));

      if(temp == 0.0) {
        Yhat_individual[j, i] = 0;
      }
      else {
        Yhat_individual[j, i] = Yhat[j,i]/temp*100.0;
      }
    }
    Yhat_gene[1, j] = inv_logit(alpha_gene[j]+beta_gene[j]*1.0)*100.0;
    Yhat_gene[2, j] = inv_logit(alpha_gene[j]+beta_gene[j]*(-1.0))*100.0;
  }
}
