functions {
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

  int zibb_rng(int y, int trials, real a, real b, real zi) {
    if (bernoulli_rng(zi) == 1) {
      return (0);
    } else {
      return (beta_binomial_rng(trials, a, b));
    }
  }
}


data {
  int <lower = 0> N_sample; // number of samples
  int <lower = 0> N_gene; // number of genes
  int Y [N_gene, N_sample]; // number of cells in patient x gene
  int N [N_sample]; // number of total cells
  int <lower = -1, upper = 1> X[N_sample]; // condition
}


parameters {
  real <lower = 0> beta_sigma;
  vector [N_gene] beta_raw [N_sample];
  real <lower = 0> phi;
  real<lower = 0> tau;

  vector [N_gene] alpha_gene;
  vector [N_gene] beta_gene;

  //inflation param
  vector <lower = 0, upper = 1> [N_gene] z_gene;
  real <lower = 0> z_phi;
  real <lower = 0> z_tau;
  real <lower = 0, upper = 1> z_mu;
}



transformed parameters {
  vector [N_gene] beta [N_sample];
  vector <lower = 0> [N_gene] a [N_sample];
  vector <lower = 0> [N_gene] b [N_sample];

  real <lower = 0> z_a;
  real <lower = 0> z_b;

  for(i in 1:N_sample) {
    beta[i] = beta_gene + beta_sigma * beta_raw[i];
    a[i] = inv_logit(alpha_gene + beta[i]*X[i]) * phi;
    b[i] = phi - a[i];
  }


  //inflation
  z_a = z_mu * z_phi;
  z_b = z_phi - z_a;
}


model {
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      target += zibb_lpmf(Y[j, i] | N[i], a[i][j], b[i][j], z_gene[j]);
    }
  }

  alpha_gene ~ normal(0, 20);
  beta_gene ~ normal(0, 5);

  beta_sigma ~ cauchy(0, 1);

  for(i in 1:N_sample) {
    beta_raw[i] ~ normal(0, 1);
  }

  phi ~ exponential(tau);
  tau ~ gamma(3, 0.1);

  //inflation
  // z ~ beta(1, 1);
  z_gene ~ beta(z_a, z_b);
  z_mu ~ beta(1, 3);
  z_phi ~ exponential(z_tau);
  z_tau ~ gamma(3, 0.1);
}


generated quantities {
  int Yhat [N_gene, N_sample];
  real Yhat_individual [N_gene, N_sample];
  matrix [2, N_gene] Yhat_gene;
  matrix [N_gene, N_sample] log_lik;
  real temp;

  for(j in 1:N_gene) {
    for(i in 1:N_sample) {
      log_lik[j,i] = zibb_lpmf(Y[j, i] | N[i], a[i][j], b[i][j], z_gene[j]);
      temp = N[i];
      Yhat[j, i] = zibb_rng(Y[j, i], N[i], a[i][j], b[i][j], z_gene[j]);


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
