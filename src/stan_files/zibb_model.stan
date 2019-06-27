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
  real zero_inflated_beta_binomial_lpmf(int y, int trials,
                                        real mu, real phi,
                                        real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         beta_binomial_lpmf(0 | trials, mu*phi, (1-mu)*phi));
    } else {
      return bernoulli_lpmf(0 | zi) +
             beta_binomial_lpmf(y | trials, mu*phi, (1-mu)*phi);
    }
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
  real <lower = 0> beta_gene_sigma;
  vector [N_gene] alpha_raw;
  vector [N_gene] beta_raw [N_sample];
  vector [N_gene] beta_gene_raw;
  real <lower = 0> phi;
  vector <lower=0,upper=1> [N_gene] z;  // zero-inflation probability
}



transformed parameters {
  vector [N_gene] beta [N_sample];
  vector <lower = 0, upper = 1> [N_gene] mu [N_sample];
  vector [N_gene] alpha_gene;
  vector [N_gene] beta_gene;

  // non-centered params
  alpha_gene = alpha_grand + alpha_sigma * alpha_raw;
  beta_gene = beta_grand + beta_gene_sigma * beta_gene_raw;
  for(i in 1:N_sample) {
    beta[i] = beta_gene + beta_sigma * beta_raw[i];
  }

  for(i in 1:N_sample) {
    mu[i] = inv_logit(alpha_gene + beta[i]*X[i]);
  }
}


model {
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      target += zero_inflated_beta_binomial_lpmf(Y[j, i]|N[i], mu[i][j] * phi, (1 - mu[i][j]) * phi, z[j]);
    }
    // Y[, i] ~ beta_binomial(N[i], mu[i] * phi, (1 - mu[i]) * phi);
  }

  alpha_grand ~ normal(0, 20);
  beta_grand ~ normal(0, 5);

  alpha_sigma ~ cauchy(0, 3);
  beta_sigma ~ cauchy(0, 3);
  beta_gene_sigma ~ cauchy(0, 3);

  alpha_raw ~ normal(0, 1);
  for(i in 1:N_sample) {
    beta_raw[i] ~ normal(0, 1);
  }
  beta_gene_raw ~ normal(0, 1);
  phi ~ gamma(0.01, 0.01);
  z ~ beta(1, 1);
}

//
// generated quantities {
//   int Yhat [N_gene, N_sample];
//   real Yhat_individual [N_gene, N_sample];
//   matrix [2, N_gene] Yhat_gene;
//   matrix [N_gene, N_sample] log_lik;
//   real temp;
//
//   for(j in 1:N_gene) {
//     for(i in 1:N_sample) {
//       temp = N[i];
//
//       log_lik[j,i] = beta_binomial_lpmf(Y[j,i] | N[i], mu[i][j] * phi, (1 - mu[i][j]) * phi);
//       Yhat[j, i] = beta_binomial_rng(N[i], mu[i][j] * phi, (1 - mu[i][j]) * phi);
//
//
//       if(temp == 0.0) {
//         Yhat_individual[j, i] = 0;
//       }
//       else {
//         Yhat_individual[j, i] = Yhat[j,i]/temp*100.0;
//       }
//     }
//     Yhat_gene[1, j] = inv_logit(alpha_gene[j]+beta_gene[j]*1.0)*100.0;
//     Yhat_gene[2, j] = inv_logit(alpha_gene[j]+beta_gene[j]*(-1.0))*100.0;
//   }
// }
