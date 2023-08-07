functions {
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
  int <lower=0> N_sample; // number of samples (repertoires)
  int <lower=0> N_gene; // number of genes
  int <lower=0> N_group; // number of groups
  int Y [N_gene, N_sample]; // number of successes (cells) in samples x gene
  int N [N_sample]; // number of total tries (repertoire size)
  int <lower=1> group_id [N_sample]; // group index for each condition
}

transformed data {
  // convert int N -> real N fo convenient division
  // in generated quantities block
  real Nr [N_sample];
  Nr = N;
}


parameters {
  // pop intercept mean
  real alpha_pop_mu;
  
  // scales
  real <lower = 0> alpha_pop_sigma;
  real <lower = 0> beta_pop_sigma;
  real <lower = 0> alpha_gene_sigma;
  real <lower = 0> beta_gene_sigma;
  
  // aux variables
  // vector [N_gene] alpha_z [N_sample];
  vector [N_gene] beta_z [N_sample];
  vector [N_gene] alpha_gene_z;
  vector [N_gene] beta_gene_z [N_group];
  
  // overdispersion
  real <lower = 0> phi;
  real<lower = 0> tau;
  
  // zero-inflation probability
  vector <lower = 0, upper = 1> [N_gene] z;
}



transformed parameters {
  vector [N_gene] beta [N_sample];
  vector [N_gene] alpha_gene_mu;
  vector [N_gene] beta_gene_mu [N_group];
  vector <lower = 0> [N_gene] a [N_sample];
  vector <lower = 0> [N_gene] b [N_sample];
  
  // non-centered params (at repertoire level)
  alpha_gene_mu = alpha_pop_mu + alpha_pop_sigma * alpha_gene_z;
  for(j in 1:N_group) {
    beta_gene_mu[j] = 0 + beta_pop_sigma * beta_gene_z[j];
  }

  // non-centered params (at gene level)
  for(i in 1:N_sample) {
    beta[i] = beta_gene_mu[group_id[i]] + beta_gene_sigma * beta_z[i];
    a[i] = inv_logit(alpha_gene_mu + beta[i]) * phi;
    b[i] = phi - a[i];
  }
}


model {
  // priors
  target += normal_lpdf(alpha_pop_mu | 0.0, 5.0);
  
  target += cauchy_lpdf(alpha_pop_sigma | 0.0, 1.0);
  target += cauchy_lpdf(beta_pop_sigma | 0.0, 1.0);
  target += cauchy_lpdf(alpha_gene_sigma | 0.0, 1.0);
  target += cauchy_lpdf(beta_gene_sigma | 0.0, 1.0);
  
  // zero-inflation
  target += beta_lpdf(z | 1.0, 10.0);
  
  //pareto 2 for overdispersion
  target += gamma_lpdf(tau | 3.0, 0.1);
  target += exponential_lpdf(phi | tau);
  
  // dummy
  for(i in 1:N_sample) {
    // target += std_normal_lpdf(alpha_z[i]);
    target += std_normal_lpdf(beta_z[i]);
  }
  target += std_normal_lpdf(alpha_gene_z);
  for(j in 1:N_group) {
    target += std_normal_lpdf(beta_gene_z[j]);
  }
  
  // likelihood
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      target += zibb_lpmf(Y[j, i] | N[i], a[i][j], b[i][j], z[j]);
    }
  }
}


generated quantities {
  // PPC: count usage
  int Yhat [N_gene, N_sample];

  // PPC: proportion usage
  real Yhat_rep [N_gene, N_sample];

  // PPC: proportion usage at a gene level in condition
  vector [N_gene] Yhat_condition [N_group];

  // LOG-LIK
  vector [N_gene] log_lik [N_sample];

  //TODO: speedup, run in C++ not big factor on performance
  for(j in 1:N_gene) {
    for(i in 1:N_sample) {
      log_lik[i][j] = zibb_lpmf(Y[j, i] | N[i], a[i][j], b[i][j], z[j]);
      Yhat[j, i] = zibb_rng(Y[j, i], N[i], a[i][j], b[i][j], z[j]);

      if(Nr[i] == 0.0) {
        Yhat_rep[j, i] = 0;
      }
      else {
        Yhat_rep[j, i] = Yhat[j,i]/Nr[i];
      }
    }
  }
  for(g in 1:N_group) {
    Yhat_condition[g] = inv_logit(alpha_gene_mu+beta_gene_mu[g]);
  }
}
