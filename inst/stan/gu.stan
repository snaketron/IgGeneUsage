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
  
  real z_rng(real a, real b, real zi) {
    if (bernoulli_rng(zi) == 1) {
      return (0);
    } else {
      return(inv_logit(a+b));
    }
  }
}

data {
  int <lower = 0> N_sample; // number of samples (repertoires)
  int <lower = 0> N_gene; // number of genes
  int Y [N_gene, N_sample]; // number of successes (cells) in samples x gene
  int N [N_sample]; // number of total tries (repertoire size)
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
  
  // aux variables
  vector [N_gene] alpha_gene_z;
  
  // overdispersion
  real <lower = 0> phi;
  // real<lower = 0> tau;
  
  // zero-inflation probability
  vector <lower = 0, upper = 1> [N_gene] kappa;
}

transformed parameters {
  vector [N_gene] alpha_gene_mu;
  vector <lower = 0> [N_gene] a [N_sample];
  vector <lower = 0> [N_gene] b [N_sample];
  
  // non-centered params (at repertoire level)
  alpha_gene_mu = alpha_pop_mu + alpha_pop_sigma * alpha_gene_z;
  // non-centered params (at gene level)
  for(i in 1:N_sample) {
    a[i] = inv_logit(alpha_gene_mu) * phi;
    b[i] = phi - a[i];
  }
}


model {
  // priors
  target += normal_lpdf(alpha_pop_mu | 0.0, 5.0);
  
  target += cauchy_lpdf(alpha_pop_sigma | 0.0, 1.0);
  
  // zero-inflation
  target += beta_lpdf(kappa | 0.1, 1.0);
  
  //pareto 2 for overdispersion
  // target += gamma_lpdf(tau | 3.0, 0.1);
  // target += exponential_lpdf(phi | tau);
  target += exponential_lpdf(phi | 0.01);
  target += std_normal_lpdf(alpha_gene_z);
  
  // likelihood
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      target += zibb_lpmf(Y[j, i] | N[i], a[i][j], b[i][j], kappa[j]);
    }
  }
}

generated quantities {
  // PPC: count usage
  int Yhat [N_gene, N_sample];

  // PPC: proportion usage
  real Yhat_rep [N_gene, N_sample];

  // PPC: proportion usage at a gene level in condition
  vector [N_gene] Yhat_condition;

  // LOG-LIK
  vector [N_gene] log_lik [N_sample];
  
  // probability
  vector [N_gene] theta_condition;

  //TODO: speedup, run in C++ not big factor on performance
  for(j in 1:N_gene) {
    for(i in 1:N_sample) {
      log_lik[i][j] = zibb_lpmf(Y[j, i] | N[i], a[i][j], b[i][j], kappa[j]);
      Yhat[j, i] = zibb_rng(Y[j, i], N[i], a[i][j], b[i][j], kappa[j]);

      if(Nr[i] == 0.0) {
        Yhat_rep[j, i] = 0;
      }
      else {
        Yhat_rep[j, i] = Yhat[j,i]/Nr[i];
      }
    }
    
    // if kappa=0 ->0, else -> inverse_logit(a)
    Yhat_condition[j] = z_rng(alpha_gene_mu[j], 0, kappa[j]);
    theta_condition[j] = z_rng(alpha_gene_mu[j], 0, kappa[j]);
  }
}
