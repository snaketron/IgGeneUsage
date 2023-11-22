functions {
  real zibb_lpmf(int y, int n, real theta, real phi, real kappa) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | kappa),
                         bernoulli_lpmf(0 | kappa) +
                         beta_binomial_lpmf(0 | n, theta * phi, (1 - theta) * phi));
    } else {
      return bernoulli_lpmf(0 | kappa) +
             beta_binomial_lpmf(y | n, theta * phi, (1 - theta) * phi);
    }
  }
  
  int zibb_rng(int y, int n, real theta, real phi, real kappa) {
    if (bernoulli_rng(kappa) == 1) {
      return (0);
    } else {
      return (beta_binomial_rng(n, theta * phi, (1 - theta) * phi));
    }
  }
  
  real z_rng(real a, real b, real kappa) {
    if (bernoulli_rng(kappa) == 1) {
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
  // dispersion + zero-inflation
  real <lower = 0> phi;
  real <lower = 0, upper = 1> kappa;
  
  // gene
  vector [N_gene] alpha_gene_mu;
  real <lower = 0> alpha_gene_sigma;
  vector [N_gene] alpha_z [N_sample];
}

transformed parameters {
  vector [N_gene] alpha [N_sample];
  vector <lower = 0, upper=1> [N_gene] theta [N_sample];
  
  for(i in 1:N_sample) {
    alpha[i] = alpha_gene_mu + alpha_gene_sigma * alpha_z[i];
    theta[i] = inv_logit(alpha[i]);
  }
}

model {
  // priors
  target += normal_lpdf(alpha_gene_mu | -5.0, 5.0);
  target += cauchy_lpdf(alpha_gene_sigma | 0.0, 1.0);
  target += beta_lpdf(kappa | 0.1, 1.0);
  target += exponential_lpdf(phi | 0.01);
  for(i in 1:N_sample) {
    target += std_normal_lpdf(alpha_z[i]);
  }

  // likelihood
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      target += zibb_lpmf(Y[j, i] | N[i], theta[i][j], phi, kappa);
    }
  }
}

generated quantities {
  // PPC: count usage (repertoire-level)
  int Yhat_rep [N_gene, N_sample];
  
  // PPC: proportion usage (repertoire-level)
  real Yhat_rep_prop [N_gene, N_sample];
  
  // PPC: proportion usage at a gene level in condition
  vector [N_gene] Yhat_condition_prop;
  
  // LOG-LIK
  vector [N_gene] log_lik [N_sample];
  
  
  //TODO: speedup, run in C++ not big factor on performance
  for(j in 1:N_gene) {
    for(i in 1:N_sample) {
      Yhat_rep[j, i] = zibb_rng(Y[j, i], N[i], theta[i][j], phi, kappa);
      log_lik[i][j] = zibb_lpmf(Y[j, i] | N[i], theta[i][j], phi, kappa);
      
      if(Nr[i] == 0.0) {
        Yhat_rep_prop[j, i] = 0;
      }
      else {
        Yhat_rep_prop[j, i] = Yhat_rep[j,i]/Nr[i];
      }
    }
    Yhat_condition_prop[j] = z_rng(alpha_gene_mu[j], 0, 0);
  }
}
