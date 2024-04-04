functions {
  real zibb_lpmf(int y, int n, real mu, real phi, real kappa) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | kappa),
      bernoulli_lpmf(0 | kappa) +
      beta_binomial_lpmf(0 | n, mu * phi, (1 - mu) * phi));
    } else {
      return bernoulli_lpmf(0 | kappa) +
      beta_binomial_lpmf(y | n, mu * phi, (1 - mu) * phi);
    }
  }
  
  int zibb_rng(int y, int n, real mu, real phi, real kappa) {
    if (bernoulli_rng(kappa) == 1) {
      return (0);
    } else {
      return (beta_binomial_rng(n, mu * phi, (1 - mu) * phi));
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
  int<lower=0> N_sample;                    // number of repertoires
  int<lower=0> N_gene;                      // gene
  int<lower=0> N_individual;                // number of individuals
  array [N_sample] int N;                   // number of tries (repertoire size)
  array [N_gene, N_sample] int Y;           // number of heads for each coin
  array [N_sample] int individual_id;       // id of replicate
}

transformed data {
  // convert int to real N for convenient division in generated quantities block
  array [N_sample] real Nr;
  Nr = N;
}

parameters {
  real <lower=0> phi;
  real <lower=0, upper=1> kappa;
  
  real <lower=0> sigma_individual;
  real <lower=0> sigma_beta_rep;
  
  array [N_sample] vector [N_gene] z_beta_sample;
  array [N_individual] vector [N_gene] z_beta_individual;
  
  vector [N_gene] beta_condition;
}

transformed parameters {
  array [N_sample] vector <lower=0, upper=1> [N_gene] theta;
  array [N_sample] vector [N_gene] beta_sample;
  array [N_individual] vector [N_gene] beta_individual;
  
  for(i in 1:N_individual) {
    beta_individual[i] = beta_condition + sigma_individual * z_beta_individual[i];
  }
  
  for(i in 1:N_sample) {
    beta_sample[i]  = beta_individual[individual_id[i]] + sigma_beta_rep * z_beta_sample[i];
    theta[i] = inv_logit(beta_sample[i]);
  }
}

model {
  target += beta_lpdf(kappa | 1.0, 5.0);
  target += exponential_lpdf(phi | 0.01);
  target += normal_lpdf(beta_condition | -3.0, 3.0);
  
  for(i in 1:N_individual) {
    target += std_normal_lpdf(z_beta_individual[i]);
  }
  for(i in 1:N_sample) {
    target += std_normal_lpdf(z_beta_sample[i]);
  }
  
  target += cauchy_lpdf(sigma_beta_rep | 0.0, 1.0);
  target += cauchy_lpdf(sigma_individual | 0.0, 1.0);
  
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      target += zibb_lpmf(Y[j,i] | N[i], theta[i][j], phi, kappa);
    }
  }
}

generated quantities {
  // PPC: count usage (repertoire-level)
  array [N_gene, N_sample] int Yhat_rep;
  
  // PPC: proportion usage (repertoire-level)
  array [N_gene, N_sample] real Yhat_rep_prop;
  
  // PPC: proportion usage at a gene level in condition
  vector [N_gene] Yhat_condition_prop;
  
  // LOG-LIK
  array [N_sample] vector [N_gene] log_lik;

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
    Yhat_condition_prop[j] = z_rng(0, beta_condition[j], 0);
  }
}
