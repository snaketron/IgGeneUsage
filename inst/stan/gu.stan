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
  // here N_individual = N_sample as there are not replicates
  int<lower=0> N_sample;                   // number of repertoires
  int<lower=0> N_gene;                     // gene
  int<lower=0> N_individual;               // number of individuals
  int<lower=0> N_condition;                // number of conditions
  array [N_individual] int N;               // number of tries (repertoire size)
  array [N_gene, N_individual] int Y;       // number of heads for each coin
  array [N_individual] int condition_id;    // id of conditions
  array [N_individual] int individual_id;   // id of replicate
}

transformed data {
  // convert int N -> real N fo convenient division
  // in generated quantities block
  array [N_individual] real Nr;
  Nr = N;
}

parameters {
  real <lower=0> phi;
  real <lower=0, upper=1> kappa;
  real <lower=0> sigma_individual;
  array [N_individual] vector [N_gene] z_beta_individual;
  vector [N_gene] beta_condition;
}

transformed parameters {
  array [N_individual] vector <lower=0, upper=1> [N_gene] theta;
  array [N_individual] vector [N_gene] beta_sample;
  array [N_individual] vector [N_gene] beta_individual;
  
  for(i in 1:N_individual) {
    beta_individual[i] = beta_condition + sigma_individual * z_beta_individual[i];
    theta[i] = inv_logit(beta_individual[i]);
  }
}

model {
  target += beta_lpdf(kappa | 1.0, 5.0);
  target += exponential_lpdf(phi | 0.01);
  target += normal_lpdf(beta_condition | -3.0, 3.0);
  
  for(i in 1:N_individual) {
    target += std_normal_lpdf(z_beta_individual[i]);
  }
  
  target += cauchy_lpdf(sigma_individual | 0.0, 1.0);
  
  for(i in 1:N_individual) {
    for(j in 1:N_gene) {
      target += zibb_lpmf(Y[j,i] | N[i], theta[i][j], phi, kappa);
    }
  }
}

generated quantities {
  // PPC: count usage (repertoire-level)
  array [N_gene, N_individual] int Yhat_rep;
  
  // PPC: proportion usage (repertoire-level)
  array [N_gene, N_individual] real Yhat_rep_prop;
  
  // PPC: proportion usage at a gene level in condition
  vector [N_gene] Yhat_condition_prop;
  
  // LOG-LIK
  array [N_individual] vector [N_gene] log_lik;

  //TODO: speedup, run in C++ not big factor on performance
  for(j in 1:N_gene) {
    for(i in 1:N_individual) {
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
