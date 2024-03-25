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
  int<lower=0> N_sample;                   // number of repertoires 
  int<lower=0> N_gene;                     // gene
  int<lower=0> N_individual;               // number of individuals
  int<lower=0> N_condition;                // number of conditions
  int<lower=0> N_replicate;                // number of replicates
  array [N_individual] int N;              // number of tries
  array [N_gene, N_individual] int Y;      // number of heads for each coin
  array [N_individual] int condition_id;   // id of conditions
  array [N_sample] int individual_id;      // id of individual
  array [N_sample] int replicate_id;       // id of replicate
}

transformed data {
  // convert int to real N for easier division in generated quantities block
  array [N_individual] real Nr;
  Nr = N;
}

parameters {
  real <lower=0> phi;
  real <lower=0, upper=1> kappa;
  
  vector [N_gene] alpha;
  
  vector <lower=0> [N_condition] sigma_condition;
  vector <lower=0> [N_condition] sigma_individual;
  real <lower=0> sigma_alpha;
  real <lower=0> sigma_alpha_rep;
  real <lower=0> sigma_beta_rep;
  
  array [N_individual] vector [N_gene] z_alpha_individual;
  array [N_individual] vector [N_gene] z_beta_individual;
  array [N_condition] vector [N_gene] z_beta_condition;
  array [N_individual, N_replicate] vector [N_gene] z_alpha_sample;
  array [N_individual, N_replicate] vector [N_gene] z_beta_sample;
}

transformed parameters {
  array [N_condition] vector [N_gene] beta_condition;
  array [N_individual] vector [N_gene] alpha_individual;
  array [N_individual] vector [N_gene] beta_individual;
  array [N_individual, N_replicate] vector [N_gene] alpha_sample;
  array [N_individual, N_replicate] vector [N_gene] beta_sample;
  array [N_individual] vector <lower=0, upper=1> [N_gene] theta;
  
  for(i in 1:N_condition) {
    beta_condition[i] = 0 + sigma_condition[i] * z_beta_condition[i];
  }
  
  for(i in 1:N_individual) {
    alpha_individual[i] = alpha + sigma_alpha * z_alpha_individual[i];
    beta_individual[i]  = beta_condition[condition_id[i]] + sigma_individual[condition_id[i]] * z_beta_individual[i];
  }
  
  for(i in 1:N_sample) {
    alpha_sample[individual_id[i], replicate_id[i]] = alpha_individual[individual_id[i]] + sigma_alpha_rep * z_alpha_sample[individual_id[i], replicate_id[i]];
    beta_sample[individual_id[i], replicate_id[i]] = beta_individual[individual_id[i]] + sigma_beta_rep * z_beta_sample[individual_id[i], replicate_id[i]];
    theta[i] = inv_logit(alpha_sample[individual_id[i], replicate_id[i]] + beta_sample[individual_id[i], replicate_id[i]]);
  }
}

model {
  target += beta_lpdf(kappa | 1.0, 5.0);
  target += exponential_lpdf(phi | 0.01);
  target += normal_lpdf(alpha | -3.0, 3.0);
  
  for(i in 1:N_condition) {
    target += std_normal_lpdf(z_beta_condition[i]);
  }
  for(i in 1:N_individual) {
    target += std_normal_lpdf(z_beta_individual[i]);
  }
  for(i in 1:N_sample) {
    target += std_normal_lpdf(z_alpha_sample[individual_id[i], replicate_id[i]]);
    target += std_normal_lpdf(z_beta_sample[individual_id[i], replicate_id[i]]);
  }
  
  target += cauchy_lpdf(sigma_individual | 0.0, 1.0);
  target += cauchy_lpdf(sigma_condition | 0.0, 1.0);
  target += cauchy_lpdf(sigma_alpha | 0.0, 1.0);
  target += cauchy_lpdf(sigma_alpha_rep | 0.0, 1.0);
  target += cauchy_lpdf(sigma_beta_rep | 0.0, 1.0);
  
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
  array [N_condition] vector [N_gene] Yhat_condition_prop;
  
  // LOG-LIK
  array [N_individual] vector [N_gene] log_lik;
  
  // DGU matrix
  matrix [N_gene, N_condition*(N_condition-1)/2] dgu;
  matrix [N_gene, N_condition*(N_condition-1)/2] dgu_prob;
  int c = 1;
  
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
    for(g in 1:N_condition) {
      Yhat_condition_prop[g][j] = z_rng(alpha[j], beta_condition[g][j], 0);
    }
  }
  
  // DGU analysis
  for(i in 1:(N_condition-1)) {
    for(j in (i+1):N_condition) {
      dgu[,c] = beta_condition[i]-beta_condition[j];
      dgu_prob[,c]=to_vector(Yhat_condition_prop[i])-to_vector(Yhat_condition_prop[j]);
      c = c + 1;
    }
  }
}
