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
  int<lower=0> N_sample;                        // # samples 
  int<lower=0> N_gene;                          // # genes
  int<lower=0> N_individual;                    // # individuals
  int<lower=0> N_condition;                     // # conditions
  array [N_sample] int N;                       // rep size
  array [N_gene, N_sample] int Y;               // gene usage
  array [N_sample] int condition_id_of_sample;  // condition of sample
  array [N_sample] int individual_id;           // individual id
}

transformed data {
  // convert int to real N for easier division in generated quantities block
  array [N_sample] real Nr;
  Nr = N;
}

parameters {
  real <lower=0> phi;
  real <lower=0, upper=1> kappa;
  vector [N_gene] alpha;
  real <lower=0> sigma_individual;
  vector <lower=0> [N_condition] sigma_condition;
  array [N_condition] vector [N_gene] z_beta_condition;
  array [N_individual] vector [N_gene] z_beta_individual;
  vector [N_individual] mu_individual;
}

transformed parameters {
  array [N_sample] vector <lower=0, upper=1> [N_gene] theta;
  array [N_individual] vector [N_gene] beta_individual;
  array [N_condition] vector [N_gene] beta_condition;
  
  for(j in 1:N_condition) {
    beta_condition[j] = 0 + sigma_condition[j] * z_beta_condition[j];
  }
  
  
  for(i in 1:N_individual) {
    beta_individual[i] = mu_individual[i] + sigma_individual * z_beta_individual[i];
  }

  for(i in 1:N_sample) {
    theta[i] = inv_logit(alpha + beta_individual[individual_id[i]] + beta_condition[condition_id_of_sample[i]]);
  }
}

model {
  target += beta_lpdf(kappa | 1.0, 5.0);
  target += exponential_lpdf(phi | 0.01);
  target += normal_lpdf(alpha | -3.0, 3.0);
  for(i in 1:N_individual) {
    target += normal_lpdf(mu_individual[i] | 0.0, 1.0);
    target += std_normal_lpdf(z_beta_individual[i]);
  }
  for(j in 1:N_condition) {
    target += std_normal_lpdf(z_beta_condition[j]);
  }
  target += cauchy_lpdf(sigma_condition | 0.0, 1.0);
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
  array [N_condition] vector [N_gene] Yhat_condition_prop;
  
  // LOG-LIK
  array [N_sample] vector [N_gene] log_lik;
  
  // DGU matrix
  matrix [N_gene, N_condition*(N_condition-1)/2] dgu;
  matrix [N_gene, N_condition*(N_condition-1)/2] dgu_prob;
  int c = 1;
  
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
