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
    if(y==0) {
      if (bernoulli_rng(kappa) == 1) {
        return (0);
      } else {
        return (beta_binomial_rng(n, mu * phi, (1 - mu) * phi));
      }
    }
    return (beta_binomial_rng(n, mu * phi, (1 - mu) * phi));
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
  int<lower=0> N_sample;           // number of repertoires
  int<lower=0> N_gene;             // gene
  array[N_sample] int N;           // number of total tries (repertoire size)
  array[N_gene, N_sample] int Y;   // number of heads for each coin
  array[N_sample] int group_id;    // number of groups
}


transformed data {
  // convert int N -> real N fo convenient division
  // in generated quantities block
  array[N_sample] real Nr;
  Nr = N;
}


parameters {
  real <lower=0> phi;
  real <lower=0, upper=1> kappa;
  
  // gene
  vector [N_gene] alpha_gene_mu;
  vector <lower=0> [max(group_id)] beta_gene_sigma;
  array[N_sample] vector [N_gene] beta_z;
  
  // pop
  vector <lower=0> [max(group_id)] beta_pop_sigma;
  array[max(group_id)] vector [N_gene] beta_gene_mu_z;
}


transformed parameters {
  array[N_sample] vector <lower=0, upper=1> [N_gene] theta;
  array[N_sample] vector [N_gene] beta;
  array[max(group_id)] vector [N_gene] beta_gene_mu;
  
  
  for(i in 1:max(group_id)) {
    beta_gene_mu[i] = 0 + beta_pop_sigma[i] * beta_gene_mu_z[i];
  }
  
  for(i in 1:N_sample) {
    beta[i]  = beta_gene_mu[group_id[i]] + beta_gene_sigma[group_id[i]] * beta_z[i];
    theta[i] = inv_logit(alpha_gene_mu + beta[i]);
  }
}

model {
  target += beta_lpdf(kappa | 0.1, 1.0);
  target += exponential_lpdf(phi | 0.01);
  target += normal_lpdf(alpha_gene_mu | -5, 5);
  
  for(i in 1:max(group_id)) {
    target += std_normal_lpdf(beta_gene_mu_z[i]);
  }
  for(i in 1:N_sample) {
    target += std_normal_lpdf(beta_z[i]);
  }
  target += cauchy_lpdf(beta_gene_sigma | 0, 1);
  target += cauchy_lpdf(beta_pop_sigma | 0, 1);
  
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      target += zibb_lpmf(Y[j,i] | N[i], theta[i][j], phi, kappa);
    }
  }
}

generated quantities {
  // PPC: count usage (repertoire-level)
  array[N_gene, N_sample] int Yhat_rep;
  
  // PPC: proportion usage (repertoire-level)
  array[N_gene, N_sample] real Yhat_rep_prop;
  
  // PPC: proportion usage at a gene level in condition
  array[max(group_id)] vector [N_gene] Yhat_condition_prop;
  
  // LOG-LIK
  array[N_sample] vector [N_gene] log_lik;
  
  // DGU matrix
  matrix [N_gene, max(group_id)*(max(group_id)-1)/2] dgu;
  matrix [N_gene, max(group_id)*(max(group_id)-1)/2] dgu_prob;
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
    
    for(g in 1:max(group_id)) {
      Yhat_condition_prop[g][j] = z_rng(alpha_gene_mu[j], beta_gene_mu[g][j], 0);
    }
  }
  
  // DGU analysis
  for(i in 1:(max(group_id)-1)) {
    for(j in (i+1):max(group_id)) {
      dgu[,c] = beta_gene_mu[i]-beta_gene_mu[j];
      dgu_prob[,c]=to_vector(Yhat_condition_prop[i])-to_vector(Yhat_condition_prop[j]);
      c = c + 1;
    }
  }
}
