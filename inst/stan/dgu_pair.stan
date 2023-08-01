functions {

  real zibb_lpmf(int y, int trials, real mu, real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         beta_binomial_lpmf(0 | trials,
                                            mu * phi,
                                            (1 - mu) * phi));
    } else {
      return bernoulli_lpmf(0 | zi) +
             beta_binomial_lpmf(y | trials, mu * phi, (1 - mu) * phi);
    }
  }

  int zibb_rng(int y, int trials, real mu, real phi, real zi) {
    if (bernoulli_rng(zi) == 1) {
      return (0);
    } else {
      return (beta_binomial_rng(trials, mu*phi, (1-mu)*phi));
    }
  }
}


data {
  int <lower = 0> N_sample; // number of samples (repertoires)
  int <lower = 0> N_gene; // number of genes
  int Y_1 [N_gene, N_sample]; // number of successes (cells) in samples x gene
  int Y_2 [N_gene, N_sample]; // number of successes (cells) in samples x gene
  int N_1 [N_sample]; // number of total tries 
  int N_2 [N_sample]; // number of total tries
}

transformed data {
  real Nr_1 [N_sample];
  real Nr_2 [N_sample];
  Nr_1 = N_1;
  Nr_2 = N_2;
}

parameters {
  // pop intercept mean
  real alpha_pop_mu;
  
  // scales
  real <lower = 0> alpha_gene_sigma;
  real <lower = 0> alpha_pop_sigma;
  real <lower = 0> beta_gene_sigma;
  real <lower = 0> beta_pop_sigma;
  
  // aux variables
  vector [N_gene] alpha_gene_z;
  vector [N_gene] beta_gene_z;
  vector [N_gene] alpha_z [N_sample];
  vector [N_gene] beta_z [N_sample];
  
  // overdispersion
  real <lower = 0> phi;
  real<lower = 0> tau;
  
  // zero-inflation probability
  vector <lower = 0, upper = 1> [N_gene] z;
  real<lower=0, upper=1> z_mu;
  real<lower=0> z_phi;
}


transformed parameters {
  vector [N_gene] alpha [N_sample];
  vector [N_gene] beta [N_sample];
  vector [N_gene] alpha_gene_mu;
  vector [N_gene] beta_gene_mu;
  
  // non-centered params (at gene pop. level)
  alpha_gene_mu = alpha_pop_mu + alpha_pop_sigma * alpha_gene_z;
  beta_gene_mu = 0 + beta_pop_sigma * beta_gene_z;

  // non-centered params (at gene level)
  for(i in 1:N_sample) {
    alpha[i] = alpha_gene_mu + alpha_gene_sigma * alpha_z[i];
    beta[i] = beta_gene_mu + beta_gene_sigma * beta_z[i];
  }
}


model {
  // prior pop. intercept
  target += normal_lpdf(alpha_pop_mu | 0.0, 5.0);
  
  // prior on scales
  target += cauchy_lpdf(alpha_pop_sigma | 0.0, 1.0);
  target += cauchy_lpdf(beta_pop_sigma | 0.0, 1.0);
  target += cauchy_lpdf(alpha_gene_sigma | 0.0, 1.0);
  target += cauchy_lpdf(beta_gene_sigma | 0.0, 1.0);
  
  // zero-inflation
  target += exponential_lpdf(z_phi | 0.05);
  target += beta_lpdf(z_mu | 1.0, 10.0);
  target += beta_proportion_lpdf(z | z_mu, z_phi);
  
  //pareto 2 for overdispersion
  target += gamma_lpdf(tau | 3.0, 0.1);
  target += exponential_lpdf(phi | tau);
  
  // prior on aux. parameters
  for(i in 1:N_sample) {
    target += std_normal_lpdf(alpha_z[i]);
    target += std_normal_lpdf(beta_z[i]);
  }
  target += std_normal_lpdf(alpha_gene_z);
  target += std_normal_lpdf(beta_gene_z);
  
  // likelihood: TODO speedup
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      target += zibb_lpmf(Y_1[j,i] | N_1[i], inv_logit(alpha[i][j]-beta[i][j]), phi, z[j]);
      target += zibb_lpmf(Y_2[j,i] | N_2[i], inv_logit(alpha[i][j]+beta[i][j]), phi, z[j]);
    }
  }
}

generated quantities {
  // PPC: count usage
  int Yhat_1 [N_gene,N_sample];
  int Yhat_2 [N_gene,N_sample];
  
  real Yhat_rep_1 [N_gene,N_sample];
  real Yhat_rep_2 [N_gene,N_sample];

  // LOG-LIK
  real a [N_gene, N_sample, 2];
  real log_lik [N_sample, 2];
  
  // PPC: proportion usage at a gene level in condition
  vector [N_gene] Yhat_condition [2];
  
  real mu [2];

  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      mu[1] = inv_logit(alpha[i][j]-beta[i][j]);
      mu[2] = inv_logit(alpha[i][j]+beta[i][j]);
      a[j,i,1] = zibb_lpmf(Y_1[j,i] | N_1[i], mu[1], phi, z[j]);
      a[j,i,2] = zibb_lpmf(Y_2[j,i] | N_2[i], mu[2], phi, z[j]);
      Yhat_1[j,i] = zibb_rng(Y_1[j,i],N_1[i], mu[1], phi, z[j]);
      Yhat_2[j,i] = zibb_rng(Y_2[j,i],N_2[i], mu[2], phi, z[j]);
      
      if(Nr_1[i] == 0.0) {
        Yhat_rep_1[j,i] = 0;
      }
      else {
        Yhat_rep_1[j,i] = Yhat_1[j,i]/Nr_1[i];
      }
      if(Nr_2[i] == 0.0) {
        Yhat_rep_2[j,i] = 0;
      }
      else {
        Yhat_rep_2[j,i] = Yhat_2[j,i]/Nr_2[i];
      }
      
      // j loop only (i==1)
      if(i == 1) {
        Yhat_condition[1][j] = inv_logit(alpha_gene_mu[j]-beta_gene_mu[j]);
        Yhat_condition[2][j] = inv_logit(alpha_gene_mu[j]+beta_gene_mu[j]);
      }
    }
    log_lik[i,1] = sum(a[,i,1]);
    log_lik[i,2] = sum(a[,i,2]);
  }
}
