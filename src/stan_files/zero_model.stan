data {
  int <lower = 0> N_sample; // number of samples
  int <lower = 0> N_gene; // number of genes
  int Y [N_gene, N_sample]; // number of cells in patient x gene
  int N [N_sample]; // number of total cells
  int <lower = -1, upper = 1> X[N_sample]; //condition
}

parameters {
  real alpha_grand;
  real beta_grand;
  real <lower = 0> alpha_gene_sigma;
  real <lower = 0> beta_gene_sigma;
  vector [N_gene] alpha_gene_raw;
  vector [N_gene] beta_gene_raw;
}


transformed parameters {
  vector [N_gene] alpha_gene;
  vector [N_gene] beta_gene;
  vector [N_gene] mu [N_sample];

  // non-centered params
  alpha_gene = alpha_grand + alpha_gene_sigma * alpha_gene_raw;
  beta_gene = beta_grand + beta_gene_sigma * beta_gene_raw;
  for(i in 1:N_sample) {
    mu[i] = alpha_gene + beta_gene*X[i];
  }
}


model {
  for(i in 1:N_sample) {
    Y[, i] ~ binomial_logit(N[i], mu[i]);
  }

  alpha_grand ~ normal(0, 20);
  beta_grand ~ normal(0, 5);

  alpha_gene_sigma ~ cauchy(0, 3);
  beta_gene_sigma ~ cauchy(0, 3);

  alpha_gene_raw ~ normal(0, 1);
  beta_gene_raw ~ normal(0, 1);
}


generated quantities {
  int Yhat [N_gene, N_sample];
  real Yhat_individual [N_gene, N_sample];
  matrix [2, N_gene] Yhat_gene;
  matrix [N_gene, N_sample] log_lik;
  real temp;

  for(j in 1:N_gene) {
    for(i in 1:N_sample) {
      temp = N[i];

      log_lik[j,i] = binomial_logit_lpmf(Y[j,i]|N[i],alpha_gene[j]+beta_gene[j]*X[i]);
      Yhat[j, i] = binomial_rng(N[i], inv_logit(alpha_gene[j]+beta_gene[j]*X[i]));

      if(temp == 0.0) {
        Yhat_individual[j, i] = 0;
      }
      else {
        Yhat_individual[j, i] = Yhat[j,i]/temp*100.0;
      }
    }
    Yhat_gene[1, j] = inv_logit(alpha_gene[j]+beta_gene[j]*1.0)*100.0;
    Yhat_gene[2, j] = inv_logit(alpha_gene[j]+beta_gene[j]*(-1.0))*100.0;
  }
}