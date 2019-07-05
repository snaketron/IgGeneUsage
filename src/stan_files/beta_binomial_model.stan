data {
  int <lower = 0> N_sample; // number of samples
  int <lower = 0> N_gene; // number of genes
  int Y [N_gene, N_sample]; // number of cells in patient x gene
  int N [N_sample]; // number of total cells
  int <lower = -1, upper = 1> X[N_sample]; // condition
}


parameters {
  real alpha_grand;
  real beta_grand;
  real <lower = 0> alpha_sigma;
  real <lower = 0> beta_sigma;
  real <lower = 0> beta_gene_sigma;
  vector [N_gene] alpha_raw;
  vector [N_gene] beta_raw [N_sample];
  vector [N_gene] beta_gene_raw;
  real <lower = 0> phi;
  real<lower = 0> tau;
}


transformed parameters {
  vector [N_gene] beta [N_sample];
  vector [N_gene] alpha_gene;
  vector [N_gene] beta_gene;
  vector <lower = 0> [N_gene] a [N_sample];
  vector <lower = 0> [N_gene] b [N_sample];

  // non-centered params
  alpha_gene = alpha_grand + alpha_sigma * alpha_raw;
  beta_gene = beta_grand + beta_gene_sigma * beta_gene_raw;

  for(i in 1:N_sample) {
    beta[i] = beta_gene + beta_sigma * beta_raw[i];
    a[i] = inv_logit(alpha_gene + beta[i]*X[i]) * phi;
    b[i] = phi - a[i];
  }
}


model {
  int Ytemp [N_gene];

  for(i in 1:N_sample) {
    Ytemp = Y[, i];
    Ytemp ~ beta_binomial(N[i], a[i], b[i]);
  }

  alpha_grand ~ normal(0, 20);
  beta_grand ~ normal(0, 5);

  alpha_sigma ~ cauchy(0, 1);
  beta_sigma ~ cauchy(0, 1);
  beta_gene_sigma ~ cauchy(0, 1);

  alpha_raw ~ normal(0, 1);
  for(i in 1:N_sample) {
    beta_raw[i] ~ normal(0, 1);
  }
  beta_gene_raw ~ normal(0, 1);

  phi ~ exponential(tau);
  tau ~ gamma(3, 0.1);
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

      log_lik[j,i] = beta_binomial_lpmf(Y[j,i] | N[i], a[i][j], b[i][j]);
      Yhat[j, i] = beta_binomial_rng(N[i], a[i][j], b[i][j]);


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
