data {
  int <lower = 0> N_sample;
  int <lower = 0> N_gene;
  real Yp [N_gene, N_sample];
  vector <lower = -1, upper = 1> [N_sample] X;
}

parameters {
  real alpha_grand;
  real beta_grand;

  real <lower = 0> alpha_sigma;
  real <lower = 0> beta_sigma;
  real <lower = 0> beta_gene_sigma;

  vector [N_gene] alpha_raw;
  vector [N_gene] beta_gene_raw;
  vector [N_gene] beta_raw [N_sample];

  real <lower = 0> phi;
  real<lower = 0> tau;
}


transformed parameters {
  vector [N_gene] alpha_gene;
  vector [N_gene] beta_gene;
  vector <lower = 0> [N_gene] a [N_sample];
  vector <lower = 0> [N_gene] b [N_sample];
  vector [N_gene] beta_sample [N_sample];

  // non-centered params
  alpha_gene = alpha_grand + alpha_sigma * alpha_raw;
  beta_gene = beta_grand + beta_gene_sigma * beta_gene_raw;
  for(i in 1:N_sample) {
    beta_sample[i] = beta_gene + beta_sigma * beta_raw[i];
    a[i] = inv_logit(alpha_gene + beta_sample[i]*X[i]) * phi;
    b[i] = phi - a[i];
  }
}


model {
  for(i in 1:N_sample) {
    Yp[, i] ~ beta(a[i], b[i]);
  }

  alpha_grand ~ normal(0, 10);
  beta_grand ~ normal(0, 5);

  beta_sigma ~ cauchy(0, 1);
  for(i in 1:N_sample) {
    beta_raw[i] ~ normal(0, 1);
  }
  alpha_raw ~ normal(0, 1);
  beta_gene_raw ~ normal(0, 1);


  phi ~ exponential(tau);
  tau ~ gamma(3, 0.1);
}

//
// generated quantities {
//   real Yhat [N_sample];
//   real Yhat_group [2];
//
//   for(i in 1:N_sample) {
//     Yhat[i] = beta_rng(a[i], b[i]);
//   }
//
//   Yhat_group[1] = inv_logit(alpha_gene + beta_gene * 1.0);
//   Yhat_group[2] = inv_logit(alpha_gene + beta_gene * (-1.0));
// }
