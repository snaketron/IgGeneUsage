data {
  int <lower = 0> N_sample;
  int Y [N_sample];
  int N [N_sample];
  vector <lower = -1, upper = 1> [N_sample] X;
}


transformed data {
  real Nreal [N_sample];
  Nreal = N;
}


parameters {
  real <lower = 0> phi;
  real<lower = 0> tau;
  real alpha_gene;
  real beta_gene;
  vector <lower = 0, upper = 1> [N_sample] beta_raw;
  real <lower = 0> beta_sigma;
}


transformed parameters {
  vector <lower = 0> [N_sample] a;
  vector <lower = 0> [N_sample] b;
  vector [N_sample] beta;

  for(i in 1:N_sample) {
    beta[i] = beta_gene + beta_sigma * beta_raw[i];
    a[i] = inv_logit(alpha_gene + beta[i]*X[i]) * phi;
    b[i] = phi - a[i];
  }
}


model {
  for(i in 1:N_sample) {
    Y[i] ~ beta_binomial(N[i], a[i], b[i]);
  }

  alpha_gene ~ normal(0, 10);
  beta_gene ~ normal(0, 5);

  beta_sigma ~ cauchy(0, 1);
  beta_raw ~ normal(0, 1);

  phi ~ exponential(tau);
  tau ~ gamma(3, 0.1);
}


generated quantities {
  int Yhat [N_sample];
  real Yhat_individual [N_sample];
  real Yhat_group [2];

  for(i in 1:N_sample) {
    Yhat[i] = beta_binomial_rng(N[i], a[i], b[i]);

    if(Nreal[i] == 0.0) {
      Yhat_individual[i] = 0;
    } else {
      Yhat_individual[i] = Yhat[i]/Nreal[i]*100.0;
    }
  }

  Yhat_group[1] = 1/(1 + exp(-(alpha_gene + beta_gene)));
  Yhat_group[2] = 1/(1 + exp(-(alpha_gene + (-1)*beta_gene)));
}
