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
  real a;
  real b;
}


parameters {
  real <lower = 0> beta_sigma;
  vector [N_gene] beta_raw [N_sample];
  real <lower = 0> phi;
  real<lower = 0> tau;

  vector [N_gene] alpha_gene;
  vector [N_gene] beta_gene;
}


transformed parameters {
  vector [N_gene] beta [N_sample];
  vector <lower = 0> [N_gene] a [N_sample];
  vector <lower = 0> [N_gene] b [N_sample];

  for(i in 1:N_sample) {
    beta[i] = beta_gene + beta_sigma * beta_raw[i];
    a[i] = inv_logit(alpha_gene + beta[i]*X[i]) * phi;
    b[i] = phi - a[i];
  }
}


model {
  Y ~ beta_binomial(N, a, b);
}


generated quantities {
  int Yhat [N_sample];
  real Yhat_individual [N_sample];

  for(i in 1:N_sample) {
    Yhat[i] = binomial_rng(N[i], 1/(1 + exp(-(a+b*X[i]))));

    if(Nreal[i] == 0.0) {
      Yhat_individual[i] = 0;
    } else {
      Yhat_individual[i] = Yhat[i]/Nreal[i]*100.0;
    }
  }
}
