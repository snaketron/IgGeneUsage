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
  real alpha_gene;
  vector [N_sample] b;
  real beta_gene;
  real <lower = 0> b_sigma;
}


model {
  for(i in 1:N_sample) {
    Y[i] ~ binomial_logit(N[i], alpha_gene + b[i] * X[i]);
  }
  alpha_gene ~ normal(0, 10);
  b ~ normal(beta_gene, b_sigma);
  beta_gene ~ normal(0, 5);
  b_sigma ~ cauchy(0, 1);
}

generated quantities {
  int Yhat [N_sample];
  real Yhat_individual [N_sample];
  real Yhat_group [2];

  for(i in 1:N_sample) {
    Yhat[i] = binomial_rng(N[i], 1/(1 + exp(-(alpha_gene+b[i]*X[i]))));

    if(Nreal[i] == 0.0) {
      Yhat_individual[i] = 0;
    } else {
      Yhat_individual[i] = Yhat[i]/Nreal[i]*100.0;
    }
  }

  Yhat_group[1] = 1/(1 + exp(-(alpha_gene + beta_gene)));
  Yhat_group[2] = 1/(1 + exp(-(alpha_gene + (-1)*beta_gene)));
}
