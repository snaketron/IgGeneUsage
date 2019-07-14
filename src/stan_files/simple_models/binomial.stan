data {
  int <lower = 0> N_sample; // number of samples
  int Y [N_sample]; // number of cells in patient x gene
  int N [N_sample]; // number of total cells
  int <lower = -1, upper = 1> X[N_sample]; // condition
}

transformed data {
  // convert int N -> real N fo convenient division
  // in generated quantities block
  real Nreal [N_sample];
  Nreal = N;
}

parameters {
  real <lower = 0> beta_sigma;
  vector [N_gene] beta_raw [N_sample];
  real <lower = 0> phi;
  real<lower = 0> tau;

  vector [N_gene] alpha_gene;
  vector [N_gene] beta_gene;
}




model {
    Y ~ binomial_logit(N, a+b*X);

  alpha_gene ~ normal(0, 20);
  beta_gene ~ normal(0, 5);

  beta_sigma ~ cauchy(0, 1);

  for(i in 1:N_sample) {
    beta_raw[i] ~ normal(0, 1);
  }

  phi ~ exponential(tau);
  tau ~ gamma(3, 0.1);
}


generated quantities {
  int Yhat [N_sample];
  real Yhat_individual [N_sample];

  for(i in 1:N_sample) {
      Yhat[j, i] = binomial_logit_rng(N[i], a+b);


      if(Nreal[i] == 0.0) {
        Yhat_individual[j, i] = 0;
      }
      else {
        Yhat_individual[j, i] = Yhat[j,i]/Nreal[i]*100.0;
      }
    }
    Yhat_gene[1, j] = inv_logit(alpha_gene[j]+beta_gene[j]*1.0)*100.0;
    Yhat_gene[2, j] = inv_logit(alpha_gene[j]+beta_gene[j]*(-1.0))*100.0;

}
