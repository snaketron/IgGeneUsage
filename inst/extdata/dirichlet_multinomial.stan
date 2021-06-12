data {
  int <lower = 0> N_sample; // number of samples (repertoires)
  int <lower = 0> N_gene; // number of genes
  int Y [N_gene, N_sample]; // number of successes (cells) in samples x gene
  int <lower = -1, upper = 1> X[N_sample]; // condition
}

parameters {
  simplex [N_gene] theta [N_sample];
  real<lower=0> phi;
  real<lower = 0> tau;
  
  real alpha_grand;
  real beta_grand;
  
  real <lower = 0> alpha_sigma;
  real <lower = 0> beta_sigma;
  real <lower = 0> alpha_gene_sigma;
  real <lower = 0> beta_gene_sigma;
  vector [N_gene] alpha_raw [N_sample];
  vector [N_gene] beta_raw [N_sample];
  vector [N_gene] alpha_gene_raw;
  vector [N_gene] beta_gene_raw;
}

transformed parameters {
  vector [N_gene] alpha [N_sample];
  vector [N_gene] beta [N_sample];
  vector [N_gene] alpha_gene;
  vector [N_gene] beta_gene;
  vector <lower=0> [N_gene] a [N_sample];
  
  // non-centered params (at repertoire level)
  alpha_gene = alpha_grand + alpha_gene_sigma * alpha_gene_raw;
  beta_gene = beta_grand + beta_gene_sigma * beta_gene_raw;
  
  for(i in 1:N_sample) {
    beta[i] = beta_gene + beta_sigma * beta_raw[i];
    alpha[i] = alpha_gene + alpha_sigma * alpha_raw[i];
    a[i] = phi * softmax(alpha_gene + beta[i]*X[i]);
  }
}

model {
  
  for(i in 1:N_sample) {
    Y[, i] ~ multinomial(theta[i]);
    theta[i] ~ dirichlet(a[i]);
  }
  
  // priors
  alpha_grand ~ normal(0, 10);
  beta_grand ~ normal(0, 5);

  alpha_sigma ~ cauchy(0, 1);
  beta_sigma ~ cauchy(0, 1);
  
  alpha_gene_sigma ~ cauchy(0, 1);
  beta_gene_sigma ~ cauchy(0, 1);

  for(i in 1:N_sample) {
    alpha_raw[i] ~ normal(0, 1);
    beta_raw[i] ~ normal(0, 1);
  }
  alpha_gene_raw ~ normal(0, 1);
  beta_gene_raw ~ normal(0, 1);

  phi ~ exponential(tau); //pareto 2
  tau ~ gamma(3, 0.1);
}

