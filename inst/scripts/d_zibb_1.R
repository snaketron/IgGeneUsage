require(rstan)

# Stan generative model
sim_stan <- "
functions {
  int zibb_rng(int n, real mu, real phi, real kappa) {
    if (bernoulli_rng(kappa) == 1) {
      return (0);
    } else {
      return (beta_binomial_rng(n, mu * phi, (1 - mu) * phi));
    }
  }
}

data {
  int<lower=0> N_individual;
  int<lower=0> N_gene;
  int<lower=0> N;
  real sigma;
  real phi;
  real kappa;
}

generated quantities {
  array [N_gene, N_individual] int Y;
  real mu;
  real a;
  
  for(g in 1:N_gene) {
    a = normal_rng(-3, 1);
    for(s in 1:N_individual) {
      mu = normal_rng(a, sigma);
      Y[g,s] = zibb_rng(N, 1/(1+exp(-(mu))), phi, kappa);
    }
  }
}
" 

# compile model
m <- rstan::stan_model(model_code = sim_stan)

# generate data based on following fixed parameters
set.seed(123456)
N_individual <- 5
N_gene <- 15
N <- 10^3
sigma <- 0.5
phi <- 200
kappa <- 0.03

l <- list(N_individual = N_individual, 
          N_gene = N_gene, 
          N = N,
          sigma = sigma,
          phi = phi,
          kappa = kappa)

# simulate
sim <- rstan::sampling(object = m,
                       data = l, 
                       iter = 1, 
                       chains = 1, 
                       algorithm = "Fixed_param",
                       seed = 123456)

# extract simulation and convert to be used as input of IgGeneUsage
ysim <- rstan::extract(object = sim, par = "Y")$Y
ysim <- ysim[1,,]

ysim_df <- reshape2::melt(ysim)
colnames(ysim_df) <- c("gene_name", "individual_id", "gene_usage_count")
ysim_df$condition <- "C_1"
ysim_df <- ysim_df[, c("individual_id", "condition", "gene_name", "gene_usage_count")]
ysim_df$individual_id <- paste0("I_", as.character(ysim_df$individual_id))
ysim_df$gene_name <- paste0("G_", as.character(ysim_df$gene_name))
d_zibb_1 <- ysim_df

# save
save(d_zibb_1, file = "data/d_zibb_1.RData", compress = TRUE)
