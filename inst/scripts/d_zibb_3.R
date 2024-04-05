set.seed(seed = 12345)
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
  int<lower=0> N_gene;                         // gene
  int<lower=0> N_individual;                   // number of individuals
  int<lower=0> N_condition;                    // number of conditions
  int<lower=0> N;                              // repertoire size
  array [N_individual] int condition_id;       // id of conditions
  real <lower=0> phi;
  real <lower=0, upper=1> kappa;
  vector <lower=0> [N_condition] sigma_individual;
  vector <lower=0> [N_condition] sigma_condition;
}

generated quantities {
  vector [N_gene] alpha;
  array [N_individual] vector <lower=0, upper=1> [N_gene] theta;
  array [N_individual] vector [N_gene] beta_individual;
  array [N_condition] vector [N_gene] beta_condition;
  // generate usage
  array [N_gene, N_individual] int Y;
  
  for(i in 1:N_condition) {
    for(j in 1:N_gene) {
      beta_condition[i][j] = normal_rng(0, sigma_condition[i]);
    }
  }
  
  for(i in 1:N_gene) {
    alpha[i] = normal_rng(-3.0, 1.0);
    
    for(j in 1:N_individual) {
      beta_individual[j][i] = normal_rng(beta_condition[condition_id[j]][i], sigma_individual[condition_id[j]]);
      theta[j][i] = inv_logit(alpha[i] + beta_individual[j][i]);
      Y[i,j] = zibb_rng(N, theta[j][i], phi, kappa);
    }
  }
}
"


# compile model
m <- rstan::stan_model(model_code = sim_stan)

# generate data based on following fixed parameters
set.seed(123456)
N_condition <- 3
N_individual <- 5
N_gene <- 8
N <- 10^3
sigma_individual <- runif(n = N_condition, min = 0.1, max = 0.2)
sigma_condition <- runif(n = N_condition, min = 0.2, max = 0.6)
phi <- 200
kappa <- 0.015

condition_id <- rep(x = 1:N_condition, each = N_individual)

l <- list(N_individual = N_individual*N_condition, 
          N_gene = N_gene, 
          N_condition = N_condition,
          N = N,
          condition_id = condition_id,
          sigma_individual = sigma_individual,
          sigma_condition = sigma_condition,
          phi = phi,
          kappa = kappa)

# simulate
sim <- rstan::sampling(object = m,
                       data = l, 
                       iter = 1, 
                       chains = 1, 
                       algorithm = "Fixed_param",
                       seed = 12346)


# extract simulation and convert into data frame which can 
# be used as input of IgGeneUsage
ysim <- rstan::extract(object = sim, par = "Y")$Y
ysim <- ysim[1,,]

ysim_df <- reshape2::melt(ysim)
colnames(ysim_df) <- c("gene_name", "individual_id", "gene_usage_count")

m <- data.frame(individual_id = 1:l$N_individual, condition = l$condition_id)

ysim_df <- merge(x = ysim_df, y = m, by = "individual_id", all.x = TRUE)

ysim_df$condition <- paste0("C_", ysim_df$condition)
ysim_df$gene_name <- paste0("G_", ysim_df$gene_name)
ysim_df$individual_id <- paste0("I_", ysim_df$individual_id)

d_zibb_3 <- ysim_df

# save
save(d_zibb_3, file = "data/d_zibb_3.RData", compress = TRUE)
