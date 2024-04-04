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
  int<lower=0> N_sample;                   // number of repertoires
  int<lower=0> N_gene;                     // gene
  int<lower=0> N_individual;               // number of individuals
  int<lower=0> N_condition;                // number of conditions
  int N;                                   // repertoire size
  array [N_individual] int condition_id;   // id of conditions
  array [N_sample] int individual_id;      // id of replicate
  real <lower=0> phi;
  real <lower=0, upper=1> kappa;
  vector <lower=0> [N_condition] sigma_condition;
  vector <lower=0> [N_condition] sigma_individual;
  real <lower=0> sigma_replicate;
}

generated quantities {
  vector [N_gene] alpha;
  array [N_sample] vector <lower=0, upper=1> [N_gene] theta;
  array [N_sample] vector [N_gene] beta_sample;
  array [N_individual] vector [N_gene] beta_individual;
  array [N_condition] vector [N_gene] beta_condition;
  // generate usage
  array [N_gene, N_sample] int Y;
  
  for(i in 1:N_gene) {
    alpha[i] = normal_rng(-3, 2);
  }
  
  for(i in 1:N_condition) {
    for(j in 1:N_gene) {
      beta_condition[i][j] = normal_rng(0, sigma_condition[i]);
    }
  }
  
  for(i in 1:N_individual) {
    for(j in 1:N_gene) {
      beta_individual[i][j] = normal_rng(beta_condition[condition_id[i]][j], sigma_individual[condition_id[i]]);
    }
  }
  
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      beta_sample[i][j]  = normal_rng(beta_individual[individual_id[i]][j], sigma_replicate);
      theta[i][j] = inv_logit(alpha[j] + beta_sample[i][j]);
      Y[j, i] = zibb_rng(N, theta[i][j], phi, kappa);
    }
  }
}
"

# compile model
m <- rstan::stan_model(model_code = sim_stan)


# generate data based on the following parameters parameters
set.seed(1021)
N_gene <- 10
N_replicates <- 4
N_condition <- 3
N_individual_per_condition <- 7
N_individual <- N_individual_per_condition * N_condition
N_sample <- N_individual * N_replicates
condition_id <- rep(1:N_condition, each = N_individual_per_condition)
individual_id <- rep(1:N_individual, each = N_replicates)

N <- 1000
phi <- 200
kappa <- 0.02
sigma_condition <- runif(n = N_condition, min = 0.5, max = 0.75)
sigma_individual <- runif(n = N_condition, min = 0.25, max = 0.4)
sigma_replicate <- 0.2


l <- list(N_sample = N_sample, 
          N_gene = N_gene, 
          N_individual = N_individual,
          N_condition = N_condition,
          N = N,
          condition_id = condition_id,
          individual_id = individual_id,
          phi = phi,
          kappa = kappa,
          sigma_condition = sigma_condition,
          sigma_individual = sigma_individual,
          sigma_replicate = sigma_replicate)


# simulate
sim <- rstan::sampling(object = m,
                       data = l, 
                       iter = 1, 
                       chains = 1, 
                       algorithm="Fixed_param")

# extract simulation and convert into data frame which can 
# be used as input of IgGeneUsage
ysim <- rstan::extract(object = sim, par = "Y")$Y
ysim <- ysim[1,,]

ysim_df <- reshape2::melt(ysim)
colnames(ysim_df) <- c("gene_name", "sample_id", "gene_usage_count")

m <- data.frame(sample_id = 1:l$N_sample,
                individual_id = l$individual_id)

ysim_df <- merge(x = ysim_df, y = m, by = "sample_id", all.x = TRUE)

m <- data.frame(individual_id = 1:l$N_individual,
                condition_id = l$condition_id)

ysim_df <- merge(x = ysim_df, y = m, by = "individual_id", all.x = TRUE)


ysim_df$condition <- paste0("C_", ysim_df$condition_id)
ysim_df$gene_name <- paste0("G_", ysim_df$gene_name)
ysim_df$individual_id <- paste0("I_", ysim_df$individual_id)

ysim_df$replicate <- rep(rep(x = paste0("R_", 1:N_replicates), 
                             each = N_gene), times = N_individual)
ysim_df$condition_id <- NULL
ysim_df$sample_id <- NULL
ysim_df <- ysim_df[, c("individual_id", "condition", "gene_name", 
                       "replicate", "gene_usage_count")]

d_zibb_4 <- ysim_df

# save
save(d_zibb_4, file = "data/d_zibb_4.RData", compress = TRUE)

