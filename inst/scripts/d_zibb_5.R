require(rstan)

# Stan generative model
sim_stan <- "
functions {
  int zibb_rng(int y, int n, real mu, real phi, real kappa) {
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
  array [N_sample] int N;                  // repertoire size
  array [N_sample] int condition_id;   // id of conditions
  array [N_sample] int individual_id;      // id of replicate
  real <lower=0> phi;
  real <lower=0, upper=1> kappa;
  array [N_condition] vector [N_gene] beta_condition;
  vector <lower=0> [N_condition] sigma_condition;
  real <lower=0> sigma_alpha;
}
generated quantities {
  vector [N_gene] alpha;
  array [N_individual] vector [N_gene] alpha_individual;
  array [N_sample] vector [N_gene] beta_individual;
  // generate usage
  array [N_sample] vector <lower=0, upper=1> [N_gene] theta;
  array [N_gene, N_sample] int Y;
  
  for(i in 1:N_gene) {
    alpha[i] = normal_rng(-3, 0.5);
  }
  
  for(i in 1:N_sample) {
    for(j in 1:N_gene) {
      
      alpha_individual[individual_id[i]][j] = normal_rng(alpha[j], sigma_alpha);
      
      beta_individual[i][j] = normal_rng(beta_condition[condition_id[i]][j], sigma_condition[condition_id[i]]);
      
      theta[i][j] = inv_logit(alpha_individual[individual_id[i]][j] + beta_individual[i][j]);
      Y[j, i] = zibb_rng(Y[j, i], N[i], theta[i][j], phi, kappa);
    }
  }
}
"

# compile model
model <- rstan::stan_model(model_code = sim_stan)



# generate data based on the following parameters parameters
set.seed(11132)
N_gene <- 10
N_individual <- 6
N_condition <- 3
N_sample <- N_individual*N_condition

condition_id <- rep(x = 1:N_condition, each = N_individual) 

N <- rep(x = 1000, times = N_sample)

individual_id <- rep(x = 1:N_individual, times = N_condition)

phi <- 200
kappa <- 0.02

beta_condition <- array(data = 0, dim = c(N_condition, N_gene))
for(c in 1:N_condition) {
  for(g in 1:N_gene) {
    u <- runif(n = 1, min = 0, max = 1)
    if(u <= 0.8) {
      beta_condition[c,g] <- rnorm(n = 1, mean = 0, sd = 0.1)
    } else {
      beta_condition[c,g] <- rnorm(n = 1, mean = 0, sd = 2)
    }
  }
}

sigma_condition <- rep(x = 0.5, times = N_condition)
sigma_alpha <- 0.25


l <- list(N_sample = N_sample, 
          N_gene = N_gene, 
          N_individual = N_individual,
          N_condition = N_condition,
          N = N,
          condition_id = condition_id,
          individual_id = individual_id,
          phi = phi,
          kappa = kappa,
          beta_condition = beta_condition,
          sigma_condition = sigma_condition,
          sigma_alpha = sigma_alpha)

# simulate
sim <- rstan::sampling(object = model,
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
                individual_id = l$individual_id,
                condition_id = l$condition_id)

ysim_df <- merge(x = ysim_df, y = m, by = "sample_id", all.x = TRUE)

ysim_df$condition <- paste0("C_", ysim_df$condition_id)
ysim_df$gene_name <- paste0("G_", ysim_df$gene_name)
ysim_df$individual_id <- paste0("I_", ysim_df$individual_id)

ysim_df$condition_id <- NULL
ysim_df$sample_id <- NULL
ysim_df <- ysim_df[, c("individual_id", "condition", "gene_name",
                       "gene_usage_count")]

d_zibb_5 <- ysim_df

# save
save(d_zibb_5, file = "data/d_zibb_5.RData", compress = TRUE)
