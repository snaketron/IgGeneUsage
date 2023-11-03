set.seed(seed = 12345)
require(rstan)

# Stan generative model
sim_stan <- "
functions {
  int zibb_rng(int trials, real a, real b, real zi) {
    if (bernoulli_rng(zi) == 1) {
      return (0);
    } else {
      return (beta_binomial_rng(trials, a, b));
    }
  }
}

data {
  int<lower=0> N_rep;
  int<lower=0> N_gene;
  int<lower=0> N;
  vector [N_gene] as;
  vector [N_gene] zs;
  real sigma_a;
  real phi;
}

generated quantities {
  real a_sim [N_gene, N_rep];
  int ysim [N_gene, N_rep];
  real a;
  real b;
  
  for(s in 1:N_rep) {
    for(g in 1:N_gene) {
      a_sim[g,s] = normal_rng(as[g], sigma_a);
      a = inv_logit(a_sim[g,s]) * phi;
      b = phi - a;
      ysim[g,s] = zibb_rng(N, a, b, zs[g]);
    }
  }
}
" 

# compile model
m <- rstan::stan_model(model_code = sim_stan)

# generate data based on following fixed parameters
N_rep <- 5
N_gene <- 15
Y_max <- 10^3
as <- rnorm(n = N_gene, mean = -2, sd = 2)
zs <- c(runif(n = N_gene-3, min = 0, max = 0),
        runif(n = 3, min = 0, max = 0.05))
d <- list(N_rep = N_rep, 
          N_gene = N_gene, 
          N = Y_max,
          as = as,
          zs = zs,
          sigma_a = 0.1,
          phi = 50)

# simulate
sim <- rstan::sampling(object = m,
                       data = d, 
                       iter = 1, 
                       chains = 1, 
                       algorithm="Fixed_param")

# extract simulation and convert into data frame which can 
# be used as input of IgGeneUsage
ysim <- rstan::extract(object = sim, par = "ysim")$ysim
ysim <- ysim[1, ,]

ysim_df <- reshape2::melt(ysim)
colnames(ysim_df) <- c("gene_name", "sample_id", "gene_usage_count")
ysim_df$condition <- "C1"
ysim_df <- ysim_df[, c("sample_id", "condition", "gene_name", "gene_usage_count")]
ysim_df$sample_id <- paste0("S", as.character(ysim_df$sample_id))
ysim_df$gene_name <- paste0("G", as.character(ysim_df$gene_name))
d_zibb_1 <- ysim_df

# save
save(d_zibb_1, file = "data/d_zibb_1.RData")