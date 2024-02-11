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
  vector [N_gene] bs;
  vector [N_gene] zs;
  real sigma_a;
  real sigma_b;
  real phi;
}

generated quantities {
  real a_sim [N_gene];
  real b_k_sim [N_gene, 3];
  real b_sim [N_gene, N_rep, 3];
  int ysim [N_gene, N_rep, 3];
  real a;
  real b;
  
  for(g in 1:N_gene) {
    a_sim[g] = normal_rng(as[g], sigma_a);
    for(k in 1:3) {
      b_k_sim[g,k] = normal_rng(bs[g], sigma_b);
      for(s in 1:N_rep) {
        b_sim[g,s,k] = normal_rng(b_k_sim[g,k], 0.1);
        a = inv_logit(a_sim[g] + b_sim[g,s,k]) * phi;
        b = phi - a;
        ysim[g,s,k] = zibb_rng(N, a, b, zs[g]);
      }
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
as <- rnorm(n = N_gene, mean = -2, sd = 1.5)
bs <- rnorm(n = N_gene, mean = 0, sd = 0.35)
zs <- c(runif(n = N_gene-3, min = 0, max = 0),
        runif(n = 3, min = 0, max = 0.05))
d <- list(N_rep = N_rep, 
          N_gene = N_gene, 
          N = Y_max,
          as = as,
          bs = bs,
          zs = zs,
          sigma_a = 0.1,
          sigma_b = 0.2,
          phi = 20)

# simulate
sim <- rstan::sampling(object = m,
                       data = d, 
                       iter = 1, 
                       chains = 1, 
                       algorithm="Fixed_param")

# extract simulation and convert into data frame which can 
# be used as input of IgGeneUsage
ysim <- rstan::extract(object = sim, par = "ysim")$ysim
ysim <- ysim[1,,,]

ysim_df <- reshape2::melt(ysim)
colnames(ysim_df) <- c("gene_name", "individual_id", "condition", "gene_usage_count")
ysim_df$condition <- paste0("C", ysim_df$condition)
ysim_df <- ysim_df[, c("individual_id", "condition", "gene_name", "gene_usage_count")]
ysim_df$individual_id <- paste0("pt_", as.character(ysim_df$individual_id))
ysim_df$gene_name <- paste0("gene_", as.character(ysim_df$gene_name))
d_zibb_3 <- ysim_df

# save
save(d_zibb_3, file = "data/d_zibb_3.RData", compress = T)
