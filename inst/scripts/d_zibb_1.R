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
      if(bernoulli_rng(0.05) == 1) {
        ysim[g,s] = zibb_rng(N, a, b, 1);
      } else {
      ysim[g,s] = zibb_rng(N, a, b, 0);
      }
    }
  }
}
" 

# compile model
m <- rstan::stan_model(model_code = sim_stan)

# generate data based on following fixed parameters
set.seed(123456)
N_rep <- 10
N_gene <- 30
Y_max <- 10^3
as <- rnorm(n = N_gene, mean = -5, sd = 2)
zs <- runif(n = N_gene, min = 0, max = 0)
d <- list(N_rep = N_rep, 
          N_gene = N_gene, 
          N = Y_max,
          as = as,
          zs = zs,
          sigma_a = 1,
          phi = 100)
# zs[order(as, decreasing = T)[1:5]] <- runif(n = 5, min = 0, max = 0.3)

# simulate
sim <- rstan::sampling(object = m,
                       data = d, 
                       iter = 1, 
                       chains = 1, 
                       algorithm="Fixed_param",
                       seed = 123456)

# extract simulation and convert into data frame which can 
# be used as input of IgGeneUsage
ysim <- rstan::extract(object = sim, par = "ysim")$ysim
ysim <- ysim[1, ,]

ysim_df <- reshape2::melt(ysim)
colnames(ysim_df) <- c("gene_name", "individual_id", "gene_usage_count")
ysim_df$condition <- "C1"
ysim_df <- ysim_df[, c("individual_id", "condition", "gene_name", "gene_usage_count")]
ysim_df$individual_id <- paste0("pt_", as.character(ysim_df$individual_id))
ysim_df$gene_name <- paste0("gene_", as.character(ysim_df$gene_name))
d_zibb_1 <- ysim_df

# save
save(d_zibb_1, file = "data/d_zibb_1.RData", compress T)
