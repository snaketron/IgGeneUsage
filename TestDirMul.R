require(IgGeneUsage)
data("IGHV_HCV")
IGHV_HCV

require(rstan)
source("R/Dgu.R")
source("R/Loo.R")
source("R/Util.R")
source("R/UtilInput.R")

model.dir <- rstan::stan_model(file = "inst/extdata/dirichlet_multinomial.stan")
model.flex <- rstan::stan_model(file = "inst/extdata/zibb_flex.stan")


# int <lower = 0> N_sample; // number of samples (repertoires)
# int <lower = 0> N_gene; // number of genes
# int Y [N_gene, N_sample]; // number of successes (cells) in samples x gene
# int <lower = -1, upper = 1> X[N_sample]; // condition

data("Ig_SE")
l <- getUsageData(usage = IGHV_HCV)
# l <- getUsageData(usage = Ig)

# Ig$gene_usage_count
ggplot(data = Ig)+
  geom_point(aes(x = gene_name, y = gene_usage_count, col = condition))

fit.dm <- rstan::sampling(object = model.dir, 
                          data = l,
                          chains = 4,
                          cores = 4,
                          control = list(adapt_delta = 0.99))


fit.flex <- rstan::sampling(object = model.flex, 
                          data = l,
                          chains = 4,
                          cores = 4,
                          control = list(adapt_delta = 0.99))


summary(fit.dm, par = "beta_gene")$summary
summary(fit.flex, par = "beta_gene")$summary


rstan::check_hmc_diagnostics(object = fit.dm)
rstan::check_hmc_diagnostics(object = fit.flex)
