\name{LOO}



\alias{LOO}




\title{
Leave-one-out analysis for quantitative evaluation of the probability of DGU
}



\description{
IgGeneUsage detects differential gene usage in immune repertoires that belong
to two biological conditions with its function DGU. To assert quantitatively
the robustness of the estimated probability of DGU (pi), IgGeneUsage has a 
built-in procedure for a fully Bayesian leave-one-out (LOO) analysis. During 
each step of LOO, we discard the data of one of the repertoires, and use the 
remaining data to analyze for DGU with IgGeneUsage. In each step we recorded 
pi for all genes. Therefore, by evaluating the variability of pi for a given 
gene, we can we assert quantitatively its robustness.

Notice, however, that for datasets that include many repertoires (e.g. 100) LOO 
can be computationally costly.}



\usage{
LOO(usage.data, mcmc.warmup, mcmc.steps,
    mcmc.chains, mcmc.cores, hdi.level,
    adapt.delta, max.treedepth)
}


\arguments{
\item{usage.data}{Data.frame with 4 columns: 'sample_id' = character identifier
of each repertoire, 'condition' = character key representing each of the two
biological conditions, 'gene_name' = character name of each gene to be tested
for differential usage, 'gene_usage_count' = number of rearrangements belonging
to a specific sample_id x condition x gene_name. Alternatively, usage.data can
be a SummarizedExperiment object. See examplary data 'data(Ig_SE)' for more 
information.}
\item{mcmc.chains, mcmc.warmup, mcmc.steps, mcmc.cores}{Number of MCMC chains
(default = 4), number of cores to use (default = 1), length of MCMC chains
(default = 1,500), length of adaptive part of MCMC chains (default = 500).}
\item{hdi.level}{Highest density interval (HDI) (default = 0.95).}
\item{adapt.delta}{MCMC setting (default = 0.95).}
\item{max.treedepth}{MCMC setting (default = 12).}
}






\details{
IgGeneUsage invokes the function DGU in each LOO step. For more details
see help for DGU or vignette 'User Manual: IgGeneUsage'.
}

\value{
\item{loo.summary}{differential gene usage statistics for each gene of a given 
LOO step. 1) es = effect size on differential gene usage (mean, median standard 
error (se), standard deviation (sd), L (low boundary of HDI), H (high boundary 
of HDI); 2) contrast = direction of the effect; 3) pmax = probability of 
differential gene usage; 4) loo.id (LOO step ID); 5 Neff (effective sample 
size), Rhat (potential scale reduction factor)}
}

\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}


\examples{
# input data:
# data(Ig)
# head(Ig)

# Alternative:
# use SummarizedExperiment input data
# data(Ig_SE)


# run leave-one-out (LOO)
# L <- LOO(usage.data = Ig,
#          mcmc.warmup = 250,
#          mcmc.steps = 1000,
#          mcmc.chains = 2,
#          mcmc.cores = 1,
#          hdi.level = 0.95,
#          adapt.delta = 0.95,
#          max.treedepth = 13)
}


\seealso{DGU, Ig, IGHV_Epitopes, IGHV_HCV, Ig_SE}
