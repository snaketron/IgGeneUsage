\name{DGU}

\alias{DGU}

\title{
Differential immunoglobulin (Ig) gene usage in immune repertoires
}

\description{
IgGeneUsage detects differential gene usage (DGU) in immune repertoires that 
belong to two biological conditions.}

\usage{
DGU(ud, 
    mcmc_warmup, 
    mcmc_steps,
    mcmc_chains, 
    mcmc_cores, 
    hdi_lvl,
    adapt_delta, 
    max_treedepth)
}

\arguments{
\item{ud}{Data.frame with 4 columns: 
\itemize{
\item 'sample_id' = character, repertoire name (e.g. R1)
\item 'condition' = character, name of biological conditions (e.g. tumor)
\item 'gene_name' = character, Ig gene name (e.g. IGHV1-69)
\item 'gene_usage_count' = number, frequency (=usage) of rearrangements from 
       sample_id x condition x gene_name
}
ud can also be be a SummarizedExperiment object. See examplary data 
'data(Ig_SE)' for more information.
}
\item{mcmc_chains, mcmc_warmup, mcmc_steps, mcmc_cores}{Number of MCMC chains
(default = 4), number of cores to use (default = 1), length of MCMC chains
(default = 1,500), length of adaptive part of MCMC chains (default = 500).}
\item{hdi_lvl}{Highest density interval (HDI) (default = 0.95).}
\item{adapt_delta}{MCMC setting (default = 0.95).}
\item{max_treedepth}{MCMC setting (default = 12).}
}

\details{
The main input of IgGeneUsage is a table with Ig gene usage frequencies for 
a set of repertoires that belong to one of two biological condition. For the 
DGU analysis between two biological conditions, IgGeneUsage employs a Bayesian 
hierarchical model for zero-inflated beta-binomial (ZIBB) regression (see 
vignette 'User Manual: IgGeneUsage').
}

\value{
\item{glm_summary}{DGU statistics for each gene: 1) es = effect size on 
DGU (mean, median standard error (se), standard deviation (sd), L (low boundary 
of HDI), H (high boundary of HDI); 2) contrast = direction of the effect; 3) 
pmax = probability of DGU}
\item{test_summary}{DGU statistics computed with the Welch's t-test (columns 
start with 't'), and Wilcoxon signed-rank test (columns start with 'u'). For 
both test report P-values, FDR-corrected P-values, Bonferroni-corrected 
P-values. Additionally, we report t-value and 95\% CI (from the t-test) and 
U-value (from the Wilcoxon signed-rank test).}
\item{glm}{stanfit object}
\item{ppc}{two types of posterior predictive checks: 1) repertoire-
specific, 2) gene-specific}
\item{ud}{processed gene usage data used for the model}
}


\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}


\examples{
# input data
data(d_zibb)
head(d_zibb)

# run differential gene usage (DGU)
M <- DGU(ud = d_zibb,
         mcmc_warmup = 350,
         mcmc_steps = 1500,
         mcmc_chains = 2,
         mcmc_cores = 1,
         hdi_lvl = 0.95,
         adapt_delta = 0.8,
         max_treedepth = 10)

# look at DGU results
head(M$glm_summary)

# look at DGU results (by frequentist methods)
head(M$test_summary)

# look at posterior predictive checks (PPC)
head(M$ppc)
}

\seealso{LOO, GU, d_zibb, Ig, IGHV_Epitopes, IGHV_HCV, Ig_SE}
