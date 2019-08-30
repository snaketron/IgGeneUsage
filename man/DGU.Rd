\name{DGU}



\alias{DGU}




\title{
Differential gene usage in immune repertoires
}



\description{
IgGeneUsage detects differential gene usage in immune repertoires that belong
to two biological conditions.}



\usage{
DGU(usage.data, mcmc.warmup, mcmc.steps,
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
The input to IgGeneUsage is a table with usage frequencies for each gene of a
repertoire that belongs to a particular biological condition. For the analysis 
of differential gene usage between two biological conditions, IgGeneUsage 
employs a Bayesian hierarchical model for zero-inflated beta-binomial (ZIBB)
regression (see vignette 'User Manual: IgGeneUsage').
}

\value{
\item{glm.summary}{differential gene usage statistics for each gene.
1) es = effect size on differential gene usage (mean, median standard error 
(se), standard deviation (sd), L (low boundary of HDI), H (high boundary of 
HDI); 2) contrast = direction of the effect 3) pmax = probability of 
differential gene usage}
\item{test.summary}{differential gene usage statistics computed with the 
Welch's t-test (columns start with 't'), and Wilcoxon signed-rank test (columns 
start with 'u'). For both test report P-values, FDR-corrected P-values, 
Bonferroni-corrected P-values. Additionally, we report t-value and 95\% CI 
(from the t-test) and U-value (from the Wilcoxon signed-rank test).}
\item{glm}{stanfit object}
\item{ppc.data}{two types of posterior predictive checks: 1) repertoire-
specific, 2) gene-specific}
}

\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}


\examples{
# input data
data(Ig)
head(Ig)

# Alternative:
# use SummarizedExperiment input data
# data(Ig_SE)


# run differential gene usage (DGU)
M <- DGU(usage.data = Ig,
         mcmc.warmup = 500,
         mcmc.steps = 1500,
         mcmc.chains = 2,
         mcmc.cores = 2,
         hdi.level = 0.95,
         adapt.delta = 0.95,
         max.treedepth = 13)
}


\seealso{Ig, IGHV_Epitopes, IGHV_HCV, Ig_SE}
