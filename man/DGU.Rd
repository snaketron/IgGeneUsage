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
to a specific sample_id x condition x gene_name.}
\item{mcmc.chains, mcmc.warmup, mcmc.steps, mcmc.cores}{Number of MCMC chains
(default = 4), number of cores to use (default = 1), length of MCMC chains
(default = 1,500), length of adaptive part of MCMC chains (default = 500).}
\item{hdi.level}{Highest density interval (HDI) (default = 0.95).}
\item{adapt.delta}{MCMC setting (default = 0.95).}
\item{max_treedepth}{MCMC setting (default = 12).}
}






\details{
The input to IgGeneUsage is a table with usage frequencies for each gene of a
repertoire that belongs to a particular biological condition. Given R
repertoires, each having G genes, IgGeneUsage generates a gene usage matrix
Y^(R x G). Row sums in Y define the total usage in each repertoire (N). The
design variable X is set to X = 1 for repertoires that belong to the first
condition, and X = -1 otherwise. For the evaluation of differential gene usage
between two biological conditions, IgGeneUsage employs a Bayesian hierarchical
model for zero-inflated beta-binomial (ZIBB) regression (see vignette
'ModelOverview').
}

\value{
\item{glm}{stanfit object}
\item{glm.summary}{differential gene usage statistics for each gene.
1) effect = effect on differential gene usage (mean, median standard error (se),
standard deviation (sd), L (low boundary of HDI), H (high boundary of HDI);
2) contrast = direction of the effect
3) pmax = probability of differential gene usage}
\item{ppc.data}{two types of posterior predictive checks: 1) repertoire-
specific, 2) gene-specific}
\item{test.stats}{differential gene usage statistics computed with the Welch's
t-test (columns start with 't'), and Wilcoxon signed-rank test (columns start
with 'u'). For both test report P-values, FDR-corrected P-values, Bonferroni-
corrected P-values. Additionally, we report t-value and 95\% CI (from the
t-test) and U-value (from the Wilcoxon signed-rank test).}
}

\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}


\examples{
# input data
data(Ig)
head(data(Ig))


# run differential gene usage (DGU)
M <- DGU(usage.data = Ig,
         mcmc.warmup = 500,
         mcmc.steps = 1500,
         mcmc.chains = 4,
         mcmc.cores = 4,
         hdi.level = 0.95,
         adapt.delta = 0.95,
         max.treedepth = 13)
}


\seealso{
}
