\name{diffUsage}



\alias{diffUsage}




\title{
Differential gene usage of immune repertoires
}



\description{
Given a set of repertoires, grouped into two biological conditions, IgGeneUsage
quantifies the degree of differential usage for each gene between the two
conditions.}



\usage{
diffUsage(usage.data, mcmc.warmup, mcmc.steps,
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
(default = 1,000), length of adaptive part of MCMC chains (default = 500).}
\item{hdi.level}{Highest density interval (HDI) (default = 0.95).}
\item{...}{Optional parameters include adapt_delta: STAN configuration
(default = 0.9); max_treedepth: STAN configuration (default = 10)}
}






\details{
Given are $R$ repertoires, distributed across two biological conditions, which
are denoted by a design variable $X$ with $X = 1$ for repertoires belonging to
the first condition, and $X = -1$ otherwise. Also given is the number of all
genes $G$. IgGeneUsage generates a gene usage matrix $Y^{G x R}$ where each
cell in $Y$ contains a gene usage frequency, i.e. number of rearrangements that
belong to particular gene and repertoire. Finally, $N$ is a vector containing
the total usage (number of all rearrangements) in each repertoire. For the
evaluation of differential gene usage between two biological conditions,
IgGeneUsage uses a Bayesian hierarchical model for zero-inflated beta-binomial
(ZIBB) regression (see vignette 'ModelOverview').
}

\value{
\item{glm}{stanfit object}
\item{glm.summary}{differential gene usage statistics for each gene.
1) effect = effect on differential gene usage (se = standard error, L = low
boundary of HDI, H = high boundary of HDI); 2) pmax = probability of
differential gene usage}
\item{ppc.data}{two types of posterior predictive checks: 1) for each gene of
a particular repertoire, 2) for each gene and condition}
\item{test.stats}{differential gene usage statistics computed with the Welch's
t-test (columns start with 't'), and Wilcoxon rank sum test (columns start with
'u'). For both test we report a P-value, FDR-corrected P-value and Bonferroni-
corrected P-value. Additionally, we report t-value and 95\% CI (from the
t-test) and U-value (from the Wilcoxon rank sum test).}
}

\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}


\examples{
# input data
data(Ig)
head(data(Ig))


# run diffUsage
M <- diffUsage(usage.data = Ig,
               mcmc.warmup = 500,
               mcmc.steps = 1500,
               mcmc.chains = 4,
               mcmc.cores = 4,
               hdi.level = 0.95,
               adapt.delta = 0.95,
               max.treedepth = 13,
               dev.model = "")
}


\seealso{
}
