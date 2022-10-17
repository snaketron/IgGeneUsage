\name{LOO}



\alias{LOO}




\title{
Leave-one-out analysis for quantitative evaluation of the probability of DGU
}



\description{
IgGeneUsage detects differential gene usage (DGU) in immune repertoires that 
belong to two biological conditions.

To quantify the robustness of the estimated probability of DGU (pmax), 
IgGeneUsage has a built-in procedure for a fully Bayesian leave-one-out (LOO) 
analysis. In each LOO step we discard the data of one of the repertoires, and 
use the remaining data to analyze for DGU with IgGeneUsage. In each step we 
record pmax for all genes. Finally, we evaluate the variability of pmax for a 
given across the different LOO steps. Low variability in pmax: robust DGU; 
high variability: unrobust DGU.

For datasets that include many repertoires (e.g. 100) LOO can be 
computationally costly.}



\usage{
LOO(usage.data, 
    mcmc.warmup, 
    mcmc.steps,
    mcmc.chains, 
    mcmc.cores, 
    hdi.level,
    adapt.delta, 
    max.treedepth)
}


\arguments{
\item{usage.data}{Data.frame with 4 columns:
\itemize{
\item 'sample_id' = character, repertoire name (e.g. R1)
\item 'condition' = character, name of biological conditions (e.g. tumor)
\item 'gene_name' = character, Ig gene name (e.g. IGHV1-69)
\item 'gene_usage_count' = number, frequency (=usage) of rearrangements from 
       sample_id x condition x gene_name
}
usage.data can also be be a SummarizedExperiment object. See dataset
'data(Ig_SE)' for more information.
}
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
\item{loo.summary}{DGU statistics for each Ig gene for specific LOO step:
\itemize{
\item es = effect size statistics: mean, median, standard error (se), 
standard deviation (sd), L/H (low/high boundary of HDI)
\item contrast = direction of the effect
\item pmax = DGU probability
\item loo.id (LOO step)
\item Neff (effective sample size), Rhat (potential scale reduction factor)
}
}
}

\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}


\examples{
# input data:
data("Ig", package = "IgGeneUsage")
head(Ig)

# run leave-one-out (LOO)
L <- LOO(usage.data = Ig,
         mcmc.warmup = 500,
         mcmc.steps = 2000,
         mcmc.chains = 3,
         mcmc.cores = 1,
         hdi.level = 0.95,
         adapt.delta = 0.99,
         max.treedepth = 10)
         
head(L$loo.summary)
}


\seealso{DGU, Ig, IGHV_Epitopes, IGHV_HCV, Ig_SE, d_zibb}
