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
LOO(ud, 
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
\item 'individual_id' = character, name of the donor (e.g. Pt1)
\item 'condition' = character, name of biological conditions (e.g. tumor)
\item 'gene_name' = character, Ig gene name (e.g. IGHV1-69)
\item 'gene_usage_count' = number, frequency (=usage) of rearrangements from 
       individual_id x condition x gene_name
\item [optional] 'replicate' = character or number. Replicate id, if more than 
      one repertoire (biological replicates) is available per individual 
}
ud can also be be a SummarizedExperiment object. See dataset
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
IgGeneUsage invokes the function DGU in each LOO step. For more details
see help for DGU or vignette 'User Manual: IgGeneUsage'.
}

\value{
\item{loo}{DGU statistics for each Ig gene for specific LOO step:
\itemize{
\item es = effect size statistics: mean, median, standard error (se), 
standard deviation (sd), L/H (low/high boundary of HDI)
\item contrast = direction of the effect
\item pmax = DGU probability
\item loo_id (LOO step)
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
L <- LOO(ud = Ig,
         mcmc_warmup = 500,
         mcmc_steps = 2000,
         mcmc_chains = 3,
         mcmc_cores = 1,
         hdi_lvl = 0.95,
         adapt_delta = 0.99,
         max_treedepth = 10)

# how many LOOs?  
names(L)

# elements in first LOO, see vignette about how to extract results
names(L[[1]])
}


\seealso{DGU, Ig, IGHV_Epitopes, IGHV_HCV, Ig_SE, d_zibb_1, d_zibb_2, d_zibb_3}
