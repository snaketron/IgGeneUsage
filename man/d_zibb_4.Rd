\name{d_zibb_4}
\alias{d_zibb_4}
\docType{data}
\title{Simulated Ig gene usage data}

\description{
A small example dataset that has the following features:

  \itemize{
    \item 2 conditions
    \item 7 individuals per condition
    \item 4 replicates per individual
    \item 8 Ig genes
  }
This dataset was simulated from zero-inflated beta-binomial (ZIBB) 
distribution. Simulation code is available in inst/scripts/d_zibb_4.R
}

\usage{
data("d_zibb_4", package = "IgGeneUsage")
}

\format{
A data frame with columns:  
\itemize{
  \item "individual_id"
  \item "condition"
  \item "gene_name"
  \item "replicate"
  \item "gene_name_count"
}
This format is accepted by IgGeneUsage.
}
\source{
Simulation code is provided in inst/scripts/d_zibb_4.R
}
\examples{
data("d_zibb_4", package = "IgGeneUsage")
head(d_zibb_4)
}
\keyword{d_zibb_4}
