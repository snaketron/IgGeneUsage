\name{d_zibb}
\alias{d_zibb}
\docType{data}
\title{Simulated Ig gene usage data}

\description{
A small example dataset that has the following features:

  \itemize{
    \item 2 conditions
    \item 10 replicates (samples) per condition
    \item 15 Ig genes
  }
This dataset was simulated from zero-inflated beta-binomial (ZIBB) 
distribution. Simulation code is available in inst/scripts/d_zibb.R
}

\usage{
data("d_zibb", package = "IgGeneUsage")
}

\format{
A data frame with 4 columns:  
\itemize{
  \item "sample_id"
  \item "condition"
  \item "gene_name"
  \item "gene_name_count"
}
This format is accepted by IgGeneUsage.
}
\source{
Simulation code is provided in inst/scripts/d_zibb.R
}
\examples{
data("d_zibb", package = "IgGeneUsage")
head(d_zibb)
}
\keyword{d_zibb}
