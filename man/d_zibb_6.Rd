\name{d_zibb_6}
\alias{d_zibb_6}
\docType{data}
\title{Simulated Ig gene usage data}

\description{
A small example of paired-sample *with replicates* IRRs with these features:

  \itemize{
    \item 3 conditions
    \item 6 individuals with one IRRs per condition
    \item 10 Ig genes
    \item 4 replicates
  }
This dataset was simulated from zero-inflated beta-binomial (ZIBB) 
distribution. Simulation code is available in inst/scripts/d_zibb_6.R
}

\usage{
data("d_zibb_6", package = "IgGeneUsage")
}

\format{
A data frame with 4 columns:  
\itemize{
  \item "individual_id"
  \item "condition"
  \item "gene_name"
  \item "gene_name_count"
}
This format is accepted by IgGeneUsage.
}
\source{
Simulation code is provided in inst/scripts/d_zibb_6.R
}
\examples{
data("d_zibb_6", package = "IgGeneUsage")
head(d_zibb_6)
}
\keyword{d_zibb_6}