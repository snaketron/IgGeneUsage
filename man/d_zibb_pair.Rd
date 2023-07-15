\name{d_zibb_pair}
\alias{d_zibb_pair}
\docType{data}
\title{Simulated Ig gene usage data}

\description{
A small example dataset that has the following features:

  \itemize{
    \item 5 repertoires (samples) sampled at two timepoints (paired data)
    \item 20 Ig genes
  }
This dataset was simulated from zero-inflated beta-binomial (ZIBB) 
distribution. Simulation code is available in inst/scripts/d_zibb_pair.R
}

\usage{
data("d_zibb_pair", package = "IgGeneUsage")
}

\format{
A data frame with 4 columns:  
\itemize{
  \item "sample_id"
  \item "gene_name"
  \item "gene_name_count_1"
  \item "gene_name_count_2"
}
This format is accepted by IgGeneUsage.
}
\source{
Simulation code is provided in inst/scripts/d_zibb_pair.R
}
\examples{
data("d_zibb_pair", package = "IgGeneUsage")
head(d_zibb_pair)
}
\keyword{d_zibb_pair}
