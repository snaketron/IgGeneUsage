\name{Ig}
\alias{Ig}
\docType{data}
\title{IGHV gene family usage in vaccine-challenged B-cell repertoires}

\description{
A small example database subset from study evaluating vaccine-induced changes in
B-cell populations publicly provided by R-package alakazam (version 0.2.11). It
contains IGHV gene family usage, reported in four B-cell populations (samples
IgM, IgD, IgG and IgA) across two timepoints (conditions = -1 hour and +7 days).
}

\usage{data("Ig")}

\format{
  A data frame with 4 columns:  "sample_id", "condition", "gene_name" and
  "gene_usage_count". The format of the data is suitible to be used as input in
  IgGeneUsage
}
\source{
R package: alakazam version 0.2.11
}
\references{
Laserson U and Vigneault F, et al. High-resolution antibody dynamics of
vaccine-induced immune responses. Proc Natl Acad Sci USA. 2014 111:4928-33.
}
\examples{
data(Ig)
head(Ig)
}
\keyword{Ig}
