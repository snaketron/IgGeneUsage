\name{Ig_SE}
\alias{Ig_SE}
\docType{data}
\title{IGHV gene family usage in vaccine-challenged B-cell repertoires 
(SummarizedExperiment object)}

\description{
A small example database subset from study evaluating vaccine-induced changes in
B-cell populations publicly provided by R-package alakazam (version 0.2.11). It
contains IGHV gene family usage, reported in four B-cell populations (samples
IgM, IgD, IgG and IgA) across two timepoints (conditions = -1 hour and +7 days).
}

\usage{data("Ig_SE")}

\format{
  A SummarizedExperiment object with 1) assay data (rows = gene name, 
  columns = repertoires) and 2) column data.frame in which the sample names
  and the corresponding biological condition labels are noted.
}
\source{
R package: alakazam version 0.2.11
}
\references{
Laserson U and Vigneault F, et al. High-resolution antibody dynamics of
vaccine-induced immune responses. Proc Natl Acad Sci USA. 2014 111:4928-33.
}
\examples{
# inspect the data
data(Ig_SE)

# repertoire information: must have the two columns: 'condition' and 
'individual_id'

SummarizedExperiment::colData(Ig_SE)

# assay counts (gene frequency usage)
SummarizedExperiment::assay(x = Ig_SE)
}
\keyword{Ig_SE}
