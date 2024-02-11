\name{IGHV_HCV}
\alias{IGHV_HCV}
\docType{data}
\title{IGHV gene usage in HCV+ and healthy individuals}

\description{
Publicly available dataset of IGHV segment usage in memory B-cells of 22 HCV+
individuals and 7 healthy donors.
}

\usage{data("IGHV_HCV")}

\format{
  A data frame with 4 columns:  "individual_id", "condition", "gene_name" and
  "gene_usage_count". The format of the data is suitible to be used as input in
  IgGeneUsage
}
\source{
Tucci, Felicia A., et al. "Biased IGH VDJ gene repertoire and clonal
expansions in B cells of chronically hepatitis C virus–infected individuals."
Blood 131.5 (2018): 546-557.
}
\examples{
data(IGHV_HCV)
head(IGHV_HCV)
}
\keyword{IGHV_HCV}
