\name{CDR3_Epitopes}
\alias{CDR3_Epitopes}
\docType{data}
\title{Net charge usage in CDR3 sequences of T-cell receptor repertoires
disturbed by Influenza-A and CMV}

\description{
Data of CDR3 sequence from human T-cells receptors (TRB-chain) downloaded from
VDJdb. CDR3 sequences annotated to epitopes in Influenza-A and CMV were selected
from different publications, as long as the publication contains at least 100
CDR3 sequences. Each publication is considered as a repertoire (sample).

To compute the net CDR3 sequence charge, we consider the amino acids K, R and
H as +1 charged, while D and E as -1 charged. Thus, we computed the net charge
of a CDR3 sequence by adding up the individual residue charges.
}

\usage{data("CDR3_Epitopes")}

\format{
  A data frame with 4 columns:  "sample_id", "condition", "gene_name" and
  "gene_usage_count". The format of the data is suitible to be used as input in
  IgGeneUsage

  gene_name = net charge group
}
\source{
https://vdjdb.cdr3.net/
}
\examples{
data(CDR3_Epitopes)
head(CDR3_Epitopes)
}
\keyword{CDR3_Epitopes}
