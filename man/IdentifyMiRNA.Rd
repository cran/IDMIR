% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IdentifyMiRNA.R
\name{IdentifyMiRNA}
\alias{IdentifyMiRNA}
\title{IdentifyMiRNA}
\usage{
IdentifyMiRNA(GDEscore.table,nperm=1000,damping=0.90)
}
\arguments{
\item{GDEscore.table}{A matrix with one column of GDEscore.}

\item{nperm}{The Number of random permutations (default: 100).}

\item{damping}{Restart the probability of the random-walk algorithm (default: 0.9).}
}
\value{
A data frame with seven columns those are "MiRNA", "Target", "Number" (number of targets), "Score" (Centrality score), "P-value", and "FDR".
}
\description{
The function "IdentifyMiRNA" is used to identify significantly dysregulated miRNAs by calculating the eigenvector centrality of miRNAs.
}
\examples{
# Obtain the example data
GEP<-GetData_Mirna("GEP")
label<-GetData_Mirna("label")
# Run the function
GDEscore<-GetGDEscore(GEP,label)
\donttest{MiRNAScore<-IdentifyMiRNA(GDEscore,nperm=5, damping=0.90)}
}
