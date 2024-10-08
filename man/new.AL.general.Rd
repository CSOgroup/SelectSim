% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectX_create.r
\name{new.AL.general}
\alias{new.AL.general}
\title{Create an AL object}
\usage{
new.AL.general(
  am,
  feat.covariates = NULL,
  sample.covariates = NULL,
  min.freq,
  verbose = FALSE
)
}
\arguments{
\item{am}{The binary alteration matrix with row as features and column as samples.}

\item{feat.covariates}{gene/feature covariate.}

\item{sample.covariates}{sample covariate.}

\item{min.freq}{minimum frequency of genes to be mutated.}

\item{verbose}{print the time and each steps.}
}
\value{
An Alteration Landscape (AL) object with the gam.
}
\description{
Create an Alteration Landscape (AL) object which contains gams and mutation burden of samples of associated gams.
}
