% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectX_stats.r
\name{interaction.table}
\alias{interaction.table}
\title{Compute weight overlap stats}
\usage{
interaction.table(
  al,
  als,
  obs,
  wobs,
  r.obs = NULL,
  r.wobs = NULL,
  null,
  maxFDR = 0.2,
  n.cores = 1,
  estimate_pairwise = FALSE,
  n.permut = 1000
)
}
\arguments{
\item{al}{The alteration matrix}

\item{als}{The alteration matrix stats object}

\item{obs}{The observed overlap}

\item{wobs}{The weighted observed overlap}

\item{r.obs}{The random observed overlap}

\item{r.wobs}{The random weighted observed overlap}

\item{null}{The null model}

\item{maxFDR}{The maxFDR cutoff}

\item{n.cores}{The number of cores}

\item{estimate_pairwise}{Compute pairwise or not}

\item{n.permut}{The number of permutation}
}
\value{
The table of results
}
\description{
Compute weight overlap stats
}
