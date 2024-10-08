% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gam_utils.r
\name{maf2gam}
\alias{maf2gam}
\title{Generate gam from the maf file}
\usage{
maf2gam(
  maf,
  sample.col = "Tumor_Sample_Barcode",
  gene.col = "Hugo_Symbol",
  value.var = "HGVSp_Short",
  samples = NULL,
  genes = NULL,
  fun.aggregate = length,
  binarize = TRUE,
  fill = NA,
  ...
)
}
\arguments{
\item{maf}{a maf as dataframe}

\item{sample.col}{a list containing the elements to filter}

\item{gene.col}{column in maf file to filter}

\item{value.var}{column in maf file to filter}

\item{samples}{column in maf file to filter}

\item{genes}{column in maf file to filter}

\item{fun.aggregate}{column in maf file to filter}

\item{binarize}{a boolena to include or exclude the dataframe with values in list provided}

\item{fill}{a grep argument to specify if grep use the argumnet as string or not}

\item{...}{Other options}
}
\value{
filtered_maf  a GAM matrix
}
\description{
\code{maf2gam()} takes a maf file and converts into gam
}
