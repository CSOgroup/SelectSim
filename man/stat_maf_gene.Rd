% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gam_utils.r
\name{stat_maf_gene}
\alias{stat_maf_gene}
\title{Summary functions for MAF file}
\usage{
stat_maf_gene(maf, column = "Hugo_Symbol", ...)
}
\arguments{
\item{maf}{a maf as dataframe}

\item{column}{a schema of datafrane check Select::TCGA_maf_schema for example}

\item{...}{Other options}
}
\value{
filtered_maf  a filtered maf file
}
\description{
\code{stat_maf_column()} takes a maf file and filters a MAF dataframe by retaining only the rows with a column value included in the values list
}
