###
# Author  : Marco Mina , Arvind Iyer
# Project : SelectSim
# Desc    : Functions to process the maf to gam
# Version : 0.1.6
###


## General Schemas

###
#' Mutation list object
#'
#' @export
mutation_type <- list(
  "truncating" = c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "Splice_Site", "In_Frame_Ins", "In_Frame_Del"),
  "missense" = c("Missense_Mutation", "Splice_Site"),
  "ignore" = c("Silent", "lincRNA", "IGR", "3'UTR", "5'UTR", "Intron", "5'Flank", "3'Flank", "RNA", "synonymous_variant", "upstream_gene_variant", "intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "downstream_gene_variant", "5_prime_UTR_premature_start_codon_gain_variant")
)

###
#' TCGA_maf_schema: schema for TCGA maf file to process the mutations
#'
#' @export
TCGA_maf_schema <- list(
  "name" = "TCGA_maf",
  "column" = list(
    "gene" = "Hugo_Symbol",
    "gene.name" = "Hugo_Symbol",
    "sample" = "Tumor_Sample_Barcode",
    "sample.name" = "Tumor_Sample_Barcode",
    "mutation.type" = "Variant_Classification",
    "mutation" = "HGVSp_Short"
  ),
  "mutation.type" = list(
    "truncating" = unique(c(mutation_type$truncating, "Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "Splice_Site")),
    "missense" = unique(c(mutation_type$missense, "Missense_Mutation", "Splice_Site", "In_Frame_Ins", "In_Frame_Del")),
    "ignore" = unique(c(mutation_type$ignore))
  )
)

###
#' GENIE_maf_schema: schema for GENIE maf file to process the mutations
#'
#' @export
GENIE_maf_schema <- list(
  "column" = list(
    "gene" = "Hugo_Symbol",
    "gene.name" = "Hugo_Symbol",
    "sample" = "Tumor_Sample_Barcode",
    "sample.name" = "Tumor_Sample_Barcode",
    "mutation.type" = "Variant_Classification",
    "mutation" = "HGVSp_Short"
  ),
  "mutation.type" = list(
    "truncating" = unique(c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "Splice_Site", "In_Frame_Ins", "In_Frame_Del")),
    "missense" = unique(c("Missense_Mutation", "Splice_Site")),
    "ignore" = unique(c(mutation_type$ignore))
  )
)


###
#' Filter maf function
#'
#' @description
#' `filter_maf_column()` takes a maf file and filters a MAF dataframe by retaining only the rows with a column value included in the values list
#'
#' @param maf a maf as dataframe
#' @param values a list containing the elements to filter
#' @param column column in maf file to filter
#' @param inclusive a boolena to include or exclude the dataframe with values in list provided
#' @param fixed a grep argument to specify if grep use the argumnet as string or not
#' @param ... Other options
#' @return filtered_maf  a filtered maf file
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' filter_maf_column(luad_maf, values = "Missense_Mutation",
#'                   column = "Variant_Classification")
#'
#' @export
filter_maf_column <- function(maf, # the maf as dataframe
                              values, # values as list to filter
                              column, # column in the maf file to filter from
                              inclusive = TRUE, # Boolean to decide to include the values provided or not
                              fixed = TRUE, # grep argument to mention if it a regex or string
                              ...) {
  if (!column %in% colnames(maf)) stop(paste(column, "is not a valid column for the specified maf"))
  filtered_maf <- maf

  if (fixed) {
    if (inclusive) filtered_maf <- filtered_maf[which(filtered_maf[, column] %in% values), , drop = FALSE]
    if (!inclusive) filtered_maf <- filtered_maf[which(!filtered_maf[, column] %in% values), , drop = FALSE]
  } else {
    if ((length(unique(values)) * nrow(filtered_maf)) > 100000) warning(paste("Running grep on", length(unique(values)), "unique values vs", nrow(filtered_maf), "maf rows..."))
    if (inclusive) {
      temp <- lapply(unique(values), function(x) {
        filtered_maf[grep(x, filtered_maf[, column], fixed = !fixed), , drop = FALSE]
      })
    }
    if (!inclusive) {
      temp <- lapply(unique(values), function(x) {
        filtered_maf[grep(x, filtered_maf[, column], invert = TRUE, fixed = !fixed), , drop = FALSE]
      })
    }
    temp <- temp[sapply(temp, nrow) > 0]
    filtered_maf <- do.call(rbind, temp)
  }
  filtered_maf <- unique(filtered_maf)
  return(filtered_maf)
}

#' Filter a MAF dataframe by a combination of column values
#'
#' @param maf A MAF dataframe
#' @param values A dataframe of (column, value) pairs to match against
#' @param ... Additional arguments passed to merge
#' @return Filtered MAF dataframe containing only rows matching the value combinations.
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' combos <- data.frame(Hugo_Symbol = "TP53",
#'                      Variant_Classification = "Missense_Mutation")
#' filter_maf_complex(luad_maf, combos,
#'                    by.x = c("Hugo_Symbol","Variant_Classification"),
#'                    by.y = c("Hugo_Symbol","Variant_Classification"))
#'
#' @export
filter_maf_complex <- function(maf, values, ...) {
  filtered_maf <- merge(maf, values, suffixes = c("", "_to.ignore"), ...)
  filtered_maf <- unique(filtered_maf)
  return(filtered_maf)
}

#' Filter a MAF dataframe by sample ID
#'
#' @param maf A MAF dataframe
#' @param samples Vector of sample IDs to retain
#' @param sample.col Column name containing sample IDs
#' @param ... Additional arguments passed to filter_maf_column
#' @return Filtered MAF dataframe.
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' ids <- unique(luad_maf$Tumor_Sample_Barcode)[1:5]
#' filter_maf_sample(luad_maf, samples = ids)
#'
#' @export
filter_maf_sample <- function(maf, samples, sample.col = "Tumor_Sample_Barcode", ...) {
  filtered_maf <- filter_maf_column(maf, values = samples, column = sample.col, ...)
  return(filtered_maf)
}

#' Filter a MAF dataframe by gene name
#'
#' @param maf A MAF dataframe
#' @param genes Vector of gene names to retain
#' @param gene.col Column name containing gene symbols
#' @param ... Additional arguments passed to filter_maf_column
#' @return Filtered MAF dataframe.
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' filter_maf_gene.name(luad_maf, genes = c("TP53", "KRAS"))
#'
#' @export
filter_maf_gene.name <- function(maf, genes, gene.col = "Hugo_Symbol", ...) {
  filtered_maf <- filter_maf_column(maf, values = genes, column = gene.col, ...)
  return(filtered_maf)
}

#' Filter a MAF dataframe by mutation type
#'
#' @param maf A MAF dataframe
#' @param variants Vector of variant classification values to retain
#' @param variant.col Column name containing variant classifications
#' @param ... Additional arguments passed to filter_maf_column
#' @return Filtered MAF dataframe.
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' filter_maf_mutation.type(luad_maf, variants = "Missense_Mutation")
#'
#' @export
filter_maf_mutation.type <- function(maf, variants, variant.col = "Variant_Classification", ...) {
  filtered_maf <- filter_maf_column(maf, values = variants, column = variant.col, ...)
  return(filtered_maf)
}

#' Filter a MAF dataframe by specific gene-mutation combinations
#'
#' @param maf A MAF dataframe
#' @param values Dataframe of allowed (gene, mutation) combinations
#' @param maf.col Columns in maf to join on
#' @param values.col Corresponding columns in values to join on
#' @param ... Additional arguments passed to filter_maf_complex
#' @return Filtered MAF dataframe containing only rows matching the allowed combinations.
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' allowed <- data.frame(Hugo_Symbol  = c("TP53", "KRAS"),
#'                       HGVSp_Short  = c("p.R175H", "p.G12C"))
#' filter_maf_mutations(luad_maf, allowed)
#'
#' @export
filter_maf_mutations <- function(maf, values, maf.col = c("Hugo_Symbol", "HGVSp_Short"), values.col = maf.col, ...) {
  filtered_maf <- filter_maf_complex(maf, values = values, by.x = maf.col, by.y = values.col, ...)
  return(filtered_maf)
}


# Schema functions

###
#' This function filters a MAF dataframe by sample id
#'
#' @description
#' `filter_maf_schema()` takes a maf file and filters a MAF dataframe by retaining only the rows with a column value included in the values list
#'
#' @param maf a maf as dataframe
#' @param schema a schema of datafrane check Select::TCGA_maf_schema for example
#' @param column column in maf file to filter
#' @param values a list containing the elements to file
#' @param ... Other options
#' @return filtered_maf  a filtered maf file
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' filter_maf_schema(luad_maf, TCGA_maf_schema,
#'                   column = "mutation.type",
#'                   values = TCGA_maf_schema$mutation.type$truncating)
#'
#' @export
filter_maf_schema <- function(maf, schema = TCGA_maf_schema, column, values, ...) {
  variant.col <- schema[["column"]][[column]]
  filtered_maf <- filter_maf_column(maf, values = values, column = variant.col, ...)
  return(filtered_maf)
}


###
#' This function filters a MAF dataframe by retaining (or discarding) truncating mutations
#'
#' @description
#' `filter_maf_truncating()` takes a maf file and filters a MAF dataframe by retaining only the rows with a column value included in the values list
#'
#' @param maf a maf as dataframe
#' @param schema a schema of datafrane check Select::TCGA_maf_schema for example
#' @param ... Other options
#' @return filtered_maf  a filtered maf file
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' filter_maf_truncating(luad_maf)
#'
#' @export
filter_maf_truncating <- function(maf, schema = TCGA_maf_schema, ...) {
  column <- "mutation.type"
  filtered_maf <- filter_maf_schema(maf, schema = schema, column = column, values = schema[[column]][["truncating"]], ...)
  return(filtered_maf)
}


###
#' This function filters a MAF dataframe by retaining (or discarding) missense mutations
#'
#' @description
#' `filter_maf_missense()` takes a maf file and filters a MAF dataframe by retaining only the rows with a column value included in the values list
#'
#' @param maf a maf as dataframe
#' @param schema a schema of datafrane check Select::TCGA_maf_schema for example
#' @param ... Other options
#' @return filtered_maf  a filtered maf file
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' filter_maf_missense(luad_maf)
#'
#' @export
filter_maf_missense <- function(maf, schema = TCGA_maf_schema, ...) {
  column <- "mutation.type"
  filtered_maf <- filter_maf_schema(maf, schema = schema, column = column, values = schema[[column]][["missense"]], ...)
  return(filtered_maf)
}


###
#' This function filters a MAF dataframe by retaining (or discarding) ignore mutations
#'
#' @description
#' `filter_maf_ignore()` takes a maf file and filters a MAF dataframe by retaining only the rows with a column value included in the values list
#'
#' @param maf a maf as dataframe
#' @param schema a schema of datafrane check Select::TCGA_maf_schema for example
#' @param ... Other options
#' @return filtered_maf  a filtered maf file
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' filter_maf_ignore(luad_maf)
#'
#' @export
filter_maf_ignore <- function(maf, schema = TCGA_maf_schema, ...) {
  column <- "mutation.type"
  filtered_maf <- filter_maf_schema(maf, schema = schema, column = column, values = schema[[column]][["ignore"]], ...)
  return(filtered_maf)
}


## Summary functions

###
#' Summary functions for MAF file
#'
#' @description
#' `stat_maf_column()` takes a maf file and filters a MAF dataframe by retaining only the rows with a column value included in the values list
#'
#' @param maf a maf as dataframe
#' @param column a schema of datafrane check Select::TCGA_maf_schema for example
#' @param ... Other options
#' @return filtered_maf  a filtered maf file
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' stat_maf_column(luad_maf, column = "Variant_Classification")
#'
#' @export

stat_maf_column <- function(maf, column, ...) {
  if (!column %in% colnames(maf)) stop(paste(column, "is not a valid column for the specified maf"))
  v <- maf[, column]
  return(table(v))
}
###
#' Count mutations per sample in a MAF file
#'
#' @param maf A MAF dataframe
#' @param column Column name containing sample IDs
#' @param ... Additional arguments passed to stat_maf_column
#' @return Table of mutation counts per sample.
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' head(stat_maf_sample(luad_maf))
#'
#' @export
stat_maf_sample <- function(maf, column = "Tumor_Sample_Barcode", ...) {
  return(stat_maf_column(maf, column, ...))
}
###
#' Count mutations per gene in a MAF file
#'
#' @param maf A MAF dataframe
#' @param column Column name containing gene symbols
#' @param ... Additional arguments passed to stat_maf_column
#' @return Table of mutation counts per gene.
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' head(stat_maf_gene(luad_maf))
#'
#' @export
stat_maf_gene <- function(maf, column = "Hugo_Symbol", ...) {
  return(stat_maf_column(maf, column, ...))
}

## GAM generating function

###
#' Generate gam from the maf file
#'
#' @description
#' `maf2gam()` takes a maf file and converts into gam
#'
#' @param maf A MAF dataframe
#' @param sample.col Column name for sample IDs
#' @param gene.col Column name for gene symbols
#' @param value.var Column used as the value to aggregate
#' @param samples Vector of sample IDs to include; NULL keeps all samples present in the MAF
#' @param genes Vector of gene names to include; NULL keeps all genes present in the MAF
#' @param fun.aggregate Aggregation function applied per (sample, gene) cell
#' @param binarize If TRUE, convert aggregated counts to binary presence/absence
#' @param fill Value used for missing (sample, gene) combinations
#' @return Numeric matrix (samples x genes) representing the gene alteration matrix.
#'
#' @examples
#' data(luad_maf, package = "SelectSim")
#' gam <- maf2gam(luad_maf)
#' dim(gam)
#'
#' @export
maf2gam <- function(maf,
                    sample.col = "Tumor_Sample_Barcode",
                    gene.col = "Hugo_Symbol",
                    value.var = "HGVSp_Short",
                    samples = NULL,
                    genes = NULL,
                    fun.aggregate = length,
                    binarize = TRUE,
                    fill = NA) {
  mut_gam <- tapply(maf[[value.var]], list(maf[[sample.col]], maf[[gene.col]]), FUN = fun.aggregate)
  if (binarize) {
    if (mode(mut_gam) == "character") {
      mut_gam <- mut_gam != ""
    } else if (mode(mut_gam) == "numeric") {
      mut_gam <- mut_gam > 0
    } else {
      mut_gam <- mut_gam > 0
    }
  }
  if (!is.null(samples)) {
    samples <- unique(samples)
    tokeep <- intersect(samples, rownames(mut_gam))
    mut_gam <- mut_gam[tokeep, , drop = FALSE]
    toadd <- samples[which(!samples %in% rownames(mut_gam))]
    voidadd <- matrix(fill, nrow = length(toadd), ncol = ncol(mut_gam))
    rownames(voidadd) <- toadd
    colnames(voidadd) <- colnames(mut_gam)
    mut_gam <- rbind(mut_gam, voidadd)
  }
  if (!is.null(genes)) {
    genes <- unique(genes)
    tokeep <- intersect(genes, colnames(mut_gam))
    mut_gam <- mut_gam[, tokeep, drop = FALSE]
    toadd <- genes[which(!genes %in% colnames(mut_gam))]
    voidadd <- matrix(fill, ncol = length(toadd), nrow = nrow(mut_gam))
    colnames(voidadd) <- toadd
    rownames(voidadd) <- rownames(mut_gam)
    mut_gam <- cbind(mut_gam, voidadd)
  }
  return(mut_gam)
}
