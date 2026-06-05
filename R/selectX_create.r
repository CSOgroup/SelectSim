###
# Author  : Arvind Iyer, Miljan Petrovic
# Project : SelectSim
# Desc    : Implementation of SelectSim algorithm
# Version : 0.1.6
###

#' Create an Alteration Landscape (AL) object
#'
#' Builds an Alteration Landscape object from a list of genome alteration matrices
#' and their corresponding tumor mutation burdens.
#'
#' @param am A named list with two required elements: \code{M} (a named list of binary
#'   alteration matrices, each genes x samples) and \code{tmb} (a named list of data
#'   frames, one per matrix in \code{M}, each with columns \code{sample} and
#'   \code{mutation}). The names of \code{M} and \code{tmb} must match. All matrices
#'   in \code{M} must have identical row and column names.
#' @param feat.covariates Named character vector of alteration-type annotations, one
#'   entry per feature (gene). Names must match rownames of the matrices in \code{M}.
#'   If \code{NULL}, all features are labelled \code{"MUT"}.
#' @param sample.covariates Named character vector of sample-type annotations, one
#'   entry per sample. Names must match colnames of the matrices in \code{M}. If
#'   \code{NULL}, all samples are labelled \code{"sample"}.
#' @param min.freq Minimum number of samples a gene must be mutated in (strictly
#'   greater than) to be retained. Features with \code{rowSums <= min.freq} are
#'   dropped.
#' @param verbose Logical; print progress messages.
#' @return An Alteration Landscape (AL) object (list of class \code{"AL"}) containing
#'   the filtered alteration matrices, TMB vectors, and covariate assignments.
#'
#' @examples
#' \donttest{
#' data(luad_run_data, package = "SelectSim")
#' al <- new.AL.general(
#'   am               = luad_run_data$M,
#'   feat.covariates  = luad_run_data$alteration.class,
#'   sample.covariates = luad_run_data$sample.class,
#'   min.freq         = 10
#' )
#' }
#'
#' @export
new.AL.general <- function(am,
                           feat.covariates = NULL,
                           sample.covariates = NULL,
                           min.freq,
                           verbose = FALSE) {
  if (is.null(am$M)) stop("am$M is NULL: provide a named list of alteration matrices.")
  if (is.null(am$tmb)) stop("am$tmb is NULL: provide a named list of TMB data frames.")
  if (!all(names(am$M) %in% names(am$tmb))) {
    stop("All names in am$M must have a matching entry in am$tmb.")
  }
  if (!all(sapply(am$tmb, function(t) all(c("sample", "mutation") %in% colnames(t))))) {
    stop("Each data frame in am$tmb must have 'sample' and 'mutation' columns.")
  }
  if (length(am$M) > 1) {
    ref_rows <- rownames(am$M[[1]])
    ref_cols <- colnames(am$M[[1]])
    if (!all(sapply(am$M[-1], function(m) identical(rownames(m), ref_rows)))) {
      stop("All matrices in am$M must have identical rownames.")
    }
    if (!all(sapply(am$M[-1], function(m) identical(colnames(m), ref_cols)))) {
      stop("All matrices in am$M must have identical colnames.")
    }
  }

  # create the alteration landscape object
  al <- list("am" = list())
  # Create a full gam
  al$am[["full"]] <- matrix(0, nrow = nrow(am$M[[1]]), ncol = ncol(am$M[[1]]))
  # Keep the order of rownames and colnames as in first gam.
  row.order <- sort(rownames(am$M[[1]])) # fixed gene ordering issuses
  col.order <- colnames(am$M[[1]])
  for (i in names(am$M)) {
    al$am[[i]] <- as.matrix(am$M[[i]][row.order, col.order])
    al$am[["full"]] <- al$am[["full"]] + as.matrix(am$M[[i]][row.order, col.order])
  }
  al$am[["full"]][al$am[["full"]] >= 1] <- 1 # binarize

  # Filtering steps
  # Filtering the al to keep only genes mutated with min.freq
  if (verbose) {
    print(paste("Number of features:", nrow(al$am$full), sep = " "))
  }
  feat <- rownames(al$am$full[rowSums(al$am$full) > min.freq, ])
  for (i in names(al$am)) {
    al$am[[i]] <- al$am[[i]][feat, ]
  }
  if (verbose) {
    print(paste("Number of features after filtering:", nrow(al$am$full), sep = " "))
  }
  # Set the feature and sample co-varaites
  al$alterations <- list()
  al$samples <- list()
  if (is.null(feat.covariates)) {
    al$alterations$alteration.class <- rep("MUT", nrow(al$am[[1]]))
    names(al$alterations$alteration.class) <- rownames(al$am[[1]])
  } else {
    al$alterations$alteration.class <- feat.covariates
  }
  if (is.null(sample.covariates)) {
    al$samples$sample.class <- rep("sample", ncol(al$am[[1]]))
    names(al$samples$sample.class) <- colnames(al$am[[1]])
  } else {
    al$samples$sample.class <- sample.covariates
  }
  # set the tumor mutation burden vector
  al$tmb <- list()
  for (i in names(am$M)) {
    al$tmb[[i]] <- am$tmb[[i]]
  }
  al$tmb[["total"]] <- c(rep(0, ncol(am$M[[1]])))
  for (i in names(am$tmb)) {
    al$tmb[["total"]] <- al$tmb[["total"]] + am$tmb[[i]][, c("mutation")]
  }
  names(al$tmb$total) <- am$tmb[[1]]$sample

  class(al) <- "AL"
  return(al)
}

#' Get sample/alteration blocks
#'
#' @param al The alteration landscape
#' @return Classification of samples and alterations in blocks.
#'
#' @examples
#' \donttest{
#' data(luad_run_data, package = "SelectSim")
#' al <- new.AL.general(luad_run_data$M,
#'                      feat.covariates  = luad_run_data$alteration.class,
#'                      sample.covariates = luad_run_data$sample.class,
#'                      min.freq = 10)
#' get.blocks(al)
#' }
#'
#' @export
get.blocks <- function(al) {
  if (is.null(al$alterations$alteration.class)) {
    al$alterations$alteration.class <- rep("MUT", nrow(al$am[[1]]))
    names(al$alterations$alteration.class) <- rownames(al$am[[1]])
  }
  if (is.null(al$samples$sample.class)) {
    al$samples$sample.class <- rep("sample", ncol(al$am[[1]]))
    names(al$samples$sample.class) <- colnames(al$am[[1]])
  }
  alteration.class <- al$alterations$alteration.class[rownames(al$am$full)]
  sample.class <- al$samples$sample.class[colnames(al$am$full)]
  feature.blocks <- lapply(unique(alteration.class), function(x) which(alteration.class == x))
  names(feature.blocks) <- unique(alteration.class)
  sample.blocks <- lapply(unique(sample.class), function(x) which(sample.class == x))
  names(sample.blocks) <- unique(sample.class)
  missing_rows <- which(!1:nrow(al$am$full) %in% unlist(feature.blocks))
  missing_cols <- which(!1:ncol(al$am$full) %in% unlist(sample.blocks))
  return(list(
    "sample.blocks" = sample.blocks,
    "alteration.blocks" = feature.blocks,
    "missing.samples" = missing_cols,
    "missing.alterations" = missing_rows
  ))
}

#' Generate S matrix
#'
#' @param gam the gam with genes*samples
#' @param sample.weights the samples weights
#' @param upperBound clip the values greater than 1 to keep it bounded between 0 to 1
#' @return S the S matrix
#'
#' @examples
#' gam <- matrix(c(0,1,1,0,1,1), nrow = 2,
#'               dimnames = list(c("geneA","geneB"), c("s1","s2","s3")))
#' generateS(gam, sample.weights = c(s1 = 1, s2 = 1, s3 = 1))
#'
#' @export
generateS <- function(gam,
                      sample.weights,
                      upperBound = 1) {
  gene.freq <- as.matrix(rowSums(gam) / ncol(gam), ncol = 1)
  tmb.weight <- (as.matrix(sample.weights, nrow = 1))
  sim.gam <- tcrossprod(gene.freq, tmb.weight)
  colnames(sim.gam) <- colnames(gam)
  sim.gam[sim.gam > 1] <- upperBound
  return(sim.gam)
}


#' Generate the template matrix
#'
#' Computes the expected mutation probability matrices (one per GAM type, per
#' sample block) used as the simulation background in \code{null_model_parallel}.
#'
#' @param al Alteration landscape object
#' @return A list with \code{template.obj} (per-block S matrices) and
#'   \code{temp_mat} (full concatenated template matrices, one per GAM type).
#'
#' @examples
#' \donttest{
#' data(luad_run_data, package = "SelectSim")
#' al <- new.AL.general(luad_run_data$M,
#'                      feat.covariates  = luad_run_data$alteration.class,
#'                      sample.covariates = luad_run_data$sample.class,
#'                      min.freq = 10)
#' template.obj.gen(al)
#' }
#'
#' @export
template.obj.gen <- function(al) {
  temp.al <- list()
  category <- c()
  for (name in names(al$am)) {
    if (name == "full") {} else {
      temp.al[[name]] <- list()
      category <- c(category, name)
    }
  }
  template_list <- list()
  for (name in category) {
    template_list[[name]] <- NULL
  }
  blocks <- get.blocks(al)
  for (block in names(blocks$sample.blocks)) {
    for (i in c(1:length(category))) {
      name <- category[i]
      samples <- names(blocks$sample.blocks[[block]])
      sample.weights <- length(samples) * al$tmb[[name]][samples, "mutation"] / sum(al$tmb[[name]][samples, "mutation"])
      S <- generateS(al$am[[name]][, samples], sample.weights)
      temp.al[[name]][[block]] <- S
      temp.al[[name]][[paste(block, "_weight", sep = "")]] <- sample.weights
      template_list[[name]] <- cbind(template_list[[name]], S)
    }
  }
  return(list("template.obj" = temp.al, "temp_mat" = template_list))
}


#' Generate sample weight matrix from TMB values
#'
#' @description
#' Computes a per-sample weight matrix based on the ratio of each sample's TMB to
#' the expected (mean) TMB. Samples with higher-than-expected TMB receive lower
#' weights via a penalty \code{lambda}, controlled by the fold-change threshold
#' \code{tau}.
#'
#' @param tmb Numeric vector of per-sample TMB values.
#' @param mean_tmb Numeric scalar; the reference (expected) TMB used to compute
#'   fold changes.
#' @param ngenes Integer; number of genes (rows) in the output weight matrix.
#' @param lambda Numeric; weight penalty factor. Higher values penalise
#'   high-TMB samples more strongly (default 0.3).
#' @param tau Numeric; fold-change threshold below which no penalty is applied
#'   (default 1).
#' @param discrete Logical; if \code{TRUE}, fold changes are rounded up before
#'   applying the penalty (default \code{TRUE}).
#' @return Numeric matrix of sample weights (ngenes x length(tmb)).
#'
#' @examples
#' tmb      <- c(s1 = 10, s2 = 50, s3 = 20)
#' mean_tmb <- 25
#' generateW_mean_tmb(tmb, mean_tmb, ngenes = 3)
#'
#' @export
generateW_mean_tmb <- function(tmb,
                               mean_tmb,
                               ngenes,
                               lambda = 0.3,
                               tau = 1,
                               discrete = TRUE) {
  exp.tmb <- mean_tmb
  tmb.FC <- tmb / exp.tmb
  tmb.FC[tmb.FC <= tau] <- tau
  w <- 1 / (1 + lambda * (ceiling(tmb.FC) - tau))
  if (!discrete) {
    w <- 1 / (1 + lambda * (tmb.FC - tau))
  }
  W <- matrix(rep(w, ngenes), nrow = ngenes, byrow = T)
  return(W)
}


#' Generate block-aware sample weight matrix
#'
#' @description
#' Computes a sample weight matrix that accounts for sample-class covariates
#' (blocks). Within each block the reference TMB is the block median; the
#' overall reference used by \code{generateW_mean_tmb} is the mean of those
#' block medians.
#'
#' @param al Alteration landscape object (from \code{new.AL.general}).
#' @param lambda Numeric penalty factor for the weight computation.
#' @param tau Numeric fold-change threshold below which no penalty is applied.
#' @return List with \code{W_block} (per-block weight matrices), \code{W} (full
#'   concatenated weight matrix), \code{W_median} (per-block median TMBs), and
#'   \code{mean_TMB} (mean of block medians).
#'
#' @importFrom  stats median
#'
#' @examples
#' \donttest{
#' data(luad_run_data, package = "SelectSim")
#' al <- new.AL.general(luad_run_data$M,
#'                      feat.covariates  = luad_run_data$alteration.class,
#'                      sample.covariates = luad_run_data$sample.class,
#'                      min.freq = 10)
#' generateW_block(al, lambda = 0.3, tau = 1)
#' }
#'
#' @export
generateW_block <- function(al, lambda, tau) {
  blocks <- get.blocks(al)
  W <- list()
  W_mtx <- NULL
  median_TMB <- list()
  mean_TMB <- 0
  # Mean of medians of tmb
  for (block in names(blocks$sample.blocks)) {
    samples <- names(blocks$sample.blocks[[block]])
    exp.tmb <- median(al$tmb$total[samples])
    median_TMB[[block]] <- exp.tmb
    mean_TMB <- mean_TMB + exp.tmb
  }
  mean_TMB <- mean_TMB / length(names(blocks$sample.blocks))
  for (block in names(blocks$sample.blocks)) {
    samples <- names(blocks$sample.blocks[[block]])
    W[[block]] <- generateW_mean_tmb(al$tmb$total[samples], mean_TMB, nrow(al$am$full), lambda = lambda, tau = tau)
    colnames(W[[block]]) <- samples
    rownames(W[[block]]) <- rownames(al$am$full)
    W_mtx <- cbind(W_mtx, W[[block]])
  }
  return(list("W_block" = W, "W" = W_mtx, "W_median" = median_TMB, "mean_TMB" = mean_TMB))
}

#' Generating the null_simulation matrix
#'
#' @import doParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom Rfast rowSort
#' @import doRNG
#' @param al Alteration landscape object
#' @param temp_mat template matrices
#' @param W weight matrix
#' @param n.cores Number of cores
#' @param n.permut Number of simulations
#' @param seed Random seed
#' @return List of simulated binary matrices (one per permutation)
#'
#' @examples
#' \donttest{
#' data(luad_run_data, package = "SelectSim")
#' al       <- new.AL.general(luad_run_data$M,
#'                            feat.covariates  = luad_run_data$alteration.class,
#'                            sample.covariates = luad_run_data$sample.class,
#'                            min.freq = 10)
#' temp_obj <- template.obj.gen(al)
#' W        <- generateW_block(al, lambda = 0.3, tau = 1)
#' sims     <- null_model_parallel(al, temp_obj$temp_mat, W$W,
#'                                 n.cores = 1, n.permut = 10)
#' length(sims)
#' }
#'
#' @export
null_model_parallel <- function(al,
                                temp_mat,
                                W,
                                n.cores = 1,
                                n.permut,
                                seed = 42) {
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  # Set RNG kind and seed for reproducibility
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
  simulationStep <- function(template, total_mut) {
    S <- template
    nvalues <- nrow(S) * ncol(S)
    r <- matrix(runif(nvalues, min = 0, max = 1), nrow = nrow(S), ncol = ncol(S))
    test <- S - r
    rownames(r) <- rownames(S)
    colnames(r) <- colnames(S)
    return(test)
  }

  combine <- function(residuals, residuals_sort, val, index) {
    return(residuals[index, ] >= residuals_sort[index, ][val[index]])
  }

  simulationFixedOnes <- function(al, temp_mat, W) {
    exp <- matrix(0, nrow = nrow(W), ncol = ncol(W))
    gam_names <- names(al$am)
    gam_incidence <- list()
    k <- 1
    for (name in gam_names) {
      if (name == "full") {} else {
        gam_incidence[[k]] <- sum(al$am[[name]])
        k <- k + 1
      }
    }
    residuals <- list()
    for (i in 1:length(temp_mat)) {
      S <- temp_mat[[i]]
      test <- simulationStep(template = S, total_mut = gam_incidence[[i]])
      residuals[[i]] <- test
    }
    residual_mtx <- Reduce(pmax, residuals)
    residual_mtx_sort <- Rfast::rowSort(residual_mtx, descending = TRUE)
    temp <- lapply(seq_len(nrow(al$am$full)), function(i) {
      combine(residual_mtx, residual_mtx_sort, rowSums(al$am$full), i)
    })
    exp <- do.call(rbind, temp) * 1
    rownames(exp) <- rownames(al$am$full)
    return(exp)
  }


  if (n.cores > 1) {
    log_file <- tempfile(pattern = "gen.random.am_", fileext = ".log")
    cl <- parallel::makeCluster(n.cores, outfile = log_file)
    doParallel::registerDoParallel(cl)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
      unlink(log_file)
      add <- TRUE
    })
    registerDoRNG(seed)
    randomMs <- foreach::foreach(i = 1:n.permut) %dopar% simulationFixedOnes(al, temp_mat, W)
    return(randomMs)
  } else {
    set.seed(seed)
    foreach::registerDoSEQ()
    registerDoRNG(seed)
    randomMs <- foreach(i = 1:n.permut) %do% simulationFixedOnes(al, temp_mat, W)
    return(randomMs)
  }
}

#' Identify outlier null-model matrices
#'
#' Flags simulations whose mean per-gene and per-sample deviation from the
#' observed counts falls in the top 10%, indicating numerical instability.
#' These are removed before computing effect sizes to prevent them from
#' inflating the null distribution.
#'
#' @param obj SelectX object (list with al and null fields)
#' @param nSim Number of simulations in the null model
#' @return Logical vector; TRUE for each null matrix flagged as an outlier.
#'
#' @examples
#' \donttest{
#' data(luad_run_data, package = "SelectSim")
#' result <- selectX(M = luad_run_data$M,
#'                   sample.class = luad_run_data$sample.class,
#'                   alteration.class = luad_run_data$alteration.class,
#'                   n.cores = 1, min.freq = 10, n.permut = 10,
#'                   verbose = FALSE)
#' outliers <- retrieveOutliers(result$obj, nSim = result$obj$nSim)
#' sum(outliers)
#' }
#'
#' @export
retrieveOutliers <- function(obj, nSim = 1000) {
  gnMut <- matrix(0, nrow = nrow(obj$al$am$full), ncol = nSim)
  snMut <- matrix(0, nrow = ncol(obj$al$am$full), ncol = nSim)
  for (i in 1:nSim) {
    gnMut[, i] <- rowSums(obj$null[[i]])
    snMut[, i] <- colSums(obj$null[[i]])
  }
  rgn <- gnMut - rowSums(obj$al$am$full)
  rsn <- snMut - colSums(obj$al$am$full)
  mean.gn <- apply(abs(rgn), 2, mean)
  mean.sn <- apply(abs(rsn), 2, mean)
  dev <- mean.gn + mean.sn
  dev2 <- sort(dev)
  maxcut <- dev2[round(0.90 * length(dev2))]
  outliers <- dev >= maxcut
  return(outliers)
}
