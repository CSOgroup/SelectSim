###
# Author  : Arvind Iyer
# Email   : arvind.iyer@unil.ch
# Project : SelectSim
# Desc    : Main file which call the selectX function to create alteration object with background model and funtion to generate the table.
# Version : 0.1.5
# Todo:
# - Better Error message and running text
# - Edge case: When sample size in less than 2 there is error in computation (need to fix a number to do this analysis)
# - parallel::makeCluster(2, setup_strategy = "sequential") a possible fix to remove the erorr of not able to connect problem (https://github.com/rstudio/rstudio/issues/6692)
# - Try to fix the parallel processing issuses
###

###
#' SelectX main function from SelectSim to create alteration object with background model
#'
#' @description
#' `selectX()` takes a list object which consist of genome alteration matrix and tumor mutation burden.
#'
#' @import foreach
#' @import doParallel
#' @importFrom  stats as.dist
#' @importFrom  stats density
#' @importFrom  stats ecdf
#' @importFrom  stats runif
#' @importFrom  stats setNames
#' @param M a list data object which consist of gams and tmbs.
#' @param sample.class sample covariates as named list.
#' @param alteration.class alteration covariates as named list.
#' @param n.cores no of cores.
#' @param min.freq number of samples for features to be atleast mutated in.
#' @param n.permut number of simulations.
#' @param lambda lambda parameter.
#' @param tau tau (fold change) parameter.
#' @param folder folder path to store the results.
#' @param save.object store the SelectX object.
#' @param verbose print the time and each steps.
#' @param estimate_pairwise Compute pairwise p-value
#' @param maxFDR FDR value
#' @param seed a random seed
#' @return result  a SelectSim object with background model and other info along with result table
#'
#' @export
selectX <- function(M,
                    sample.class,
                    alteration.class,
                    n.cores = 1,
                    min.freq = 10,
                    n.permut = 1000,
                    lambda = 0.3,
                    tau = 1,
                    save.object = FALSE,
                    folder = './',
                    verbose = TRUE,
                    estimate_pairwise = FALSE,
                    maxFDR = 0.25,
                    seed = 42) {
    set.seed(seed)
    if (verbose) { print(paste('#### Creating SelectX object ####')); tictoc::tic('Total time taken:') }

    if (verbose) { print(paste('Step1-> Parsing and Filtering GAM...')); tictoc::tic('Time:') }
    al <- new.AL.general(M, alteration.class, sample.class, min.freq, verbose)
    if (verbose) { print(paste('-> Alteration Landscape object created')); tictoc::toc() }

    if (verbose) { print(paste('Step2-> Generating Template object...')); tictoc::tic('Time:') }
    temp.data <- template.obj.gen(al)
    if (verbose) { print(paste('-> Template object created')); tictoc::toc() }

    if (verbose) { print(paste('Step3-> Generating sample weight matrix...')); tictoc::tic('Time:') }
    W <- generateW_block(al, lambda, tau)
    if (verbose) { print(paste('-> Weight Matrix created')); tictoc::toc() }

    if (verbose) { print(paste('Step4-> Generating null model...')); tictoc::tic('Time:') }
    sim <- null_model_parallel(al, temp.data$temp_mat, W$W, n.cores, n.permut)
    obj <- list('al'=al, 'W'=W, 'T'=temp.data, 'null'=sim, 'nSim'=n.permut)
    if (verbose) print(paste('-> Removing the outliers matrix from null model...'))
    outliers <- retrieveOutliers(obj=obj, nSim=n.permut)
    if (length(which(outliers)) == 0) {
        if (verbose) print(paste('Removed null-matrix:', length(which(outliers)), sep=" "))
        obj$nSim <- n.permut
    } else {
        if (verbose) print(paste(' Removed null-matrix:', length(which(outliers)), sep=" "))
        obj$null <- obj$null[which(!outliers)]
        if (verbose) print(paste(' Updated the null-model and nSim variables...'))
        obj$nSim <- length(which(!outliers))
    }
    if (verbose) { tictoc::toc(); print(paste('-> Null model generated')); print(paste('### SelectSim object created ###')) }

    if (verbose) { print(paste('#### Computing EDs on the dataset ####')); tictoc::tic('Time:') }
    al <- obj$al
    als <- al.stats(obj)
    alp <- al.pairwise.alteration.stats(obj, als, do.blocks=FALSE)
    als[['alteration.pairwise']] <- alp
    als[['alteration.pairwise']]$sample.blocks <- NULL
    for (i in names(als$sample.blocks)) als$sample.blocks[[i]][['alteration.pairwise']] <- alp$sample.blocks[[i]]
    obs.co <- as.matrix(als$alteration.pairwise$overlap)
    wobs.co <- as.matrix(als$alteration.pairwise$w_overlap)
    robs.co <- r.am.pairwise.alteration.overlap(null=obj$null, n.permut=obj$nSim, n.cores=1)
    wrobs.co <- w.r.am.pairwise.alteration.overlap(null=obj$null, W=obj$W$W, n.permut=obj$nSim, n.cores=1)
    selectX_result <- interaction.table(al, als, obs.co, wobs.co, robs.co, wrobs.co,
                                        null=obj$null, maxFDR=maxFDR, n.cores=1,
                                        estimate_pairwise=estimate_pairwise, n.permut=obj$nSim)
    obj$robs.co <- robs.co
    obj$wrobs.co <- wrobs.co
    if (verbose) { tictoc::toc(); print(paste('#### EDs computed ####')) }

    if (verbose) tictoc::toc()

    if (save.object) {
        saveRDS(obj, file=paste(folder, 'selectsim_object.rds', sep=''))
        saveRDS(selectX_result, file=paste(folder, 'selectsim_results.rds', sep=''))
    }
    return(list(obj=obj, result=selectX_result))
}

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

## usethis namespace: start
#' @useDynLib SelectSim, .registration = TRUE
## usethis namespace: end
NULL
