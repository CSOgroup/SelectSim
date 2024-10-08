###
# Author  : Arvind Iyer
# Email   : arvind.iyer@unil.ch
# Project : SelectSim
# Desc    : The alteration landscape object, template object and selectX object creating functions
# Version : 0.1
# Updates : Re wrote the whole code
# Todo:
# - Random seed to set to one value for parallel foreach?
# - c++ backend for faster computation
###

#' Create an AL object
#'
#' Create an Alteration Landscape (AL) object which contains gams and mutation burden of samples of associated gams.
#'
#' @param am The binary alteration matrix with row as features and column as samples.
#' @param feat.covariates gene/feature covariate.
#' @param sample.covariates sample covariate.
#' @param min.freq minimum frequency of genes to be mutated.
#' @param verbose print the time and each steps.
#' @return An Alteration Landscape (AL) object with the gam. 
#' @export
new.AL.general  <- function(am,
							feat.covariates=NULL,
							sample.covariates=NULL,
							min.freq,
                            verbose=FALSE) {
    
    # Assert if the AM are provided or not
    if(is.null(am$M)) 
        stop('Problem:-->  Input data is null')


    # Assert if the data structre is in correct format
    # if(is.list(am)){
    #     if(names(am) != c('M','tmb'))
    #         stop('Problem:-->  Input data (M) is not in correct structure\n Check the format which should be: 
    #             list(M=list("missesne"=matrix()), # which has features on rows and sample on columns
    #                  tmb=list("missesne"=data.frame()) # which has sample,mutation as column names and rownames as sample
    #                 )')
    # }
  
    # create the alteration landscape object
    al = list('am'=list()) 
    # Create a full gam  
    al$am[['full']] = matrix(0, nrow = nrow(am$M[[1]]), ncol = ncol(am$M[[1]]))
    # Keep the order of rownames and colnames as in first gam.
    row.order = rownames(am$M[[1]])
    col.order = colnames(am$M[[1]])
    for ( i in names(am$M)){
        al$am[[i]] <-as.matrix(am$M[[i]][row.order,col.order])
        al$am[['full']] <- al$am[['full']]+as.matrix(am$M[[i]][row.order,col.order])
     }
    al$am[['full']][al$am[['full']]>=1] <-1 #binarize             
    
    #Filtering steps
    #Filtering the al to keep only genes mutated with min.freq
    if(verbose)
        print(paste('Number of features:',nrow(al$am$full),sep=" "))
    feat <- rownames(al$am$full[rowSums(al$am$full)>min.freq,])
    for ( i in names(al$am)){
     al$am[[i]] <- al$am[[i]][feat,]
    }
    if(verbose)         
        print(paste('Number of features after filtering:',nrow(al$am$full),sep=" "))
    #Set the feature and sample co-varaites
    al$alterations = list()
    al$samples = list()
    if(is.null(feat.covariates)) { 
        al$alterations$alteration.class = rep('MUT', nrow(al$am[[1]]));
        names(al$alterations$alteration.class) = rownames(al$am[[1]]) 
    }
    else
    {
    	al$alterations$alteration.class <- feat.covariates
    }
    if(is.null(sample.covariates)) { 
        al$samples$sample.class = rep('sample', ncol(al$am[[1]]));
        names(al$samples$sample.class) = colnames(al$am[[1]]) 
    }
    else
    {
    	al$samples$sample.class <- sample.covariates
    }
    #set the tumor mutation burden vector
    al$tmb =list()
    for (i in names(am$M)){
        al$tmb[[i]]<-am$tmb[[i]]
    }
    al$tmb[['total']]<- c(rep(0,ncol(am$M[[1]])))
    for (i in names(am$tmb)){
        al$tmb[['total']] <- al$tmb[['total']] +  am$tmb[[i]][,c('mutation')]
    }
    names(al$tmb$total)<-am$tmb[[1]]$sample
    
    class(al) <- "AL";
    return(al)
}

#' Get sample/alteration blocks
#' 
#' @param al The alteration landscape
#' @return Classification of samples and alterations in blocks.
#' 
#' @export
get.blocks <- function(al) {
    if(is.null(al$alterations$alteration.class)) { 
        al$alterations$alteration.class = rep('MUT', nrow(al$am[[1]]));
        names(al$alterations$alteration.class) = rownames(al$am[[1]])
    }
    if(is.null(al$samples$sample.class)) { 
        al$samples$sample.class = rep('sample', ncol(al$am[[1]]));
        names(al$samples$sample.class) = colnames(al$am[[1]])
    }
    alteration.class = al$alterations$alteration.class[rownames(al$am$full)]
    sample.class = al$samples$sample.class[colnames(al$am$full)]
    feature.blocks = lapply(unique(alteration.class), function(x) which(alteration.class==x))
    names(feature.blocks) = unique(alteration.class)
    sample.blocks = lapply(unique(sample.class), function(x) which(sample.class==x))
    names(sample.blocks) = unique(sample.class)
    missing_rows = which(! 1:nrow(al$am$full) %in% unlist(feature.blocks))
    missing_cols = which(! 1:ncol(al$am$full) %in% unlist(sample.blocks))
    return(list('sample.blocks' = sample.blocks,
                'alteration.blocks' = feature.blocks,
                'missing.samples' = missing_cols,
                'missing.alterations' = missing_rows))
}

#' Generate S matrix
#' 
#' @param gam the gam with genes*samples
#' @param sample.weights the samples weights
#' @param upperBound clip the values greater than 1 to keep it bounded between 0 to 1
#' @return S the S matrix 
#'
#' @export
generateS = function(gam,
                     sample.weights,
                     upperBound = 1){

    gene.freq = as.matrix(rowSums(gam)/ncol(gam), ncol = 1)
    tmb.weight = t(as.matrix(sample.weights, ncol = 1))
    sim.gam = gene.freq %*% tmb.weight
    colnames(sim.gam) = colnames(gam)
    sim.gam[ sim.gam > 1 ] = upperBound
    return(sim.gam)
}


#' Generating the template matrix   
#' 
#' Computing the expected mutation background matrix using which simulation is done.
#'
#' @param al Alteration landscape object
#' @return Template matrix as list object
#' 
#' @export
templeate.obj.gen <- function(al){
    temp.al <- list()
    category <- c()
    for (name in names(al$am)){
        if (name == 'full'){}
        else{
         temp.al[[name]]<-list()       
         category <- c(category,name)
        }
    }
    T <- list()
    for (name in category){
        T[[name]]<-NULL
    }
    blocks <- get.blocks(al)    
    for(block in names(blocks$sample.blocks)){
        for(i in c(1:length(category))){
            name=category[i]
            samples= names(blocks$sample.blocks[[block]]) 
            sample.weights=length(samples)*al$tmb[[name]][samples,'mutation']/sum(al$tmb[[name]][samples,'mutation'])            
            S=generateS(al$am[[name]][,samples],sample.weights)
            temp.al[[name]][[block]]<-S
            temp.al[[name]][[paste(block,'_weight',sep="")]]<-sample.weights
            T[[name]]<-cbind(T[[name]],S)
        }
    }
    return (list('templeate.obj'=temp.al,'temp_mat'=T))                        
}


#' Generating the weight matrix
#' 
#' @param tmb TMB dataframe
#' @param mean_tmb TMB dataframe
#' @param ngenes Number of genes 
#' @param lambda 0.1 weight factor
#' @param tao 1 Fold change factor
#' @param discrete True discrete weights
#' @return weight matrix (i.e vector as matrix)
#' 
#' @export
generateW_mean_tmb = function(tmb,
                     mean_tmb,
                     ngenes,
                     lambda = 0.3,
                     tao = 1,
                     discrete = TRUE ){

    exp.tmb = mean_tmb
    tmb.FC = tmb/exp.tmb
    tmb.FC[ tmb.FC <= tao ] = tao
    w = 1 /(1+ lambda*(ceiling(tmb.FC)-tao))
    if(!discrete)
         w = 1 /(1+ lambda*(tmb.FC-tao))
    W = matrix( rep(w, ngenes), nrow = ngenes, byrow = T)
    return(W)
}


#' Generating the weight matrix taking sample covariate
#' 
#' @param al altearting landscape onject
#' @param lambda penalty factor parameter
#' @param tao fold change paramenter
#' @return Weight matrix 
#'
#' @importFrom  stats median
#' @export
generateW_block = function(al,lambda,tao) {
    blocks <- get.blocks(al)
    W <-list()
    W_mtx<-NULL
    median_TMB<-list()
    mean_TMB<-0
    # Mean of medians of tmb 
    for(block in names(blocks$sample.blocks)){
        samples= names(blocks$sample.blocks[[block]])
        exp.tmb = median(al$tmb$total[samples])
        median_TMB[[block]]<- exp.tmb
        mean_TMB=mean_TMB+exp.tmb
    }
    mean_TMB<-mean_TMB/length(names(blocks$sample.blocks))
    #print(mean_TMB)
    for(block in names(blocks$sample.blocks)){
        samples= names(blocks$sample.blocks[[block]])
        W[[block]]<-generateW_mean_tmb(al$tmb$total[samples],mean_TMB,nrow(al$am$full),lambda=lambda,tao=tao)
        colnames(W[[block]])<-samples
        rownames(W[[block]])<-rownames(al$am$full)
        W_mtx <- cbind(W_mtx,W[[block]])
    }
    return(list('W_block'=W,'W'=W_mtx,'W_median'=median_TMB,'mean_TMB'=mean_TMB))

}

#' Generating the null_simulation matrix  
#' 
#' @import doParallel
#' @import parallel
#' @importFrom Rfast rowSort
#' @import doRNG
#' @param al Alteration landscape object
#' @param temp_mat template matrixes
#' @param W weight matrix 
#' @param n.cores Number of cores
#' @param n.permut Number of simulation
#' @return Template matrix as list object
#'
#' Ordering of genes impacts the results as residual subraction is not correct.
#' @export
null_model_parallel <-function(al,
                               temp_mat,
                               W,
                               n.cores=1,
                               n.permut) {

    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`


    simulationStep = function(template, total_mut){
        S = template
        nvalues = nrow(S)*ncol(S)
        r = matrix( runif(nvalues, min = 0, max = 1), nrow = nrow(S), ncol = ncol(S))
        test = S - r
        return(test)
    }

    combine <- function(residuals,residuals_sort,val,index){
        return(residuals[index,]>=residuals_sort[index,][val[index]])
    }

    simulationFixedOnes = function(al,temp_mat, W,CORRECT_T = 0.05){
        exp = matrix(0, nrow = nrow(W), ncol = ncol(W))
        gam_names = names(al$am)
        gam_incidence = list()
        k=1
        for(name in gam_names){
            if (name=='full'){
            }
            else{
                gam_incidence[[k]]=sum(al$am[[name]])
                k=k+1
            }
        }
        residuals = list()
        for(i in 1:length(temp_mat)){    
            S = temp_mat[[i]]
            test = simulationStep(template = S, total_mut = gam_incidence[[i]])
            residuals[[i]]<-test
        }
        total_mut = sum(al$am$full)
        gene_freq = rowSums(al$am$full)
        if(length(gam_names)>2){
            residual_mtx <- pmax(residuals[[1]],residuals[[2]])
            residual_mtx_sort <- Rfast::rowSort(residual_mtx, descending = TRUE)
            temp<-list()
            for(i in c(1:nrow(al$am$full))){
                temp[[i]]<-combine(residual_mtx,residual_mtx_sort,rowSums(al$am$full),i)
            }
            exp<-do.call(rbind,temp)*1
            rownames(exp) <- rownames(al$am$full)
            tot = sum(exp)
            exp.freq = rowSums(exp)
            diff = abs(exp.freq - rowSums(al$am$full))/ncol(W)
            return(exp)
        }
        else{
            residual_mtx <- residuals[[1]]
            residual_mtx_sort <- Rfast::rowSort(residual_mtx, descending = TRUE)
            #print(dim(residual_mtx))
            #print(dim(residual_mtx_sort))
            temp<-list()
            for(i in c(1:nrow(al$am$full))){
                temp[[i]]<-combine(residual_mtx,residual_mtx_sort,rowSums(al$am$full),i)
            }
            exp<-do.call(rbind,temp)*1
            rownames(exp) <- rownames(al$am$full)
            tot = sum(exp)
            #print(tot)
            #print(gam_incidence[[1]])
            exp.freq = rowSums(exp)
            diff = abs(exp.freq - rowSums(al$am$full))/ncol(W)
            return(exp)
        }   
    }

    
    if(n.cores>1){
        cl <-  parallel::makeCluster(n.cores, outfile=paste("gen.random.am.log", sep=''))
        doParallel::registerDoParallel(cl)  
        if( foreach::getDoParRegistered()) {
            registerDoRNG(42)
            randomMs <- foreach::foreach(i=1:n.permut) %dopar% simulationFixedOnes(al,temp_mat,W,CORRECT_T=0.05)
            parallel::stopCluster(cl)
            return (randomMs)
        } else {
            stop('Error in registering a parallel cluster for randomization.')
        }
    }
    else{
        foreach::registerDoSEQ()
        randomMs <- foreach(i=1:n.permut) %do% simulationFixedOnes(al,temp_mat,W,CORRECT_T=0.05)
        return (randomMs)
    }
}



#' Generating the null_simulation matrix  
#' 
#' @import doParallel
#' @import parallel
#' @importFrom Rfast rowSort
#' @import doRNG
#' @param al Alteration landscape object
#' @param temp_mat template matrixes
#' @param W weight matrix 
#' @param n.cores Number of cores
#' @param n.permut Number of simulation
#' @return Template matrix as list object
#'
#' Ordering of genes impacts the results as residual subraction is not correct.
#' @export
null_model_parallel_debug <-function(al,
                               temp_mat,
                               W,
                               n.cores=1,
                               n.permut) {

    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`


    simulationStep = function(template, total_mut){
        S = template
        nvalues = nrow(S)*ncol(S)
        r = matrix( runif(nvalues, min = 0, max = 1), nrow = nrow(S), ncol = ncol(S))
        rownames(r)<-rownames(S)
        colnames(r)<-colnames(S)
        test = S - r[rownames(S),colnames(S)]
        return(test)
    }

    combine <- function(residuals,residuals_sort,val,index){
        return(residuals[index,]>=residuals_sort[index,][val[index]])
    }

    simulationFixedOnes = function(al,temp_mat, W,CORRECT_T = 0.05){
        exp = matrix(0, nrow = nrow(W), ncol = ncol(W))
        gam_names = names(al$am)
        gam_incidence = list()
        k=1
        for(name in gam_names){
            if (name=='full'){
            }
            else{
                gam_incidence[[k]]=sum(al$am[[name]])
                k=k+1
            }
        }
        residuals = list()
        for(i in 1:length(temp_mat)){    
            S = temp_mat[[i]]
            test = simulationStep(template = S, total_mut = gam_incidence[[i]])
            residuals[[i]]<-test
        }
        total_mut = sum(al$am$full)
        gene_freq = rowSums(al$am$full)
        if(length(gam_names)>2){
            residual_mtx <- pmax(residuals[[1]],residuals[[2]])
            residual_mtx_sort <- Rfast::rowSort(residual_mtx, descending = TRUE)
            temp<-list()
            for(i in c(1:nrow(al$am$full))){
                temp[[i]]<-combine(residual_mtx,residual_mtx_sort,rowSums(al$am$full),i)
            }
            exp<-do.call(rbind,temp)*1
            rownames(exp) <- rownames(al$am$full)
            tot = sum(exp)
            exp.freq = rowSums(exp)
            diff = abs(exp.freq - rowSums(al$am$full))/ncol(W)
            return(exp)
        }
        else{
            residual_mtx <- residuals[[1]]
            residual_mtx_sort <- Rfast::rowSort(residual_mtx, descending = TRUE)
            #print(dim(residual_mtx))
            #print(dim(residual_mtx_sort))
            temp<-list()
            for(i in c(1:nrow(al$am$full))){
                temp[[i]]<-combine(residual_mtx,residual_mtx_sort,rowSums(al$am$full),i)
            }
            exp<-do.call(rbind,temp)*1
            rownames(exp) <- rownames(al$am$full)
            tot = sum(exp)
            #print(tot)
            #print(gam_incidence[[1]])
            exp.freq = rowSums(exp)
            diff = abs(exp.freq - rowSums(al$am$full))/ncol(W)
            return(exp)
        }   
    }

    
    if(n.cores>1){
        cl <-  parallel::makeCluster(n.cores, outfile=paste("gen.random.am.log", sep=''))
        doParallel::registerDoParallel(cl)  
        if( foreach::getDoParRegistered()) {
            registerDoRNG(42)
            randomMs <- foreach::foreach(i=1:n.permut) %dopar% simulationFixedOnes(al,temp_mat,W,CORRECT_T=0.05)
            parallel::stopCluster(cl)
            return (randomMs)
        } else {
            stop('Error in registering a parallel cluster for randomization.')
        }
    }
    else{
        foreach::registerDoSEQ()
        randomMs <- foreach(i=1:n.permut) %do% simulationFixedOnes(al,temp_mat,W,CORRECT_T=0.05)
        return (randomMs)
    }
}

#' Removing the Outliers
#' 
#' @param obj selectX object
#' @param nSim Number of simulation
#' @return outliers
#'
#' @export
retrieveOutliers = function(obj, nSim=1000){
    gnMut = matrix(0, nrow = nrow(obj$al$am$full), ncol = nSim)
    snMut = matrix(0, nrow = ncol(obj$al$am$full), ncol = nSim)
    for(i in 1:nSim){
       gnMut[,i] = rowSums(obj$null[[i]])
       snMut[,i] = colSums(obj$null[[i]])
    }
    rgn = gnMut
    rsn = snMut
    ogn = rowSums(obj$al$am$full)
    osn = colSums(obj$al$am$full)
    rgn = rgn-ogn
    rsn = rsn-osn
    mean.gn = apply(abs(rgn), 2, mean)
    mean.sn = apply(abs(rsn), 2, mean)
    dev = mean.gn + mean.sn    
    dev2 = sort(dev)
    mincut = dev2[ round(0.05*length(dev2)) ]
    maxcut = dev2[ round(0.95*length(dev2)) ]
    outliers = dev < mincut | dev >= maxcut
    
    return(outliers)
}
