###
# Author  : Arvind Iyer
# Email   : arvind.iyer@unil.ch
# Project : SelectSim
# Desc    : The file which contains the function to generate the stats and table
# Version : 0.1.5
# Updates : Re wrote the whole code
# Todo:
###

#' Create an alteration landscape object
#' 
#' @param al The alteration landscape
#' @return AlS object.
#' 
#' @export
new.ALS <- function(al) {
    if(is.null(al)) stop('Input al is NULL')
    als = list()
    class(als) <- "ALS";
    return(als)
}
#' Create an alteration matrix stats object
#' 
#' @param am The alteration matrix stats object
#' @return AMS object.
#'
#' @export
new.AMS <- function(am) {
    if(is.null(am)) stop('Input am is NULL')
    ams = list()
    class(ams) <- "AMS";
    return(ams)
}
#' Create an alteration matrix stats object
#' 
#' @param am The alteration matrix stat object
#' @return AMS object with basic stats for the matrix
#'
#' @export
am.stats <- function(am) {
    ams = new.AMS(am)
    am = am * 1.0
    sample.alt.n = colSums(am, na.rm=TRUE)
    feat.alt.n = rowSums(am, na.rm=TRUE)
    numOfEdges <- sum(abs(am), na.rm=TRUE)
    temp = list('n.samples' = ncol(am),
                'n.alterations' = nrow(am),
                'n.occurrences' = numOfEdges,
                'alterations.per.sample'=sample.alt.n,
                'alteration.count'=feat.alt.n)
    ams[names(temp)] = temp
    return(ams)
}
#' Create an alteration stats object
#' 
#' @param al The alteration landscape object
#' @return als alteration landscape stats object
#'
#' @export
al.stats <- function(al) {
    als = new.ALS(al)
    temp = am.stats(al$al$am$full)
    ams = temp
    als[names(ams)] = ams
    # get stats for each single block
    blocks = get.blocks(al$al)
    als$sample.blocks = list()
    #print(length(blocks$sample.blocks))
    for(ib in 1:length(blocks$sample.blocks)) { 
    	fib = blocks$sample.blocks[[ib]]; 
    	if(length(fib)>0) {
        subM <- al$al$am$full[,fib, drop=FALSE]
        temp = am.stats(subM)
        als$sample.blocks[[ib]] = temp
    }}
    names(als$sample.blocks) = names(blocks$sample.blocks)
    als$alteration.blocks = list()
    for(ib in 1:length(blocks$alteration.blocks)) { 
    	fib = blocks$alteration.blocks[[ib]]; if(length(fib)>0) {
        subM <- al$al$am$full[fib, , drop=FALSE]
        temp = am.stats(subM)
        als$alteration.blocks[[ib]] = temp
    }}
    names(als$alteration.blocks) = names(blocks$alteration.blocks)
    
    return(als)
}

###
#' Compute overlap stats
#' 
#' @importFrom  Matrix sparseMatrix
#' @importFrom  Matrix t
#' @param am The alteration matrix
#' @return overlap the overlap between the pairs
#'
#' @export
am.pairwise.alteration.overlap <- function(am) {
    A = am * 1
    overlap =A %*% t(A)
    return(overlap)
}
#' Compute weight overlap stats
#' 
#' @importFrom  Matrix sparseMatrix
#' @importFrom  Matrix t
#' @param am The alteration matrix
#' @param W The weight matrix
#' @return overlap the weight overlap between the pairs
#'
#' @export
am.weight.pairwise.alteration.overlap <- function(am,W) {

     A = am * 1
     col_order= colnames(A)
     overlap = (W[,col_order]*A) %*% t(A)
     return(overlap)
}
#' Compute weight overlap stats
#'
#' @importFrom  Matrix Matrix 
#' @param overlap_M The overlap matrix
#' @param M.stats The am matrix statis
#' @param w_overlap_M The weighted overlap matrix
#' @return stats of the am matrix
#'
#' @export
am.pairwise.alteration.coverage <- function(overlap_M, M.stats,w_overlap_M) {
    stats = list()
    overlap_M = as.matrix(overlap_M)
    w_overlap_M=as.matrix(w_overlap_M)
    marginal_f1 = replicate(ncol(overlap_M),M.stats$alteration.count[rownames(overlap_M)])
    colnames(marginal_f1) = rownames(marginal_f1)
    marginal_f2 = t(marginal_f1)
    diag(overlap_M) = diag(marginal_f1)
    # me = marginal_f1 + marginal_f2 - 2*overlap_M # not used. Nov 2018
    # coverage = marginal_f2 + marginal_f1 - overlap_M # not used. Nov 2018
    stats[['overlap']] = Matrix(overlap_M)
    stats[['w_overlap']] = Matrix(w_overlap_M)
    # stats[['coverage']] = coverage
    # stats[['me']] = me # not used. Nov 2018
    return(stats)
}
#' Compute weight overlap stats
#' 
#' @param al The alteration landscape object
#' @param als The alteration stats object
#' @param do.blocks blockwise comptutation
#' @return overlap the weight overlap between the pairs
#'
#' @export
# TODO do.block is error prone.
al.pairwise.alteration.stats <- function(al, als=NULL, do.blocks=FALSE) {
    if(is.null(als)) als = al.stats(al)
    M.overlap = am.pairwise.alteration.overlap(as.matrix(al$al$am$full))
    M.woverlap = am.weight.pairwise.alteration.overlap(al$al$am$full,al$W$W[,colnames(al$al$am$full)])
    M.pairwise = am.pairwise.alteration.coverage(M.overlap, als, M.woverlap)

    if(do.blocks) {
        blocks = get.blocks(al$al)
        pairwise.blocks = list()
        for(ib in 1:length(blocks$sample.blocks)) { 
        	fib = blocks$sample.blocks[[ib]]; 
        	if(length(fib)>0) {
            subM <- al$al$am$full[,fib, drop=FALSE]
            subW <- al$W$W[,fib,drop=FALSE]
            M.overlap.block  = am.pairwise.alteration.overlap(subM)
            M.woverlap.block = am.weight.pairwise.alteration.overlap(subM[,order(fib)],subW[,order(fib)])
            M.pairwise.block = am.pairwise.alteration.coverage(M.overlap.block, als$sample.blocks[[ib]],M.woverlap.block)
            pairwise.blocks[[ib]] = M.pairwise.block
        }}
        names(pairwise.blocks) = names(blocks$sample.blocks)
        M.pairwise[['sample.blocks']] = pairwise.blocks
    }
     return(M.pairwise)
}
#' Compute weight overlap stats
#' 
#' @import doParallel
#' @import parallel
#' @param null The null model
#' @param n.permut The number of permutation steps
#' @param n.cores The number of cores
#' @return overlap the weight overlap between the pairs
#'
#' @export
r.am.pairwise.alteration.overlap <- function(null,n.permut,n.cores=1) {
    return(rcpp_overlap(v=null,t=1))
}
#' Compute weight overlap stats
#' 
#' @import doParallel
#' @import parallel
#' @param null The null model
#' @param W The weight matrix
#' @param n.permut The number of permutation steps
#' @param n.cores The number of cores
#' @return overlap the weight overlap between the pairs
#'
#' @export
w.r.am.pairwise.alteration.overlap <- function(null,W,n.permut,n.cores=1) {	  
    return(rcpp_w_overlap(v=null,w=W,t=1))
}
#' Compute weight overlap stats
#' 
#' @param x The alteration matrix
#' @return average
#'
#' @export
add <- function(x) Reduce("+", x)
#' Compute weight overlap stats
#' 
#' @param obs The observed overlap
#' @param exp The expected overlap
#' @return effect size
#'
#' @export
effectSize = function(obs, exp){
    es = (obs - exp)*sin(pi/4)
    return(es)
}
#' Compute weight overlap stats
#' 
#' @param null_overlap The null model overlap
#' @param mean_mat The mean effect size vector
#' @param n.permut The number of permutation steps
#' @param n.cores The number of cores
#' @return overlap the weight overlap between the pairs
#'
#' @export
r.effectSize<- function(null_overlap,mean_mat,n.permut=1000,n.cores=1){

    foreach::registerDoSEQ()
    i<-NULL
    r.effect <- foreach(i=1:n.permut) %do% as.vector(as.dist((null_overlap[[i]]-mean_mat)*sin(pi/4)))
    return(r.effect)

}
#' Compute yule coefficent
#' 
#' @param overlap The overlap matrixx
#' @param mat The gam
#' @return yule coefficent
#'
#' @export
binary.yule <- function(overlap,mat){
    overlap = as.matrix(overlap)
    marginal_f1 = replicate(ncol(overlap), diag(overlap))
    colnames(marginal_f1) = rownames(marginal_f1)
    marginal_f2 = t(marginal_f1)
    v_11 = overlap
    v_10 = marginal_f1 - v_11
    v_01 = marginal_f2 - v_11
    v_00 = nrow(mat) - v_11 - v_01 - v_10
    
    # Calculate Odds Ratio
    OR = (v_00 * v_11) / (v_10 * v_01)
    Yule = (sqrt(OR)-1) / (sqrt(OR)+1)
    return(Yule)
}

#' Compute FDR
#' 
#' @param obs The observed values
#' @param exp The expected values aka null model
#' @param nSim The weight matrix
#' @param maxFDR The maxFDR to keep
#' @return overlap the weight overlap between the pairs
#'
#' @export
estimateFDR2 = function(obs, exp, nSim, maxFDR = 0.25){

    all.fdr = rep(1, length(obs))
    orig.obs = obs
    obs = sort(obs, decreasing = TRUE)
    exp = sort(exp, decreasing = TRUE)
    obs.pos = 1
    exp.pos = 1
    scan = TRUE
    while(scan){
        if(obs.pos>length(obs)){
            scan=FALSE
        }
        else{
            value = obs[obs.pos]
            while(exp[exp.pos] >= value){
                exp.pos = exp.pos + 1
                if(exp.pos==length(exp))
                    break
            }
            fp = (exp.pos-1)/nSim
            fdr = min(fp/sum(obs >= value), 1)
            #print(paste('FP',exp.pos,'Test_value',value,sep=":"))
            #print(sum(obs >= value))
            all.fdr[ orig.obs == value ] = fdr
            if(fdr >= maxFDR)
                scan = FALSE
            obs.pos = obs.pos + 1
        }
    }
    #print('------------')
    return(all.fdr)
}

#' Compute pairwise p-value
#' 
#' @param robs_co The observed values
#' @param obs.co The expected values aka null model
#' @param gene1 The weight matrix
#' @param gene2 The maxFDR to keep
#' @return overlap the weight overlap between the pairs
#'
#' @export

estimate_p_val <- function(robs_co,obs.co,gene1,gene2){
	background <- c()
	for(obj in robs_co){
    	background<-c(background,obj[c(gene1,gene2),c(gene1,gene2)][1,2])
	}
	actual_ratio <- obs.co[c(gene1,gene2),c(gene1,gene2)][1,2]
	P <- ecdf(unlist(background))
    p_val = 2*min(P(actual_ratio) , 1 - P(actual_ratio))
    return(p_val)
}

#' Compute pairwise p-value
#' 
#' @param obs The observed values
#' @param exp The expected values aka null model
#' @param results The weight matrix
#' @param nSim The maxFDR to keep
#' @return overlap the weight overlap between the pairs
#'
#' @export
estimate_pairwise_p = function(obs, exp,results, nSim){     
    robs_co <-list()
	k=1
	for(obj in exp){
	    robs_co[[k]]<-obj
	    rownames(robs_co[[k]])<-rownames(obs)
	    colnames(robs_co[[k]])<-rownames(obs)
	    k=k+1
	}
	p_val <- c()
	for( row in rownames(results)){
		gene1 = results[row,'SFE_1']	
		gene2 = results[row,'SFE_2']
		p_val<-c(p_val,estimate_p_val(robs_co,obs,gene1,gene2))	
	}
	return(p_val)
}


#' Compute weight overlap stats
#' 
#' @param al The alteration matrix
#' @param als The alteration matrix stats object
#' @param obs The observed overlap
#' @param wobs The weighted observed overlap
#' @param r.obs The random observed overlap
#' @param r.wobs The random weighted observed overlap
#' @param null The null model
#' @param maxFDR The maxFDR cutoff
#' @param estimate_pairwise Compute pairwise or not
#' @param n.cores The number of cores
#' @param n.permut The number of permutation
#' @return The table of results
#' @export 
#'
interaction.table <- function(al,
							  als,
							  obs,
							  wobs,
							  r.obs=NULL,
							  r.wobs=NULL,
							  null,
							  maxFDR=0.2,
							  n.cores=1,
							  estimate_pairwise=FALSE,
							  n.permut=1000){
    features = rownames(al$am$full)
    num_features = length(features)
    result_cols = c('SFE_1', 'SFE_2')
    results = data.frame(matrix(NA, nrow=num_features*(num_features-1)/2, ncol=length(result_cols)))
    colnames(results) = result_cols
    temp = rep(1:length(features), length(features))
    dim(temp) = c(length(features), length(features))
    diag(temp) = NA
    temp1 = t(temp)
    temp1[which(temp1 > temp)] = NA
    temp[which(temp < temp1)] = NA
    temp = cbind(as.vector(temp1), as.vector(temp))
    temp = temp[which(!is.na(temp[,1]) & !is.na(temp[,2])),,drop=FALSE]
    results[, 'SFE_1'] = features[temp[,1]]
    results[, 'SFE_2'] = features[temp[,2]]
    rownames(results) = paste(results[,'SFE_1'], results[,'SFE_2'], sep= ' - ')
    results[, 'name'] = rownames(results)
    results[,'support_1'] = als$alteration.count[results[,'SFE_1']]
    results[,'support_2'] = als$alteration.count[results[,'SFE_2']]
    results[,'freq_1'] = results[,'support_1']/als$n.samples
    results[,'freq_2'] = results[,'support_2']/als$n.samples
    results[,'overlap'] = as.vector(as.dist(obs))
    results[,'w_overlap'] = as.vector(as.dist(wobs))
    # Get max possible overlap
    if(!is.null(als$sample.blocks)) {
        A = lapply(als$sample.blocks, function(x) {
            a = rep(x$alteration.count, length(x$alteration.count))
            dim(a) = c(length(x$alteration.count),length(x$alteration.count))
            rownames(a) = names(x$alteration.count)
            colnames(a) = names(x$alteration.count)
            b = t(a)
            c = pmin(a,b)
        })
        A = Reduce('+', A)
        results$max_overlap = as.vector(as.dist(A))
    }
    results[,'freq_overlap'] = results[,'overlap'] / (results[,'max_overlap'])       
    results[,'r_overlap'] = as.vector(as.dist(add(r.obs)/length(r.obs)))
    results[,'w_r_overlap'] = as.vector(as.dist(add(r.wobs)/length(r.wobs)))
	if(estimate_pairwise){
	        obs.co = as.matrix(als$alteration.pairwise$overlap)
	        results[,'pairwise_p']<-rep(1,nrow(results)) 
	        results[,'pairwise_p'] = estimate_pairwise_p(obs=obs.co,
	        											 exp = r.obs,
	        											 results=results,
	        											 nSim=n.permut
	    												)
	}
	exp.ES = (results[,'w_overlap'] - results[,'w_r_overlap'])*sin(pi/4)
    exp.ES = as.numeric(exp.ES) 
    results[,'wES']<-exp.ES
    exp.rES = r.effectSize (r.wobs,(add(r.wobs)/length(r.wobs)),n.permut=n.permut,n.cores=n.cores)
    #tic('compute FDR on wES')
        exp.ES = unlist(exp.rES)
        exp.ES = as.numeric(exp.ES)
        results[,'wFDR']<-rep(1,nrow(results)) 
        results[,'wFDR'] = estimateFDR2(obs = abs(results[,'wES']),
                                       exp = abs(exp.ES),
                                       maxFDR = maxFDR,
                                       nSim = n.permut)
    #toc()
    ES_mtx <- matrix(abs(unlist(exp.rES)), ncol = length(exp.rES), nrow = length(exp.rES[[1]]))
    mean_ES <- abs(rowSums(ES_mtx)/n.permut)
    n.ES = (abs(results[,'wES']) - abs(mean_ES))
    n.ES[ n.ES < 0 ] = 0
    results[,'nES']<- sign(results[,'wES']) * n.ES
    results[,'mean_r_nES']<- sign(results[,'wES']) * abs(mean_ES)    
    #tic('compute FDR on nES')
        exp.rES <- (abs(ES_mtx)-abs(mean_ES))
        exp.rES[ exp.rES < 0 ] = 0
        exp.rES<- sign(ES_mtx) * exp.rES
        exp.ES = unlist(exp.rES)
        exp.ES = as.numeric(exp.ES)
        results[,'nFDR']<-rep(1,nrow(results)) 
        results[,'nFDR'] = estimateFDR2(obs = abs(results[,'nES']),
                                       exp = abs(exp.ES),
                                       maxFDR = maxFDR,
                                       nSim = n.permut)
    #toc()
    #tic('compute FDR on nES (frequency specific)')
        results[,'cum_freq'] = results$support_1 + results$support_2
        samples=length(colnames(al$am$full))
        freq.cats = c(0, 0.02, 0.05, 0.10, 1)
        results$nFDR2 = rep(1, nrow(results))
        for(i in 2:length(freq.cats)){
            select = results[,'cum_freq'] > freq.cats[i-1]*2*samples & results[,'cum_freq'] < freq.cats[i]*2*samples
            if (sum(select)>1){
                results$nFDR2[select] = estimateFDR2(obs = abs(results$nES[select]),
                                                     exp = abs(as.numeric(unlist(exp.rES[select,]))), 
                                                     maxFDR = maxFDR,
                                                     nSim = n.permut)

            }
        }
    #toc()
    results = results[order(abs(results$nES), decreasing = T),]
    results$type = rep('ME', nrow(results))
    results$type[ results$nES > 0 ] = 'CO'
    results$FDR<-results$nFDR2<=maxFDR
    return (results)
}
