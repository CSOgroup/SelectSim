###
# Author  : Arvind Iyer
# Email   : arvind.iyer@unil.ch
# Project : SelectSim
# Desc    : The file which contains the function to generate the stats and table
# Version : 0.1.5
# Updates : Re wrote the whole code
# Todo:
###

#' Initialize an Alteration Landscape Stats (ALS) container
#'
#' @param al The alteration landscape object (checked for NULL)
#' @return Empty ALS list object.
#'
#' @export
new.ALS <- function(al) {
    if(is.null(al)) stop('Input al is NULL')
    als = list()
    class(als) <- "ALS";
    return(als)
}
#' Initialize an Alteration Matrix Stats (AMS) container
#'
#' @param am The alteration matrix (checked for NULL)
#' @return Empty AMS list object.
#'
#' @export
new.AMS <- function(am) {
    if(is.null(am)) stop('Input am is NULL')
    ams = list()
    class(ams) <- "AMS";
    return(ams)
}
#' Compute summary statistics for a binary alteration matrix
#'
#' @param am Binary alteration matrix (features x samples)
#' @return AMS object with basic counts: n.samples, n.alterations, n.occurrences, per-sample and per-feature counts.
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
#' Compute alteration landscape statistics
#'
#' @param al SelectX object (list containing al, W, etc. as returned by selectX)
#' @return ALS object with overall and per-block alteration statistics.
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
#' Compute TMB-weighted pairwise alteration overlap
#'
#' @importFrom  Matrix sparseMatrix
#' @importFrom  Matrix t
#' @param am Binary alteration matrix (features x samples)
#' @param W Weight matrix (features x samples) with per-sample TMB weights
#' @return Weighted pairwise overlap matrix (features x features).
#'
#' @export
am.weight.pairwise.alteration.overlap <- function(am,W) {

     A = am * 1
     col_order= colnames(A)
     overlap = (W[,col_order]*A) %*% t(A)
     return(overlap)
}
#' Compute pairwise alteration coverage statistics
#'
#' @importFrom  Matrix Matrix
#' @param overlap_M The pairwise overlap matrix
#' @param M.stats The alteration matrix stats (from am.stats)
#' @param w_overlap_M The weighted pairwise overlap matrix
#' @return List with overlap and w_overlap sparse matrices.
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
    stats[['overlap']] = Matrix(overlap_M)
    stats[['w_overlap']] = Matrix(w_overlap_M)
    return(stats)
}
#' Compute pairwise alteration statistics for an alteration landscape
#'
#' @param al SelectX object (list containing al, W, etc.)
#' @param als Alteration landscape stats (from al.stats); computed internally if NULL
#' @param do.blocks Whether to also compute block-level pairwise stats
#' @return List with overlap and w_overlap matrices, plus optional sample.blocks entries.
#'
#' @export
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
#' Compute null overlap matrix
#'
#' @param null The null model (list of simulated binary matrices)
#' @param n.permut The number of permutation steps
#' @param n.cores The number of cores
#' @return overlap the overlap summed across permutations
#'
#' @export
r.am.pairwise.alteration.overlap <- function(null,n.permut,n.cores=1) {
    return(rcpp_overlap(v=null,t=1))
}
#' Compute null weighted overlap matrix
#'
#' @param null The null model (list of simulated binary matrices)
#' @param W The weight matrix
#' @param n.permut The number of permutation steps
#' @param n.cores The number of cores
#' @return weighted overlap summed across permutations
#'
#' @export
w.r.am.pairwise.alteration.overlap <- function(null,W,n.permut,n.cores=1) {	  
    return(rcpp_w_overlap(v=null,w=W,t=1))
}
#' Sum a list of matrices element-wise
#'
#' @param x List of matrices of identical dimensions
#' @return Single matrix that is the element-wise sum of all matrices in x
#'
#' @export
add <- function(x) Reduce("+", x)
#' Compute effect size between observed and expected overlap
#'
#' @param obs The observed overlap values
#' @param exp The expected (null model mean) overlap values
#' @return Effect size value(s)
#'
#' @export
effectSize = function(obs, exp){
    es = (obs - exp)*sin(pi/4)
    return(es)
}
#' Compute effect sizes for null model permutations
#'
#' @param null_overlap List of null model overlap matrices (one per permutation)
#' @param mean_mat Mean overlap matrix across all permutations
#' @param n.permut Number of permutations
#' @param n.cores Number of cores (currently unused; sequential only)
#' @return List of effect size vectors, one per permutation.
#'
#' @export
r.effectSize<- function(null_overlap,mean_mat,n.permut=1000,n.cores=1){

    foreach::registerDoSEQ()
    i<-NULL
    r.effect <- foreach(i=1:n.permut) %do% as.vector(as.dist((null_overlap[[i]]-mean_mat)*sin(pi/4)))
    return(r.effect)

}
#' Compute Yule Q coefficient for all gene pairs
#'
#' @param overlap The pairwise overlap matrix
#' @param mat The binary GAM (features x samples)
#' @return Matrix of Yule Q coefficients
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
    v_00 = ncol(mat) - v_11 - v_01 - v_10
    
    # Calculate Odds Ratio
    OR = (v_00 * v_11) / (v_10 * v_01)
    Yule = (sqrt(OR)-1) / (sqrt(OR)+1)
    return(Yule)
}

#' Estimate FDR by scanning observed vs null effect sizes
#'
#' @param obs Vector of observed effect sizes
#' @param exp Vector of null model effect sizes (all permutations concatenated)
#' @param nSim Number of permutations used to generate exp
#' @param maxFDR FDR cutoff; scanning stops once FDR exceeds this value
#' @return Vector of FDR values, one per element of obs.
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
            all.fdr[ orig.obs == value ] = fdr
            if(fdr >= maxFDR)
                scan = FALSE
            obs.pos = obs.pos + 1
        }
    }
    #print('------------')
    return(all.fdr)
}

#' Compute empirical two-sided p-value for a gene pair
#'
#' @param robs_co List of null model overlap matrices (one per permutation)
#' @param obs.co Observed pairwise overlap matrix
#' @param gene1 Name of the first gene/alteration
#' @param gene2 Name of the second gene/alteration
#' @return Two-sided empirical p-value
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

#' Compute p-values for all gene pairs in a results table
#'
#' @param obs Observed pairwise overlap matrix
#' @param exp List of null model overlap matrices (one per permutation)
#' @param results Results data frame with SFE_1 and SFE_2 columns
#' @param nSim Number of permutations
#' @return Vector of p-values, one per row in results.
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


#' Build the full interaction results table from selectX outputs
#'
#' @param al Alteration landscape object (al$am, etc.)
#' @param als Alteration landscape stats (from al.stats)
#' @param obs Observed pairwise overlap matrix
#' @param wobs Weighted observed pairwise overlap matrix
#' @param r.obs List of null model overlap matrices
#' @param r.wobs List of null model weighted overlap matrices
#' @param null Null model list (simulated binary matrices)
#' @param maxFDR FDR cutoff for calling significant results
#' @param estimate_pairwise Whether to compute per-pair empirical p-values
#' @param n.cores Number of cores
#' @param n.permut Number of permutations
#' @return Data frame with one row per gene pair and columns for overlap, effect sizes, FDR, and interaction type.
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
        exp.ES = unlist(exp.rES)
        exp.ES = as.numeric(exp.ES)
        results[,'wFDR']<-rep(1,nrow(results)) 
        results[,'wFDR'] = estimateFDR2(obs = abs(results[,'wES']),
                                       exp = abs(exp.ES),
                                       maxFDR = maxFDR,
                                       nSim = n.permut)
    ES_mtx <- matrix(abs(unlist(exp.rES)), ncol = length(exp.rES), nrow = length(exp.rES[[1]]))
    mean_ES <- abs(rowSums(ES_mtx)/n.permut)
    n.ES = (abs(results[,'wES']) - abs(mean_ES))
    n.ES[ n.ES < 0 ] = 0
    results[,'nES']<- sign(results[,'wES']) * n.ES
    results[,'mean_r_nES']<- sign(results[,'wES']) * abs(mean_ES)    
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
    results = results[order(abs(results$nES), decreasing = T),]
    results$type = rep('ME', nrow(results))
    results$type[ results$nES > 0 ] = 'CO'
    results$FDR<-results$nFDR2<=maxFDR
    return (results)
}
