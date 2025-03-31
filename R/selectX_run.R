###
# Author  : Arvind Iyer
# Email   : arvind.iyer@unil.ch
# Project : SelectSim
# Desc    : Main file which call the selectX function to create alteration object with background model and funtion to generate the table.
# Version : 0.1.4
# Updates : Re wrote the whole code
# Todo:
# - Same ordering of colnames in input gam as the one used in null. 
# - Better Error message and running text
# - Edge case: When sample size in less than 2 there is error in computation (need to fix a number to do this analysis)
# - parallel::makeCluster(2, setup_strategy = "sequential") a possible fix to remove the erorr of not able to connect problem (https://github.com/rstudio/rstudio/issues/6692)
# - Try to fix the parallel processing issuses
###

###
#' SelectX main function to create alteration object with background model
#' 
#' @description
#' `selectX()` takes a list object which consist of genome alteration matrix and tumor mutation burden.
#'
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom  tictoc tic
#' @importFrom  tictoc toc
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
#' @param tao tao parameter.
#' @param folder folder path to store the results.
#' @param save.object store the SelectX object.
#' @param verbose print the time and each steps.
#' @param estimate_pairwise Compute pairwise p-value
#' @param maxFDR FDR value
#' @param seed a random seed
#' @return result  a SelectX object with background model and other info along with result table
#'  
#' @export                          
selectX <- function(M, #A list object consiting of GAM and TMB data
				    sample.class, # Sample Covariates as named list
				    alteration.class, #Alteration Covariates as named list
					n.cores = 1, #number of cores
				    min.freq=10,  # Min freq of gene to be mutated to do the analysis
				    n.permut = 1000, # number of simulations
    	            lambda=0.3, #wieght factor
    	            tao=1, #the fold change factor
    	            save.object = FALSE, #save the object
	                folder='./', # folder
				    verbose = TRUE, #verbose to print the steps
				    estimate_pairwise=FALSE, #compute pairwise p-val
				    maxFDR=0.25, #max fdr
				    seed=42 # a random seed,
				    ){
	## Set a global seed for computation
	set.seed(seed)
	if(verbose){
		tic('Total time taken:')
			print(paste('#### Creating SelectX object ####'))
			tic('Time:')
				print(paste('Step1-> Parsing and Filtering GAM...'))
					al = new.AL.general(M,alteration.class,sample.class,min.freq,verbose)
				print(paste('-> Alteration Landscape object created'))
			toc()	
			print(paste('Step2-> Generating Templetate object...'))
			tic('Time:')
					temp.data<-templeate.obj.gen(al)
				print(paste('-> Templetate object created'))
			toc()
			print(paste('Step3-> Generating penalty matrix...'))
			tic('Time:')
					W<-generateW_block(al,lambda,tao)
				print(paste('-> Penalty Matrix created'))
			toc()
			print(paste('Step4-> Generating null model...'))
			tic('Time:')
				sim <- null_model_parallel(al,temp.data$temp_mat,W$W,n.cores,n.permut)
		 		obj<- list('al'=al,'W'=W,'T'=temp.data,'null'=sim,'nSim'=n.permut)
		    	print(paste('-> Removing the outliers matrix from null model...'))
		    		outliers <- retrieveOutliers(obj = obj,nSim = n.permut)
				    if(length(which(outliers))==0){
				            print(paste('Removed null-matrix:',length(which(outliers)),sep=" "))
				            obj$nSim  <- n.permut
				    }
				    else{
				            print(paste(' Removed null-matrix:',length(which(outliers)),sep=" "))
				            obj$null <- obj$null[which(!outliers)]
				            print(paste(' Updated the null-model and nSim variables...'))
				            obj$nSim <- length(which(!outliers))
				    }
			toc()
			print(paste('-> Null model generated'))
			print(paste('### SelectX object created ###'))

			print(paste('#### Computing EDs on the dataset ####'))
			tic('Time:')
				al<-obj$al
			    als = al.stats(obj) # Univariate stats
			    alp = al.pairwise.alteration.stats(obj, als, do.blocks=FALSE)
			    als[['alteration.pairwise']] = alp
			    als[['alteration.pairwise']]$sample.blocks = NULL
			    for(i in names(als$sample.blocks)) als$sample.blocks[[i]][['alteration.pairwise']] = alp$sample.blocks[[i]]
				obs.co = as.matrix(als$alteration.pairwise$overlap)
	    		wobs.co = as.matrix(als$alteration.pairwise$w_overlap)
	    		robs.co<-r.am.pairwise.alteration.overlap(null = obj$null,
								 						  n.permut = obj$nSim,
                  		                              	  n.cores = n.cores)
    			wrobs.co<-w.r.am.pairwise.alteration.overlap(null = obj$null,
    											 			W= obj$W$W,
                                                 			n.permut = obj$nSim,
                                                 			n.cores = n.cores)

			   selectX_result <- interaction.table(al,
			                                        als,
			                                        obs.co,
			                                        wobs.co,
			                                        robs.co,
			                                        wrobs.co,
			                                        null=obj$null,
			                                        maxFDR=maxFDR,
			                                        n.cores=n.cores,
			                                        estimate_pairwise=estimate_pairwise,
			                                        n.permut=obj$nSim)
			   obj$robs.co <- robs.co
			   obj$wrobs.co <- wrobs.co
	    	toc()
	    	print(paste('#### EDs computed ####'))
		toc()
		if(save.object) {
			saveRDS(obj, file=paste(folder,'selectX_object.rds',sep=''))
			saveRDS(selectX_result, file=paste(folder,'selectX_results.rds',sep=''))
		}
		return(list(obj=obj,result=selectX_result))		
	}
	else{
		tic('Total Time')
			al = new.AL.general(M,alteration.class,sample.class,min.freq,verbose)
			temp.data<-templeate.obj.gen(al)
			W<-generateW_block(al,lambda,tao)
			sim <- null_model_parallel(al,temp.data$temp_mat,W$W,n.cores,n.permut)
		 	obj<- list('al'=al,'W'=W,'T'=temp.data,'null'=sim,'nSim'=n.permut)
		 	outliers <- retrieveOutliers(obj = obj,nSim = n.permut)
		 	if(length(which(outliers))==0){
		 		obj$nSim  <- n.permut
		 	}
		 	else{
		 		obj$null <- obj$null[which(!outliers)]
		 		obj$nSim <- length(which(!outliers))
		 	}
			al<-obj$al
		    als = al.stats(obj) # Univariate stats
		    alp = al.pairwise.alteration.stats(obj, als, do.blocks=FALSE)
		    als[['alteration.pairwise']] = alp
		    als[['alteration.pairwise']]$sample.blocks = NULL
		    for(i in names(als$sample.blocks)) als$sample.blocks[[i]][['alteration.pairwise']] = alp$sample.blocks[[i]]
			obs.co = as.matrix(als$alteration.pairwise$overlap)
    		wobs.co = as.matrix(als$alteration.pairwise$w_overlap)
    		robs.co<-r.am.pairwise.alteration.overlap(null = obj$null,
									 					    	   n.permut = obj$nSim,
                  		                           n.cores = 1)
			wrobs.co<-w.r.am.pairwise.alteration.overlap(null = obj$null,
						  											   W= obj$W$W,
                                             			n.permut = obj$nSim,
                                             			n.cores = 1)

		   selectX_result <- interaction.table(al,
		                                        als,
		                                        obs.co,
		                                        wobs.co,
		                                        robs.co,
		                                        wrobs.co,
		                                        null=obj$null,
		                                        maxFDR=maxFDR,
		                                        n.cores=1,
		                                        estimate_pairwise=estimate_pairwise,
		                                        n.permut=obj$nSim)
		   obj$robs.co <- robs.co
		   obj$wrobs.co <- wrobs.co
		toc()

		if(save.object) {
			saveRDS(obj, file=paste(folder,'selectX_object.rds',sep=''))
			saveRDS(selectX_result, file=paste(folder,'selectX_results.rds',sep=''))
		}
		return(list(obj=obj,result=selectX_result))
	}        
}



###
#' SelectX main function to create alteration object with background model
#' 
#' @description
#' `selectX()` takes a list object which consist of genome alteration matrix and tumor mutation burden.
#'
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom  tictoc tic
#' @importFrom  tictoc toc
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
#' @param tao tao parameter.
#' @param folder folder path to store the results.
#' @param save.object store the SelectX object.
#' @param verbose print the time and each steps.
#' @param estimate_pairwise Compute pairwise p-value
#' @param maxFDR FDR value
#' @param seed a random seed
#' @return result  a SelectX object with background model and other info along with result table
#'  
#' @export                          
selectX_debug <- function(M, #A list object consiting of GAM and TMB data
				    sample.class, # Sample Covariates as named list
				    alteration.class, #Alteration Covariates as named list
					 n.cores = 1, #number of cores
				    min.freq=10,  # Min freq of gene to be mutated to do the analysis
				    n.permut = 1000, # number of simulations
					lambda=0.3, #wieght factor
					tao=1, #the fold change factor
					save.object = FALSE, #save the object
					folder='./', # folder
				    verbose = TRUE, #verbose to print the steps
				    estimate_pairwise=FALSE, #compute pairwise p-val
				    maxFDR=0.25, #max fdr
				    seed=42 # a random seed,
				    ){
	## Set a global seed for computation
	set.seed(seed)
	if(verbose){
		tic('Total time taken:')
			print(paste('#### Creating SelectX object ####'))
			tic('Time:')
				print(paste('Step1-> Parsing and Filtering GAM...'))
					al = new.AL.general(M,alteration.class,sample.class,min.freq,verbose)
				print(paste('-> Alteration Landscape object created'))
			toc()	
			print(paste('Step2-> Generating Templetate object...'))
			tic('Time:')
					temp.data<-templeate.obj.gen(al)
				print(paste('-> Templetate object created'))
			toc()
			print(paste('Step3-> Generating penalty matrix...'))
			tic('Time:')
					W<-generateW_block(al,lambda,tao)
				print(paste('-> Penalty Matrix created'))
			toc()
			print(paste('Step4-> Generating null model...'))
			tic('Time:')
				sim <- null_model_parallel_debug(al,temp.data$temp_mat,W$W,n.cores,n.permut)
		 		obj<- list('al'=al,'W'=W,'T'=temp.data,'null'=sim,'nSim'=n.permut)
		    	print(paste('-> Removing the outliers matrix from null model...'))
		    		outliers <- retrieveOutliers(obj = obj,nSim = n.permut)
				    if(length(which(outliers))==0){
				            print(paste('Removed null-matrix:',length(which(outliers)),sep=" "))
				            obj$nSim  <- n.permut
				    }
				    else{
				            print(paste(' Removed null-matrix:',length(which(outliers)),sep=" "))
				            obj$null <- obj$null[which(!outliers)]
				            print(paste(' Updated the null-model and nSim variables...'))
				            obj$nSim <- length(which(!outliers))
				    }
			toc()
			print(paste('-> Null model generated'))
			print(paste('### SelectX object created ###'))

			print(paste('#### Computing EDs on the dataset ####'))
			tic('Time:')
				al<-obj$al
			    als = al.stats(obj) # Univariate stats
			    alp = al.pairwise.alteration.stats(obj, als, do.blocks=FALSE)
			    als[['alteration.pairwise']] = alp
			    als[['alteration.pairwise']]$sample.blocks = NULL
			    for(i in names(als$sample.blocks)) als$sample.blocks[[i]][['alteration.pairwise']] = alp$sample.blocks[[i]]
				obs.co = as.matrix(als$alteration.pairwise$overlap)
	    		wobs.co = as.matrix(als$alteration.pairwise$w_overlap)
	    		robs.co<-r.am.pairwise.alteration.overlap(null = obj$null,
								 							    	   n.permut = obj$nSim,
                  		                              n.cores = 1)
    			wrobs.co<-w.r.am.pairwise.alteration.overlap(null = obj$null,
    											 				 W= obj$W$W,
                                                 n.permut = obj$nSim,
                                                 n.cores = 1)

			   selectX_result <- interaction.table(al,
			                                        als,
			                                        obs.co,
			                                        wobs.co,
			                                        robs.co,
			                                        wrobs.co,
			                                        null=obj$null,
			                                        maxFDR=maxFDR,
			                                        n.cores=1,
			                                        estimate_pairwise=estimate_pairwise,
			                                        n.permut=obj$nSim)
			   obj$robs.co <- robs.co
			   obj$wrobs.co <- wrobs.co
	    	toc()
	    	print(paste('#### EDs computed ####'))
		toc()
		if(save.object) {
			saveRDS(obj, file=paste(folder,'selectX_object.rds',sep=''))
			saveRDS(selectX_result, file=paste(folder,'selectX_results.rds',sep=''))
		}
		return(list(obj=obj,result=selectX_result))		
	}
	else{
		tic('Total Time')
			al = new.AL.general(M,alteration.class,sample.class,min.freq,verbose)
			temp.data<-templeate.obj.gen(al)
			W<-generateW_block(al,lambda,tao)
			sim <- null_model_parallel_debug(al,temp.data$temp_mat,W$W,n.cores,n.permut)
		 	obj<- list('al'=al,'W'=W,'T'=temp.data,'null'=sim,'nSim'=n.permut)
		 	outliers <- retrieveOutliers(obj = obj,nSim = n.permut)
		 	if(length(which(outliers))==0){
		 		obj$nSim  <- n.permut
		 	}
		 	else{
		 		obj$null <- obj$null[which(!outliers)]
		 		obj$nSim <- length(which(!outliers))
		 	}
			al<-obj$al
		    als = al.stats(obj) # Univariate stats
		    alp = al.pairwise.alteration.stats(obj, als, do.blocks=FALSE)
		    als[['alteration.pairwise']] = alp
		    als[['alteration.pairwise']]$sample.blocks = NULL
		    for(i in names(als$sample.blocks)) als$sample.blocks[[i]][['alteration.pairwise']] = alp$sample.blocks[[i]]
			obs.co = as.matrix(als$alteration.pairwise$overlap)
    		wobs.co = as.matrix(als$alteration.pairwise$w_overlap)
    		robs.co<-r.am.pairwise.alteration.overlap(null = obj$null,
									 					    	   n.permut = obj$nSim,
                  		                           n.cores = 1)
			wrobs.co<-w.r.am.pairwise.alteration.overlap(null = obj$null,
						  											   W= obj$W$W,
                                             			n.permut = obj$nSim,
                                             			n.cores = 1)

		   selectX_result <- interaction.table(al,
		                                        als,
		                                        obs.co,
		                                        wobs.co,
		                                        robs.co,
		                                        wrobs.co,
		                                        null=obj$null,
		                                        maxFDR=maxFDR,
		                                        n.cores=1,
		                                        estimate_pairwise=estimate_pairwise,
		                                        n.permut=obj$nSim)
		   obj$robs.co <- robs.co
		   obj$wrobs.co <- wrobs.co
		toc()

		if(save.object) {
			saveRDS(obj, file=paste(folder,'selectX_object.rds',sep=''))
			saveRDS(selectX_result, file=paste(folder,'selectX_results.rds',sep=''))
		}
		return(list(obj=obj,result=selectX_result))
	}        
}

# ## usethis namespace: start
# #' @importFrom Rcpp sourceCpp
# ## usethis namespace: end
# NULL

# ## usethis namespace: start
# #' @useDynLib SelectSim, .registration = TRUE
# ## usethis namespace: end
# NULL