###
# Author  : Arvind Iyer
# Email   : arvind.iyer@unil.ch
# Project : SelectX
# Desc    : The file contains plot related functions
# Version : 0.1
# Updates : Re wrote the whole code
# Todo:
# - Oncoprint & other usefull plot functions
###


#' Create an AL object
#'
#' Create an Alteration Landscape (AL) object which contains gams and mutation burden of samples of associated gams.
#'
#' @import ggpubr
#' @import ggplot2
#' @importFrom  dplyr %>%
#' @importFrom  dplyr mutate
#' @importFrom  dplyr case_when
#' @param result  result table of selectX obj$result
#' @param title   title of the plot
#' @return a ggplot2 object with observed vs random overlap plot. 
#' @export


obs_exp_scatter <- function(result,title){

	result$log_overlap <- log10(result$w_overlap+1)
	result$log_r_overlap <- log10(result$w_r_overlap+1)
	result <- result %>% mutate(cat=case_when(FDR==TRUE ~ type, FDR==FALSE ~ 'NS'))

	plot <- ggscatter(result, 
	          x = "log_r_overlap",
	          y = "log_overlap",
	          color='cat',
	          palette=c('forestgreen','purple','grey93'),
	          repel=TRUE,
	          alpha=1
	          )+xlab('Random weighted co-mutation (log10)')+ylab('Actual weighted co-mutation(log10)')+geom_abline(slope = 1,intercept = 0)+
	scale_y_continuous(breaks = seq(0, max(result$log_overlap), 1))+
    scale_x_continuous(breaks= seq(0, max(result$log_overlap), 1))

	plot <- plot + ggtitle(title)+theme_pubr()+
	theme(	legend.position = "top",
			legend.title = element_blank(),
			plot.title =  element_text(hjust = 0.5),
			text = element_text(size = 16))

	return(plot)
}


#' Extract the backgroun distribution
#'
#' Create an Alteration Landscape (AL) object which contains gams and mutation burden of samples of associated gams.
#'
#' @param gene1  gene1
#' @param gene2  gene2
#' @param obj    selectX obj
#' @return a backgorund distribution vector. 
#' @export

overlap_pair_extract <- function(gene1,gene2,obj){
    wes_dist <- c()
    panel_dist <- c()
    for(i in c(1:obj$nSim)){
        rownames(obj$wrobs.co[[i]])<-rownames(obj$al$am$full)
        colnames(obj$wrobs.co[[i]])<-rownames(obj$al$am$full)
        wes_dist <- c(wes_dist,obj$wrobs.co[[i]][gene1,gene2])
    }
    return(wes_dist)
}


#' Generate a pairs background plot
#'
#' Create an Alteration Landscape (AL) object which contains gams and mutation burden of samples of associated gams.
#'
#' @import ggpubr
#' @import ggplot2
#' @import ggridges
#' @import graphics
#' @importFrom  dplyr %>%
#' @importFrom  reshape2 melt
#' @importFrom  dplyr case_when
#' @importFrom  dplyr group_by
#' @importFrom  dplyr ungroup
#' @importFrom  dplyr mutate
#' @importFrom  dplyr row_number
#' @importFrom  dplyr filter
#' @importFrom  dplyr left_join
#' @importFrom  dplyr arrange
#' @importFrom  stats as.dist
#' @importFrom  stats density
#' @importFrom  stats ecdf
#' @importFrom  stats runif
#' @importFrom  stats setNames
#' @param result_df  result table of selectX obj$result
#' @param obj    selectX obj
#' @return a ggplot2 object with observed vs random overlap plot. 
#' @export

ridge_plot_ed <- function(result_df,obj){
    pair_mat <- matrix(0,nrow = nrow(result_df),ncol=obj$nSim)
    for(row in c(1:nrow(result_df))){
        gene1<-(result_df[row,'SFE_1'])
        gene2<-(result_df[row,'SFE_2'])
        pair<-(result_df[row,'name'])
        pair_mat[row,]<-overlap_pair_extract(gene1,gene2,obj)
    }
    pair_mat_plot <- t(pair_mat)
    colnames(pair_mat_plot)<-result_df$name
    df <- setNames(reshape2::melt(pair_mat_plot), c('rows', 'pairs', 'freq'))
    p<-ggplot(df, aes(x=freq, y=pairs)) +
             xlab('Co-mutated samples')+
             geom_density_ridges_gradient(quantile_lines = FALSE, quantiles = 2,color='black',fill='lightgrey') +
             theme_ridges(grid = TRUE, center = TRUE)+scale_x_continuous(breaks = seq(0, max(df$freq), by = 2))
    
    q <- ggplot_build(p)$data[[1]]
    density_lines <- q %>%
                      group_by(group) %>% 
                      dplyr::filter(density == max(density)) %>% 
                      ungroup()
    temp <- df %>% distinct(pairs) %>% mutate(number = row_number())
    density_lines_complete <- left_join(density_lines, temp, by = c("group" = "number"))
    density_lines_complete$x_actual<-result_df$w_overlap

    a <- ifelse(result_df$type == 'ME', "purple", "forestgreen")    
    plot<- ggplot(df ,
                 aes(x=freq,
                     y=pairs))+
                xlab('Weighted overlap')+ylab('')+
                geom_density_ridges_gradient(quantile_lines = FALSE, quantiles = 2,color='black',fill='lightgrey') +
                theme_ridges(grid = TRUE, center = TRUE,font_size = 18)+ scale_x_continuous(breaks = seq(0, max(df$freq), by = 5))+
                geom_segment(data = density_lines_complete, aes(x = x_actual, xend = x_actual, y = ymin, yend = ymin+density*scale*iscale,color="Actual overlap"))+
                geom_segment(data = density_lines_complete, aes(x = x, xend = x, y = ymin, yend = ymin+density*scale*iscale,color="Mean Background overlap"))+
                scale_color_manual("Legend",values = c("red","blue"))+
                theme(axis.text.y = element_text(colour = a))
    
  return(plot)
}

#' Generate a pairs background plot for two dataset comaprision
#'
#' Create an Alteration Landscape (AL) object which contains gams and mutation burden of samples of associated gams.
#'
#' @import ggpubr
#' @import ggplot2
#' @import ggridges
#' @import graphics
#' @importFrom  dplyr %>%
#' @importFrom  reshape2 melt
#' @importFrom  dplyr case_when
#' @importFrom  dplyr group_by
#' @importFrom  dplyr ungroup
#' @importFrom  dplyr mutate
#' @importFrom  dplyr row_number
#' @importFrom  dplyr filter
#' @importFrom  dplyr left_join
#' @importFrom  dplyr arrange
#' @importFrom  stats as.dist
#' @importFrom  stats density
#' @importFrom  stats ecdf
#' @importFrom  stats runif
#' @importFrom  stats setNames
#' @param result_df  a common result table of selectX obj$result with dataset1_w_overlap and dataset2_w_overlap columns
#' @param obj1    selectX obj of dataset1
#' @param obj2    selectX obj of dataset2
#' @param name1    dataset1 name
#' @param name2    dataset2 name
#' @return a ggplot2 object with observed vs random overlap plot. 
#' @export

ridge_plot_ed_compare <- function(result_df,obj1,obj2,name1,name2){
    pair_mat <- matrix(0,nrow = nrow(result_df),ncol=obj1$nSim)
    for(row in c(1:nrow(result_df))){
        gene1<-(result_df[row,'SFE_1'])
        gene2<-(result_df[row,'SFE_2'])
        pair<-(result_df[row,'name'])
        pair_mat[row,]<-overlap_pair_extract(gene1,gene2,obj1)
    }
    pair_mat_plot <- t(pair_mat)
    colnames(pair_mat_plot)<-result_df$name
    df1 <- setNames(reshape2::melt(pair_mat_plot), c('rows', 'pairs', 'freq'))
    
    pair_mat <- matrix(0,nrow = nrow(result_df),ncol=obj2$nSim)
    for(row in c(1:nrow(result_df))){
        gene1<-(result_df[row,'SFE_1'])
        gene2<-(result_df[row,'SFE_2'])
        pair<-(result_df[row,'name'])
        pair_mat[row,]<-overlap_pair_extract(gene1,gene2,obj2)
    }
    pair_mat_plot <- t(pair_mat)
    colnames(pair_mat_plot)<-result_df$name
    df2 <- setNames(reshape2::melt(pair_mat_plot), c('rows', 'pairs', 'freq'))
    
    df1$name<-rep(name1,nrow(df1))
    df2$name<-rep(name2,nrow(df2))
    
    df <- rbind(df1,df2)
    
    p<-ggplot(df, aes(x=freq, y=pairs,fill=name)) +
             xlab('Co-mutated samples')+
             geom_density_ridges_gradient(quantile_lines = FALSE, quantiles = 2,color='black') + scale_fill_manual(values = c("#E69F00", "#56B4E9"))+
             theme_ridges(grid = TRUE, center = TRUE)+scale_x_continuous(breaks = seq(0, max(df$freq), by = 5))
    
    q <- ggplot_build(p)$data[[1]]
    density_lines <- q %>% group_by(fill,y) %>% 
                      dplyr::filter(density == max(density)) %>% arrange(fill) %>%
                      ungroup()
    
    density_lines$pairs <- c(result_df$name,result_df$name)
    density_lines$x_actual <- c(result_df$dataset1_w_overlap,result_df$dataset2_w_overlap)
    density_lines$name <- c(rep(name2,nrow(result_df)),rep(name1,nrow(result_df)))
    a <- ifelse(result_df$type == 'ME', "purple", "forestgreen")    
    plot<- ggplot(df,
                 aes(x=freq,
                     y=pairs,
                     fill=name),alpha = 0.2)+
                xlab('Weighted overlap')+ylab('')+
                geom_density_ridges_gradient(quantile_lines = FALSE, quantiles = 2) +
                theme_ridges(grid = TRUE, center = TRUE,font_size = 18)+ scale_x_continuous(breaks = seq(0, max(df$freq), by = 5))+
                geom_segment(data = density_lines, aes(x = x_actual, xend = x_actual, y = ymin, yend = ymin+density*scale*iscale,color=name,linetype='Actual overlap'))+
                geom_segment(data = density_lines, aes(x = x, xend = x, y = ymin, yend = ymin+density*scale*iscale,color=name,linetype='Mean Background overlap'))+
                scale_color_manual("Line Color",values = c("red", "blue"))+
                scale_linetype_manual("Linetype",values = c("Actual overlap"=1,"Mean Background overlap"=2))+
                scale_fill_manual("Dataset",values = c("#E69F00", "#56B4E9"))+
                theme(axis.text.y = element_text(colour = a))
    
  return(plot)

}