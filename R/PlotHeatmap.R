#'@title PlotHeatmap
#'@description The function "PlotHeatmap" is used to plot a heat map of miRNA targets expression.
#'@param MiRNA_CRData A list includes a data frame with seven parts those are "sample", "status", "time", "target genes expression", "risk score", "group", and a data frame with five columns those are "Gene", "HR", "HR.95L", "HR.95H", "beta", and "P-value".
#'@param yaxis The upper and lower limits of this heat map.
#'@param scale character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none".
#'@param cluster.rows A logical value that represents whether row clustering is used.
#'@param cluster.cols A logical value that represents whether col clustering is used.
#'@param show.colnames This parameter controls whether column names are displayed.
#'@param ann_colors Vector of colors used to define groups.
#'@param col Vector of colors used in the heatmap.
#'@return A heat map of miRNA targets expression.
#'@importFrom pheatmap pheatmap
#'@importFrom stats wilcox.test
#'@importFrom grDevices colorRampPalette
#'@usage PlotHeatmap(MiRNA_CRData,yaxis=c(-2,2),scale="row",
#'cluster.rows=FALSE,cluster.cols=FALSE,show.colnames=FALSE,
#'ann_colors=c("#ef6d6d","#5470c6"),col=c("#ef6d6d","#5470c6"))
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-GetData_Mirna("GEP")
#' survival<-GetData_Mirna("survival")
#' MiRNAs<-c("hsa-miR-21-5p","hsa-miR-26a-5p","hsa-miR-369-5p","hsa-miR-1238-3p","hsa-miR-10b-5p")
#' # Run the function
#' SingleMiRNA_CRData<-SingleMiRNA_CRModel(GEP,
#' "hsa-miR-21-5p",survival,cutoff.point=NULL)
#' PlotHeatmap(SingleMiRNA_CRData)
#' MutiMiRNA_CRData<-MutiMiRNA_CRModel(GEP,
#' MiRNAs,survival,cutoff.point=NULL)
#' PlotHeatmap(MutiMiRNA_CRData)
PlotHeatmap<-function(MiRNA_CRData,yaxis=c(-2,2),scale="row",cluster.rows=FALSE,cluster.cols=FALSE,show.colnames=FALSE,ann_colors=c("#ef6d6d","#5470c6"),col=c("#ef6d6d","#5470c6")){
  riskresult<-data.frame(MiRNA_CRData[[1]])

  cox.result = data.frame(MiRNA_CRData[[2]])
  cox.result<-cox.result[order(cox.result$HR,decreasing=T),]
  rownames(cox.result)

  riskresult<-riskresult[order(riskresult[,"score"],decreasing=T),]
  t_riskresult<-t(riskresult[,3:dim(riskresult)[2]])
  t_riskresult<-cbind(t_riskresult[,which(t_riskresult[dim(t_riskresult)[1],]=="low")],t_riskresult[,which(t_riskresult[dim(t_riskresult)[1],]=="high")])
  group<-data.frame(t_riskresult[dim(t_riskresult)[1],])
  colnames(group)<-"group"
  t_riskresult<-t_riskresult[1:(dim(t_riskresult)[1]-2),]
  t_riskresult1<-apply(t_riskresult,2,as.numeric)
  rownames(t_riskresult1)<-rownames(t_riskresult)
  t_riskresult1<-t_riskresult1[rownames(cox.result),]
  rowns<-c()
  names<-c()

  for(i in 1:length(t_riskresult[,1])){
    test<-cbind(group,data.frame(t_riskresult1[i,]))
    colnames(test)[2]<-"value"
    test[,2]<-as.numeric(test[,2])
    t<-wilcox.test(value~group,var.equal=T,data=test)

    if(t$p.value<0.00001){
      rown<-paste0(rownames(t_riskresult1)[i]," (P<0.00001)")
    }else{
      rown<-paste0(rownames(t_riskresult1)[i]," (P= ",round(t$p.value,5),")")
    }
    rowns<-c(rowns,rown)
    names<-c(names,rownames(t_riskresult1)[i])


  }

  bk <- c(seq(yaxis[1],mean(yaxis),by=0.01),seq(mean(yaxis),yaxis[2],by=0.01))
  ann_colors=list(group=c(high=ann_colors[1],low=ann_colors[2]))
  ga<-pheatmap(t_riskresult1, cluster_rows = cluster.rows, cluster_cols = cluster.cols,annotation_col = group, color = c(colorRampPalette(colors = c(col[2],"white"))(length(bk)/2),colorRampPalette(colors = c("white",col[1]))(length(bk)/2)),
               breaks=seq(yaxis[1],yaxis[2],0.01),scale=scale,show_colnames = show.colnames,annotation_colors=ann_colors, labels_row = rowns )


}

