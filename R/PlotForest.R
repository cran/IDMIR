#'@title PlotForest
#'@description Function "PlotForest" can visualize the result of Cox regression analysis through forest plot.
#'@param MiRNA_CRData A list includes a data frame with seven parts those are "sample", "status", "time", "target genes expression", "risk score", "group", and a data frame with five columns those are "Gene", "HR", "HR.95L", "HR.95H", "beta", and "P-value".
#'@param g.pos The position of the graph element within the table of text. The position can be 1-(ncol(labeltext) + 1). You can also choose set the position to "left" or "right".
#'@param b.size Override the default box size based on precision.
#'@param col Set the colors for all the elements in the plot.
#'@param lwd.zero lwd for the vertical line that gives the no-effect line, see gpar.
#'@param lwd.ci lwd for the confidence bands, see gpar.
#'@param x.lab x-axis label.
#'@return Forest maps associated with the Cox risk model.
#'@importFrom forestplot forestplot
#'@importFrom forestplot fpColors
#'@importFrom forestplot fpTxtGp
#'@importFrom grid gpar
#'@importFrom stats qnorm
#'@usage PlotForest(MiRNA_CRData,g.pos = 2,b.size = 3,col = c("#FE0101", "#1C61B6", "#A4A4A4"),
#'lwd.zero = 2,lwd.ci = 3,x.lab = "Hazard Ratio Plot")
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-GetData_Mirna("GEP")
#' survival<-GetData_Mirna("survival")
#' MiRNAs<-c("hsa-miR-21-5p","hsa-miR-26a-5p","hsa-miR-369-5p","hsa-miR-1238-3p","hsa-miR-10b-5p")
#' # Run the function
#' SingleMiRNA_CRData<-SingleMiRNA_CRModel(GEP,
#' "hsa-miR-21-5p",survival,cutoff.point=NULL)
#' PlotForest(SingleMiRNA_CRData)
#' MutiMiRNA_CRData<-MutiMiRNA_CRModel(GEP,
#' MiRNAs,survival,cutoff.point=NULL)
#' PlotForest(MutiMiRNA_CRData)

PlotForest<-function (MiRNA_CRData, g.pos = 2, b.size = 3,
                       col = c("#FE0101", "#1C61B6", "#A4A4A4"),
                       lwd.zero = 2, lwd.ci = 3, x.lab = "Hazard Ratio Plot") {

  cox.result = data.frame(MiRNA_CRData[[2]])
  cox.result<-cox.result[order(cox.result$HR,decreasing=T),]
  HR = cox.result[, 1:3]
    hr = round(cox.result[, "HR"],3)
    hrLow = round(cox.result[, "HR.95L"],3)
    hrHigh = round(cox.result[, "HR.95H"],3)
    beta <- round(cox.result[, "beta"],3)
    pVal = round(cox.result[, "P.value"],3)
    tabletext <- list(c("ID", rownames(cox.result)),
                      append("HR",hr),
                      append("95% CI for HR", paste(hrLow, "-", hrHigh,sep = "")),
                      append("Cox.Beta", beta) )
    boxsize <- c(NA, (as.numeric(cox.result$HR))/b.size)
    forestplot(tabletext, rbind(rep(NA, 3), HR), graph.pos = g.pos,
               col = fpColors(box = col[1], lines = col[2], zero = col[3]),
               xlog = F, zero = 1, lwd.zero = lwd.zero, lwd.ci = lwd.ci,
               boxsize = boxsize, xlab = x.lab, ci.vertices = T, ci.vertices.height = 0.2,
               hrzl_lines = list(`2` = gpar(lwd = 2, col = "gray")),
               txt_gp = fpTxtGp(label = gpar(cex = 1), ticks = gpar(cex = 0.8),
                                xlab = gpar(cex = 1.5), ))
  }


