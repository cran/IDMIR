#'@title PlotSurvival
#'@description Function "PlotSurvival" is used to draw a Kaplan-Meier curve.
#'@param MiRNA_CRData A list includes a data frame with seven parts those are "sample", "status", "time", "target genes expression", "risk score", "group", and a data frame with five columns those are "Gene", "HR", "HR.95L", "HR.95H", "beta", and "P-value".
#'@param colors Vector of colors used to define groups.
#'@return A survival curve of a data set.
#'@importFrom survival Surv
#'@importFrom survival survfit
#'@importFrom survminer ggsurvplot
#'@usage PlotSurvival(MiRNA_CRData,colors=c("#ef6d6d","#5470c6"))
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-GetData_Mirna("GEP")
#' survival<-GetData_Mirna("survival")
#' MiRNAs<-c("hsa-miR-21-5p","hsa-miR-26a-5p","hsa-miR-369-5p","hsa-miR-1238-3p","hsa-miR-10b-5p")
#' # Run the function
#' SingleMiRNA_CRData<-SingleMiRNA_CRModel(GEP,
#' "hsa-miR-21-5p",survival,cutoff.point=NULL)
#' PlotSurvival(SingleMiRNA_CRData)
#' MutiMiRNA_CRData<-MutiMiRNA_CRModel(GEP,
#' MiRNAs,survival,cutoff.point=NULL)
#' PlotSurvival(MutiMiRNA_CRData)
PlotSurvival<-function(MiRNA_CRData,colors=c("#ef6d6d","#5470c6")){
  data<-data.frame(MiRNA_CRData[[1]])
  kmfit1<-survfit(Surv(time,status)~group,data=data)

  ggsurvplot(kmfit1, data=data,
             pval =TRUE,
             xlab="Time",
             title="Survival Curve",
             palette=colors
  )

}
