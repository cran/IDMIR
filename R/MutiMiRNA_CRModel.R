#'@title MutiMiRNA_CRModel
#'@description Function "MutiMiRNA_CRModel" can build a multivariate Cox model through integrating the models constructed separately based on different mirna targets.
#'@param ExpData A gene expression profile of interest (rows are genes, columns are samples).
#'@param MiRNAs An interest miRNA vector.
#'@param SurvivalData Survival data (the column names are: "sample", "status", and "time") corresponding to the samples in gene expression profile of interest.
#'@param cutoff.point A numeric value used to divide high-risk and low-risk groups.
#'@return A list includes a data frame with seven parts those are "sample", "status", "time", "target gene expression", "risk score", "group", and a data frame with five columns those are "Gene", "HR", "HR.95L", "HR.95H", "beta", and "P-value".
#'@importFrom survival Surv
#'@importFrom survival coxph
#'@importFrom stats qnorm
#'@usage MutiMiRNA_CRModel(ExpData, MiRNAs,SurvivalData,cutoff.point=NULL)
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-GetData_Mirna("GEP")
#' survival<-GetData_Mirna("survival")
#' MiRNAs<-c("hsa-miR-21-5p","hsa-miR-26a-5p","hsa-miR-369-5p","hsa-miR-1238-3p","hsa-miR-10b-5p")
#' # Run the function
#' MutiMiRNA_CRData<-MutiMiRNA_CRModel(GEP,
#' MiRNAs,survival,cutoff.point=NULL)


MutiMiRNA_CRModel<-function(ExpData,MiRNAs,SurvivalData,cutoff.point=NULL){


  MiRNA_Target<- GetData_Mirna("MiRNA_Target")

  riskresult_Intergrate<-data.frame(patient=intersect(colnames(ExpData),SurvivalData[,1]))

  for(i in 1:length(MiRNAs)){
    data_MiRNA_genes<-MiRNA_Target[which(MiRNA_Target[,1]==MiRNAs[i]),]
    data_MiRNA_gene<-strsplit(unlist(data_MiRNA_genes[2]),",")
    data_MiRNA_gene_exp<-ExpData[which(as.vector(rownames(ExpData))%in%as.vector(unlist(data_MiRNA_gene))),]
    index<-length(which(as.vector(rownames(ExpData))%in%as.vector(unlist(data_MiRNA_gene))))
    if(index>1){

      t_data_MiRNA_gene_exp<-data.frame(t(data_MiRNA_gene_exp))
      t_data_MiRNA_gene_exp<-cbind(rownames(t_data_MiRNA_gene_exp),t_data_MiRNA_gene_exp)
      colnames(t_data_MiRNA_gene_exp)[1]<-"sample"
      survival_Intergrate<-merge(SurvivalData,t_data_MiRNA_gene_exp,by.x=colnames(SurvivalData)[1],by.y="sample")

      rownames(survival_Intergrate)<-survival_Intergrate[,1]
      survival_Intergrate<-survival_Intergrate[,-1]

      cox=coxph(Surv(time,status)~.,survival_Intergrate)
      coxsummary = summary(cox)
      xishu<-coxsummary$coefficients

      #  gene_exp<-a[,which(colnames(a)%in%"TGFBR2")]
      gene_exp<-survival_Intergrate[,intersect(colnames(survival_Intergrate),rownames(xishu))]
      xishu<-xishu[intersect(colnames(survival_Intergrate),rownames(xishu)),]
      gene_exp<-as.matrix(gene_exp)
      #gene_beat<-xishu1[2]
      gene_beat<-as.matrix(xishu[,1])
      colnames(gene_beat)<-"score"
      gene_exp_beta<-gene_exp%*%gene_beat
      #gene_exp_beta<-gene_exp*gene_beat
      riskresult_single = cbind(patient = rownames(survival_Intergrate),gene_exp_beta)

      colnames(riskresult_single)[dim(riskresult_single)[2]]<-MiRNAs[i]
      riskresult_Intergrate<-merge(riskresult_Intergrate,riskresult_single,by.x="patient",by.y="patient")

    }else if(index==1){
      gene<-rownames(ExpData)[which(as.vector(rownames(ExpData))%in%as.vector(unlist(data_MiRNA_gene)))]
      t_data_MiRNA_gene_exp<-data.frame(t(data_MiRNA_gene_exp))
      t_data_MiRNA_gene_exp<-cbind(rownames(t_data_MiRNA_gene_exp),t_data_MiRNA_gene_exp)
      colnames(t_data_MiRNA_gene_exp)[1]<-"sample"
      survival_Intergrate<-merge(SurvivalData,t_data_MiRNA_gene_exp,by.x=colnames(SurvivalData)[1],by.y="sample")

      rownames(survival_Intergrate)<-survival_Intergrate[,1]
      survival_Intergrate<-survival_Intergrate[,-1]

      cox=coxph(Surv(time,status)~.,survival_Intergrate)
      coxsummary = summary(cox)
      xishu<-coxsummary$coefficients

      gene_exp<-survival_Intergrate[,which(colnames(survival_Intergrate)%in%gene)]
      xishu<-xishu[intersect(colnames(survival_Intergrate),rownames(xishu)),]

      gene_beat<-xishu[1]



      gene_exp_beta<-gene_exp*gene_beat
      riskresult_single = cbind(patient = rownames(survival_Intergrate),gene_exp_beta)

      colnames(riskresult_single)[dim(riskresult_single)[2]]<-MiRNAs[i]
      riskresult_Intergrate<-merge(riskresult_Intergrate,riskresult_single,by.x="patient",by.y="patient")

    }else if(index==0){
print(paste0("No corresponding target gene of ",MiRNAs[i]," mirna was found in the difference list"))
    }
      }
  TimeStatus<-cbind(rownames(survival_Intergrate),survival_Intergrate[,1:2])
  colnames(TimeStatus)[1]<-"patient"
  riskresult_Intergrate<-merge(TimeStatus,riskresult_Intergrate,by.x="patient",by.y="patient")

  rownames(riskresult_Intergrate)<-riskresult_Intergrate[,1]
  riskresult_Intergrate<-riskresult_Intergrate[,-1]
  result<-apply(riskresult_Intergrate,2,as.numeric)
  rownames(result)<-rownames(riskresult_Intergrate)
  result<- data.frame(result)



  cox=coxph(Surv(time,status)~.,result)
  coxsummary = summary(cox)
  xishu<-coxsummary$coefficients

  outTab<-cbind(
    HR=coxsummary$coefficients[,"exp(coef)"],
    HR.95L=coxsummary$conf.int[,"lower .95"],
    HR.95H=coxsummary$conf.int[,"upper .95"],
    beta=coxsummary$coefficients[,"coef"],
    pvalue=coxsummary$coefficients[,"Pr(>|z|)"]


  )
  colnames(outTab)<-c("HR","HR.95L","HR.95H","beta","P-value")


  #  gene_exp<-a[,which(colnames(a)%in%"TGFBR2")]
  gene_exp<-result[,which(colnames(result)%in%rownames(xishu))]
  gene_exp<-as.matrix(gene_exp)
  #gene_beat<-xishu1[2]
  gene_beat<-as.matrix(xishu[,1])
  colnames(gene_beat)<-"score"
  gene_exp_beta<-gene_exp%*%gene_beat
  #gene_exp_beta<-gene_exp*gene_beat
  if(is.null(cutoff.point)){
    cutoff.point<-median(gene_exp_beta)
  }
  risk = ifelse(gene_exp_beta>cutoff.point,"high","low")

  riskresult = cbind(patient = rownames(survival_Intergrate),survival_Intergrate[,1:2],gene_exp,gene_exp_beta, risk)
  riskresult<-riskresult[,-1]
  colnames(riskresult)[dim(riskresult)[2]]<-"group"
  resultlist<-c(list(riskresult),list(data.frame(outTab)))
  names(resultlist)<-c("riskresult","outTab")



  return(resultlist)
}
