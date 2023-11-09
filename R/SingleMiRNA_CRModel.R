#'@title SingleMiRNA_CRModel
#'@description Function "SingleMiRNA_CRModel" uses survival data to build a multivariate Cox model using the targets of a single miRNA.
#'@param ExpData A gene expression profile of interest (rows are genes, columns are samples).
#'@param MiRNA A miRNA ID.
#'@param SurvivalData Survival data (the column names are: "sample", "status", "time") corresponding to samples in the gene expression profile of interest.
#'@param cutoff.point A numeric value is used to divide high-risk and low-risk groups.
#'@return A list includes a data frame with seven parts those are "sample", "status", "time", "target genes expression", "risk score", "group", and a dataframe with five columns those are "Gene", "HR", "HR.95L", "HR.95H", "beta", and "P-value".
#'@importFrom survival Surv
#'@importFrom survival coxph
#'@importFrom stats qnorm
#'@usage SingleMiRNA_CRModel(ExpData,MiRNA,cutoff.point=NULL,SurvivalData)
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-GetData_Mirna("GEP")
#' survival<-GetData_Mirna("survival")
#' # Run the function
#' SingleMiRNA_CRData<-SingleMiRNA_CRModel(GEP,
#' "hsa-miR-21-5p",cutoff.point=NULL,survival)



SingleMiRNA_CRModel<-function(ExpData,MiRNA,cutoff.point=NULL,SurvivalData){
  MiRNA_Target<- GetData_Mirna("MiRNA_Target")

    data_MiRNA_genes<-MiRNA_Target[which(MiRNA_Target[,1]==MiRNA),]
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

    outTab<-cbind(

      HR=coxsummary$coefficients[,"exp(coef)"],
      HR.95L=coxsummary$conf.int[,"lower .95"],
      HR.95H=coxsummary$conf.int[,"upper .95"],
      beta=coxsummary$coefficients[,"coef"],
      pvalue=coxsummary$coefficients[,"Pr(>|z|)"]


    )

    colnames(outTab)<-c("HR","HR.95L","HR.95H","beta","P-value")

    #  gene_exp<-a[,which(colnames(a)%in%"TGFBR2")]
    gene_exp<-survival_Intergrate[,which(colnames(survival_Intergrate)%in%rownames(xishu))]
    gene_exp<-as.matrix(gene_exp)
    #gene_beat<-xishu1[2]
    gene_beat<-as.matrix(xishu[,1])
    colnames(gene_beat)<-"score"
    gene_exp_beta<-gene_exp%*%gene_beat
    #gene_exp_beta<-gene_exp*gene_beat
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

      outTab<-cbind(

        HR=coxsummary$coefficients[,"exp(coef)"],
        HR.95L=coxsummary$conf.int[,"lower .95"],
        HR.95H=coxsummary$conf.int[,"upper .95"],
        beta=coxsummary$coefficients[,"coef"],
        pvalue=coxsummary$coefficients[,"Pr(>|z|)"]


      )
colnames(outTab)<-c("HR","HR.95L","HR.95H","beta","P-value")


     gene_exp<-survival_Intergrate[,which(colnames(survival_Intergrate)%in%gene)]
     xishu<-xishu[intersect(colnames(survival_Intergrate),rownames(xishu)),]
      gene_beat<-xishu[1]


      gene_exp_beta<-gene_exp*gene_beat
    }else if(index==0){
      print(paste0("No corresponding target gene of '",MiRNA,"' mirna was found in the difference list"))

    }
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
