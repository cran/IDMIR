## ----include = FALSE----------------------------------------------------------
library(IDMIR)

## ----eval=TRUE----------------------------------------------------------------
# Obtain the example data
GEP<-GetData_Mirna("GEP")
label<-GetData_Mirna("label")
# Calculate the zscore
GDEscore<-GetGDEscore(GEP,label)
head(GDEscore)

## ----eval=FALSE---------------------------------------------------------------
#  # load depend package
#  require(igraph)
#  # Calculate the centrality score of miRNAs
#  MiRNAScore<-IdentifyMiRNA(de_kich,nperm = 10,damping = 0.9)
#  ###view the result
#  MiRNAScore[1:5,-2]
#  

## ----echo=FALSE---------------------------------------------------------------
###Get the result of this function
MiRNAScore<-GetData_Mirna("MiRNAScore")
MiRNAScore[1:5,-2]

## ----eval=TRUE----------------------------------------------------------------
# load depend package
require(survival)
###Get the survival and gene expression data
GEP<-GetData_Mirna("GEP")
survival<-GetData_Mirna("survival")
# Build a multivariate cox model
SingleMiRNA_CRData<-SingleMiRNA_CRModel(GEP,"hsa-miR-21-5p",cutoff.point=NULL,survival)
###view the result
data.frame(SingleMiRNA_CRData[[1]])[1:5,]
###view the result
data.frame(SingleMiRNA_CRData[[2]])[1:5,]


## ----eval=TRUE----------------------------------------------------------------
# load depend package
require(survival)
###Get the survival and gene expression data
GEP<-GetData_Mirna("GEP")
survival<-GetData_Mirna("survival")
MiRNAs<-c("hsa-miR-21-5p","hsa-miR-26a-5p","hsa-miR-369-5p","hsa-miR-1238-3p","hsa-miR-10b-5p")
# Build a multivariate cox model
MutiMiRNA_CRData<-MutiMiRNA_CRModel(GEP,MiRNAs,cutoff.point=NULL,survival)
###view the result
data.frame(MutiMiRNA_CRData[[1]])[1:5,]
###view the result
data.frame(MutiMiRNA_CRData[[2]])[1:5,]


## ----message=FALSE------------------------------------------------------------
# load depend package
require(survminer)
require(survival)

###Get the survival and gene expression data
GEP<-GetData_Mirna("GEP")
survival<-GetData_Mirna("survival")
# plot the Kaplan-Meier curve
SingleMiRNA_CRData<-SingleMiRNA_CRModel(GEP,"hsa-miR-21-5p",cutoff.point=NULL,survival)
PlotSurvival(SingleMiRNA_CRData)

## ----message=FALSE------------------------------------------------------------
# load depend package
require(pheatmap)
require(survival)
###Get the survival and gene expression data
GEP<-GetData_Mirna("GEP")
survival<-GetData_Mirna("survival")
# plot the heatmap of miRNA targets expression
SingleMiRNA_CRData<-SingleMiRNA_CRModel(GEP,"hsa-miR-21-5p",cutoff.point=NULL,survival)
PlotHeatmap(SingleMiRNA_CRData)

## ----message=FALSE------------------------------------------------------------
# load depend package
require(forestplot)
require(survival)
###Get the survival and gene expression data
GEP<-GetData_Mirna("GEP")
survival<-GetData_Mirna("survival")
# plot the forest diagram
SingleMiRNA_CRData<-SingleMiRNA_CRModel(GEP,"hsa-miR-21-5p",cutoff.point=NULL,survival)
PlotForest(SingleMiRNA_CRData)

## ----message=FALSE------------------------------------------------------------
# load depend package
require(egg)
require(survival)
require(ggplot2)
###Get the survival and gene expression data
GEP<-GetData_Mirna("GEP")
survival<-GetData_Mirna("survival")
# plot the scatter diagram
SingleMiRNA_CRData<-SingleMiRNA_CRModel(GEP,"hsa-miR-21-5p",cutoff.point=NULL,survival)
PlotScatter(SingleMiRNA_CRData)

