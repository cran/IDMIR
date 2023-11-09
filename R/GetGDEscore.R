#'@title GetGDEscore
#'@description Function "GetGDEscore" is used to calculate gene differential expression levels.
#'@param ExpData A gene expression profile of interest (rows are genes, columns are samples).
#'@param Label A character vector consists of "0" and "1" which represent sample class in the gene expression profile. "0" means normal sample and "1" means disease sample.
#'@return A matrix with one column of GDEscore.
#'@importFrom stats pt
#'@importFrom stats qnorm
#'@usage GetGDEscore(ExpData,Label)
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-GetData_Mirna("GEP")
#' label<-GetData_Mirna("label")
#' # Run the function
#' GDEscore<-GetGDEscore(GEP,label)

GetGDEscore <- function(ExpData, Label){



  Label<-as.character(Label)

  Ttest<-matrix(nrow=nrow(ExpData),ncol=1)
  rownames(Ttest)<-rownames(ExpData)


  zscore.vector<-apply(ExpData, 1, function(x){
    ind1 <- which(Label == 1)
    ind2 <- which(Label == 0)
    m <- length(ind1 <- which(Label == 1))
    n <- length(ind2 <- which(Label == 0))
    expdata1 <- x[ind1]
    expdata2 <- x[ind2]
    rmean1 <- mean(expdata1)
    rmean2 <- mean(expdata2)
    ss1 <- sum((expdata1 - rmean1)^2)
    ss2 <- sum((expdata2 - rmean2)^2)
    Tscore <- (m + n - 2)^0.5*(rmean1 - rmean2)/((1/m + 1/n)*(ss1 + ss2))^0.5
    Pvalue <- pt(abs(Tscore),lower.tail=FALSE,df=m+n-2)
    Zvalue <- qnorm(Pvalue, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
    return(Zvalue)
  })
  Ttest[,1]<-zscore.vector
  Ttest <- as.matrix(Ttest[!is.infinite(Ttest[,1]),])
  Ttest <- as.matrix(Ttest[!is.na(Ttest[,1]),])
  Ttest <- as.matrix(Ttest[!is.nan(Ttest[,1]),])
  colnames(Ttest)<-c("GDEscore")
  return(Ttest)
}

