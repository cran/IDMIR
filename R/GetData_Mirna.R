#' @title GetData_Mirna
#' @description Get the example data
#' @param Data  A character should be one of"survival", "GEP", "MF_Target", "MiRNA_Target", "matrix_mirna_go_inter", "matrix_mirna_go_jaccard"

#' @return data
#' @export

GetData_Mirna<-function(Data){
  if(!exists("MirnaData")) {
    utils::data("MirnaData",package="miRNA404")
  }
  if (Data=="survival")
  {
    dataset<- get("survival",envir=MirnaData)
    return(dataset)
  }
  if (Data=="GEP")
  {
    dataset<- get("GEP",envir=MirnaData)
    return(dataset)
  }
  if (Data=="MF_Target")
  {
    dataset<- get("MF_Target",envir=MirnaData)
    return(dataset)
  }
  if (Data=="MiRNA_Target")
  {
    dataset<- get("MiRNA_Target",envir=MirnaData)
    return(dataset)
  }
  if (Data=="matrix_mirna_go_inter")
  {
    dataset<- get("matrix_mirna_go_inter",envir=MirnaData)
    return(dataset)
  }
  if (Data=="matrix_mirna_go_jaccard")
  {
    dataset<- get("matrix_mirna_go_jaccard",envir=MirnaData)
    return(dataset)
  }
  if (Data=="label")
  {
    dataset<- get("label",envir=MirnaData)
    return(dataset)
  }
  if (Data=="MiRNAScore")
  {
    dataset<- get("MiRNAScore",envir=MirnaData)
    return(dataset)
  }
}
