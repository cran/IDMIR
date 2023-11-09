#'@title IdentifyMiRNA
#' @description The function "IdentifyMiRNA" is used to identify significantly dysregulated miRNAs by calculating the eigenvector centrality of miRNAs.
#'@param GDEscore.table A matrix with one column of GDEscore.
#'@param nperm The Number of random permutations (default: 100).
#'@param damping Restart the probability of the random-walk algorithm (default: 0.9).
#'@return A data frame with seven columns those are "MiRNA", "Target", "Number" (number of targets), "Score" (Centrality score), "P-value", and "FDR".
#' @importFrom igraph graph.adjacency
#' @importFrom igraph V
#' @importFrom igraph page.rank
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @importFrom stats pt
#' @importFrom stats sd
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @importFrom stats na.omit
#' @importFrom fastmatch fmatch
#'@usage IdentifyMiRNA(GDEscore.table,nperm=1000,damping=0.90)
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-GetData_Mirna("GEP")
#' label<-GetData_Mirna("label")
#' # Run the function
#' GDEscore<-GetGDEscore(GEP,label)
#' \donttest{MiRNAScore<-IdentifyMiRNA(GDEscore,nperm=5, damping=0.90)}



IdentifyMiRNA<-function(GDEscore.table,nperm=1000,damping=0.90){
  random_network<-function(random,damping){
    adj.final1<-as.matrix(random)
    graph1 = graph.adjacency(adj.final1,mode=c("undirected"), weighted=TRUE,add.rownames=T)
    temp1 = page.rank(graph1, vids=V(graph1), directed=FALSE, damping=damping, weights=NULL)
    rank1 = temp1$vector
    rank2 = as.matrix(rank1)
    return(rank2)
  }
  go_mirna_score_row<-function(a,table,del,c){
    score_row<-rep(0,c)

    for(j in 1:length(a)){
      if(a[j]!=""){
        gene<-unlist(strsplit(a[j], split = ","))
        location<-fastmatch::fmatch(gene, table)
        dell<- na.omit(del[location])
        de_score1<-median(as.numeric(dell))
        if (!is.na(de_score1)) {
          score_row[j]<-de_score1
        }

      }

    }
    return(score_row)
  }
  matrix_mirna_go_inter<- GetData_Mirna("matrix_mirna_go_inter")
  matrix_mirna_go_jaccard<-GetData_Mirna("matrix_mirna_go_jaccard")
  GDEscore<-cbind(rownames(GDEscore.table),abs(GDEscore.table[,"GDEscore"]))
  median_inter<-function(GDEscore){

    table <- GDEscore[,1]
    del <- GDEscore[, 2]

    median_score<-matrix(0,nrow=length(rownames(matrix_mirna_go_inter)),ncol=length(colnames(matrix_mirna_go_inter)))
    for(k in 1:length(rownames(matrix_mirna_go_inter))){
      Genes_vector<-matrix_mirna_go_inter[k,]
      row<-go_mirna_score_row(Genes_vector,table,del,length(colnames(matrix_mirna_go_inter)))
      median_score[k,]<-row

    }
    matrix_median_genes<-median_score*matrix_mirna_go_jaccard
    colnames(matrix_median_genes)<-colnames(matrix_mirna_go_inter)
    rownames(matrix_median_genes)<-rownames(matrix_mirna_go_inter)

    matrix_m_m_score<-t(matrix_median_genes)%*%matrix_median_genes
    matrix_m_m_score[is.na(matrix_m_m_score)]<-0
    diag(matrix_m_m_score)<-0
    return(matrix_m_m_score)
  }
  matrix_m_m_score<-median_inter(GDEscore)
  adj.final<-as.matrix(matrix_m_m_score)
  graph = graph.adjacency(adj.final,mode=c("undirected"), weighted=TRUE,add.rownames=T)
  temp = page.rank(graph, vids=V(graph), directed=FALSE, damping=damping, weights=NULL)
  rank = temp$vector
  rank1 = as.matrix(rank)

  iter<-nperm
  Centrality_Scores<-matrix(nrow=dim(matrix_m_m_score)[1],ncol=iter+1)
  real.centra<-rank1
  Centrality_Scores[,1]<-real.centra
  real.subname<-rownames(real.centra)
  rownames(Centrality_Scores)<-rownames(real.centra)

  for(i in 1:iter){
    rownames(GDEscore.table)<-sample(rownames(GDEscore.table),replace = F)
    GDEscore<-cbind(rownames(GDEscore.table),abs(GDEscore.table[,"GDEscore"]))
    per.re<-median_inter(GDEscore)
    per.adj.final<-as.matrix(per.re)
    per.graph = graph.adjacency(per.adj.final,mode=c("undirected"), weighted=TRUE,add.rownames=T)
    per.temp = page.rank(per.graph, vids=V(per.graph), directed=FALSE, damping=0.90, weights=NULL)
    per.rank = per.temp$vector
    per.rank = as.matrix(per.rank)
    Centrality_Scores[,i+1]<-per.rank[real.subname,]

    print(i)
  }


  adj = as.matrix(Centrality_Scores)
  allresult<-c()
  for(i in 1:length(Centrality_Scores[,1])){
    perm_rank = adj[i,2:(iter+1)]
    orig_rank = adj[i,1]
    pval=pnorm(orig_rank,mean = mean(perm_rank),sd=sd(perm_rank),lower.tail = FALSE)
    ## p_value<-cbind(kegg1[,1],)
    ## colnames(p_value)<-c("subpathway","pvalue")
    p_padjust<-p.adjust(pval,method = "fdr")
    pa<-as.numeric(p_padjust)
    fdr<-round(pa,3)
    re<-cbind(orig_rank,pval,fdr)
    allresult<-rbind(allresult,re)
  }
  allresult<-cbind(rownames(allresult),allresult)
  allresult<-data.frame(allresult)
  colnames(allresult)<-c("miRNA","Score","P-value","FDR")
  allresult[,2]<-as.numeric(round(as.numeric(as.character(allresult[,"Score"])),5))
  allresult[,3]<-as.numeric(round(as.numeric(as.character(allresult[,"P-value"])),3))
  allresult[,4]<-as.numeric(round(as.numeric(as.character(allresult[,"FDR"])),3))

  MiRNA_Target<- GetData_Mirna("MiRNA_Target")

  result<-merge(MiRNA_Target,allresult,by.x="miRNA",by.y="miRNA")
  result<-result[order(result[,"FDR"],decreasing=F),]
  return(result)
}
