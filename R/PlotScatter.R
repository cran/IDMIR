#'@title PlotScatter
#'@description Function "PlotScatter" is used to plot a scatter diagram.
#' @param MiRNA_CRData A list includes a data frame with seven parts those are "sample", "status", "time", "target genes expression", "risk score", "group", and a data frame with five columns those are "Gene", "HR", "HR.95L", "HR.95H", "beta", and "P-value".
#' @param status.0 string. Code for event 0. Default is 'Alive'
#' @param status.1 string. Code for event 1. Default is 'Dead'
#' @param TitleYlab_A string, y-lab title for figure A. Default is 'Riskscore'
#' @param TitleYlab_B string, y-lab title for figure B. Default is 'Survival Time'
#' @param TitleXlab string, x-lab title for figure B. Default is 'Rank'
#' @param TitleLegend_A string, legend title for figure A. Default is 'Risk Group'
#' @param TitleLegend_B string, legend title for figure B. Default is 'Status'
#' @param color.A color for figure A. Default is low = 'blue', high = 'red'
#' @param color.B color for figure B. Default is status.0 = 'blue', status.1 = 'red'
#' @return A riskscore picture
#' @importFrom egg ggarrange
#' @importFrom ggplot2 aes aes_string geom_point geom_vline theme element_blank element_text scale_colour_hue coord_trans
#' @importFrom ggplot2 ylab geom_tile unit scale_fill_gradient2 scale_x_continuous geom_raster theme_classic annotate
#' @importFrom ggplot2 scale_color_manual element_line scale_fill_manual ggplot scale_fill_manual xlab
#' @usage PlotScatter(MiRNA_CRData,status.0='Alive',status.1='Dead',
#' TitleYlab_A='Risk Score',TitleYlab_B='Survival Time',TitleXlab='Rank',
#' TitleLegend_A='Risk Group',TitleLegend_B='Status',
#' color.A=c(low='blue',high='red'),color.B=c(status.0='blue',status.1='red'))
#' @export
#' @examples
#' # Obtain the example data
#' GEP<-GetData_Mirna("GEP")
#' survival<-GetData_Mirna("survival")
#' MiRNAs<-c("hsa-miR-21-5p","hsa-miR-26a-5p","hsa-miR-369-5p","hsa-miR-1238-3p","hsa-miR-10b-5p")
#' # Run the function
#' SingleMiRNA_CRData<-SingleMiRNA_CRModel(GEP,
#' "hsa-miR-21-5p",survival,cutoff.point=NULL)
#' PlotScatter(SingleMiRNA_CRData)
#' MutiMiRNA_CRData<-MutiMiRNA_CRModel(GEP,
#' MiRNAs,survival,cutoff.point=NULL)
#' PlotScatter(MutiMiRNA_CRData)

PlotScatter<-function(MiRNA_CRData,
                      status.0='Alive',
                      status.1='Dead',
                      TitleYlab_A='Risk Score',
                      TitleYlab_B='Survival Time',
                      TitleXlab='Rank',
                      TitleLegend_A='Risk Group',
                      TitleLegend_B='Status',
                      color.A=c(low='blue',high='red'),
                      color.B=c(status.0='blue',status.1='red')
                      ){
  data<-data.frame(MiRNA_CRData[[1]])
  data=data[order(data[,"score"],decreasing = F),]
  cutoff.point.x<-length(which(data$group=="low"))
  cutoff.point.y<-median(data$score)

  data$score=round(data$score,1)

  `Group` = data$group
  #figure A risk plot
  #rearange colorA
  color.A=c(color.A['low'],color.A['high'])
  names(color.A)=c("low","high")
  fA = ggplot(data = data,
              aes_string(
                x = 1:nrow(data),
                y = data$score,
                color=factor(`Group`)
              )
  ) +
    geom_point(size = 2) +
    scale_color_manual(name=TitleLegend_A,values = color.A) +
    geom_vline(
      xintercept = cutoff.point.x,
      linetype = 'dotted',
      size = 1
    ) +
    #bg
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank())+
    #x-axis
    theme(
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    #y-axis
    theme(
      axis.title.y = element_text(
        size = 14,vjust = 1,angle = 90,family="sans"),
      axis.text.y = element_text(size=11,family = "sans"),
      axis.line.y = element_line(size=0.5,colour = "black"),
      axis.ticks.y = element_line(size = 0.5,colour = "black"))+
    #legend
    theme(legend.title = element_text(size = 13,family = "sans"),
          legend.text = element_text(size=12,family = "sans"))+
    coord_trans()+
    ylab(TitleYlab_A)+
    scale_x_continuous(expand = c(0,3))

    cutoff.label=paste0('cutoff: ',round(cutoff.point.y,2))

    fA=fA+ annotate("text",
                    x=cutoff.point.x,
                    y=cutoff.point.y,
                    label=cutoff.label,
                    family="sans",
                    size=5,
                    fontface="plain",
                    colour="black")


  #fB
  color.B=c(color.B['status.0'],color.B['status.1'])
  names(color.B)=c(status.0,status.1)
  fB=ggplot(data = data,
            aes_string(
              x = 1:nrow(data),
              y = data[, 2],
              color=factor(ifelse(data[,"status"]==1,status.0,status.1)))
  ) +
    geom_point(size=2)+
    scale_color_manual(name=TitleLegend_B,values = color.B) +
    geom_vline(
      xintercept = cutoff.point.x,
      linetype = 'dotted',
      size = 1
    )  +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank())+
    #x a-xis
    theme(
      axis.line.x = element_line(size=0.5,colour = "black"),
      axis.ticks.x = element_line(size=0.5,colour = "black"),
      axis.text.x = element_text(size=11,family = "sans"),
      axis.title.x = element_text(size = 14,family="sans")
    ) +
    #y-axis
    theme(
      axis.title.y = element_text(
        size = 14,vjust = 2,angle = 90,family="sans"),
      axis.text.y = element_text(size=11,family = "sans"),
      axis.ticks.y = element_line(size = 0.5),
      axis.line.y = element_line(size=0.5,colour = "black")
    )+
    theme(legend.title = element_text(size = 13,family = "sans"),
          legend.text = element_text(size=12,family = "sans"))+
    ylab(TitleYlab_B)+xlab(TitleXlab)+
    coord_trans()+
    scale_x_continuous(expand = c(0,3))


  ggarrange(
    fA,
    fB,
    ncol = 1,
    labels = c('A', 'B'),
    label.args = list(gp = grid::gpar(font = 2, cex =1.5,
                                      family="sans"))
  )



}
