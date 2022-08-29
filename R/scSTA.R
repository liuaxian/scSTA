#' scSTA main algorithm
#'
#' scSTA  main algorithm: takes a preprocessed gene expression matrix as input. Then applies standard Seurat pipeline matrix over multiple resolution parameters.
#' @param sc A Seurat object with normalized count
#' @importFrom stats ave
#' @return  a figure as ggplot2 object of stability score system
#' @importFrom ggplot2 ggplot aes geom_line labs scale_y_continuous scale_x_continuous geom_text geom_hline scale_linetype_manual guides guide_legend geom_point scale_colour_continuous theme element_text element_blank element_line rel arrow unit margin
#' @export


#usethis::use_package(mclust, type = "Imports", min_version = NULL)

scSTA <- function(sc) {
  #suppressPackageStartupMessages(library(Seurat))
  #suppressPackageStartupMessages(library(clustree))
  suppressPackageStartupMessages(library(mclust))
  #suppressPackageStartupMessages(library(ggplot2))
  #suppressPackageStartupMessages(library(grid))
  # normalizing the data

  #sc <- NormalizeData(object = sc, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)
  # sc <- FindClusters(object = sc,resolution = c(seq(0.1,1.6,.1)))
  #clustree(sc@meta.data,prefix="RNA_snn_res.")

  # Calculate Cluster number
  k1 = max(as.numeric(sc@meta.data[["RNA_snn_res.0.1"]]))
  k2 = max(as.numeric(sc@meta.data[["RNA_snn_res.0.2"]]))
  k3 = max(as.numeric(sc@meta.data[["RNA_snn_res.0.3"]]))
  k4 = max(as.numeric(sc@meta.data[["RNA_snn_res.0.4"]]))
  k5 = max(as.numeric(sc@meta.data[["RNA_snn_res.0.5"]]))
  k6 - max(as.numeric(sc@meta.data[["RNA_snn_res.0.6"]]))
  k7 = max(as.numeric(sc@meta.data[["RNA_snn_res.0.7"]]))
  k8 = max(as.numeric(sc@meta.data[["RNA_snn_res.0.8"]]))
  k9 = max(as.numeric(sc@meta.data[["RNA_snn_res.0.9"]]))
  k10 = max(as.numeric(sc@meta.data[["RNA_snn_res.1"]]))
  k11 = max(as.numeric(sc@meta.data[["RNA_snn_res.1.1"]]))
  k12 = max(as.numeric(sc@meta.data[["RNA_snn_res.1.2"]]))
  k13 = max(as.numeric(sc@meta.data[["RNA_snn_res.1.3"]]))
  k14 = max(as.numeric(sc@meta.data[["RNA_snn_res.1.4"]]))
  k15 = max(as.numeric(sc@meta.data[["RNA_snn_res.1.5"]]))

  # Calculate Cluster stability score
  Res0.1 = adjustedRandIndex(sc$RNA_snn_res.0.1,sc$RNA_snn_res.0.2)
  Res0.2 = adjustedRandIndex(sc$RNA_snn_res.0.2,sc$RNA_snn_res.0.3)
  Res0.3 = adjustedRandIndex(sc$RNA_snn_res.0.3,sc$RNA_snn_res.0.4)
  Res0.4 = adjustedRandIndex(sc$RNA_snn_res.0.4,sc$RNA_snn_res.0.5)
  Res0.5 = adjustedRandIndex(sc$RNA_snn_res.0.5,sc$RNA_snn_res.0.6)
  Res0.6 = adjustedRandIndex(sc$RNA_snn_res.0.6,sc$RNA_snn_res.0.7)
  Res0.7 = adjustedRandIndex(sc$RNA_snn_res.0.7,sc$RNA_snn_res.0.8)
  Res0.8 = adjustedRandIndex(sc$RNA_snn_res.0.8,sc$RNA_snn_res.0.9)
  Res0.9 = adjustedRandIndex(sc$RNA_snn_res.0.9,sc$RNA_snn_res.1)
  Res1 = adjustedRandIndex(sc$RNA_snn_res.1,sc$RNA_snn_res.1.1)
  Res1.1 = adjustedRandIndex(sc$RNA_snn_res.1.1,sc$RNA_snn_res.1.2)
  Res1.2 = adjustedRandIndex(sc$RNA_snn_res.1.2,sc$RNA_snn_res.1.3)
  Res1.3 = adjustedRandIndex(sc$RNA_snn_res.1.3,sc$RNA_snn_res.1.4)
  Res1.4 = adjustedRandIndex(sc$RNA_snn_res.1.4,sc$RNA_snn_res.1.5)
  Res1.5 = adjustedRandIndex(sc$RNA_snn_res.1.5,sc$RNA_snn_res.1.6)

  k_all <- as.vector(NULL)
  
  for(i in seq(0.1,1.5,.1)){
  k_all <- c(k_all, get(paste("Res", i, sep="")))}
  
  Stability_score <- k_all

  #Calulate the mean stability score
  a = ave(x = k_all)

  z_all = as.vector(NULL)
  for(i in seq(1,15,1)){
  z_all = c(z_all, get(paste("k", i, sep="")))}
  
  Cluster_number <- z_all
  b <- max(Cluster_number)
  df1 <- data.frame(x=c(seq(0.1,1.5,0.1)),y=Stability_score)

  sb <- ggplot(data=df1,aes(x = x, y=Stability_score),main="Stability")+
    geom_line(size = 0.5,color = "black")+
    labs(x = "Resolution",y = "ARI")+
    labs(title = "Stability")+
    scale_y_continuous(breaks = seq(0,1,0.01))+
    scale_x_continuous(limits=c(0,1.5),breaks = seq(0,1.5,0.1))+
    geom_text(aes(label = round(Cluster_number,4),vjust=-1), color="black",size = 5)+
    geom_text(aes(label = round(Stability_score,4),vjust=2), color="red",size = 4)+
    geom_hline(aes(yintercept = a,linetype="Mean score"),color= "red")+
    scale_linetype_manual("Mean score",values = 2) +
    guides(fill = guide_legend(override.aes = list(linetype = 0)))+
    geom_point(aes(size=Stability_score,colour =Cluster_number))+
    scale_colour_continuous(limits=c(1,b+5),low = "gray", high = "blue")+
    theme(legend.title= element_text(colour = "black"))+
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))),
          axis.title = element_text(size = rel(1.2)),
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(size = rel(1.5),hjust = 0.5),
          axis.text.y = element_text(hjust = 1,size = rel(1.5)),
          axis.ticks = element_line(size = rel(1.5)),
          plot.margin = margin(15,9,9,9))

  print(sb)

}


