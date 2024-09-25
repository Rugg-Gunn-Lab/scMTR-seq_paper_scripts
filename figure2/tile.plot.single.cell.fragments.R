library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(purrr)
library(data.table)
library(stringr)
H3K27me3.color="#E6332A"
H3K27ac.color="#365FAA"
H3K4me1.color="#AFCB31"
H3K4me3.color="#FFE609"
H3K36me3.color="#7CC291"
IgG.color="grey60"


in.file.dir<-"/Users/yang/data/scMTR_method.dev/figure3.data"
signac<-readRDS(paste0(in.file.dir,"/seurat.rna.hist_bin5k.rds"))


################################################################################
######################### tile plot single cell fragments ######################
################################################################################


# tile plot

in.file.dir<-"/Users/yang/data/scMTR_method.dev/figure3.data"
signac<-readRDS(paste0(in.file.dir,"/seurat.rna.hist_bin5k.processed.rds"))

for (i in c("H3K27me3", "H3K27ac","H3K4me1","H3K4me3","H3K36me3", "IgG")){
  if (i =="H3K27me3") {
    plot.color="#E6332A"
  } 
  
  if (i =="H3K27ac"){
    plot.color="#365FAA"
  }
  
  if (i =="H3K4me1"){
    plot.color="#AFCB31"
  }
  
  if (i =="H3K4me3"){
    plot.color="#FFE609"
  }
  
  if (i =="H3K36me3"){
    plot.color="#7CC291"
  }
  
  if (i =="IgG"){
    plot.color="grey60"
  }
  plot.color="black"
  
  DefaultAssay(signac)<- i
  p1<-TilePlot(
    object = signac,
    region = "chr1-35000000-38500000",
    tile.size = 20000
  )
  to.plot<-p1$data%>% as.data.table()
  order<- to.plot[,sum(value),name]%>% .$name
  to.plot$name<-factor(to.plot$name,levels = rev(order))
  pdf(paste0(in.file.dir,"/",i,".tile.plot.pdf"),height=10,width=6)
  
  p2<- ggplot(to.plot,aes(bin,name,fill=value))+
    geom_tile(color = NA)+
    scale_fill_gradient(low = "white", high = plot.color,space = "Lab")+
    labs(x="",y="",title = i)+
    theme_minimal()+
    theme(
      panel.grid = element_blank()
    )
  print(p2)
  dev.off()
  
}

gene_plot <- AnnotationPlot(
  object = signac,
  region = "chr1-35000000-38500000"
)
gene_plot

# increase tile size could apparently increase the visibility, 200 tiles
#chr1:35,000,000-38,500,000
signac@meta.data$celltype<-"H9"
saveRDS(signac, paste0(in.file.dir,"/seurat.rna.hist_bin5k.processed.rds"))


