
library(alluvial)
library(ggalluvial)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(Seurat)
library(Signac)

################################################################################
########################### alluvial plot clusters #############################
################################################################################

in.file.dir <- "/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments"

signac<-readRDS(paste0(in.file.dir,"/seurat.rna.hist_bin5k.def.endo.rds"))

to.plot<-signac@meta.data%>% as.data.table()
to.plot<-to.plot[,.(sample_bio.id,RNA_snn_res.0.2,H3K27me3_snn_res.0.2, H3K4me1_snn_res.0.2,H3K27ac_snn_res.0.2, H3K4me3_snn_res.0.2, H3K36me3_snn_res.0.2, IgG_snn_res.0.2,wknn_res.0.2)]

 # sub.to.plot.tmp<-to.plot[D0 %in% c(i,j) |D1 %in% c(i,j)|D2 %in% c(i,j)|D3 %in% c(i,j)]
  to.plot.tmp<-as.data.frame(to.plot[,.(freq=.N),.( sample_bio.id,RNA_snn_res.0.2,H3K27me3_snn_res.0.2, H3K4me1_snn_res.0.2,H3K27ac_snn_res.0.2, H3K4me3_snn_res.0.2, H3K36me3_snn_res.0.2, IgG_snn_res.0.2,wknn_res.0.2)] )
  pdf(paste0(in.file.dir,"/alluvial.plot.cluster.pdf"),width = 8,height = 4)
  #  pdf(paste0(inputdir,"/alluvial.plot.states.promoters.pdf"),width = 8,height = 4)
  p1<-ggplot(to.plot.tmp,
             aes(y=freq,
                 axis1= sample_bio.id, 
                 axis2= RNA_snn_res.0.2,
                 axis3=H3K27me3_snn_res.0.2,
                 axis4=H3K27ac_snn_res.0.2,
                 axis5=H3K4me1_snn_res.0.2,
                 axis6= H3K4me3_snn_res.0.2,
                 axis7=H3K36me3_snn_res.0.2,
                 axis8=IgG_snn_res.0.2,
                 axis9=wknn_res.0.2
                  ) )+
    geom_alluvium(aes(fill =as.factor (sample_bio.id)), width = 1/3) +
    geom_stratum(aes(fill = as.factor (sample_bio.id)), width = 1/3) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
    #geom_label(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
    theme_void()
  print(p1)
  
  dev.off()
  fwrite(to.plot.tmp,paste0(in.file.dir,"/alluvial.cluster.tsv"),sep="\t")

  # manual adjust to match 
  to.plot<-signac@meta.data%>% as.data.table()
  to.plot<-to.plot[,.(sample_bio.id,RNA_snn_res.0.2,H3K27me3_snn_res.0.2, H3K4me1_snn_res.0.2,H3K27ac_snn_res.0.2, H3K4me3_snn_res.0.2, H3K36me3_snn_res.0.2, IgG_snn_res.0.2,wknn_res.0.2)]
  
  to.plot[RNA_snn_res.0.2=="2",RNA_snn_res.0.2:="a"]
  to.plot[RNA_snn_res.0.2=="0",RNA_snn_res.0.2:="b"]
  to.plot[RNA_snn_res.0.2=="1",RNA_snn_res.0.2:="c"]
  to.plot[RNA_snn_res.0.2=="4",RNA_snn_res.0.2:="d"]
  to.plot[RNA_snn_res.0.2=="3",RNA_snn_res.0.2:="e"]
  
  to.plot[H3K27me3_snn_res.0.2=="1",H3K27me3_snn_res.0.2:="a"]
  to.plot[H3K27me3_snn_res.0.2=="0",H3K27me3_snn_res.0.2:="b"]
  to.plot[H3K27me3_snn_res.0.2=="2",H3K27me3_snn_res.0.2:="c"]
  
  to.plot[H3K27ac_snn_res.0.2=="2",H3K27ac_snn_res.0.2:="a"]
  to.plot[H3K27ac_snn_res.0.2=="1",H3K27ac_snn_res.0.2:="b"]
  to.plot[H3K27ac_snn_res.0.2=="3",H3K27ac_snn_res.0.2:="c"]
  to.plot[H3K27ac_snn_res.0.2=="0",H3K27ac_snn_res.0.2:="d"]
  
  to.plot[H3K4me1_snn_res.0.2=="2",H3K4me1_snn_res.0.2:="a"]
  to.plot[H3K4me1_snn_res.0.2=="0",H3K4me1_snn_res.0.2:="b"]
  to.plot[H3K4me1_snn_res.0.2=="3",H3K4me1_snn_res.0.2:="c"]
  to.plot[H3K4me1_snn_res.0.2=="4",H3K4me1_snn_res.0.2:="d"]
  to.plot[H3K4me1_snn_res.0.2=="1",H3K4me1_snn_res.0.2:="e"]
  
  to.plot[H3K4me3_snn_res.0.2=="2",H3K4me3_snn_res.0.2:="a"]
  to.plot[H3K4me3_snn_res.0.2=="0",H3K4me3_snn_res.0.2:="b"]
  to.plot[H3K4me3_snn_res.0.2=="1",H3K4me3_snn_res.0.2:="c"]
  to.plot[H3K4me3_snn_res.0.2=="3",H3K4me3_snn_res.0.2:="d"]
  
  
  to.plot[H3K36me3_snn_res.0.2=="2",H3K36me3_snn_res.0.2:="a"]
  to.plot[H3K36me3_snn_res.0.2=="0",H3K36me3_snn_res.0.2:="b"]
  to.plot[H3K36me3_snn_res.0.2=="1",H3K36me3_snn_res.0.2:="c"]
  to.plot[H3K36me3_snn_res.0.2=="3",H3K36me3_snn_res.0.2:="d"]
  
  to.plot[IgG_snn_res.0.2=="0",IgG_snn_res.0.2:="a"]
  
  to.plot[wknn_res.0.2=="2",wknn_res.0.2:="a"]
  to.plot[wknn_res.0.2=="0",wknn_res.0.2:="b"]
  to.plot[wknn_res.0.2=="1",wknn_res.0.2:="c"]
  to.plot[wknn_res.0.2=="3",wknn_res.0.2:="d"]
  
  to.plot.tmp<-as.data.frame(to.plot[,.(freq=.N),.( sample_bio.id,RNA_snn_res.0.2,H3K27me3_snn_res.0.2, H3K4me1_snn_res.0.2,H3K27ac_snn_res.0.2, H3K4me3_snn_res.0.2, H3K36me3_snn_res.0.2, IgG_snn_res.0.2,wknn_res.0.2)] )
  pdf(paste0(in.file.dir,"/alluvial.plot.cluster.manual.ordering.pdf"),width = 8,height = 4)
  p1<-ggplot(to.plot.tmp,
             aes(y=freq,
                 axis1= sample_bio.id, 
                 axis2= RNA_snn_res.0.2,
                 axis3=H3K27me3_snn_res.0.2,
                 axis4=H3K27ac_snn_res.0.2,
                 axis5=H3K4me1_snn_res.0.2,
                 axis6= H3K4me3_snn_res.0.2,
                 axis7=H3K36me3_snn_res.0.2,
                 axis8=IgG_snn_res.0.2,
                 axis9=wknn_res.0.2
             ) )+
    geom_alluvium(aes(fill =as.factor (sample_bio.id)), width = 1/3) +
    geom_stratum(aes(fill = as.factor (sample_bio.id)), width = 1/3) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
    theme_void()
  print(p1)
  
  dev.off()
  fwrite(to.plot.tmp,paste0(in.file.dir,"/alluvial.cluster.manual.ordering.tsv"),sep="\t")
  
  

