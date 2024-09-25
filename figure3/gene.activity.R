
## generate signac files for definitive endoderm 
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
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


signac<-readRDS("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/fragments/seurat.rna.hist_bin5k.def.endo.rds")

##############################################################
######## calculate gene activity with different hist #########
##############################################################


for (i in c ( "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "IgG" )) {
  #  i<-"H3K27me3"
  # i<-"H3K4me1"
 # DefaultAssay(signac)<-i
  
  DefaultAssay(signac)<-i
  gene.activities<-GeneActivity(signac)
  signac[[paste0(i,"_genes")]]<-CreateAssayObject(counts=gene.activities)
  signac<-NormalizeData(object=signac, assay=paste0(i,"_genes"), normalization.method="LogNormalize", scale.factor= signac@meta.data %>% .[,paste0("nCount_",i,"_genes")] %>% median() )
  DefaultAssay(signac)<-paste0(i,"_genes")
}

saveRDS(signac,paste0(in.file.dir,"/seurat.rna.hist_bin5k.def.endo.rds"))




  
  
  
  
  
  
  
  
  
  
  


