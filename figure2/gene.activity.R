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


##############################################################
######## calculate gene activity with different hist #########
##############################################################

DefaultAssay(signac)<-"H3K27ac"
gene.activities.H3K27ac<-GeneActivity(signac)
signac[["H3K27ac_genes"]]<-CreateAssayObject(counts=gene.activities.H3K27ac)
signac<-NormalizeData(object=signac, assay="H3K27ac_genes", normalization.method="LogNormalize", scale.factor=median(signac$nCount_H3K27ac_genes))
DefaultAssay(signac)<-"H3K27ac_genes"

DefaultAssay(signac)<-"H3K27me3"
gene.activities.H3K27me3<-GeneActivity(signac)
signac[["H3K27me3_genes"]]<-CreateAssayObject(counts=gene.activities.H3K27me3)
signac<-NormalizeData(object=signac, assay="H3K27me3_genes", normalization.method="LogNormalize", scale.factor=median(signac$nCount_H3K27me3_genes))
DefaultAssay(signac)<-"H3K27me3_genes"

DefaultAssay(signac)<-"H3K4me3"
gene.activities.H3K4me3<-GeneActivity(signac)
signac[["H3K4me3_genes"]]<-CreateAssayObject(counts=gene.activities.H3K4me3)
signac<-NormalizeData(object=signac, assay="H3K4me3_genes", normalization.method="LogNormalize", scale.factor=median(signac$nCount_H3K4me3_genes))
DefaultAssay(signac)<-"H3K4me3_genes"

DefaultAssay(signac)<-"H3K4me1"
gene.activities.H3K4me1<-GeneActivity(signac)
signac[["H3K4me1_genes"]]<-CreateAssayObject(counts=gene.activities.H3K4me1)
signac<-NormalizeData(object=signac, assay="H3K4me1_genes", normalization.method="LogNormalize", scale.factor=median(signac$nCount_H3K4me1_genes))
DefaultAssay(signac)<-"H3K4me1_genes"

DefaultAssay(signac)<-"H3K36me3"
gene.activities.H3K36me3<-GeneActivity(signac)
signac[["H3K36me3_genes"]]<-CreateAssayObject(counts=gene.activities.H3K36me3)
signac<-NormalizeData(object=signac, assay="H3K36me3_genes", normalization.method="LogNormalize", scale.factor=median(signac$nCount_H3K36me3_genes))
DefaultAssay(signac)<-"H3K36me3_genes"

in.file.dir<-"/Users/yang/data/scMTR_method.dev/figure3.data"
saveRDS(signac,paste0(in.file.dir,"/seurat.rna.hist_bin5k.processed.rds"))
#signac<-readRDS(paste0(in.file.dir,"/seurat.rna.hist_bin5k.processed.rds"))

