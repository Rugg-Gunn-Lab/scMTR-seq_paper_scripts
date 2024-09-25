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
############# perform dimensionality reduction ###############
##############################################################

############################ RNA #############################

DefaultAssay(signac)<-"RNA"
signac <- NormalizeData(signac, normalization.method = "LogNormalize", scale.factor = 10000)
signac <- FindVariableFeatures(signac, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(signac), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(signac)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(signac)
signac <- ScaleData(signac, features = all.genes)
signac <- RunPCA(signac, features = VariableFeatures(object = signac))
print(signac[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(signac, dims = 1:2, reduction = "pca")
DimPlot(signac, reduction = "pca",dims = c(1,2)) #,group.by= "seurat_clusters"
DimPlot(signac, reduction = "pca",dims = c(1,2),group.by= "orig.ident")
DimHeatmap(signac, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(signac, dims = 1:20, balanced = TRUE) #cells = 500, 
signac <- JackStraw(signac, dims = 20, num.replicate = 100)
signac <- ScoreJackStraw(signac, dims = 1:20)
dev.off()
pdf(paste0(in.file.dir,"/hE_jackstrawplot.pdf"),width = 5,height = 4)
JackStrawPlot(signac, dims = 1:20)
dev.off()
pdf(paste0(in.file.dir,"/hE_ElbowPlot.pdf"),width = 5,height = 4)
ElbowPlot(signac,ndims = 20)
dev.off()
signac <- FindNeighbors(signac, dims = 1:10)
i<-0.5
signac <- FindClusters(signac, resolution = i)
head(Idents(signac), 5)
signac <- RunUMAP(signac, dims = 1:10)

pdf(paste0(in.file.dir,"/RNA_res.",i,".pdf"),width = 8,height = 6)
p1<-DimPlot(signac, reduction = "umap", label = TRUE)#+NoLegend()
print(p1)
dev.off()

############################ H3K27ac #############################

DefaultAssay(signac)<-"H3K27ac"
signac <- RunTFIDF(signac, scale.factor = median(signac$nCount_H3K27ac))
VariableFeatures(signac) <- rownames(signac)
signac <- RunSVD(signac, reduction.name="lsi.H3K27ac")
signac <- RunUMAP(signac, reduction = "lsi.H3K27ac", reduction.name="umap.H3K27ac", dims = 2:8)
signac <- FindNeighbors(signac, reduction = "lsi.H3K27ac", dims = 2:8)
signac <- FindClusters(signac, algorithm = 3)

DimPlot(signac)


############################ H3K27me3 #############################

DefaultAssay(signac)<-"H3K27me3"
signac <- RunTFIDF(signac, scale.factor = median(signac$nCount_H3K27me3))
VariableFeatures(signac) <- rownames(signac)
signac <- RunSVD(signac, reduction.name="lsi.H3K27me3")
signac <- RunUMAP(signac, reduction = "lsi.H3K27me3", reduction.name="umap.H3K27me3", dims = 2:8)
signac <- FindNeighbors(signac, reduction = "lsi.H3K27me3", dims = 2:8)
signac <- FindClusters(signac, algorithm = 3)

DimPlot(signac)

############################ H3K4me3 #############################

DefaultAssay(signac)<-"H3K4me3"
signac <- RunTFIDF(signac, scale.factor = median(signac$nCount_H3K4me3))
VariableFeatures(signac) <- rownames(signac)
signac <- RunSVD(signac, reduction.name="lsi.H3K4me3")
signac <- RunUMAP(signac, reduction = "lsi.H3K4me3", reduction.name="umap.H3K4me3", dims = 2:8)
signac <- FindNeighbors(signac, reduction = "lsi.H3K4me3", dims = 2:8)
signac <- FindClusters(signac, algorithm = 3)

DimPlot(signac)

############################ H3K4me1 #############################

DefaultAssay(signac)<-"H3K4me1"
signac <- RunTFIDF(signac, scale.factor = median(signac$nCount_H3K4me1))
VariableFeatures(signac) <- rownames(signac)
signac <- RunSVD(signac, reduction.name="lsi.H3K4me1")
signac <- RunUMAP(signac, reduction = "lsi.H3K4me1", reduction.name="umap.H3K4me1", dims = 2:8)
signac <- FindNeighbors(signac, reduction = "lsi.H3K4me1", dims = 2:8)
signac <- FindClusters(signac, algorithm = 3)

DimPlot(signac)

############################ H3K36me3 #############################

DefaultAssay(signac)<-"H3K36me3"
signac <- RunTFIDF(signac, scale.factor = median(signac$nCount_H3K36me3))
VariableFeatures(signac) <- rownames(signac)
signac <- RunSVD(signac, reduction.name="lsi.H3K36me3")
signac <- RunUMAP(signac, reduction = "lsi.H3K36me3", reduction.name="umap.H3K36me3", dims = 2:8)
signac <- FindNeighbors(signac, reduction = "lsi.H3K36me3", dims = 2:8)
signac <- FindClusters(signac, algorithm = 3)

DimPlot(signac)

############################ IgG #############################

DefaultAssay(signac)<-"IgG"
signac <- RunTFIDF(signac, scale.factor = median(signac$nCount_IgG))
VariableFeatures(signac) <- rownames(signac)
signac <- RunSVD(signac, reduction.name="lsi.IgG")
signac <- RunUMAP(signac, reduction = "lsi.IgG", reduction.name="umap.IgG", dims = 2:8)
signac <- FindNeighbors(signac, reduction = "lsi.IgG", dims = 2:8)
signac <- FindClusters(signac, algorithm = 3)

DimPlot(signac)


############## Integrate multimodal neighbors ################
signac<-FindMultiModalNeighbors(object=signac, 
                                reduction.list=list("pca","lsi.H3K27ac", "lsi.H3K27me3", "lsi.H3K4me1", "lsi.H3K4me3", "lsi.H3K36me3" ), 
                                dims.list=list(1:10, 2:8, 2:8, 2:8,2:8, 2:8), 
                                modality.weight.name=list("RNA.weight","H3K27ac.weight", "H3K27me3.weight","H3K4me1.weight", "H3K4me3.weight","H3K36me3.weight"), 
                                verbose=TRUE)

### Run UMAP on WNN graph and find clusters ###
signac<-RunUMAP(signac, nn.name="weighted.nn", reduction.name="wnn.umap", reduction.key = "wnnUMAP_")
signac <- FindClusters(signac, graph.name="wsnn", algorithm = 3)

DimPlot(signac,reduction = "wnn.umap" ,label = T)

saveRDS(signac,paste0(in.file.dir,"/seurat.rna.hist_bin5k.processed.rds"))



##############################################################
##################### plot connected UMAP ####################
##############################################################

signac@meta.data$celltype<-"H9"

to.plot.coord<-list()
for (i in c("umap.H3K27ac", "umap.H3K27me3","umap.H3K4me3", "umap.H3K4me1","umap.H3K36me3","umap","wnn.umap" )){
  tmp.coord<-as.data.frame(signac[[i]]@cell.embeddings)
  colnames(tmp.coord)<-c("UMAP_1", "UMAP_2")
  tmp.coord$Target<-i 
  tmp.coord$Cell<-rownames(tmp.coord)
  tmp.coord$Celltype<-signac@meta.data$celltype
  tmp.coord.mod<-tmp.coord
  
  tmp.coord.mod$xend<-tmp.coord.mod$UMAP_1
  tmp.coord.mod$yend<-tmp.coord.mod$UMAP_2
  to.plot.coord[[i]]<-tmp.coord.mod
  
}
to.plot.coord<-rbindlist(to.plot.coord)
to.plot.coord[Target=="umap",Target:="RNA"]
to.plot.coord[Target=="wnn.umap",Target:="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA"]

xyend<-to.plot.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA",.(Cell,xend,yend)]
to.plot.coord<-merge(to.plot.coord[,-c("xend","yend")],xyend,by="Cell")
library(ggplot2)
to.plot.coord[Target=="umap.H3K27ac",UMAP_1:=UMAP_1-16]
to.plot.coord[Target=="umap.H3K27ac",UMAP_2:=UMAP_2-8]
to.plot.coord[Target=="umap.H3K4me1",UMAP_1:=UMAP_1-16]
to.plot.coord[Target=="umap.H3K4me1",UMAP_2:=UMAP_2+8]
to.plot.coord[Target=="umap.H3K27me3",UMAP_1:=UMAP_1]
to.plot.coord[Target=="umap.H3K27me3",UMAP_2:=UMAP_2-16]
to.plot.coord[Target=="umap.H3K4me3",UMAP_1:=UMAP_1+16]
to.plot.coord[Target=="umap.H3K4me3",UMAP_2:=UMAP_2-8]
to.plot.coord[Target=="umap.H3K36me3",UMAP_1:=UMAP_1]
to.plot.coord[Target=="umap.H3K36me3",UMAP_2:=UMAP_2+16]
to.plot.coord[Target=="RNA",UMAP_1:=UMAP_1+16]
to.plot.coord[Target=="RNA",UMAP_2:=UMAP_2+8]


for (i in c("umap.H3K27ac", "umap.H3K27me3", "umap.H3K4me3", "umap.H3K4me1", "umap.H3K36me3", "umap", "wnn.umap")){
  pdf(paste0(in.file.dir,"/dimplot.reduction.",i,".pdf"),width=5,height=4)
  p1<-DimPlot(signac,reduction = i )
  p1<-ggplot(p1$data,aes(p1$data[,1],p1$data[,2],color=p1$data[,3]))+
    geom_point()+
    labs(x="UMAP 1",y="UMAP 2",title = i)+
    theme_classic()+
    theme(
      axis.text = element_text(size=12)
    )
  print(p1)
  dev.off()
}

pdf(paste0(in.file.dir,"/joint.projection.pdf"),width=16,height=10)
p1<-ggplot(to.plot.coord) + 
  geom_segment(aes(x=UMAP_1, y=UMAP_2, color=Target, xend=xend, yend=yend), alpha=0.01) + 
  geom_point(aes(x=UMAP_1, y=UMAP_2, color=Target), cex=0.5,size=0.05,alpha=0.8) + 
  theme_light() +
  xlab("UMAP 1") + 
  ylab("UMAP2") + 
  scale_color_manual(values=c(`H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA`="#F8766D", RNA="#DB72FB", umap.H3K27ac="#93AA00",umap.H3K4me1= "#00BA38" ,umap.H3K27me3="#D39200",umap.H3K4me3="#00C19F" ,umap.H3K36me3="#00B9E3")) + 
  # scale_color_viridis(discrete=TRUE) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=14))
print(p1)
dev.off()

fwrite(to.plot.coord,paste0(in.file.dir,"/joint.projection.tsv"),sep="\t")



