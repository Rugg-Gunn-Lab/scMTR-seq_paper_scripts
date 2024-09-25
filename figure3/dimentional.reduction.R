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


##########################
##########################
##########################
rm(scRNA)

signac<-readRDS("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/fragments/seurat.rna.hist_bin5k.def.endo.rds")


##############################################################
############# perform dimensionality reduction ###############
##############################################################

############################ RNA #############################

DefaultAssay(signac)<-"RNA"

counts <- GetAssayData(signac,slot = "counts",assay = "RNA")   
genes.percent.expression <- rowMeans(counts>0 )*100 
counts[counts>0]<-1
genes.percent.expression <- rowSums(counts) /31827 *100
table(names(genes.percent.expression[genes.percent.expression>0.1]) %in%proteincodesgenes )
hist(genes.percent.expression)
presence 
table(rowSums(counts>0) >30000)
## only subtract protein coding genes
proteincodesgenes<-annotation %>% as.data.table() %>% .[gene_biotype=="protein_coding",unique(gene_name)]
scRNA2<-subset(signac,features= proteincodesgenes)


signac <- NormalizeData(signac)
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
DimPlot(signac, reduction = "pca",dims = c(3,2)) #,group.by= "seurat_clusters"

DimPlot(signac, reduction = "pca",dims = c(1,2),group.by= "sample_bio.id",raster = T,raster.dpi = c(261,261))
dev.off()

dir.create(paste0(in.file.dir,"/plot"))

pdf(paste0(in.file.dir,"/plot/seurat.rna.PCA.pdf"),width = 8,height=6)
DimPlot(signac, reduction = "pca",dims = c(1,2),group.by= "sample_bio.id")
dev.off()
pdf(paste0(in.file.dir,"/plot/seurat.rna.PCA.split.pdf"),width = 12,height=6)
DimPlot(signac, reduction = "pca",dims = c(1,2),group.by= "sample_bio.id",split.by = "sample_bio.id")
dev.off()

DimHeatmap(signac, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(signac, dims = 1:20, cells = 500, balanced = TRUE) #cells = 500, 
# could be time consuming 
#signac <- JackStraw(signac, dims = 10, num.replicate = 100)
#signac <- ScoreJackStraw(signac, dims = 1:20)
dev.off()

#pdf(paste0(in.file.dir,"/plot/RNA_jackstrawplot.pdf"),width = 5,height = 4)
#JackStrawPlot(signac, dims = 1:20)
#dev.off()

pdf(paste0(in.file.dir,"/RNA_ElbowPlot.pdf"),width = 5,height = 4)
ElbowPlot(signac,ndims = 20)
dev.off()

signac <- FindNeighbors(signac, dims = 1:10)
i<-0.2
signac <- FindClusters(signac, resolution = i)
head(Idents(signac), 5)
signac <- RunUMAP(signac,
                  # n.neighbors=50,
                  dims = 1:10)

pdf(paste0(in.file.dir,"/plot/RNA_res.",i,".pdf"),width = 8,height = 6)
p1<-DimPlot(signac, reduction = "umap", label = TRUE)#+NoLegend()
p2<-DimPlot(signac, reduction = "umap",group.by ="sample_bio.id",label = TRUE)#+NoLegend()
#p3<-DimPlot(signac, reduction = "umap",group.by ="sample.anno",label = TRUE)#+NoLegend()
print(p1)
print(p2)
# print(p3)
dev.off()

saveRDS(signac,paste0(in.file.dir,"/seurat.rna.hist_bin5k.def.endo.rds"))

signac<-readRDS(paste0(in.file.dir,"/seurat.rna.hist_bin5k.def.endo.rds"))



############################ hist marks #############################
for (i in c ("H3K27me3", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "IgG" )) {
  #  i<-"H3K27me3"
  # i<-"H3K4me1"
  DefaultAssay(signac)<-i
  signac <- RunTFIDF(signac, scale.factor = signac@meta.data %>% .[,paste0("nCount_",i)] %>% median() )
  VariableFeatures(signac) <- rownames(signac)
  signac <- RunSVD(signac, reduction.name= paste0("lsi.",i),scale.embeddings= T) # default is use scale.embeddings, use or not doesn't seem to affect the separation, in NTT not use, MULTI-Tag used 
  signac <- RunUMAP(signac, reduction = paste0("lsi.",i), reduction.name=paste0("umap.",i), dims = 2:10)
  signac <- FindNeighbors(signac, reduction = paste0("lsi.",i), dims = 2:10)
  signac <- FindClusters(signac, algorithm = 3,resolution=0.2)
  
  pdf(paste0(in.file.dir,"/plot/UMAP_",i,".pdf"),width = 8,height = 6)
  p1<-DimPlot(signac,group.by = "sample_bio.id")
  p2<-DimPlot(signac,group.by = paste0(i,"_snn_res.0.2"))
  print(p1)
  print(p2)
  dev.off()
}


######################## Integrate multimodal neighbors ####################

signac<-FindMultiModalNeighbors(object=signac, 
                                reduction.list=list("pca","lsi.H3K27ac", "lsi.H3K27me3", "lsi.H3K4me1", "lsi.H3K4me3", "lsi.H3K36me3" ), 
                                dims.list=list(1:10, 2:10, 2:10, 2:10,2:10, 2:10), 
                                modality.weight.name=list("RNA.weight","H3K27ac.weight", "H3K27me3.weight","H3K4me1.weight", "H3K4me3.weight","H3K36me3.weight"), 
                                verbose=TRUE)
### Integrate multimodal neighbors ###

signac<-FindMultiModalNeighbors(object=signac, 
                                reduction.list=list("lsi.H3K27ac", "lsi.H3K27me3", "lsi.H3K4me1", "lsi.H3K4me3", "lsi.H3K36me3" ), 
                                dims.list=list(2:10, 2:10, 2:10,2:10, 2:10), 
                                modality.weight.name=list("H3K27ac.weight", "H3K27me3.weight","H3K4me1.weight", "H3K4me3.weight","H3K36me3.weight"), 
                                verbose=TRUE)



### Run UMAP on WNN graph and find clusters ###

signac<-RunUMAP(signac, nn.name="weighted.nn", reduction.name="wnn.umap", reduction.key = "wnnUMAP_")
signac <- FindClusters(signac, graph.name="wsnn", algorithm = 3, resolution=0.1 )
signac <- FindClusters(signac, graph.name="wknn", algorithm = 3, resolution=0.2 )
#367 singletons identified. 42 final clusters.



pdf(paste0(in.file.dir,"/plot/UMAP_joint.wnn.pdf"),width = 8,height = 6)

DimPlot(signac,reduction = "wnn.umap" ,group.by = "sample_bio.id",label = T)
DimPlot(signac,reduction = "wnn.umap" ,group.by = "wknn_res.0.2",label = T)
dev.off()

pdf(paste0(in.file.dir,"/plot/UMAP_joint.wnn.no.rna.pdf"),width = 8,height = 6)

DimPlot(signac,reduction = "wnn.umap" ,group.by = "sample_bio.id",label = T)
DimPlot(signac,reduction = "wnn.umap" ,group.by = "wknn_res.0.2",label = T)
dev.off()
to.plot<-signac@meta.data %>% as.data.table() %>% .[,.N,.(wknn_res.0.2,sample_bio.id)] #%>% dcast(.,sample_bio.id~wknn_res.0.2,value.var = "N")

pdf(paste0(in.file.dir,"/plot/UMAP_joint.wnn.cluster.vs.sample.pdf"),width = 5,height = 4)
ggplot(to.plot,aes(wknn_res.0.2 ,N,fill=sample_bio.id))+
  geom_col()+
  theme_bw()
dev.off()
fwrite(to.plot, paste0(in.file.dir,"/plot/UMAP_joint.wnn.cluster.vs.sample.tsv"), sep="\t")

pdf(paste0(in.file.dir,"/plot/UMAP_joint.wnn.no.rna.cluster.vs.sample.pdf"),width = 5,height = 4)
ggplot(to.plot,aes(wknn_res.0.2 ,N,fill=sample_bio.id))+
  geom_col()+
  theme_bw()
dev.off()
fwrite(to.plot, paste0(in.file.dir,"/plot/UMAP_joint.wnn.no.rna.cluster.vs.sample.tsv"), sep="\t")


###  ###

############################ PLot connected UMAP ##############################


to.plot.coord<-list()
for (i in c("umap.H3K27ac", "umap.H3K27me3","umap.H3K4me3", "umap.H3K4me1","umap.H3K36me3","umap","wnn.umap" )){
  # i<-"umap.H3K27ac"
  tmp.coord<-as.data.frame(signac[[i]]@cell.embeddings)
  colnames(tmp.coord)<-c("UMAP_1", "UMAP_2")
  tmp.coord$Target<-i 
  tmp.coord$Cell<-rownames(tmp.coord)
  tmp.coord$Celltype<-signac@meta.data$sample_bio.id
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

to.plot.coord[Target=="umap.H3K27ac",UMAP_1:=UMAP_1-24]
to.plot.coord[Target=="umap.H3K27ac",UMAP_2:=UMAP_2-16]

to.plot.coord[Target=="umap.H3K4me1",UMAP_1:=UMAP_1-24]
to.plot.coord[Target=="umap.H3K4me1",UMAP_2:=UMAP_2+16]

to.plot.coord[Target=="umap.H3K27me3",UMAP_1:=UMAP_1]
to.plot.coord[Target=="umap.H3K27me3",UMAP_2:=UMAP_2-24]

to.plot.coord[Target=="umap.H3K4me3",UMAP_1:=UMAP_1+24]
to.plot.coord[Target=="umap.H3K4me3",UMAP_2:=UMAP_2-16]

to.plot.coord[Target=="umap.H3K36me3",UMAP_1:=UMAP_1]
to.plot.coord[Target=="umap.H3K36me3",UMAP_2:=UMAP_2+24]

to.plot.coord[Target=="RNA",UMAP_1:=UMAP_1+24]
to.plot.coord[Target=="RNA",UMAP_2:=UMAP_2+16]

for (i in c("umap.H3K27ac", "umap.H3K27me3", "umap.H3K4me3", "umap.H3K4me1", "umap.H3K36me3", "umap", "wnn.umap")){
  pdf(paste0(in.file.dir,"/plot/dimplot.reduction.",i,".pdf"),width=5,height=4)
  p1<-DimPlot(signac,reduction = i ,group.by = "sample_bio.id")
  p1<-ggplot(p1$data,aes(p1$data[,1],p1$data[,2],color=p1$data[,3]))+
    geom_point(size=0.2)+
    labs(x="UMAP 1",y="UMAP 2",title = i)+
    theme_classic()+
    theme(
      axis.text = element_text(size=12)
    )
  print(p1)
  dev.off()
}

pdf(paste0(in.file.dir,"/plot/joint.projection.pdf"),width=16,height=15)
png(paste0(in.file.dir,"/plot/joint.projection.png"),width=3200,height=2000,res = 300)
p1<-ggplot(to.plot.coord) + 
  geom_segment(aes(x=UMAP_1, y=UMAP_2, color=Celltype, xend=xend, yend=yend), alpha=0.005) + 
  geom_point(aes(x=UMAP_1, y=UMAP_2, color=Celltype), cex=0.5,alpha=0.8) + 
  theme_light() +
  xlab("UMAP 1") + 
  ylab("UMAP2") + 
  #  scale_color_manual(values=c(`H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA`="#D39200", RNA="#DB72FB", umap.H3K27ac="#365FAA",umap.H3K4me1= "#AFCB31" ,umap.H3K27me3="#E6332A",umap.H3K4me3="#FFE609" ,umap.H3K36me3="#7CC291")) + 
  # scale_color_viridis(discrete=TRUE) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=14))
print(p1)
dev.off()


fwrite(to.plot.coord,paste0(in.file.dir,"/plot/joint.projection.tsv"),sep="\t")


####################### plot clusters in integrated umap ########################

pdf(paste0(in.file.dir,"/plot/barplot.num.of.clusters.pdf"),width = 8,height = 6)
ggplot(to.plot[,.N,.(sample_bio.id,RNA_snn_res.0.2)],aes(RNA_snn_res.0.2,N,fill=sample_bio.id )) +
  geom_col()+
  theme_bw()

ggplot(to.plot[,.N,.(sample_bio.id,wknn_res.0.2)],aes(wknn_res.0.2,N,fill=sample_bio.id )) +
  geom_col()+
  theme_bw()

ggplot(to.plot[,.N,.(RNA_snn_res.0.2,wknn_res.0.2)],aes(wknn_res.0.2,N,fill=RNA_snn_res.0.2 )) +
  geom_col()+
  theme_bw()

dev.off()

####################### plot marker genes in integrated umap #####################

DefaultAssay(signac)<-"RNA"

FeaturePlot(signac,reduction = "wnn.umap",features = c("POU5F1",'SOX2',"TBXT","GATA6","SOX17","PDGFRA","PRDM1","FOXA2"))

DefaultAssay(signac)<-"RNA"
DefaultAssay(signac)<-"H3K4me3_genes"
#pluripotency markers
pluri.genes<-c('SOX2',"POU5F1","PRDM14","NANOG")
# primititve streak markers
# but RNA name is TBXT
ps.genes<-c("TBXT","PDGFRA","EOMES","WNT3") 
#endoderm markers
endo.genes<-c("SOX17","GATA6","PRDM1","CER1")#"SHISA2","OTX2",'FOXA2',"HHEX","FZD8"

plot.genes<-ps.genes
"H3K4me3_genes"

#for (assays in c("RNA")){
#  for (plot.genes in c(pluri.genes,ps.genes,endo.genes)){
assays<-"RNA"
DefaultAssay(signac)<-assays
#   pdf(paste0(in.file.dir,"/plot/Gene.plot_",assays,"_",plot.genes,".pdf"),width = 6,height = 5)
plot.genes<-pluri.genes
pdf(paste0(in.file.dir,"/plot/Gene.plot_",assays,"_pluri.genes.pdf"),width = 8,height = 7)
plot.genes<-ps.genes
pdf(paste0(in.file.dir,"/plot/Gene.plot_",assays,"_ps.genes.pdf"),width = 8,height = 7)
plot.genes<-endo.genes
pdf(paste0(in.file.dir,"/plot/Gene.plot_",assays,"_endo.genes.pdf"),width = 8,height = 7)

p1<-VlnPlot(signac,group.by = "sample_bio.id",features = plot.genes,ncol  = 2)
p2<-FeaturePlot(signac,reduction = "wnn.umap",features = plot.genes,ncol = 2)
print(p1)
print(p2)
#RidgePlot(signac, features = plot.genes, ncol = 2)
dev.off()

