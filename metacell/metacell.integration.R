
## generate signac files for definitive endoderm 
#BiocManager::install("EnsDb.Hsapiens.v86")
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

#BiocManager::install('biovizBase')
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
#genome(annotation) <- "hg38"
seqlengths(annotation)

in.file.dir <- "/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments"
file.pattern <- "sorted.bed.gz.tbi"
fragment.list<-list.files(in.file.dir,pattern=file.pattern,full.names = T) %>% sub(".tbi","",.)
fragment.list <-fragment.list[grep("D3",fragment.list)]


read.summary <-fread("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/rna.dna.kept.summary.tsv")
read.summary<-read.summary[lib.type=="RNA"]
#read.summary[,num:=.N,new.cellid]
#read.summary[num>1]
#read.summary[,cellid:=new.cellid]

# from 202306.scMTR-step5.RNA.seurat
scRNA<-readRDS("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/RNA/scRNA.data.seurat.rds")
#scRNA<-readRDS("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/RNA/from.gavin/seurat_rna_complete.rds")
#scRNA<-subset(scRNA,subset= sample_bio.id=="Def.endo.D3")
#scRNA<-subset(scRNA,subset= stage=="D0")
Assays(scRNA)


# RNA was not filter before to remove the ratio lower than 0.5, as different time point RNA could match hist' s new.cellid
scRNA<-subset(scRNA,subset=CB %in% read.summary[,paste0(sub.lib,"_",cellid)] )

#load metacell info
metacell<-fread("/Users/wangy/Desktop/metacell/cell2metacell_assignment.txt.gz") %>% setnames(c("cellid","metacell"))

tmp<-metacell[,.(unique(metacell))] %>% .[,metacell.id:=paste0("metacell",.I)]
metacell<-merge(metacell,tmp[,.(metacell=V1,metacell.id)],by="metacell")

scRNA<-subset(scRNA,subset=CB %in% read.summary[new.cellid%in% metacell$cellid,paste0(sub.lib,"_",cellid)] )
cell.anno<-merge(read.summary[,.(CB=paste0(sub.lib,"_",cellid),new.cellid)],metacell[,.(new.cellid=cellid,metacell,metacell.id)],by="new.cellid") %>% as.data.frame() %>% tibble::column_to_rownames("CB")
scRNA$metacell<-cell.anno[colnames(scRNA), "metacell"]
scRNA$metacell.id<-cell.anno[colnames(scRNA), "metacell.id"]
scRNA@meta.data

## metacell consistent of cells from different time points
ov<-merge(read.summary,metacell,by="cellid")
ov.plot<-ov[,.N,.(metacell.id,sample.type)] %>% .[N>31] %>% dcast(metacell.id~sample.type,value.var = "N")
setkey(ov.plot,Def.endo.D0, Def.endo.D1 ,Def.endo.D2, Def.endo.D3)
ov.plot<-ov.plot%>% as.data.frame() %>% tibble::column_to_rownames("metacell.id")
ov.plot[is.na(ov.plot)]<-0
# order from high to low for each time point, ignor the values less than 31
row.order<-rev(rownames(ov.plot))

ov.plot<-ov[,.N,.(metacell.id,sample.type)] %>% dcast(metacell.id~sample.type,value.var = "N")
ov.plot<-ov.plot%>% as.data.frame() %>% tibble::column_to_rownames("metacell.id")
ov.plot[is.na(ov.plot)]<-0
pdf("/Users/wangy/Desktop/metacell/metacell.heatmap.cellnum.pdf",width = 5,height = 10)
pheatmap::pheatmap(ov.plot[row.order,],cluster_rows = F,cluster_cols = F)
dev.off()

ov.plot<-ov[,.N,.(metacell.id,sample.type)] %>% .[,ratio:=N/sum(N),metacell.id] %>% dcast(metacell.id~sample.type,value.var = "ratio")
ov.plot<-ov.plot%>% as.data.frame() %>% tibble::column_to_rownames("metacell.id")
ov.plot[is.na(ov.plot)]<-0
pdf("/Users/wangy/Desktop/metacell/metacell.heatmap.cellnum.ratio.pdf",width = 5,height = 10)
pheatmap::pheatmap(ov.plot[row.order,],cluster_rows = F,cluster_cols = F)
dev.off()

#scRNA$orig.ident<-scRNA$metacell
#Idents(scRNA)<-"metacell"
# aggregate single cell rna data to metacell data
aggregate.rna <- AggregateExpression(scRNA, group.by = c("metacell.id"), return.seurat = TRUE)
#aggregate.rna@assays$RNA$scale.data

cell.anno<-metacell[,.(metacell,metacell.id)] %>% unique() %>% as.data.frame() %>% tibble::column_to_rownames("metacell.id")

aggregate.rna$metacell<-cell.anno[colnames(aggregate.rna), "metacell"]
Cells(aggregate.rna )
#aggregate.rna$metacell.info<-cell.anno[colnames(aggregate.rna), "metacell"]
aggregate.rna@meta.data
colnames(aggregate.rna)<-aggregate.rna$metacell
aggregate.rna@meta.data



RNAcounts<-aggregate.rna@assays$RNA$counts


fragment.list<-list.files("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments",pattern=".tbi",full.names = T) %>% sub(".tbi","",.)

  frag.path<-fragment.list[1]
  kept.cells<-unique(fread(frag.path)%>% .$V4)

  for (binsizes in c( 5000,20000)) {
   # binsizes=20000 5000,
#AggregateTiles:Quantifies fragment counts per cell in fixed-size genome bins across the whole genome, then removes bins with less than a desired minimum number of counts in the bin, then merges adjacent tiles into a single region.
#GenomeBinMatrix:This function bins the genome and calls FeatureMatrix to construct a bin x cell matrix. 
    
    for (frag.path in fragment.list) {
      #frag.path<-fragment.list[1]
      fragments <- CreateFragmentObject(frag.path,cells = kept.cells)
      counts <- GenomeBinMatrix(
        fragments,
        genome= seqlengths(annotation),
        cells = kept.cells,
        #   min_counts = 0,
        binsize = binsizes,
        verbose = TRUE
      )
      
      saveRDS(counts,paste0(frag.path,"_bin",binsizes/1000,"k.rds"))
    }
    
  }


####
# for generate integrated seurat object, generate each modality with data from different time points, integrate all modalities into one
# creates a Seurat object based on the scRNA-seq data

scRNA <- aggregate.rna


in.file.dir <- "/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments"
file.pattern <- "sorted.bed.gz.tbi"
fragment.list<-list.files(in.file.dir,pattern=file.pattern,full.names = T) %>% sub(".tbi","",.)

file.pattern <- "bin5k.rds"
count.list<-list.files(in.file.dir,pattern=file.pattern,recursive = F,full.names = T) 
i="H3K27me3"
count.list.sub <-count.list[grep(i,count.list)]
kept.cells<-lapply(count.list.sub, function(f) data.table(CB=colnames(readRDS(f)), times=basename(f))) %>% rbindlist() #%>% .[,CB]
kept.cells[duplicated(kept.cells)]
kept.cells<-kept.cells$CB

for (i in c( i="H3K27me3","H3K27ac","H3K4me1","H3K4me3","H3K36me3", "IgG")){ #"H3K27me3",
 # i="H3K27me3"
  signac<-c()
  
  fragment.list.sub <-fragment.list[grep(i,fragment.list)]
  count.list.sub <-count.list[grep(i,count.list)]
    counts <-readRDS(count.list.sub)
    signac <- CreateChromatinAssay(
    counts = counts,
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = fragment.list.sub,
    min.cells = 0,
    min.features = 0
  )
  
  
  ## filter out bins with no detection
  counts <- GetAssayData(signac)   
  genes.percent.expression <- rowMeans(counts>0 )*100   
  genes.filter <- names(genes.percent.expression[genes.percent.expression>0])  
  counts.sub <- counts[genes.filter,]
  
  scRNA[[i]] <- CreateChromatinAssay(
    counts = counts.sub[,colnames(scRNA)],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    min.cells = 0,
    min.features = 0
  )
## add fragments list
  Fragments(scRNA[[i]])<-   signac@fragments
  
}



### Compute H3K27ac gene scores and map to KNN projection ###

for (i in c ( "H3K27me3","H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "IgG" )) {
  #  i<-"H3K27me3"
  # i<-"H3K4me1"
  # DefaultAssay(signac)<-i
  
  DefaultAssay(signac)<-i
  gene.activities<-GeneActivity(signac)
  signac[[paste0(i,"_genes")]]<-CreateAssayObject(counts=gene.activities)
  signac<-NormalizeData(object=signac, assay=paste0(i,"_genes"), normalization.method="LogNormalize", scale.factor= signac@meta.data %>% .[,paste0("nCount_",i,"_genes")] %>% mean() ) # use mean instead median to normalise data
  DefaultAssay(signac)<-paste0(i,"_genes")
}


saveRDS(signac,paste0(in.file.dir,"/seurat.rna.hist_bin5k.metacell.rds"))

##########################
##########################
##########################
rm(scRNA)


signac<-readRDS("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/seurat.rna.hist_bin5k.metacell.rds")

DefaultAssay(signac)<-"RNA"


annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
#genome(annotation) <- "hg38"
seqlengths(annotation)

## perform dimensionality reduction on RNA assay

## only subtract protein coding genes
proteincodesgenes<-annotation %>% as.data.table() %>% .[gene_biotype=="protein_coding",unique(gene_name)]
proteincodesgenes [proteincodesgenes%in% "T"]

signac<-subset(signac,features= c(proteincodesgenes,"TBXT"))

#proteincodesgenes [proteincodesgenes%in% "T"]


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
DimPlot(signac, reduction = "pca",dims = c(1,3)) #,group.by= "seurat_clusters"
FeaturePlot(signac,reduction = "pca",dims = c(1,3), features = c("POU5F1","TBXT","WNT3","EOMES","GATA6","CXCR4","GATA4"))
FeaturePlot(signac,reduction = "pca",dims = c(1,2), features = c("POU5F1","TBXT","WNT3","EOMES","GATA6","CXCR4","GATA4"))
DimPlot(signac, reduction = "pca",dims = c(1,2),group.by= "RNA_snn_res.2")
DimPlot(signac, reduction = "pca",dims = c(1,3),group.by= "RNA_snn_res.2")
dev.off()

dir.create(paste0(in.file.dir,"/plot"))

pdf(paste0(in.file.dir,"/plot/seurat.rna.PCA.pdf"),width = 8,height=6)
DimPlot(signac, reduction = "pca",dims = c(1,2))#,group.by= "sample_bio.id")
dev.off()

DimHeatmap(signac, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(signac, dims = 1:20, cells = 500, balanced = TRUE) #cells = 500, 
# could be time consuming 
signac <- JackStraw(signac, dims = 20, num.replicate = 100)
signac <- ScoreJackStraw(signac, dims = 1:20)
dev.off()

#pdf(paste0(in.file.dir,"/plot/RNA_jackstrawplot.pdf"),width = 5,height = 4)
#JackStrawPlot(signac, dims = 1:20)
#dev.off()

pdf(paste0(in.file.dir,"/RNA_ElbowPlot.pdf"),width = 5,height = 4)
ElbowPlot(signac,ndims = 20)
dev.off()

signac <- FindNeighbors(signac, dims = 1:20)
i<-2
signac <- FindClusters(signac, resolution = i)
head(Idents(signac), 5)
signac <- RunUMAP(signac,
                 # n.neighbors=50,
                  dims = 1:20)

pdf(paste0(in.file.dir,"/plot/RNA_res.metacell.",i,".pdf"),width = 8,height = 6)
p1<-DimPlot(signac, reduction = "umap", label = TRUE)#+NoLegend()
p2<-DimPlot(signac, reduction = "pca",dims = c(1,2), label = TRUE)#+NoLegend()
p3<-DimPlot(signac, reduction = "pca",dims = c(1,3), label = TRUE)#+NoLegend()
#p2<-DimPlot(signac, reduction = "umap",group.by ="sample_bio.id",label = TRUE)#+NoLegend()
#p3<-DimPlot(signac, reduction = "umap",group.by ="sample.anno",label = TRUE)#+NoLegend()
print(p1)
print(p2)
 print(p3)
dev.off()

###add metacell conponent
signac@meta.data<-cbind(signac@meta.data,ov.plot[signac$metacell.id,])


saveRDS(signac,paste0(in.file.dir,"/seurat.rna.hist_bin5k.def.endo.metacell.rds"))

signac<-readRDS(paste0(in.file.dir,"/seurat.rna.hist_bin5k.def.endo.metacell.rds"))




signac<-readRDS(paste0(in.file.dir,"/seurat.rna.hist_bin5k.def.endo.metacell.rds"))

DefaultAssay(signac)<-"RNA"
FeaturePlot(signac,features = c("POU5F1","NANOG","SOX2","TBXT","LEF1","EOMES","GATA4","GATA6"))
VlnPlot(signac,features = c("POU5F1","NANOG","SOX2","TBXT","LEF1","EOMES","GATA4","GATA6"))
DefaultAssay(signac)<-"H3K27ac_genes"
DefaultAssay(signac)<-"H3K4me3_genes"
DefaultAssay(signac)<-"H3K27me3_genes"

DefaultAssay(signac)<-"IgG_genes"
FeaturePlot(signac,reduction = "pca",features = c("POU5F1","NANOG","SOX2","T","LEF1","EOMES","GATA4","GATA6"))
FeaturePlot(signac,features = c("POU5F1","NANOG","SOX2","T","LEF1","EOMES","GATA4","GATA6"))
VlnPlot(signac,features = c("POU5F1","NANOG","SOX2","T","LEF1","EOMES","GATA4","GATA6"))


#signac2<-readRDS("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/fragments/seurat.rna.hist_bin5k.def.endo.rds")
#rm(signac2)
signac@assays$H3K27me3_genes$counts

rna<-melt(as.data.frame(signac@assays$RNA$data )%>% tibble::rownames_to_column("gene") ) %>% setnames(c("gene","cell","rna"))
igg<-melt(as.data.frame(signac@assays$IgG_genes$data +1)%>% tibble::rownames_to_column("gene") ) %>% setnames(c("gene","cell","igg"))
k27me3<-melt(as.data.frame(signac@assays$H3K27me3_genes$data)%>% tibble::rownames_to_column("gene") )%>% setnames(c("gene","cell","h3k27me3"))
k27ac<-melt(as.data.frame(signac@assays$H3K27ac_genes$data)%>% tibble::rownames_to_column("gene") )%>% setnames(c("gene","cell","h3k27ac"))
k4me3<-melt(as.data.frame(signac@assays$H3K4me3_genes$data)%>% tibble::rownames_to_column("gene") )%>% setnames(c("gene","cell","h3k4me3"))
k4me1<-melt(as.data.frame(signac@assays$H3K4me1_genes$data)%>% tibble::rownames_to_column("gene") )%>% setnames(c("gene","cell","h3k4me1"))
k36me3<-melt(as.data.frame(signac@assays$H3K36me3_genes$data)%>% tibble::rownames_to_column("gene") )%>% setnames(c("gene","cell","h3k36me3"))
ov<-merge(igg,k27me3,by=c("gene","cell")) %>% merge(k27ac,by=c("gene","cell"))%>% merge(k4me3,by=c("gene","cell"))%>% merge(k4me1,by=c("gene","cell"))%>% merge(k36me3,by=c("gene","cell")) #%>% merge(rna,by=c("gene","cell"))

ov<-as.data.table(ov)

## pick the highest value of time point as for the meta cell time point
timepoints<-data.frame(cell=colnames(signac),signac$Def.endo.D0, signac$Def.endo.D1,signac$Def.endo.D2,signac$Def.endo.D3)%>% as.data.table() %>% melt() %>% .[,.SD[which.max(value)],cell]

ov<-merge(ov,timepoints[,.(cell,timepoints=variable)])

pc1<-data.frame(colnames(signac),signac@reductions$pca@cell.embeddings[,"PC_1"]) %>% as.data.table() %>% setnames(c("cell","pc1"))
pc1[,order.pca:=rank(pc1)]
ov<-merge(ov,pc1,by="cell")
#ov<-merge(ov,as.data.table(signac2@meta.data[,c("metacell","RNA_snn_res.2",'wnn.pseudotime')]) %>% setnames(c("cell","cluster")), by="cell")
#ov<-merge(ov,as.data.table(data.frame(colnames(signac2),signac2$sample_bio.id ,signac2$wknn_res.0.2,signac2$wnn.pseudotime)) %>% setnames(c("cell","sample_bio.id","cluster","pseudotime")), by="cell")

#ov[,norm.h3k27me3:=h3k27me3/igg]
#ov[,norm.h3k27ac:=h3k27ac/igg]
#ov[,norm.h3k4me3:=h3k4me3/igg]
#ov[,norm.h3k4me1:=h3k4me1/igg]
#ov[,norm.h3k36me3:=h3k36me3/igg]
#ov[cluster=="0",cellanno:="PS"]
#ov[cluster=="1",cellanno:="DE1"]
#ov[cluster=="2",cellanno:="DE2"]
#ov[cluster=="3",cellanno:="PSC"]
#ov$cellanno<-factor(ov$cellanno,levels = c("PSC","PS","DE1","DE2"))

#cell.order<-ov[,.(cell,pseudotime)] %>% unique() %>% setkey(pseudotime)
#cell.order[,"order":=rank(pseudotime)]#sample.type,sub.lib,
#ov[cell=="11:47:03:S6" &gene=="GATA4"]
#ov<-merge(ov,cell.order[,.(cell,order)],by="cell")

#to.plot.hist<-melt(ov,id.vars = c(1:2,9:12) ) 
#to.plot.hist<-melt(ov,id.vars = c(1:2) ) 
#to.plot.hist[,scale:=(value-min(value))/(max(value)-min(value)), .(gene,variable)]


umap.coord<-fread("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/joint.projection.tsv")

to.plot<-merge( umap.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA",.(cell=Cell,UMAP_1,UMAP_2)], ov%>% unique(), by="cell")

#to.plot[is.na(order),sizes:=0.1]
#to.plot[!is.na(order),sizes:=0.5]
#to.plot$cellanno<-factor(to.plot$cellanno,levels = c("NA","DE2","DE1","PS","PSC"))
#to.plot$order<-factor(to.plot$order,levels = c("NA",1:60))
#setkey(to.plot,order)

to.plot[,h3k27me3:=(h3k27me3-min(h3k27me3))/(max(h3k27me3)-min(h3k27me3)),gene]
to.plot[,h3k27ac:=(h3k27ac-min(h3k27ac))/(max(h3k27ac)-min(h3k27ac)),gene]
to.plot[,h3k4me3:=(h3k4me3-min(h3k4me3))/(max(h3k4me3)-min(h3k4me3)),gene]
to.plot[,h3k4me1:=(h3k4me1-min(h3k4me1))/(max(h3k4me1)-min(h3k4me1)),gene]
to.plot[,h3k36me3:=(h3k36me3-min(h3k36me3))/(max(h3k36me3)-min(h3k36me3)),gene]
to.plot[,igg:=(igg-min(igg))/(max(igg)-min(igg)),gene]

to.plot
rna<-as.data.table(rna)
rna[gene=="TBXT"]
rna[gene=="TBXT",gene:="T"]
to.plot<-merge(to.plot,rna,by=c("cell","gene"))
to.plot[,rna:=(rna-min(rna))/(max(rna)-min(rna)),gene]

fwrite(to.plot,"/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.gene.activity.tsv",sep="\t")
fwrite(to.plot,"/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.gene.activity.not.scaled.tsv",sep="\t")
to.plot<-fread("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.gene.activity.tsv")



pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.k27me3.pdf",width = 10,height = 12)

ggplot(umap.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA"],aes(UMAP_1,UMAP_2) )+
  geom_point(color="grey95",size=0.1) +
  geom_point(data = to.plot[gene%in% c("NANOG","POU5F1","T","WNT3","GATA6","PRDM1","SOX17","PDGFRA","FZD8","CXCR4")],aes(color=h3k27me3),size=2)+
  #scale_color_gradientn(colours  = c("grey90",H3K27me3.color))+
  scale_colour_viridis_c(option = "plasma")+
  # scale_colour_viridis()+
 # scale_color_manual(values = c("#A080B9","#2AB6BE","#7CAF2A","#EE756E"), na.value ="grey95" )+
 # geom_text_repel(data=to.plot[!is.na(order),],min.segment.length = 0 ,max.overlaps =100)+
  facet_wrap(.~gene,scales = "free",ncol = 3)+
#  geom_point(cell)
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.k27ac.pdf",width = 10,height = 12)
ggplot(umap.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA"],aes(UMAP_1,UMAP_2) )+
  geom_point(color="grey95",size=0.1) +
  geom_point(data = to.plot[gene%in% c("NANOG","POU5F1","T","WNT3","GATA6","PRDM1","SOX17","PDGFRA","FZD8","CXCR4")],aes(color=h3k27ac),size=2)+
 # scale_color_gradientn(colours  = c("grey90",H3K27ac.color))+
  scale_colour_viridis_c(option = "plasma")+
  # scale_colour_viridis()+
  # scale_color_manual(values = c("#A080B9","#2AB6BE","#7CAF2A","#EE756E"), na.value ="grey95" )+
  # geom_text_repel(data=to.plot[!is.na(order),],min.segment.length = 0 ,max.overlaps =100)+
  facet_wrap(.~gene,scales = "free",ncol = 3)+
  #  geom_point(cell)
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.k4me3.pdf",width = 10,height = 12)

ggplot(umap.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA"],aes(UMAP_1,UMAP_2) )+
  geom_point(color="grey95",size=0.1) +
  geom_point(data = to.plot[gene%in% c("NANOG","POU5F1","T","WNT3","GATA6","PRDM1","SOX17","PDGFRA","FZD8","CXCR4")],aes(color=h3k4me3),size=2)+
  #scale_color_gradientn(colours  = c("grey90",H3K4me3.color))+
  scale_colour_viridis_c(option = "plasma")+
  # scale_color_manual(values = c("#A080B9","#2AB6BE","#7CAF2A","#EE756E"), na.value ="grey95" )+
  # geom_text_repel(data=to.plot[!is.na(order),],min.segment.length = 0 ,max.overlaps =100)+
  facet_wrap(.~gene,scales = "free",ncol = 3)+
  #  geom_point(cell)
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.k4me1.pdf",width = 10,height = 12)

ggplot(umap.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA"],aes(UMAP_1,UMAP_2) )+
  geom_point(color="grey95",size=0.1) +
  geom_point(data = to.plot[gene%in% c("NANOG","POU5F1","T","WNT3","GATA6","PRDM1","SOX17","PDGFRA","FZD8","CXCR4")],aes(color=h3k4me1),size=2)+
  #scale_color_gradientn(colours  = c("grey90",H3K4me1.color))+
  scale_colour_viridis_c(option = "plasma")+
  # scale_color_manual(values = c("#A080B9","#2AB6BE","#7CAF2A","#EE756E"), na.value ="grey95" )+
  # geom_text_repel(data=to.plot[!is.na(order),],min.segment.length = 0 ,max.overlaps =100)+
  facet_wrap(.~gene,scales = "free",ncol = 3)+
  #  geom_point(cell)
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.k36me3.pdf",width = 10,height = 12)

ggplot(umap.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA"],aes(UMAP_1,UMAP_2) )+
  geom_point(color="grey95",size=0.1) +
  geom_point(data = to.plot[gene%in% c("NANOG","POU5F1","T","WNT3","GATA6","PRDM1","SOX17","PDGFRA","FZD8","CXCR4")],aes(color=h3k36me3),size=2)+
  #scale_color_gradientn(colours  = c("grey90",H3K36me3.color))+
   scale_colour_viridis_c(option = "plasma")+
  # scale_color_manual(values = c("#A080B9","#2AB6BE","#7CAF2A","#EE756E"), na.value ="grey95" )+
  # geom_text_repel(data=to.plot[!is.na(order),],min.segment.length = 0 ,max.overlaps =100)+
  facet_wrap(.~gene,scales = "free",ncol = 3)+
  #  geom_point(cell)
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )  
dev.off()
  

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.igg.pdf",width = 10,height = 12)

ggplot(umap.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA"],aes(UMAP_1,UMAP_2) )+
  geom_point(color="grey95",size=0.1) +
  geom_point(data = to.plot[gene%in% c("NANOG","POU5F1","T","WNT3","GATA6","PRDM1","SOX17","PDGFRA","FZD8","CXCR4")],aes(color=igg),size=2)+
  #scale_color_gradientn(colours  = c("grey90",H3K36me3.color))+
  scale_colour_viridis_c(option = "plasma")+
  # scale_color_manual(values = c("#A080B9","#2AB6BE","#7CAF2A","#EE756E"), na.value ="grey95" )+
  # geom_text_repel(data=to.plot[!is.na(order),],min.segment.length = 0 ,max.overlaps =100)+
  facet_wrap(.~gene,scales = "free",ncol = 3)+
  #  geom_point(cell)
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )  
dev.off()

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.rna.pdf",width = 10,height = 12)

ggplot(umap.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA"],aes(UMAP_1,UMAP_2) )+
  geom_point(color="grey95",size=0.1) +
  geom_point(data = to.plot[gene%in% c("NANOG","POU5F1","T","WNT3","GATA6","PRDM1","SOX17","PDGFRA","FZD8","CXCR4")],aes(color=rna),size=2)+
  #scale_color_gradientn(colours  = c("grey90",H3K36me3.color))+
  scale_colour_viridis_c(option = "plasma")+
  # scale_color_manual(values = c("#A080B9","#2AB6BE","#7CAF2A","#EE756E"), na.value ="grey95" )+
  # geom_text_repel(data=to.plot[!is.na(order),],min.segment.length = 0 ,max.overlaps =100)+
  facet_wrap(.~gene,scales = "free",ncol = 3)+
  #  geom_point(cell)
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )  
dev.off()


pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.wnn.umap.pdf",width = 5,height = 4)
ggplot(umap.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA"],aes(UMAP_1,UMAP_2) )+
  geom_point(color="grey95",size=0.1) +
  geom_point(data = to.plot[gene%in% c("NANOG")],aes(color=order),size=2)+
  #scale_color_gradientn(colours  = c("grey90",H3K36me3.color))+
  scale_colour_viridis_c(option = "plasma")+
  # scale_color_manual(values = c("#A080B9","#2AB6BE","#7CAF2A","#EE756E"), na.value ="grey95" )+
   geom_text_repel(data=to.plot[!is.na(order)&gene%in% c("NANOG"),],aes(label=order),max.overlaps =100)+ #min.segment.length = 0 ,
  facet_wrap(.~gene,scales = "free",ncol = 3)+
  #  geom_point(cell)
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.rna.pca.pdf",width = 5,height = 4)
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.rna.pca.no.label.pdf",width = 5,height = 4)
ggplot(umap.coord[Target=="H3K27ac-H3K27me3-H3K4me1-H3K4me3-H3K36me3-RNA"],aes(UMAP_1,UMAP_2) )+
  geom_point(color="grey95",size=0.1) +
  geom_point(data = to.plot[gene%in% c("NANOG")],aes(color=order.pca),size=2)+
  scale_colour_viridis_c(option = "plasma")+
 # geom_text_repel(data=to.plot[!is.na(order)&gene%in% c("NANOG"),],aes(label=order.pca),max.overlaps =100)+ #min.segment.length = 0 ,
  facet_wrap(.~gene,scales = "free",ncol = 3)+
  #  geom_point(cell)
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )  
dev.off()



to.plot.hist<-melt(to.plot,id.vars = c(1:4,11:13) ) 

fwrite(to.plot.hist,"/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.order.pca.tsv" ,sep="\t")
to.plot.hist<-fread("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.order.pca.tsv")

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.sample.anno.pdf",width = 5,height = 2)
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.sample.anno.pca.new.pdf",width = 5,height = 2)
ggplot(to.plot.hist[gene=="POU5F1"& variable=="rna"], aes(order.pca,variable, fill=timepoints) )+
  geom_tile ()+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
dev.off()


plotcolors<-c("h3k27me3"= "#E6332A", 
              "h3k27ac"="#365FAA", 
              "h3k4me1"="#AFCB31",
              "h3k4me3"="#FFE609",
              "h3k36me3"="#7CC291",
              "rna"= "#232b58")


pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.pou5f1.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="POU5F1"& !variable=="igg"],aes( order.pca,value,color=variable))+
 # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+
  scale_color_manual(values = plotcolors)+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.tbxt.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="T"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+

  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.gata6.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="GATA6"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+

  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.cxcr4.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="CXCR4"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+

  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()


pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.wnt3.pdf",width = 6,height = 4)

ggplot(to.plot.hist[gene=="WNT3"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+

  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()


pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.dkk1.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="DKK1"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+

  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.eomes.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="EOMES"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+

  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.gata4.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="GATA4"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+

  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()


pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.sox17.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="SOX17"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+

  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.GSC.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="GSC"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+

  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.prdm14.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="PRDM14"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+
  
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell.pseudotime.sox2.pdf",width = 6,height = 4)
ggplot(to.plot.hist[gene=="SOX2"& !variable=="igg"],aes( order.pca,value,color=variable))+
  # geom_point(aes(order.pca,cluster))+
  geom_smooth(
    stat = "smooth",#method = "loess",
    position = "identity",se = F
  )+  scale_color_manual(values = plotcolors)+
  
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

str(p1$scales$scales)

str(p1)
p1$layers[[1]]$aes_params


DimPlot(signac,reduction = "pca",dims = c(1,3))
DimPlot(signac,reduction = "umap")

tmp<-data.frame(signac@meta.data$RNA_snn_res.2, signac@reductions$pca@cell.embeddings[,"PC_1"]) %>% as.data.table() %>% setnames(c("cluster","pc1"))

ggplot(tmp,aes(pc1,cluster))+
  geom_point()

DefaultAssay(signac)<-"RNA"
signac@active.ident<-signac$RNA_snn_res.2
scRNA2.markers <- FindAllMarkers(signac,   logfc.threshold = 0.1)
fwrite(scRNA2.markers %>% as.data.table(),"/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell_res2.markers.tsv",sep = "\t")
# overlap with protein coding genes?
#head(VariableFeatures(signac),n=1000)

scRNA2.markers<-scRNA2.markers%>% as.data.table()
scRNA2.markers[,.N,cluster]
scRNA2.markers[gene%in%"EOMES"]
scRNA2.markers[gene%in%"POU5F1"]
#top.genes
scRNA2.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = 300, order_by = abs(avg_log2FC)) ->top.genes

top.genes<-as.data.table(top.genes)
to.plot.heat["GATA4",]


to.plot[is.na(h3k27me3),unique(gene)]
to.plot[is.na(rna),unique(gene)]
###################
###################cor test for each gene
###################
###################
out.cor<-data.table()
cor.tmp<-data.table()

"h3k27ac","h3K27me3","h3k4me1","h3k4me3","h3k36me3"
for (i in to.plot[!gene %in% c(to.plot[is.na(h3k36me3),unique(gene)],
                  to.plot[is.na(rna),unique(gene)])
                  ,unique(gene)]) {
cor.hist.rna<-cor.test(to.plot[gene==i,h3k36me3],to.plot[gene==i,rna])
#cor.hist.rna$estimate
#cor.hist.rna$p.value
cor.tmp$gene<- i
cor.tmp[gene== i,c("cor","pvalue","anno"):=list(cor.hist.rna$estimate,cor.hist.rna$p.value, "h3k36me3")]
out.cor<-rbind(out.cor,cor.tmp)
}

fwrite(out.cor,"/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell_cor.test.tsv",sep = "\t")

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell_cor.plot.genes.pdf",width = 5,height = 4)
ggplot(to.plot[gene=="GATA4"],aes(h3k27me3,rna,color=order.pca))+
  geom_point(size=1)+
  geom_smooth()+
  scale_color_viridis_c()+  labs(title = "GATA4")+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
ggplot(to.plot[gene=="EOMES"],aes(h3k27me3,rna,color=order.pca))+
  geom_point(size=1)+
  geom_smooth()+
  scale_color_viridis_c()+  labs(title = "EOMES")+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
ggplot(to.plot[gene=="GATA6"],aes(h3k27me3,rna,color=order.pca))+
  geom_point(size=1)+
  geom_smooth()+
  scale_color_viridis_c()+
  labs(title = "GATA6")+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
ggplot(to.plot[gene=="DKK1"],aes(h3k27me3,rna,color=order.pca))+
  geom_point(size=1)+
  geom_smooth()+
  scale_color_viridis_c()+  labs(title = "DKK1")+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
ggplot(to.plot[gene=="WNT3"],aes(h3k27me3,rna,color=order.pca))+
  geom_point(size=1)+
  geom_smooth()+
  scale_color_viridis_c()+  labs(title = "WNT3")+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
ggplot(to.plot[gene=="SOX17"],aes(h3k27me3,rna,color=order.pca))+
  geom_point(size=1)+
  geom_smooth()+
  scale_color_viridis_c()+  labs(title = "SOX17")+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
ggplot(to.plot[gene=="GSC"],aes(h3k27me3,rna,color=order.pca))+
  geom_point(size=1)+
  geom_smooth()+
  scale_color_viridis_c()+  labs(title = "GSC")+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()


cor.test(to.plot$h3k27ac,to.plot$rna)

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell_cor.plot.pdf",width = 5,height = 20)
ggplot(out.cor,aes(cor,-log10(pvalue),label=gene))+
  geom_point()+facet_grid(anno~.,scales = "free")+
  geom_point(data=out.cor[ gene %in% c("EOMES","GATA4","GATA6","SOX17","GSC","FOXA2","CXCR4") ] ,aes(cor,-log10(pvalue), color="red"))+
  geom_text_repel(data=out.cor[gene %in% c("EOMES","GATA4","GATA6","SOX17","GSC","FOXA2","CXCR4") ], min.segment.length=0 )+
  
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell_cor.h3k27me3.plot.pdf",width = 5,height = 4)
ggplot(out.cor[anno=="h3k27me3"],aes(cor,-log10(pvalue) ,label=gene))+
  geom_point()+
  geom_point(data=out.cor[anno=="h3k27me3"& gene %in% c("EOMES","GATA4","GATA6","SOX17","GSC","WNT3","DKK1","MIXL1","FOXA2","CXCR4") ] ,aes(cor,-log10(pvalue), color="red"))+
  geom_text_repel(data=out.cor[anno=="h3k27me3"& gene %in% c("EOMES","GATA4","GATA6","SOX17","GSC","WNT3","DKK1","MIXL1","FOXA2","CXCR4") ], min.segment.length=0,  )+
    theme_bw()+
    theme(
      panel.grid = element_blank()
    )
  

dev.off()

out.cor[anno=="h3k27me3" &(cor)<(-0.7)]

out.cor[anno=="h3k27me3" &(-log10(pvalue))>10]
out.cor[gene=="T"]
out.cor[gene=="GATA4"]
out.cor[gene %in% rownames(to.plot.heat) &anno=="h3k27me3" & cor < (-0.5)]

cor.plot<-out.cor[gene %in% rownames(to.plot.heat)] %>% dcast(gene~anno,value.var = "cor") %>% as.data.frame() %>% tibble::column_to_rownames("gene")
cor.plot[abs(cor.plot)<0.5]<-0
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/metacell_cor.plot.heatmap.pdf",width = 5,height = 10)
pheatmap::pheatmap(cor.plot[p1$tree_row$labels[p1$tree_row$order], ],
                   cluster_cols = FALSE,
                   cluster_rows=F,
                   cutree_rows = 8,
                   gaps_row = c(116,198,229,307,331,508,541,653),
                   #   gaps_row = c(117,419,454,229,307,331,508,541,653),
                   # clustering_method = "ward.D2",
                  # annotation_row=ann[p1$tree_row$labels[p1$tree_row$order],], 
                   color = colorRampPalette(c("black", "grey95", "orangered"))(3),
                   fontsize = 1)
dev.off()







to.plot.heat<-to.plot.hist[gene%in%c(top.genes$gene,"T") &variable=="rna"] %>% dcast(gene~order.pca,value.var = "value") %>% as.data.frame() %>% tibble::column_to_rownames("gene")

#rownames((to.plot.heat)[rowMeans(to.plot.heat)>0.1,] )
#rownames((to.plot.heat)[row(to.plot.heat)<0.1,] )

to.plot.heat["FGF4",]
to.plot.heat[is.na(to.plot.heat)]<-0
#[rownames((to.plot.heat)[rowMeans(to.plot.heat)<0.05,] ),]

p1<-pheatmap::pheatmap(to.plot.heat, cluster_cols = F,clustering_method = "ward.D2",cutree_rows = 8,show_rownames = T,fontsize = 6)
to.plot.heat <- smoother_aggregate_nearest_nb(mat=as.matrix(to.plot.heat), D=c(1:60), k=5)

p1$tree_row$labels[p1$tree_row$order]

cl<-cutree(p1$tree_row,8)
ann = data.frame(cl)
table(ann$cl)
#1   2   3   4   5   6   7   8 
#116 112 177  78  82  31  33  24 
# without rna scale 
#1   2   3   4   5   6   7   8 
#81 117 302  35  18  65  29   6 

## final used due to the markers got different with default, top 300 per cluster got diferent results. but the pattern are overall similar, 

#. 1   2   3   4   5   6   7   8 
# 148 185 124 182  96  30  41  29 
rownames(ann) = rownames(to.plot.heat)
#2, 1, 6,8,5,4, 7,3
to.plot.heat <- smoother_aggregate_nearest_nb(mat=as.matrix(to.plot.heat), D=c(1:60), k=5)
pheatmap::pheatmap(to.plot.heat[p1$tree_row$labels[p1$tree_row$order], ],
                   cluster_cols = FALSE,
                   cluster_rows=F,
                   cutree_rows = 8,
                   gaps_row = c(185,333,363,392,488,670,711, 835),
                #   gaps_row = c(117,419,454,229,307,331,508,541,653),
                   # clustering_method = "ward.D2",
                   annotation_row=ann, 
                   color = colorRampPalette(c("black", "darkblue", "yellow"))(100),
                   border_color = NA ,
                   show_rownames = T,
                   na_col = "white",
                   fontsize = 6)
#nrow(ann)

names(cl[cl==1])
p1$tree_row$merge

H3K27me3.color="#E6332A"
H3K27ac.color="#365FAA"
H3K4me1.color="#AFCB31"
H3K4me3.color="#FFE609"
H3K36me3.color="#7CC291"
IgG.color="grey60"

#re.odrder
#. 1   2   3   4   5   6   7   8 
# 148 185 124 182  96  30  41  29 
rownames(ann) = rownames(to.plot.heat)
#2, 1, 6,8,5,4, 7,3

#4, 7, 3, 2,6,8,5,1
#182, 41, 124, 185, 30, 29, 96, 148
#182, 223, 347, 532, 562, 591, 687, 835
## as H3K27me3 looks different, we reorder the rows by the H3K27me3 for each cluster 
to.plot.heat<-to.plot.hist[gene%in%c(top.genes$gene,"T")  &variable=="h3k27me3"] %>% dcast(gene~order.pca,value.var = "value") %>% as.data.frame() %>% tibble::column_to_rownames("gene")
to.plot.heat[is.na(to.plot.heat)]<-0
to.plot.heat <- smoother_aggregate_nearest_nb(mat=as.matrix(to.plot.heat), D=c(1:60), k=5)
pheatmap::pheatmap(to.plot.heat[p1$tree_row$labels[p1$tree_row$order],], cluster_cols = F,cluster_rows = F, color = colorRampPalette(c("black", "darkblue", "yellow"))(100),border_color = NA ,show_rownames = F,na_col = "white")
#pheatmap::pheatmap(to.plot.heat[c("POU5F1","T","CXCR4"),], cluster_cols = F,cluster_rows = F, color = colorRampPalette(c("black", "darkblue", "yellow"))(100),border_color = NA ,show_rownames = F,na_col = "white")
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/h3k27me3.heatmaps.pdf",height=15,width=6)

h1<-pheatmap::pheatmap(to.plot.heat[names(cl[cl==1]),], cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("grey80","white", H3K27me3.color))(100),border_color = NA ,show_rownames = T,na_col = "white",fontsize = 6)
h2<-pheatmap::pheatmap(to.plot.heat[names(cl[cl==2]),], cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("grey80","white", H3K27me3.color))(100),border_color = NA ,show_rownames = T,na_col = "white",fontsize = 6)
h3<-pheatmap::pheatmap(to.plot.heat[names(cl[cl==3]),], cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("grey80","white", H3K27me3.color))(100),border_color = NA ,show_rownames = T,na_col = "white",fontsize = 6)
h4<-pheatmap::pheatmap(to.plot.heat[names(cl[cl==4]),], cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("grey80","white", H3K27me3.color))(100),border_color = NA ,show_rownames = T,na_col = "white",fontsize = 6)
h5<-pheatmap::pheatmap(to.plot.heat[names(cl[cl==5]),], cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("grey80","white", H3K27me3.color))(100),border_color = NA ,show_rownames = T,na_col = "white",fontsize = 6)
h6<-pheatmap::pheatmap(to.plot.heat[names(cl[cl==6]),], cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("grey80","white", H3K27me3.color))(100),border_color = NA ,show_rownames = T,na_col = "white",fontsize = 6)
h7<-pheatmap::pheatmap(to.plot.heat[names(cl[cl==7]),], cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("grey80","white", H3K27me3.color))(100),border_color = NA ,show_rownames = T,na_col = "white",fontsize = 6)
h8<-pheatmap::pheatmap(to.plot.heat[names(cl[cl==8]),], cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("grey80","white", H3K27me3.color))(100),border_color = NA ,show_rownames = T,na_col = "white",fontsize = 6)

h9<-pheatmap::pheatmap(to.plot.heat[c(h1$tree_row$labels[h1$tree_row$order],
                                  h5$tree_row$labels[h5$tree_row$order],
                                  h6$tree_row$labels[h6$tree_row$order],
                                  h4$tree_row$labels[h4$tree_row$order],
                                  h8$tree_row$labels[h8$tree_row$order],
                                  h3$tree_row$labels[h3$tree_row$order],
                                  h7$tree_row$labels[h7$tree_row$order],
                                  h2$tree_row$labels[h2$tree_row$order]
                                  ), ],
                   cluster_cols = FALSE,
                   cluster_rows=F,
                   cutree_rows = 8,
                   gaps_row = c(185,333,363,392,488,670,711, 835),
                   # clustering_method = "ward.D2",
                   color = colorRampPalette(c("black", "darkblue", "yellow"))(100),
                   border_color = NA ,
                   show_rownames = T,
                   na_col = "white",
                   fontsize = 6)
dev.off()
top.genes$gene %in%rownames(to.plot.heat)
pheatmap::pheatmap(to.plot.heat[c(top.genes$gene[top.genes$gene %in%k27me3.marked] ,"T"), ],
                   cluster_cols = FALSE,
                   cluster_rows=F,color = colorRampPalette(c("black", "darkblue", "yellow"))(100),
                   border_color = NA ,
                   show_rownames = T,
                   na_col = "white",
                   fontsize = 6)
                   

#4, 7, 3, 2,6,8,5,1
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/h3k27me3.heatmap.order.pdf",height=8,width=6)
#pheatmap::pheatmap(to.plot.heat[c(names(cl[cl==3]),names(cl[cl==7]),names(cl[cl==8]),names(cl[cl==2]),names(cl[cl==6]),names(cl[cl==5]),names(cl[cl==4]),names(cl[cl==1])), ],
pheatmap::pheatmap(to.plot.heat[c(names(cl[cl==4]),names(cl[cl==7]),names(cl[cl==3]),names(cl[cl==2]),names(cl[cl==6]),names(cl[cl==8]),names(cl[cl==5]),names(cl[cl==1])), ],
                                      cluster_cols = FALSE,
                   cluster_rows=F,
                   cutree_rows = 8,
                   #gaps_row = c(177,210,234,346,377,459,537),
                   gaps_row = c(182, 223, 347, 532, 562, 591, 687, 835),
# clustering_method = "ward.D2",
color = colorRampPalette(c("black", "darkblue", "yellow"))(100),
border_color = NA ,
show_rownames = F,
na_col = "white",
fontsize = 6
)
dev.off()

pheatmap::pheatmap(to.plot.heat[c(rownames(to.plot.heat)[rownames(to.plot.heat) %in%k27me3.marked ],"T"), ],
                   cluster_cols = FALSE,
                   cluster_rows=T,
                   # clustering_method = "ward.D2",
                   color = colorRampPalette(c("black", "darkblue", "yellow"))(100),
                   border_color = NA ,
                   show_rownames = F,
                   na_col = "white",
                   fontsize = 6
)
rownames(to.plot.heat)[rownames(to.plot.heat) %in%k27me3.marked ]



#pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/h3k27ac.heatmap.pdf",height=8,width=6)
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/h3k27ac.heatmap.order.pdf",height=8,width=6)
to.plot.heat<-to.plot.hist[gene%in%c(top.genes$gene,"T")  &variable=="h3k27ac"] %>% dcast(gene~order.pca,value.var = "value") %>% as.data.frame() %>% tibble::column_to_rownames("gene")
to.plot.heat[is.na(to.plot.heat)]<-0
to.plot.heat <- smoother_aggregate_nearest_nb(mat=as.matrix(to.plot.heat), D=c(1:60), k=5)
pheatmap::pheatmap(to.plot.heat[c(names(cl[cl==4]),names(cl[cl==7]),names(cl[cl==3]),names(cl[cl==2]),names(cl[cl==6]),names(cl[cl==8]),names(cl[cl==5]),names(cl[cl==1])), ],
                   cluster_cols = FALSE,
                   cluster_rows=F,
                   cutree_rows = 8,
                   #gaps_row = c(177,210,234,346,377,459,537),
                   gaps_row = c(182, 223, 347, 532, 562, 591, 687, 835),
                   # clustering_method = "ward.D2",
                   color = colorRampPalette(c("black", "darkblue", "yellow"))(100),
                   border_color = NA ,
                   show_rownames = F,
                   na_col = "white",
                   fontsize = 6
)
dev.off()

#pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/h3k4me3.heatmap.pdf",height=8,width=6)
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/h3k4me3.order.heatmap.pdf",height=8,width=6)
to.plot.heat<-to.plot.hist[gene%in%c(top.genes$gene,"T")  &variable=="h3k4me3"] %>% dcast(gene~order.pca,value.var = "value") %>% as.data.frame() %>% tibble::column_to_rownames("gene")
to.plot.heat[is.na(to.plot.heat)]<-0
to.plot.heat <- smoother_aggregate_nearest_nb(mat=as.matrix(to.plot.heat), D=c(1:60), k=5)
pheatmap::pheatmap(to.plot.heat[c(names(cl[cl==4]),names(cl[cl==7]),names(cl[cl==3]),names(cl[cl==2]),names(cl[cl==6]),names(cl[cl==8]),names(cl[cl==5]),names(cl[cl==1])), ],
                   cluster_cols = FALSE,
                   cluster_rows=F,
                   cutree_rows = 8,
                   #gaps_row = c(177,210,234,346,377,459,537),
                   gaps_row = c(182, 223, 347, 532, 562, 591, 687, 835),
                   # clustering_method = "ward.D2",
                   color = colorRampPalette(c("black", "darkblue", "yellow"))(100),
                   border_color = NA ,
                   show_rownames = F,
                   na_col = "white",
                   fontsize = 6
)
dev.off()

#pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/h3k4me1.heatmap.pdf",height=8,width=6)
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/h3k4me1.heatmap.order.pdf",height=8,width=6)
to.plot.heat<-to.plot.hist[gene%in%c(top.genes$gene,"T")  &variable=="h3k4me1"] %>% dcast(gene~order.pca,value.var = "value") %>% as.data.frame() %>% tibble::column_to_rownames("gene")
to.plot.heat[is.na(to.plot.heat)]<-0
to.plot.heat <- smoother_aggregate_nearest_nb(mat=as.matrix(to.plot.heat), D=c(1:60), k=5)
pheatmap::pheatmap(to.plot.heat[c(names(cl[cl==4]),names(cl[cl==7]),names(cl[cl==3]),names(cl[cl==2]),names(cl[cl==6]),names(cl[cl==8]),names(cl[cl==5]),names(cl[cl==1])), ],
                   cluster_cols = FALSE,
                   cluster_rows=F,
                   cutree_rows = 8,
                   #gaps_row = c(177,210,234,346,377,459,537),
                   gaps_row = c(182, 223, 347, 532, 562, 591, 687, 835),
                   # clustering_method = "ward.D2",
                   color = colorRampPalette(c("black", "darkblue", "yellow"))(100),
                   border_color = NA ,
                   show_rownames = F,
                   na_col = "white",
                   fontsize = 6
)
dev.off()

#pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/h3k36me3.heatmap.pdf",height=8,width=6)
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/h3k36me3.heatmap.order.pdf",height=8,width=6)
to.plot.heat<-to.plot.hist[gene%in%c(top.genes$gene,"T")  &variable=="h3k36me3"] %>% dcast(gene~order.pca,value.var = "value") %>% as.data.frame() %>% tibble::column_to_rownames("gene")
to.plot.heat[is.na(to.plot.heat)]<-0
to.plot.heat <- smoother_aggregate_nearest_nb(mat=as.matrix(to.plot.heat), D=c(1:60), k=5)
pheatmap::pheatmap(to.plot.heat[c(names(cl[cl==4]),names(cl[cl==7]),names(cl[cl==3]),names(cl[cl==2]),names(cl[cl==6]),names(cl[cl==8]),names(cl[cl==5]),names(cl[cl==1])), ],
                   cluster_cols = FALSE,
                   cluster_rows=F,
                   cutree_rows = 8,
                   #gaps_row = c(177,210,234,346,377,459,537),
                   gaps_row = c(182, 223, 347, 532, 562, 591, 687, 835),                   # clustering_method = "ward.D2",
                   color = colorRampPalette(c("black", "darkblue", "yellow"))(100),
                   border_color = NA ,
                   show_rownames = F,
                   na_col = "white",
                   fontsize = 6
)
dev.off()



#pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/rna.heatmap.pdf",height=8,width=6)
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/rna.heatmap.withrownames.pdf",height=8,width=6)
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/rna.heatmap.order.pdf",height=8,width=6)
to.plot.heat<-to.plot.hist[gene%in%c(top.genes$gene,"T")  &variable=="rna"] %>% dcast(gene~order.pca,value.var = "value") %>% as.data.frame() %>% tibble::column_to_rownames("gene")
to.plot.heat[is.na(to.plot.heat)]<-0
to.plot.heat <- smoother_aggregate_nearest_nb(mat=as.matrix(to.plot.heat), D=c(1:60), k=5)
pheatmap::pheatmap(to.plot.heat[c(names(cl[cl==4]),names(cl[cl==7]),names(cl[cl==3]),names(cl[cl==2]),names(cl[cl==6]),names(cl[cl==8]),names(cl[cl==5]),names(cl[cl==1])), ],
                   cluster_cols = FALSE,
                   cluster_rows=F,
                   cutree_rows = 8,
                   #gaps_row = c(177,210,234,346,377,459,537),
                   gaps_row = c(182, 223, 347, 532, 562, 591, 687, 835),
                   # clustering_method = "ward.D2",
                   color = colorRampPalette(c("black", "darkblue", "yellow"))(100),
                   border_color = NA ,
                   show_rownames = T,
                   na_col = "white",
                   fontsize = 1
)
dev.off()



#pheatmap::pheatmap(to.plot.heat[names(cl[cl==4]),], cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("black", "darkblue", "yellow"))(100),border_color = NA ,na_col = "white")
#pheatmap::pheatmap(to.plot.heat, cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("black", "darkblue", "yellow"))(100),border_color = NA ,na_col = "white")
#pheatmap::pheatmap(to.plot.heat[c("POU5F1","T","CXCR4"),], cluster_cols = F,cluster_rows = F)





rna.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(to.plot.heat[c("POU5F1","T","CXCR4"),]), D=c(1:60), k=5)
pheatmap::pheatmap(rna.mtx[c("POU5F1","T","CXCR4"),], cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("black", "darkblue", "yellow"))(100),border_color = NA ,na_col = "white")


########### smooting function from gavin_LI modified nb_cid to use linear pseudotime
smoother_aggregate_nearest_nb <- function(mat, D, k){
  # @param mat A matrix in a shape of #genes x #samples.
  # @param D A predefined distance matrix in a shape of #samples x #samples.
  # @param k An integer to choose \code{k} nearest samples (self-inclusive) to
  #  aggregate based on the distance matrix \code{D}.
  denoised_mat <- sapply(seq_len(ncol(mat)), function(cid){
#    nb_cid <- head(order(D[cid, ]), k)
    nb_cid<- head(D[D>=cid],k)  
    closest_mat <- mat[, nb_cid, drop=FALSE]
    # return(Matrix::rowSums(closest_mat))
    return(Matrix::rowMeans(closest_mat))
  })
  dimnames(denoised_mat) <- dimnames(mat)
  return(denoised_mat)
}




inputdir<-"/Users/wangy/data/scMTR_method.dev/figure4.data/chromhmm_adj/states_15"

files<-list.files(inputdir,pattern = "segments.bed",full.names = T)

foo<-lapply(files,function(i) fread(i) %>% setnames(c("chr","start","end","states")) %>% .[,peaks:=paste0("peaks_",.I)] %>% .[,anno:=basename(i) %>% sub("_15_segments.bed","",.)] ) %>% rbindlist()
setkey(foo,chr,start,end)
kept.k27me3<-foo[states%in% c("E14","E12","E13")] %>% .[,.(chr,start,end)] %>% unique() %>% .[,peaks:=paste0("peaks_",.I)]
#foo[states %in% c(1,3), state.type:="Transcription"]
#foo[states %in% c( 2,4,8,11,15), state.type:="Unmarked"]
#foo[states %in% c(5), state.type:="Primed enhancer"]
#foo[states %in% c(6,7), state.type:="Active enhancer"]
#foo[states %in% c(9,10), state.type:="Active promoter"]
#foo[states %in% c(12,13), state.type:="Bivalent promoter"]
#foo[states %in% c(14), state.type:="Repressed Polycomb"] #Bivalent enhancer

fwrite(kept.k27me3,paste0(inputdir,"/all.bivalent.polycomb.chromtin.states.tsv"),sep="\t",col.names = F)
fread("/Users/wangy/data/scMTR_method.dev/figure4.data/chromhmm_adj/states_15/D0_15_segments.bed")


k27me3.marked<-fread("/Users/wangy/data/scMTR_method.dev/figure4.data/chromhmm_adj/states_15/all.bivalent.polycomb.chromtin.states.tsv.annotat.xls",
                     skip ="PeakID",select = c(2:5,8,10,15,16))%>% .[substring(Chr,1,3)=="chr"] %>% .[abs(`Distance to TSS`)<1000] %>% .[,unique(`Gene Name`)]

rownames(to.plot.heat)[rownames(to.plot.heat) %in%k27me3.marked ]










  
  
  
  
  
  
  
  
  
  


