
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
################################# calculate FRIP ###############################
################################################################################

################################ de nova MACS peaks ############################

for (i in c("H3K27me3", "H3K27ac","H3K4me1","H3K4me3","H3K36me3")){
  # i<- "IgG"   i<-"H3K36me3"
  DefaultAssay(signac)<- i
  #call peaks  
  peaks <- CallPeaks(
    object = signac,
    macs2.path="/opt/miniconda3/envs/signac/bin/macs2",
    group.by = "celltype"
  )
  #save peaks
  saveRDS(peaks,paste0(in.file.dir,"/",i,".macs2.peaks.rds"))
  #readRDS(paste0(in.file.dir,"/",i,".macs2.peaks.rds"))
  
  # feature count generate matrix 
  fragment.list<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments",
                            pattern="sorted.bed.gz.tbi",full.names = T) %>% sub(".tbi","",.)
  
  ab.type<-fragment.list[grep(i,fragment.list)]
  
  fragments <- CreateFragmentObject(ab.type,cells = kept.cells)
  counts<-FeatureMatrix(
    fragments=fragments,
    features=peaks,
    cells = kept.cells,
    process_n = 2000,
    sep = c("-", "-"),
    verbose = TRUE
  )
  signac[[paste0("Peaks_",i)]] <- CreateChromatinAssay(
    counts = counts[,rownames(signac@meta.data)],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 0,
    min.features = 0
  )
  # manual cal FRiP
  signac@meta.data[,paste0("FRiP_",i)] <- round(signac@meta.data[,paste0("nCount_Peaks_",i)] /signac@meta.data[,paste0("nCount_",i)] *10000)/100
  png(paste0(in.file.dir,"/",i,".macs2.peaks.FRiP.png"),res = 150, width = 800, height = 800)
  p1<-DensityScatter(signac, x = paste0("nCount_",i ), y = paste0('FRiP_',i), log_x = TRUE, quantiles = TRUE)
  print (p1)
  dev.off()
  
}

############################## ENCODE ChIP-seq peaks ###########################
for (i in c("H3K27me3", "H3K27ac","H3K4me1","H3K4me3","H3K36me3")){
  # i<- "IgG"
  DefaultAssay(signac)<- i
  #call peaks   load peaks from encode
  peaks <- list.files("/Users/yang/data/scMTR_method.dev/ref.peaks/ref.peaks", pattern=i,full.name=T)
  peaks<-fread(peaks[1])
  peaks<-GRanges(seqnames = peaks$V1,
                 ranges = IRanges(start =peaks$V2,end =peaks$V3),
                 strand = "*")
  # feature count generate matrix 
  
  fragment.list<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments",
                            pattern="sorted.bed.gz.tbi",full.names = T) %>% sub(".tbi","",.)
  ab.type<-fragment.list[grep(i,fragment.list)]
  fragments <- CreateFragmentObject(ab.type,cells = kept.cells)
  counts<-FeatureMatrix(
    fragments=fragments,
    features=peaks,
    cells = kept.cells,
    process_n = 2000,
    sep = c("-", "-"),
    verbose = TRUE
  )
  signac[[paste0("Peaks_ChIP_",i)]] <- CreateChromatinAssay(
    counts = counts[,rownames(signac@meta.data)],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 0,
    min.features = 0
  )
  # manual cal FRiP
  signac@meta.data[,paste0("FRiP_ChIP_",i)] <- round(signac@meta.data[,paste0("nCount_Peaks_ChIP_",i)] /signac@meta.data[,paste0("nCount_",i)] *10000)/100
  
  png(paste0(in.file.dir,"/",i,".ChIP.peaks.FRiP.png"),res = 150, width = 800, height = 800)
  p1<-DensityScatter(signac, x = paste0("nCount_",i ), y = paste0('FRiP_ChIP_',i), log_x = TRUE, quantiles = TRUE)
  print (p1)
  dev.off()
  
}

############################# de novo macs broad peaks #########################

for (i in c("H3K27me3", "H3K27ac","H3K4me1","H3K4me3","H3K36me3","IgG")){
  # i<- ""   i<-"H3K36me3"
  
  DefaultAssay(signac)<- i
  
  #call peaks  
  peaks <- CallPeaks(
    object = signac,
    broad = T,
    macs2.path="/opt/miniconda3/envs/signac/bin/macs2",
    group.by = "celltype"
  )
  #save peaks
  saveRDS(peaks,paste0(in.file.dir,"/",i,".macs2.broad.peaks.rds"))
  #readRDS(paste0(in.file.dir,"/",i,".macs2.peaks.rds"))
  # feature count generate matrix 
  
  fragment.list<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments",
                            pattern="sorted.bed.gz.tbi",full.names = T) %>% sub(".tbi","",.)
  
  ab.type<-fragment.list[grep(i,fragment.list)]
  
  fragments <- CreateFragmentObject(ab.type,cells = kept.cells)
  
  counts<-FeatureMatrix(
    fragments=fragments,
    features=peaks,
    cells = kept.cells,
    process_n = 2000,
    sep = c("-", "-"),
    verbose = TRUE
  )
  
  signac[[paste0("Peaks.broad_",i)]] <- CreateChromatinAssay(
    counts = counts[,rownames(signac@meta.data)],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 0,
    min.features = 0
  )
  # manual cal FRiP
  signac@meta.data[,paste0("FRiP.broad_",i)] <- round(signac@meta.data[,paste0("nCount_Peaks.broad_",i)] /signac@meta.data[,paste0("nCount_",i)] *10000)/100
  
  png(paste0(in.file.dir,"/",i,".macs2.broad.peaks.FRiP.png"),res = 150, width = 800, height = 800)
  p1<-DensityScatter(signac, x = paste0("nCount_",i ), y = paste0('FRiP.broad_',i), log_x = TRUE, quantiles = TRUE)
  print (p1)
  dev.off()
  
}

saveRDS(signac, paste0(in.file.dir,"/seurat.rna.hist_bin5k.processed.rds"))

################################## CUT&Tag peaks ###############################
for (i in c("H3K27me3", "H3K27ac","H3K4me1","H3K4me3","H3K36me3")){
  # i<- "IgG" i<-"H3K27me3"
  DefaultAssay(signac)<- i
  
  #load peaks from bulk cut &tag 
  
  peaks <- list.files("/Users/yang/data/scMTR_method.dev/bulk.cuttag/fragments/SEACR", pattern=i,full.name=T)
  
  peaks<-fread(peaks[grep("0.01",peaks)])
  
  peaks<-GRanges(seqnames = peaks$V1,
                 ranges = IRanges(start =peaks$V2,end =peaks$V3),
                 strand = "*")
  # feature count generate matrix 
  
  fragment.list<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments",
                            pattern="sorted.bed.gz.tbi",full.names = T) %>% sub(".tbi","",.)
  
  ab.type<-fragment.list[grep(i,fragment.list)]
  
  fragments <- CreateFragmentObject(ab.type,cells = kept.cells)
  
  counts<-FeatureMatrix(
    fragments=fragments,
    features=peaks,
    cells = kept.cells,
    process_n = 2000,
    sep = c("-", "-"),
    verbose = TRUE
  )
  
  signac[[paste0("Peaks_CT_",i)]] <- CreateChromatinAssay(
    counts = counts[,rownames(signac@meta.data)],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 0,
    min.features = 0
  )
  # manual cal FRiP
  signac@meta.data[,paste0("FRiP_CT_",i)] <- round(signac@meta.data[,paste0("nCount_Peaks_CT_",i)] /signac@meta.data[,paste0("nCount_",i)] *10000)/100
  
  png(paste0(in.file.dir,"/",i,".CT.peaks.FRiP.png"),res = 150, width = 800, height = 800)
  p1<-DensityScatter(signac, x = paste0("nCount_",i ), y = paste0('FRiP_CT_',i), log_x = TRUE, quantiles = TRUE)
  print (p1)
  dev.off()
  
}

############################ de novo SEACR relaxed peaks ######################

for (i in c("H3K27me3","H3K4me3")){
  # i<- "IgG" i<-"H3K27me3" #IgG 
  DefaultAssay(signac)<- i
  
  #call peaks   load peaks from encode
  #
  peaks <- list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/SEACR", pattern=".relaxed.bed",full.name=T)
  
  peaks<-fread(peaks[grep(i,peaks)])
  
  peaks<-GRanges(seqnames = peaks$V1,
                 ranges = IRanges(start =peaks$V2,end =peaks$V3),
                 strand = "*")
  # feature count generate matrix 
  
  fragment.list<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments",
                            pattern="sorted.bed.gz.tbi",full.names = T) %>% sub(".tbi","",.)
  
  ab.type<-fragment.list[grep(i,fragment.list)]
  
  fragments <- CreateFragmentObject(ab.type,cells = kept.cells)
  
  counts<-FeatureMatrix(
    fragments=fragments,
    features=peaks,
    cells = kept.cells,
    process_n = 2000,
    sep = c("-", "-"),
    verbose = TRUE
  )
  
  signac[[paste0("Peaks.seacr_",i)]] <- CreateChromatinAssay(
    counts = counts[,rownames(signac@meta.data)],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 0,
    min.features = 0
  )
  # manual cal FRiP
  signac@meta.data[,paste0("FRiP_seacr_",i)] <- round(signac@meta.data[,paste0("nCount_Peaks.seacr_",i)] /signac@meta.data[,paste0("nCount_",i)] *10000)/100
  
  png(paste0(in.file.dir,"/",i,".seacr.peaks.FRiP.png"),res = 150, width = 800, height = 800)
  p1<-DensityScatter(signac, x = paste0("nCount_",i ), y = paste0('FRiP_seacr_',i), log_x = TRUE, quantiles = TRUE)
  print (p1)
  dev.off()
}

######################## de novo SEACR vs IgG stringent peaks ##################

for (i in c( "H3K27ac","H3K4me1","H3K36me3")){
  # i<- "IgG" i<-"H3K27me3" #IgG 
  DefaultAssay(signac)<- i
  #call peaks   load peaks from encode
  #
  peaks <- list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/SEACR", pattern=".vs.IgG.stringent.bed",full.name=T)
  peaks<-fread(peaks[grep(i,peaks)])
  peaks<-GRanges(seqnames = peaks$V1,
                 ranges = IRanges(start =peaks$V2,end =peaks$V3),
                 strand = "*")
  # feature count generate matrix 
  fragment.list<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments",
                            pattern="sorted.bed.gz.tbi",full.names = T) %>% sub(".tbi","",.)
  
  ab.type<-fragment.list[grep(i,fragment.list)]
  fragments <- CreateFragmentObject(ab.type,cells = kept.cells)
  counts<-FeatureMatrix(
    fragments=fragments,
    features=peaks,
    cells = kept.cells,
    process_n = 2000,
    sep = c("-", "-"),
    verbose = TRUE
  )
  
  signac[[paste0("Peaks.seacr_",i)]] <- CreateChromatinAssay(
    counts = counts[,rownames(signac@meta.data)],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 0,
    min.features = 0
  )
  
  # manual cal FRiP
  signac@meta.data[,paste0("FRiP_seacr_",i)] <- round(signac@meta.data[,paste0("nCount_Peaks.seacr_",i)] /signac@meta.data[,paste0("nCount_",i)] *10000)/100
  
  png(paste0(in.file.dir,"/",i,".seacr.peaks.FRiP.png"),res = 150, width = 800, height = 800)
  p1<-DensityScatter(signac, x = paste0("nCount_",i ), y = paste0('FRiP_seacr_',i), log_x = TRUE, quantiles = TRUE)
  print (p1)
  dev.off()
}

saveRDS(signac,paste0(in.file.dir,"/seurat.rna.hist_bin5k.rds"))


################################################################################
#################################### plot FRIP #################################
################################################################################

to.plot<-signac@meta.data %>% tibble::rownames_to_column("CB") %>% as.data.table()
to.plot<-melt(to.plot)
to.plot$variable %>% unique()

frips<-unique(to.plot$variable)
to.plot<-to.plot[variable %in% frips[grep("FRiP",frips)] ]

to.plot[, type:=variable %>% str_split(.,"[.]|_",simplify = T) %>% .[,2]]

to.plot[,unique (type)]
#to.plot[!type %in% c("ChIP","broad","CT","seacr"), type:="de.novo.narrow"]
#to.plot[type %in% c("broad"), type:="de.novo.broad"]

to.plot[type %in% c("ChIP","CT","seacr"), hist:=variable %>% str_split(.,"[.]|_",simplify = T) %>% .[,2]]
to.plot<-to.plot[!is.na(hist)]

to.plot[,unique (hist)]
to.plot[hist %in% c("ChIP","broad","CT","seacr"), hist:= variable %>% str_split(.,"[.]|_",simplify = T) %>% .[,3]]
to.plot[,unique (hist)]
to.plot$hist<-factor(to.plot$hist,levels = c(   "H3K27me3",     "H3K27ac" ,"H3K4me1",  "H3K4me3", "H3K36me3" ,  "IgG"    ))
to.plot$type<-factor(to.plot$type,levels = c(   "de.novo.narrow",     "de.novo.broad","seacr" ,"CT",  "ChIP"    ))

to.label<-to.plot[,.(mean=mean(value),median=median(value)),.(variable,type,hist)]
pdf("/Users/yang/data/scMTR_method.dev/figure3.data/frip.box.pdf",width=6,height = 6)
pdf("/Users/yang/data/scMTR_method.dev/figure3.data/frip.box1.pdf",width=6,height = 6)
pdf("/Users/yang/data/scMTR_method.dev/figure3.data/frip.box2.pdf",width=6,height = 6)

ggplot(to.plot,aes(hist,value))+
  geom_jitter(aes(color=hist) ,alpha=0.5,size=0.2)+
  #geom_boxplot(color="black",alpha=0.5,fill=NA,outlier.colour = NA,width=0.2)+
  #geom_violin(fill=NA)+
  #geom_text(data = to.label, aes(label = round(mean*100)/100, y = max(to.plot$value) *0.9),size=4)+
  scale_color_manual(values=c(H3K27me3 =H3K27me3.color,H3K27ac=H3K27ac.color, H3K4me1= H3K4me1.color ,H3K4me3=H3K4me3.color ,H3K36me3=H3K36me3.color, IgG =IgG.color)) + 
  facet_grid(type~.)+
  labs(x="",y="Fraction reads in peaks (%)")+
  theme_bw() +
  theme(
    #    strip.background = element_blank(),
    strip.text = element_text(size=16),
    axis.text.x = element_text(size=16, angle = 45, vjust = 1,hjust=1),
    axis.text.y = element_text(size=16),
    axis.title = element_text(size=16),
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  ) 
dev.off()


