## generate signac files 
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

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
#genome(annotation) <- "hg38"
seqlengths(annotation)

###############################################
############# load RNA counts  ################
###############################################

in.file.dir <- "/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments"
file.pattern <- "sorted.bed.gz.tbi"
fragment.list<-list.files(in.file.dir,pattern=file.pattern,full.names = T) %>% sub(".tbi","",.)
read.summary<-fread("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/rna.dna.kept.summary.tsv")
kept.cells<-read.summary[sample.type=="Def.endo.D0", unique(new.cellid) ]

scRNA<-readRDS("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/RNA/scRNA.data.seurat.rds")
scRNA<-subset(scRNA,subset= sample_bio.id=="Def.endo.D0")
scRNA<-subset(scRNA,subset= stage=="D0")
Assays(scRNA)
scRNA$CB<-paste0(substring(colnames(scRNA),4,12),substring(colnames(scRNA),1,2)) 
rownames(scRNA@meta.data)<-scRNA$CB
RNAcounts<-scRNA@assays$RNA@counts
colnames(RNAcounts) <-paste0(substring(colnames(RNAcounts),4,12),substring(colnames(RNAcounts),1,2)) 
#7479
frag.path<-fragment.list[1]
kept.cells<-unique(fread(frag.path)%>% .$V4)
#7484
kept.cells<-kept.cells[kept.cells %in%colnames(RNAcounts)]

# creates a Seurat object based on the scRNA-seq data
signac <- CreateSeuratObject(counts = RNAcounts[,kept.cells])

library(future)
#plan(multicore, workers = 2)
options(future.globals.maxSize = 20000 * 1024^2) # for 20 Gb RAM


###############################################
############# make tile file ##################
###############################################

for (frag.path in fragment.list) {
#frag.path<-fragment.list[1]
  fragments <- CreateFragmentObject(frag.path,cells = kept.cells)
  counts <- AggregateTiles(
    fragments,
    genome= seqlengths(annotation),
    cells = kept.cells,
    min_counts = 1,
    binsize = 5000,
    verbose = TRUE
  )
  saveRDS(counts,paste0(frag.path,"_bin5000.rds"))
}

for (frag.path in fragment.list) {
  #frag.path<-fragment.list[1]
  fragments <- CreateFragmentObject(frag.path,cells = kept.cells)
  counts <- AggregateTiles(
    fragments,
    genome= seqlengths(annotation),
    cells = kept.cells,
    min_counts = 1,
    binsize = 10000,
    verbose = TRUE
  )
  saveRDS(counts,paste0(frag.path,"_bin10000.rds"))
}

for (frag.path in fragment.list) {
  #frag.path<-fragment.list[1]
  fragments <- CreateFragmentObject(frag.path,cells = kept.cells)
  counts <- AggregateTiles(
    fragments,
    genome= seqlengths(annotation),
    cells = kept.cells,
    min_counts = 1,
    binsize = 20000,
    verbose = TRUE
  )
  saveRDS(counts,paste0(frag.path,"_bin20000.rds"))
}


########################################################
############# generate signac objects ##################
########################################################

ab.type<-fragment.list[grep("H3K27ac",fragment.list)]
# counts <- readRDS(paste0("~/data/Sample_5875_SLX-22471_200223YW_CTR.HE_MTR.E75_batch1_DNA/paired/fragments/",ab.type,"_bin5000.rds"))
 counts <-readRDS(paste0(ab.type,"_bin5000.rds"))
  ##
  signac[["H3K27ac"]] <- CreateChromatinAssay(
    counts = counts[,kept.cells],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 1,
    min.features = 1
  )

  ab.type<-fragment.list[grep("H3K27me3",fragment.list)]
  counts <-readRDS(paste0(ab.type,"_bin5000.rds"))
  signac[["H3K27me3"]]<-CreateChromatinAssay(
    counts = counts[,kept.cells],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 1,
    min.features = 1)
  
 ##
   ab.type<-fragment.list[grep("H3K4me1",fragment.list)]
   counts <-readRDS(paste0(ab.type,"_bin5000.rds"))
   signac[["H3K4me1"]]<-CreateChromatinAssay(
    counts = counts[,kept.cells],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 1,
    min.features = 1) 
##
  ab.type<-fragment.list[grep("H3K4me3",fragment.list)]
  # counts <- readRDS(paste0("~/data/Sample_5875_SLX-22471_200223YW_CTR.HE_MTR.E75_batch1_DNA/paired/fragments/",ab.type,"_bin5000.rds"))
  counts <-readRDS(paste0(ab.type,"_bin5000.rds"))
  signac[["H3K4me3"]]<-CreateChromatinAssay(
    counts = counts[,kept.cells],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 1,
    min.features = 1)
  ##
  ab.type<-fragment.list[grep("H3K36me3",fragment.list)]
  
  # counts <- readRDS(paste0("~/data/Sample_5875_SLX-22471_200223YW_CTR.HE_MTR.E75_batch1_DNA/paired/fragments/",ab.type,"_bin5000.rds"))
  counts <-readRDS(paste0(ab.type,"_bin5000.rds"))
  
  signac[["H3K36me3"]]<-CreateChromatinAssay(
    counts = counts[,kept.cells],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 1,
    min.features = 1)
  ##
  ab.type<-fragment.list[grep("IgG",fragment.list)]
  # counts <- readRDS(paste0("~/data/Sample_5875_SLX-22471_200223YW_CTR.HE_MTR.E75_batch1_DNA/paired/fragments/",ab.type,"_bin5000.rds"))
  counts <-readRDS(paste0(ab.type,"_bin5000.rds"))
  signac[["IgG"]]<-CreateChromatinAssay(
    counts = counts[,kept.cells],
    annotation = annotation,
    sep = c("-", "-"),
    genome = 'hg38',
    fragments = ab.type,
    min.cells = 1,
    min.features = 1) 
  
 Assays(signac) 
 
   # saveRDS(object = signac, file = sprintf('%s/%s_signac_%d.rds',args$outputdir,args$histone,args$bin_size))
  saveRDS(signac,paste0(in.file.dir,"/seurat.rna.hist_bin5k.rds"))

  



