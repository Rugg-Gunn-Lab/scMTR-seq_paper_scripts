library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(purrr)
library(data.table)
library(stringr)

##############################################################
################### bivalent gene dynamic ####################
##############################################################

# 1. use gene and define the promoter regions, possible 1kb +- TSS
# 2. overlap chromatin state with gene-promoters 
# 3. for each time point, calculate the ratio of chromatin states in each promoter, gene1 %polycomb? %bivalent? %active promoter? %umarked? based on length

# step1
gene.anno<-fread("/Users/wangy/data/ref/genome_ref/human/refdata-gex-GRCh38-2020-A/genes/genes.gtf.gz") 
#gene.anno<-gene.anno[V2=="ENSEMBL"]
gene.anno[,.N,V3]

#gene.anno<-gene.anno[V3=="gene",] 
gene.anno[,gene:=str_split(V9,";",simplify = T)[,6]] %>% .[,gene:=substring(gene, 13,nchar(gene)-1)]
gene.anno[,type:=str_split(V9,";",simplify = T)[,5]]%>% .[,type:=substring(type, 13,nchar(type)-1)]
gene.anno[,transcript.type:=str_split(V9,";",simplify = T)[,7]]%>% .[,transcript.type:=substring(transcript.type, 19,nchar(transcript.type)-1)]
gene.anno[gene=="TBXT" &V5>166168655]
gene.anno[V3=="transcript"&gene=="WNT3"]
gene.anno[,len:=V5-V4]

#max length transcript
protein.coding<-gene.anno[V3=="transcript"&transcript.type=="protein_coding",.SD[which.max(len)],gene]
protein.coding[gene=="TBXT"]


protein.coding<-protein.coding[,c(2,5,6,8,11,12,13,1)] %>% setnames(c("chr","start","end","strand","type","len","transcript.type","gene"))
fwrite(protein.coding,"/Users/wangy/data/ref/genome_ref/human/refdata-gex-GRCh38-2020-A/genes/geneanno.gene.protein.coding.tsv",sep="\t")
protein.coding<-protein.coding[type=="protein_coding"]
protein.coding[type=="protein_coding" &strand=="+",tss:=start]
protein.coding[type=="protein_coding" &strand=="-",tss:=end]
protein.coding[,start:=tss-1000]
protein.coding[,end:=tss+1000]
protein.coding<-protein.coding[chr%in% paste0("chr",c(1:22,"X","Y"))]

#fwrite(protein.coding,"/Users/wangy/data/ref/genome_ref/human/refdata-gex-GRCh38-2020-A/genes/geneanno.tss.updown2kb.tsv",sep="\t")
fwrite(protein.coding,"/Users/wangy/data/ref/genome_ref/human/refdata-gex-GRCh38-2020-A/genes/geneanno.tss.updown1kb.tsv",sep="\t")
protein.coding[gene=="GATA4"]

#step2
inputdir<-"/Users/wangy/data/scMTR_method.dev/figure4.data/chromhmm_adj/states_15"
files<-list.files(inputdir,pattern = "segments.bed",full.names = T)
foo<-lapply(files,function(i) fread(i) %>% setnames(c("chr","start","end","states")) %>% .[,peaks:=paste0("peaks_",.I)] %>% .[,anno:=basename(i) %>% sub("_15_segments.bed","",.)] ) %>% rbindlist()
setkey(foo,chr,start,end)
setkey(protein.coding,chr,start,end)
protein.coding[is.na(start)]
ov<-foverlaps(foo,protein.coding,nomatch = 0L)

ov[,intersect.len:=min(i.end,end)-max(i.start,start),.(gene,states,peaks,anno,chr,start,end,strand)]
ov[,.N,intersect.len]
ggplot(ov,aes(intersect.len))+
  geom_density()

ov[states %in% c("E1","E3"), state.type:="Transcription"]
ov[states %in% c( "E2","E4","E8","E11","E15"), state.type:="Unmarked"]
ov[states %in% c("E5"), state.type:="Primed enhancer"]
ov[states %in% c("E6","E7"), state.type:="Active enhancer"]
ov[states %in% c("E9","E10"), state.type:="Active promoter"]
ov[states %in% c("E12","E13"), state.type:="Bivalent promoter"]
ov[states %in% c("E14"), state.type:="Repressed Polycomb"] 
ov[gene=="GATA4"]

to.plot<-ov[,sum(intersect.len),.(state.type,gene,anno)]
to.plot[,ratio:=V1/2000,.(gene,anno)]
dcast(to.plot,anno+gene~state.type,value.var = "ratio")

to.plot[,.N,state.type]
to.plot<-to.plot[!state.type=="Unmarked",.SD[which.max(ratio)],.(gene,anno)]

kept.gene<-to.plot[state.type%in%c("Bivalent promoter") &anno=="D0",gene]%>% unique()
to.plot<-to.plot[gene%in% kept.gene]

to.plot[,.N,state.type]

to.plot[state.type=="Bivalent promoter",state.quant:=1]
to.plot[state.type%in% c("Active promoter","Transcription","Active enhancer","Primed enhancer",""),state.quant:=2]
to.plot[state.type=="Unmarked",state.quant:=0]
to.plot[state.type=="Repressed Polycomb",state.quant:=-1]

to.plot[gene=="SOX17"]
ov[gene=="SOX17"&anno=="D3"]
to.plot.heat<-to.plot %>% dcast(.,gene~anno,value.var="state.quant")%>% as.data.frame() %>% tibble::column_to_rownames("gene")
to.plot.heat[is.na(to.plot.heat)]<-"0"

to.plot.heat$D0<-as.numeric(to.plot.heat$D0)
to.plot.heat$D1<-as.numeric(to.plot.heat$D1)
to.plot.heat$D2<-as.numeric(to.plot.heat$D2)
to.plot.heat$D3<-as.numeric(to.plot.heat$D3)
setorder(to.plot.heat,D3,D2,D1)
pdf("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/chromatin.states.bivalent.dynamic.heatmap.pdf",width = 6,height = 25)
pheatmap::pheatmap(to.plot.heat,cluster_cols = F,cluster_rows = F,show_rownames = T,fontsize = 1)
dev.off()

write.table(to.plot.heat,"/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/plot/chromatin.states.bivalent.dynamic.tsv",sep="\t")
