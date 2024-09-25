library(ggplot2)
library(data.table)
library(purrr)
library(viridis)
library(stringr)


################################################################################
################ metaplot lineage spcific gene promoter activity ###############
################################################################################
#to generate TSS regions for the lineage specific genes 
# gene - Tss 
genes<-fread("~/data/ref/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt")
genes<-genes[,.(chr,start,end,strand,symbol,ens_id)]
genes[strand=="+",c("start","end"):=list(start,start+1)]
genes[strand=="-",c("start","end"):=list(end-1,end)]
tss<-genes

## load gene specific markers
markers<-fread("~/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments/plot/RNA_markers.tsv")
markers[p_val_adj<0.01 &avg_log2FC>log2(4),]
markers[p_val_adj<0.01 &avg_log2FC>log2(2),.N,cluster]

ov<-merge(genes[,.(chr,start,end,gene=symbol)],
          markers[p_val_adj<0.01 &avg_log2FC>log2(4),], by ="gene"
)
ov[cluster==0,celltype:="TE"]
ov[cluster==1,celltype:="PE"]
ov[cluster==2,celltype:="EPI"]
remove.dup.genes<-ov[,.N,(gene)] %>% .[N>1,gene]

out<-ov[!gene%in% remove.dup.genes,.(chr,start,end,gene,celltype,avg_log2FC,p_val_adj)] %>% setkey(chr,start,end) %>% split(.,by="celltype",drop=F)
walk2(out , paste0("~/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments/plot/TSS.of.lineage.specific.markers.fold2.", names(out),".tsv"), fwrite,sep="\t",col.names=F)


################################################################################
######################### metaplot enhancers activity ##########################
################################################################################

# 1. to define lineage specific enhancers, use gene specificity to classify the linked enhancers
# 2. keep clear lineage spficific enhancers, collaps the regions overlapped.
# 3. use deeptools to plot metaplot the three lineage specific enhancer in three lineages in different histones

regulon<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eRegulon.csv")
enhancer.gene.links<-regulon[,.(Region, Gene)] %>% unique() %>% 
  .[,chr:=str_split(Region,":",simplify = T) %>% .[,1]] %>% 
  .[,start:=str_split(Region,":|-",simplify = T) %>% .[,2]]%>% 
  .[,end:=str_split(Region,":|-",simplify = T) %>% .[,3]] %>% .[,.(chr,start,end,Region,Gene)] %>% setkey(chr,start,end)
markers<-fread("/Users/wangy/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments/plot/RNA_markers.tsv")
markers[avg_log2FC>=2,.N,cluster]
markers[cluster==1]
markers[cluster==2, gene]
markers[gene %in% enhancer.gene.links$Gene, .N,cluster] 
#       0   288
#       1    91
#       2   123
markers[gene %in% enhancer.gene.links$Gene &cluster==2,gene] #EPI
markers[gene %in% enhancer.gene.links$Gene &cluster==1,gene] #PE


ov<-merge(enhancer.gene.links,markers[,.(cluster,Gene=gene)], by = "Gene" )%>% .[,.(chr,start,end,Region,Gene,cluster)] %>% setkey(chr,start,end)
ov[,.N,cluster]
#    0   464
#      2   151
#      1   123
ov[cluster==0,type:="TE"]
ov[cluster==1,type:="PE"]
ov[cluster==2,type:="EPI"]
ov<-ov%>% split(by="type",drop = T)
walk2(ov, paste0("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/",names(ov),"_lineage.specific.enhancer.basedon.gene.specificity.tsv"),fwrite, sep="\t", col.names=F )
ov[["TE"]]%>%.[,.N,chr]

