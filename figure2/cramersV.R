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
############### calculate Cramer's V between different modalities ##############
################################################################################

#############################  Cramer's V function #############################

cramer.V<-function(x,y){
  temp<-table(x); 
  temp[as.character(c(0,1,2,3)[!c(0,1,2,3) %in% labels(temp)[[1]]])]<-0; 
  temp<-temp[sort(labels(temp))]; 
  sqrt((as.numeric(chisq.test(matrix(temp, nrow=2))$statistic))/nrow(y))
}

###############Calculate Cramer's V of different associations ##################

pairwise <- combn(c("H3K27me3_genes","H3K27ac_genes","H3K4me3_genes","H3K4me1_genes","H3K36me3_genes","RNA"), 2) %>%
  as.data.frame()

to.plot.cramer<-list()
#for (i in paste0("V",1:15)){
for (i in paste0("V",c(5,9,12,14,15))){
  ### Generate gene score matrices ###
  DefaultAssay(signac)<-pairwise[1,i]
  hist1<-as.data.frame(GetAssayData(signac, slot="data"))
  DefaultAssay(signac)<-pairwise[2,i]
  hist2<-as.data.frame(GetAssayData(signac, slot="data"))
  hist1<-hist1[rownames(hist1) %in% rownames(hist2), ]
  dim(hist1)
  hist2<-hist2[rownames(hist2) %in% rownames(hist1), ]
  dim(hist2)
  # rowname and colname should match
  hist1[hist1 > 0]<-1
  hist2[hist2 > 0]<-2
  hist1_2<-hist1+hist2[rownames(hist1),colnames(hist1)]
  cramer.hist1_2<-apply(hist1_2, 2, function(x) cramer.V(x, hist1_2))
  celltypes<-signac$celltype[colnames(hist1_2)]
  cramer.hist1_2<-data.frame(celltypes=celltypes, cramer=cramer.hist1_2)
  rownames(cramer.hist1_2)<-colnames(hist1_2)
  cramer.hist1_2$cramer[cramer.hist1_2$cramer > 1]<-0
  cramer.hist1_2$test<-paste0(pairwise[1,i],"-",pairwise[2,i])
  cramer.hist1_2$cellid<-rownames(cramer.hist1_2)
  to.plot.cramer[[i]]<-cramer.hist1_2
}

to.plot<-rbindlist(to.plot.cramer)
fwrite(to.plot,paste0(in.file.dir,"/cramers.v.tsv"),sep="\t")

############################## Plot Cramer's V results #########################

to.plot[,test:= sub("_genes","",test) %>% sub("_genes","",.)]

to.plot[,median(cramer),test] %>% setkey(V1) %>% .$test -> order
to.plot$test<-factor(to.plot$test,levels = order)

pdf(paste0(in.file.dir,"/cramers.v.boxplot.pdf"),height = 6,width=8)
ggplot(to.plot, aes(x=test, y=cramer, color=test)) + 
  geom_boxplot(width=0.1, position=position_dodge(width=0.7), lwd=0.75,outlier.colour = NA,color="grey") + 
  geom_violin(scale="width", width=0.7, position=position_dodge(width=0.7), lwd=0.75,fill=NA) + 
  theme_bw() + 
  labs(x="",y="Cramer's V")+
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text.x=element_text(size=12,angle = 45,hjust=1,vjust=1), 
        axis.title=element_text(size=14), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=14)) 
dev.off()



