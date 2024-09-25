library(ggplot2)
library(data.table)
library(purrr)
library(viridis)
library(stringr)

################################################################################
####################### plot eRegulon activity in heatmap ######################
################################################################################

############################# gene based activity ##############################

## heatmap plot
gene.based<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eregulon.gene.based.heatmap.correlation.jaccard.csv") %>% as.data.frame() %>% tibble::column_to_rownames("V1")
#gene.based[gene.based>0.5]<-0.5
select<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/selected.eRegulons.gene.based.csv")
select[abs(Rho)>=0.25,Cistrome]
gene.based<-gene.based[select[abs(Rho)>=0.25,Cistrome], select[abs(Rho)>=0.25,Cistrome]]
to.plot<-pheatmap::pheatmap(gene.based,cluster_rows = T,cluster_cols = F,color = colorRampPalette(c("darkblue", "purple", "orangered","yellow"))(100),border_color = NA ,na_col = "white",)
rownames(gene.based)[to.plot$tree_row$order]

pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eregulon.gene.based.heatmap.clustered.pdf",height=10,width = 6)
pheatmap::pheatmap(gene.based[,rownames(gene.based)[to.plot$tree_row$order]],cluster_rows = T,cluster_cols = F,
                   color = colorRampPalette(c(plasma(20)))(100),
                   border_color = NA ,na_col = "white",)
dev.off()

############################ region based activity #############################

region.based<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eregulon.region.based.heatmap.correlation.jaccard.csv") %>% as.data.frame() %>% tibble::column_to_rownames("V1")
#gene.based[gene.based>0.5]<-0.5
select<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/selected.eRegulons.gene.based.csv")
dict<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eRegulon.csv")%>% .[,.(Region_signature_name, Gene_signature_name,TF )] %>% unique()
kept<-dict[ Gene_signature_name %in% select[abs(Rho)>=0.25,Cistrome],Region_signature_name]
region.based<-region.based[kept,kept]
to.plot<-pheatmap::pheatmap(region.based,cluster_rows = T,cluster_cols = F,color = colorRampPalette(c("darkblue", "purple", "orangered","yellow"))(100),border_color = NA ,na_col = "white",)
rownames(region.based)[to.plot$tree_row$order]

pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eregulon.region.based.heatmap.clustered.pdf",height=10,width = 6)
pheatmap::pheatmap(region.based[,rownames(region.based)[to.plot$tree_row$order]],cluster_rows = T,cluster_cols = F,
                   color = colorRampPalette(c(plasma(20)))(100),
                   border_color = NA ,na_col = "white",)
dev.off()

gene.based[gene.based>0.5]<-0.5
pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eregulon.region.based.heatmap.cut.0.5.clustered.pdf",height=10,width = 6)
pheatmap::pheatmap(gene.based[rownames(gene.based)[to.plot$tree_row$order],rownames(gene.based)[to.plot$tree_row$order]],cluster_rows = F,cluster_cols = F,
                   color = colorRampPalette(c(plasma(20)))(100),
                   border_color = NA ,na_col = "white",)
dev.off()


################################################################################
########################### Bulb plot eRegulon  ###############################
################################################################################

#################### calculate lineage level TF expression #####################

# heatmap plot bulk level tf exp and auc region or auc gene
dict<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eRegulon.csv")%>% .[,.(Region_signature_name, Gene_signature_name,TF )] %>% unique()
select<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/selected.eRegulons.gene.based.csv")
kept<-dict[ Gene_signature_name %in% select[abs(Rho)>=0.25,Cistrome],TF]
#gene.based<-gene.based[kept,kept]

auc.gene<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/gene.based.auc.single.cell.tsv.gz")
auc.region<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/region.based.auc.single.cell.tsv.gz")
exp<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/exp.single.cell.tsv.gz")
select_cols<-c("V1", kept)
exp<-exp[,..select_cols]

auc.region<-melt(auc.region) %>% setnames(c("cell","Region_signature_name","region.AUC"))
auc.gene<-melt(auc.gene)%>% setnames(c("cell","Gene_signature_name","gene.AUC"))
exp<-melt(exp)%>% setnames(c("cell","TF","EXP"))
merged<-merge(auc.region,dict,by="Region_signature_name") %>% merge(auc.gene, by=c("Gene_signature_name","cell")) %>% merge(exp,by=c("TF","cell")) %>% .[,celltype:=str_split(cell,"___",simplify =T) %>% .[,2] ]

##save merged sc
fwrite(merged,"/Users/wangy/Desktop/scenicplus.results/plot.e4.5/single.cell.Region.AUC.Gene.AUC.TF.exp.tsv",sep='\t')

# bulk celltype level)
bulk.merged<-merged[,.(region=mean(region.AUC),gene=mean(gene.AUC),exp=mean(EXP)),.(celltype,TF,Gene_signature_name,Region_signature_name)]

# scale 0-1
bulk.merged[,scaled.region:=(region-min(region))/(max(region)-min(region)),Region_signature_name]
bulk.merged[,scaled.gene:=(gene-min(gene))/(max(gene)-min(gene)),Gene_signature_name]
bulk.merged[,scaled.exp:=(exp-min(exp))/(max(exp)-min(exp)),TF]

################ plot eregulon TF expression and regulon activity ##############

test<-dcast(bulk.merged,Region_signature_name~celltype,value.var = "scaled.exp") %>% as.data.frame()%>% tibble::column_to_rownames("Region_signature_name")
ord.tf<-pheatmap::pheatmap(test,cluster_cols = F)
dev.off()

################ region based eregulon activity
bulk.merged$Region_signature_name<-factor(bulk.merged$Region_signature_name,levels = (ord.tf$tree_row$labels[ord.tf$tree_row$order] ))
bulk.merged$celltype <-factor(bulk.merged$celltype,levels = rev(c("EPI","PE","TE")))
pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/region.based.bulb.plot.pdf",width = 8,height=7)
ggplot(data = bulk.merged, mapping = aes(x=Region_signature_name,y=celltype))+
  geom_point(aes(size= scaled.region,color=  scaled.exp))+ 
  scale_color_viridis()+
  scale_size_continuous(breaks = c(-1,0,1),range = c(0.1,8))+
  labs(y="",x="",title = "Region based regulon",color="Scaled TF expression",size="Scaled AUC score")+
  theme_minimal()+
  theme(
    axis.text.y = element_text(size=10),
    #    axis.title.x = element_text(size=16),
    axis.title = element_text(size=14),
    panel.grid = element_blank(),
    #axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=14),
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=12),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10),
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(size=10)
  )+ coord_fixed(ratio = 1.5)## coord_flip()
dev.off()

################ gene based eregulon activity
pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/gene.based.bulb.plot.pdf",width = 8,height=7)
ggplot(data = bulk.merged, mapping = aes(x=Region_signature_name,y=celltype))+
  geom_point(aes(size=scaled.gene ,color=scaled.exp))+ 
  scale_color_viridis()+
  # scale_size_continuous(breaks = c(0,10,25,50,100),range = c(0.2,4))+
  #  facet_grid(type.ndr~.)+
  labs(y="",x="",title = "Gene based regulon",color="Scaled TF expression",size="Scaled AUC score")+
  theme_minimal()+
  theme(
    axis.text.y = element_text(size=10),
    #    axis.title.x = element_text(size=16),
    axis.title = element_text(size=14),
    panel.grid = element_blank(),
    #axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=14),
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=12),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10),
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(size=10)
  )+ coord_fixed(ratio = 1.5)## coord_flip()
dev.off()

################################################################################
########## scatter plot gene-based eregulon and region-based eregulon ##########
################################################################################

dict<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eRegulon.csv")%>% .[,.(Region_signature_name, Gene_signature_name,TF )] %>% unique()
select.gene<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/selected.eRegulons.gene.based.csv")
select.region<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/selected.eRegulons.region.based.csv")
select.gene[,"Gene_signature_name":=Cistrome]
select.region[,"Region_signature_name":=Cistrome]
to.plot <-merge(select.gene,dict[,.(Region_signature_name, Gene_signature_name)],by="Gene_signature_name")  %>% merge(select.region, by="Region_signature_name") #%>% merge(cor.out, by=c("Gene_signature_name",   "Region_signature_name" ))
to.plot[,labels:=str_split(Region_signature_name, "[(]",simplify = T) %>% .[,1]]
to.plot[,labels:=substring(labels,1,nchar(labels)-1)]
to.plot[,labels:=labels %>% sub("_extended","",.)]

pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/correlations.pdf",width = 10, height=10)
ggplot(to.plot,aes(Rho.x,Rho.y,label=labels,color=(abs(Rho.x)<0.25)))+
  geom_point()+
  #scale_color_viridis(direction = -1)+ #
  labs(x="correlation coefficient of \n TF expression with target Gene AUC", 
       y="correlation coefficient of \n TF expression with target Region AUC") +
  geom_hline( yintercept =0)+
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = 0.25)+
  geom_vline(xintercept = -0.25)+
  geom_text_repel(data=to.plot[abs(Rho.x)>0.25,] ,max.overlaps =100,min.segment.length = 0 )+
  theme_minimal()+
  theme(
    panel.grid=element_blank()
  )
dev.off()


################################################################################
##################### plot Trps1 negative regulated genes ######################
################################################################################

exp<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/exp.single.cell.tsv.gz")
regulon<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eRegulon.csv")
kept<-regulon[Region_signature_name =="Trps1_extended_-_(37r)",unique(Gene)]
select_cols<-c("V1", kept,"Trps1","Gata2","Gata3")
exp<-exp[,..select_cols]%>% .[,celltype:=str_split(V1,"___",simplify =T) %>% .[,2] ]

trps1.neg.targets<-regulon[Region_signature_name =="Trps1_extended_-_(37r)",unique(Gene)]

gata3.targets<-regulon[TF=="Gata3",unique (Gene)]
gata3.targets[gata3.targets%in% trps1.neg.targets]

gata2.targets<-regulon[TF=="Gata2",unique (Gene)]
gata2.targets[gata2.targets%in% trps1.neg.targets]

gata2.3<-c(gata2.targets[gata2.targets%in% trps1.neg.targets],gata3.targets[gata3.targets%in% trps1.neg.targets]) %>% unique

to.plot<-melt(exp)

to.plot<-exp%>% setkey(celltype)%>% as.data.frame() %>% tibble::column_to_rownames("V1") 
to.plot[1:5,1:35]

################################## heat maps ###################################

rows<-data.frame(celltype =factor(to.plot[,36]))
rownames(rows)<-rownames(to.plot)
pheatmap::pheatmap(to.plot[,1:35],cluster_rows = F,#scale = "column",
                   color = colorRampPalette(c(plasma(20)))(100),show_rownames = F,
                   annotation_row = rows,
                   border_color = NA ,na_col = "white")
c1<-c('Gata3','Dstn','Efnb2','Pdzk1','Khnyn','Tspan8','Snx2','Prkcg')
c2<-c('Clic4',  'Klf5',  'Ebf4',  'Stard10',  'Spsb4',  'Enpep',  'Wipf1',  'Lrrfip1',  'Frs2',  'Wnt7b',  'Mbnl2')
c3<-c(  'Elf1',  'Rin3',  'Satb1',  'Slc22a23',  'Tmem144',  'Parva',  'P2rx3',  'Slc7a6',  'Mob3b',  'Sh3tc2')
c4<-c("Trps1","Gata2","Gata3")

pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/trps1_negative.targets.exp.heatmap.pdf",width = 6,height=6)
pheatmap::pheatmap(t(to.plot[,c(c1,c2,c3)]),cluster_rows = T,cluster_cols = F,#scale = "column",
                   color = colorRampPalette(c(plasma(20)))(100),show_colnames = F,
                   annotation_col = rows,
                   border_color = NA ,na_col = "white")
dev.off()

pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/trps1_negative.targets.gata2.3.exp.heatmap.pdf",width = 6,height=6)
pheatmap::pheatmap(t(to.plot[,gata2.3]),cluster_rows = T,cluster_cols = F,#scale = "column",
                   color = colorRampPalette(c(plasma(20)))(100),show_colnames = F,
                   annotation_col = rows,
                   border_color = NA ,na_col = "white")
dev.off()

pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/trps1.gata2.3.exp.heatmap.pdf",width = 6,height=6)
pheatmap::pheatmap(t(to.plot[,c("Trps1","Gata2","Gata3")]),cluster_rows = T,cluster_cols = F,#scale = "column",
                   color = colorRampPalette(c(plasma(20)))(100),show_colnames = F,
                   annotation_col = rows,
                   border_color = NA ,na_col = "white")
dev.off()

pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/trps1.exp.heatmap.pdf",width = 4,height=6)
pheatmap::pheatmap(to.plot[,c(c4)],cluster_rows = F,cluster_cols = T,#scale = "column",
                   color = colorRampPalette(c(plasma(20)))(100),show_rownames = F,
                   annotation_row = rows,
                   border_color = NA ,na_col = "white")
dev.off()

regulon[Region_signature_name=="Trps1_extended_-_(37r)",Region]

########################### overlap regions barplot ############################

dict<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/eRegulon.csv")%>% .[,.(Region_signature_name, Gene_signature_name,TF )] %>% unique()
select<-fread("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/selected.eRegulons.gene.based.csv")
kept<-dict[ Gene_signature_name %in% select[abs(Rho)>=0.25,Cistrome]]

to.plot<-regulon[Region_signature_name %in% kept$Region_signature_name & Region %in% regulon[Region_signature_name=="Trps1_extended_-_(37r)",Region],.( .N, TF,Gene_signature_name), Region_signature_name]%>% unique() %>% setorder(-N)
to.plot$Region_signature_name <- factor(to.plot$Region_signature_name,levels =c(to.plot$Region_signature_name))

pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/trps1.overlap.regions.num.barplot.pdf",width = 8,height=4)
ggplot(to.plot,aes(Region_signature_name,N))+
  geom_col()+
  theme_bw()+
  labs(x="",y="Number of overlapped regions")+
  theme(
    axis.text.x = element_text(angle=90, hjust=1,vjust=0.5)
  )
  dev.off()

############################ overlap genes barplot ############################
  
  to.plot<-regulon[Region_signature_name %in% to.plot$Region_signature_name & Gene %in% regulon[Region_signature_name=="Trps1_extended_-_(37r)",Gene], .(length(unique(Gene)),TF,Region_signature_name  ), Gene_signature_name]%>% unique()%>% setorder(-V1)
  to.plot$Gene_signature_name <- factor(to.plot$Gene_signature_name,levels =c(to.plot$Gene_signature_name))
  
#  pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/trps1.overlap.genes.num.barplot.pdf",width = 8,height=4)
  pdf("/Users/wangy/Desktop/scenicplus.results/plot.e4.5/trps1.overlap.genes.num.barplot.sameto.region.pdf",width = 8,height=4)
  ggplot(to.plot,aes(Gene_signature_name,V1))+
    geom_col()+theme_bw()+labs(x="",y="Number of overlapped genes")+
    theme(
      axis.text.x =  element_text(angle=90, hjust=1,vjust=0.5)
    )
  dev.off()
  
  
