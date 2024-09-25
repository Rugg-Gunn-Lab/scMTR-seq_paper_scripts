
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


##############################################################
###################### kept cell summary ####################
##############################################################

in.file.dir <- "/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments"
signac<-readRDS(paste0(in.file.dir,"/seurat.rna.hist_bin5k.def.endo.rds"))


to.plot<-signac@meta.data %>% tibble::rownames_to_column("CB") %>% as.data.table()
to.plot<-to.plot[,c(1:4,6,10,12,14,16,18,20)] %>% melt()


to.plot.all<-to.plot[!variable %in% c ("nCount_RNA", "nFeature_RNA"),.(value=sum(value)),.(CB,sample_bio.id,orig.ident)]

to.plot.all[,variable:="DNA all"]
to.plot<-rbind(to.plot,to.plot.all)
to.plot[!variable %in% c ("nCount_RNA", "nFeature_RNA"),lib:="DNA"]
to.plot[is.na(lib),lib:="RNA"]

to.plot$variable<-factor(to.plot$variable,levels = c("nCount_H3K27me3", "nCount_H3K27ac", "nCount_H3K4me1", "nCount_H3K4me3", "nCount_H3K36me3", "nCount_IgG", "DNA all", "nCount_RNA","nFeature_RNA" ))                
to.plot$lib<-factor(to.plot$lib,levels = c("DNA", "RNA"))                                                     
to.plot$sample_bio.id<-factor(to.plot$sample_bio.id,levels = c("Def.endo.D0","Def.endo.D1" ,"Def.endo.D2", "Def.endo.D3" ))

to.label<-to.plot[,.(mean=mean(value),median=median(value)),.(sample_bio.id,variable,lib)]


#to.plot$lib.type<-factor(to.plot$lib.type,levels = c("RNA","DNA"))
setkey(to.plot,variable)

pdf("/Users/yang/data/scMTR_method.dev/figure4.data/kept.reads.unique.box.pdf",width=10,height = 6)

ggplot(to.plot,aes(variable,log10(value),fill=sample_bio.id))+
  geom_violin()+
  #  geom_jitter(aes(color=sample_bio.id) ,alpha=0.5,size=0.2 )+
  geom_boxplot(color="black",outlier.colour = NA,width=0.2,position=position_dodge(0.9) )+
  geom_text(data = to.label, aes(x=variable, label = round(mean), y =6 *0.9),size=2,position = position_dodge(0.9))+
  ylim(c(0,6))+
  # scale_color_manual(values=c( H3K27me3.color, H3K27ac.color,  H3K4me1.color , H3K4me3.color , H3K36me3.color, IgG.color, "#D39200",  "#DB72FB", "#FF61C3")) + 
  facet_grid(.~lib,scales = "free",space = "free",drop = T)+#,)+
  labs(x="",y="log10(Number of Unique reads/genes)")+
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


