##### overlap peaks between ENCODE ChIP-seq peaks, single target CUT&Tag peaks, multitargets CUT&Tag peaks in different conditions.

library(data.table)
library(purrr)
library(ggplot2)
library(viridis)
library(stringr)
library(bedtoolsr)
library(ggpointdensity)
source("~/scripts/utility.R")
library(wesanderson)

####################################
##### ENCODE ChIP-seq peaksets #####
####################################

chip.peaks<-list.files("/Users/yang/data/ref/encode.h9.ren.bin.tracks.peaks",pattern=".bed.gz",full.names = T,recursive = T)
foo.chip<-lapply(chip.peaks,function(i) fread(i) %>%
         setnames(c("chr","start","end","peakid", "score","strand","signalValue","pValue","qValue","peak")) %>% 
         .[,anno:=basename(i) %>% 
             str_split(.,'[.]',simplify = T)%>% .[,1]] %>% 
           .[,assay:="ChIP-seq"]
)%>% rbindlist() %>%
  setkey(chr,start,end) %>% .[,peaktype:="ref"]
####################################
############# blacklist ############
####################################
blacklist.h<-fread("/Users/yang/data/ref/Blacklist-master/lists/hg38-blacklist.v2.bed.gz")%>% setnames(c("chr","start","end","type")) %>% setkey(chr,start,end)

####################################
############# atac-seq  ############
####################################
foo.atac<-fread("/Users/yang/data/ref/Pastor_NCB_primed.naive.hESC.ATACseq/primed.H9.ATAC-seq.hg38.tsv.gz") %>%
  #foo.atac<-fread("/Users/yang/data/scMTR_method.dev/ref.peaks/ATAC/SEACR/b10.rm.Y.M.fragments.tsv_seacr_top0.02.norm.peaks.stringent.bed",select = c(1:3)) %>% 
  setnames (c("chr","start","end")) %>% 
  .[,peakid:=paste0("Peak_",.I)]%>% 
  .[,anno:="Primed"]%>%
  .[,assay:="ATAC-seq"]
fwrite(foo.atac, "/Users/yang/data/ref/Pastor_NCB_primed.naive.hESC.ATACseq/primed.H9.ATACseq.hg38.bed",sep="\t",col.names = F )

####################################
############# CUT&Tag  ############
####################################
ct.peaks<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5742_200522YW_BMT2/Lane_8208_200522YW_BMT2/bulk.paired.mapped/fragments/SEACR",
                     pattern=".bed.gz",full.names = T) #,recursive = T) 
#ct.peaks<-ct.peaks[grep("STB",ct.peaks)]

foo<-lapply(ct.peaks,function(f)fread(f,header = F,select = c(1:3))%>% setnames(c("chr","start","end"))%>%
              .[,peakid:=paste0("Peak_",.I)]%>% 
              .[,anno:=basename(f) %>% str_split(.,"_",simplify = T) %>% .[,2] ] %>% 
              .[,assay:=paste0("CUT&Tag_",basename(f) %>% str_split(.,"[.]",simplify = T) %>% .[,2] )] %>% 
                .[,peaktype:=basename(f) %>% str_split(.,"seacr_",simplify = T)%>% .[,2] %>% substring(.,1,7) ]) %>% 
  rbindlist() %>% 
setkey(chr,start,end)



foo.tmp<-rbind(foo,foo.chip[anno %in% c("H3K27me3","H3K27ac")], foo.atac,fill=T) 

foo.tmp<- foo.tmp %>% setkey(anno,assay,peaktype,chr,start,end)

foo.tmp[,anno:=paste0(anno,"_",assay,"_",peaktype)]
foo.tmp<-split(foo.tmp,by=c("anno"))


###########################################################
################# pairwise overlap peaks  #################
########################################################### 
  pairwise <- combn(seq_along(foo.tmp), 2) %>%
    as.data.frame()
  
  overlaps <- map(pairwise, ~{
    dt <- c(.[1], .[2]) %>%
      map(~foo.tmp[[.]]) %>%
      map(setkey, chr, start, end) %>%
      purrr::reduce(foverlaps, nomatch = 0L)
  })
  
  
  venn_numbers <- map2(overlaps[names(overlaps[lapply(overlaps,nrow)>0])], pairwise[names(overlaps[lapply(overlaps,nrow)>0])], ~{
    anno1 <- foo.tmp[[.y[1]]][, .(anno = unique(anno), .N)]
    anno2 <- foo.tmp[[.y[2]]][, .(anno = unique(anno), .N)]
    list(anno1 = anno1[, anno],
         anno2 = anno2[, anno],
         cross.area = .x[, min(length(unique(peakid)), length(unique(i.peakid)))],
         cross.1=.x[, length(unique(i.peakid))],
         cross.2=.x[, length(unique(peakid))],
         area1 = anno1[, N],
         area2 = anno2[, N],
         ratio.1 = .x[, length(unique(i.peakid))]/anno1[, N]*100,
         ratio.2 = .x[, length(unique(peakid))]/anno2[, N]*100
    )
  })
  
  test<-rbindlist(venn_numbers)
  row.order<-unique(rbindlist(foo.tmp)$anno)#[c(8,4,5,1,7,3,6,2)]
  
  tmp<-test[anno1%in%row.order &anno2%in% row.order]
  tmp$anno1<-factor(tmp$anno1,levels = row.order)
  tmp$anno2<-factor(tmp$anno2,levels = row.order)
  setkey(tmp,anno1,anno2)
  tmp[,ratio.1:=round(ratio.1)]
  tmp[,ratio.2:=round(ratio.2)]
  
  sel.anno1<- tmp$anno1[grep("0.02|ATAC|ChIP",tmp$anno1)] %>% unique()
  sel.anno2<-tmp$anno2[grep("0.02|ATAC|ChIP",tmp$anno2)] %>% unique()
  
  to.plot<-tmp[anno1%in% sel.anno1 & anno2 %in% sel.anno2]
  
  fwrite(tmp, "/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.peaks.ov.tsv",sep="\t")
 # fwrite(tmp, "/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.peaks.vs.atac.ov.tsv",sep="\t")
  #fwrite(tmp, "/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.peaks.vs.ref.atac.ov.tsv",sep="\t")
  pdf(paste0("/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.ov.1.pdf"),width = 10,height=10)
#  pdf(paste0("/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.vs.atac.ov.1.pdf"),width = 10,height=10)
 # pdf(paste0("/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.vs.ref.atac.ov.1.pdf"),width = 10,height=10)
#  pdf(paste0("/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.vs.ref.atac.ov.1.sel.top02.pdf"),width = 6,height=6)
  
  p1<-ggplot(data =to.plot , aes(anno1, anno2, fill =  ratio.2))+
    geom_tile(color = "white")+
    scale_fill_viridis(discrete = F, begin = 0.1, end = 0.9, option = "magma", limit = c(0,100), space = "Lab", alpha = 0.8,name="Percentage") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1))+
    theme(axis.text.y = element_text(angle = 0, vjust = 0.5, 
                                     size = 10, hjust = 1))+
    geom_text(aes(label = ratio.2), color = "grey", size = 5) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_line(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank()#,
    ) +
    coord_fixed()
  print(p1)
  dev.off()
  

  
  pdf(paste0("/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.ov.2.pdf"),width = 10,height=10)
#  pdf(paste0("/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.vs.atac.ov.2.pdf"),width = 10,height=10)
 # pdf(paste0("/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.vs.ref.atac.ov.2.pdf"),width = 10,height=10)
  #pdf(paste0("/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.vs.ref.atac.ov.2.sel.top02.pdf"),width = 6,height=6)
  p2<-ggplot(data = to.plot, aes(anno2, anno1, fill = ratio.1))+ #ratio.umr ratio.ndr [anno1 %in% stages &anno2%in% stages ]
    geom_tile(color = "white")+
    scale_fill_viridis(discrete = F, begin = 0.1, end = 0.9, option = "magma", limit = c(0,100), space = "Lab", alpha = 0.8,name="Percentage") +
    
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1))+
    theme(axis.text.y = element_text(angle = 0, vjust = 0.5, 
                                     size = 10, hjust = 1))+
    geom_text(aes(label =ratio.1), color = "grey", size = 5) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_line(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank()#,
    ) +
    coord_fixed()
  print(p2)
  dev.off()

    pdf("/Users/yang/data/scMTR_method.dev/plot/bulk.ct.vs.chip.peaks.pdf",width = 8,height=5)
   ggplot(rbindlist(foo.tmp) %>% .[,.N,.(anno,assay,peaktype)],aes(assay,N,fill=peaktype))+
     geom_col(position = "dodge")+
     facet_grid(.~ str_split(anno,"_",simplify=T)[,1])+
     labs(x="",y="# of peaks")+
     theme_bw()
   dev.off()
 
   



