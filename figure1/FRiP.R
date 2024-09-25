# calculate counts per peaks to assess the contanimations, consistencies, among histone marks 
library(data.table)
library(purrr)
library(ggplot2)
library(viridis)
library(stringr)
library(bedtoolsr)
#BiocManager::install("ggpointdensity")
library(ggpointdensity)
source("~/scripts/utility.R")
library(wesanderson)
######################################################################################################
###### loading ENCODE ChIP-seq peak sets, ATAC-seq peak sets, ST Bulk CUT&Tag, Black list ############
######################################################################################################
# peak overlap (reproducible), peak signal, peak 
chip.peaks<-list.files("/Users/yang/data/ref/encode.h9.ren.bin.tracks.peaks",pattern=".bed.gz",full.names = T,recursive = T)
foo.chip<-lapply(chip.peaks[grep("H3K27", chip.peaks)],function(i) fread(i) %>%
         setnames(c("chr","start","end","peakid", "score","strand","signalValue","pValue","qValue","peak")) %>% 
         .[,anno:=basename(i) %>% 
             str_split(.,'[.]',simplify = T)%>% .[,1]] %>% 
           .[,assay:="ChIP-seq"]
)%>% rbindlist() %>%
  setkey(chr,start,end) %>% .[,peaktype:="ref"]

blacklist.h<-fread("/Users/yang/data/ref/Blacklist-master/lists/hg38-blacklist.v2.bed.gz")%>% setnames(c("chr","start","end","type")) %>% setkey(chr,start,end)

foo.atac<-fread("/Users/yang/data/ref/Pastor_NCB_primed.naive.hESC.ATACseq/primed.H9.ATAC-seq.hg38.tsv.gz") %>%
  #foo.atac<-fread("/Users/yang/data/scMTR_method.dev/ref.peaks/ATAC/SEACR/b10.rm.Y.M.fragments.tsv_seacr_top0.02.norm.peaks.stringent.bed",select = c(1:3)) %>% 
  setnames (c("chr","start","end")) %>% 
  .[,peakid:=paste0("Peak_",.I)]%>% 
  .[,anno:="Primed"]%>%
  .[,assay:="ATAC-seq"]
fwrite(foo.atac, "/Users/yang/data/ref/Pastor_NCB_primed.naive.hESC.ATACseq/primed.H9.ATACseq.hg38.bed",sep="\t",col.names = F )


ct.peaks<-list.files("/Users/yang/data/scMTR_method.dev/figure1.data/peaks",
                     pattern=".bed.gz",full.names = T) #,recursive = T) 
#ct.peaks<-ct.peaks[grep("STB",ct.peaks)]

foo<-lapply(ct.peaks[c(41:48)],function(f)fread(f,header = F,select = c(1:3))%>% setnames(c("chr","start","end"))%>%
              .[,peakid:=paste0("Peak_",.I)]%>% 
              .[,anno:=basename(f) %>% str_split(.,"_",simplify = T) %>% .[,2] ] %>% 
              .[,sample.anno:=basename(f) %>% str_split(.,"[.]",simplify = T) %>% .[,1] ] %>% 
              .[,assay:=paste0("CUT&Tag_",basename(f) %>% str_split(.,"[.]",simplify = T) %>% .[,2] )] %>% 
                .[,peaktype:=basename(f) %>% str_split(.,"seacr_",simplify = T)%>% .[,2] %>% substring(.,1,7) ]) %>% 
  rbindlist() %>% 
setkey(chr,start,end)



foo.tmp<-rbind(foo[,.(chr,start,end,peakid, anno,assay,peaktype)],foo.chip[,.(chr,start,end,peakid, anno,assay,peaktype)],foo.atac[,.(chr,start,end,peakid, anno="ATAC-seq",assay,peaktype="ref")],fill=T)  #foo.atac

## calculate Fraction of reads in peaks FRiPs
#load peaks

foo.tmp<-foo.tmp[,.(chr,start,end,peakid,peaktype=paste0(anno,":",assay,":",peaktype))]
  setkey(foo.tmp,chr,start,end)

  ######################################################################################################
  ############################## loading MT  CUT&Tag peaks, calculate FRIP #############################
  ######################################################################################################

fragments.files<-list.files ("/Users/yang/data/scMTR_method.dev/figure1.data/fragments",pattern = "tsv",full.names = T)
fread(fragments.files[1])
 
 frag<-lapply(fragments.files,function(f)fread(f,header = F,select = c(1:3))%>% setnames(c("chr","start","end"))%>%
               .[,anno:=basename(f) %>% str_split(.,".fragments",simplify = T) %>% .[,1] ] %>% 
               .[,hist:=basename(f) %>% str_split(.,"_",simplify = T) %>% .[,2] ] %>% 
               .[,type:=basename(f) %>% str_split(.,"[.]",simplify = T) %>% .[,2] ]%>%
                .[,totalnum:=.N,anno]%>%
                setkey(chr,start,end)%>% 
                foverlaps(foo.tmp,.,nomatch = 0L) %>%
                .[,.N,.(anno,hist,type,totalnum,peaktype)]
                ) %>% 
   rbindlist() 
 
 frag[,FRiP:=round(N/totalnum*1000)/1000]
 fwrite(frag,"/Users/yang/data/scMTR_method.dev/plot/H3K27s_CUTTag.FRiP.tsv",sep="\t")
 
 frag[,unique(peaktype)]
 
 
  to.plot<-frag[peaktype %in%c("H3K27ac:CUT&Tag_V2:top0.01","H3K27me3:CUT&Tag_V2:top0.01","H3K27ac:ChIP-seq:ref" , "H3K27me3:ChIP-seq:ref", "ATAC-seq:ATAC-seq:ref")]
 
  unique(to.plot$peaktype)
  to.plot$peaktype<-factor(to.plot$peaktype, levels = c("H3K27me3:ChIP-seq:ref", "H3K27ac:ChIP-seq:ref","H3K27me3:CUT&Tag_V2:top0.01","H3K27ac:CUT&Tag_V2:top0.01","ATAC-seq:ATAC-seq:ref"  ))
  unique(to.plot$anno)

  order<-c("STB_H3K27me3_H3K27me3.V1" ,"STB_H3K27me3_H3K27me3.V2"  ,
  "MTB_H3K27me3_H3K27s_0_5h.V2", "MTB_H3K27me3_H3K27s_1h.V2", "MTB_H3K27me3_H3K27s_2h.V2", "MTB_H3K27me3_H3K27s_4h.V2",
  "MTB_H3K27me3_H3K27s_Igg1.V2" ,  "MTB_H3K27me3_H3K27s_Igg2.V2",  "MTB_H3K27me3_H3K27s_Igg5.V2","MTB_H3K27me3_H3K27s_Igg10.V2",
  "MTB_H3K27me3_H3K27me3_Pre+D.D","MTB_H3K27me3_H3K27me3_Post+D.D"  ,  
  "STB_H3K27ac_H3K27ac.V1" ,"STB_H3K27ac_H3K27ac.V2", 
  "MTB_H3K27ac_H3K27s_0_5h.V2","MTB_H3K27ac_H3K27s_1h.V2", "MTB_H3K27ac_H3K27s_2h.V2", "MTB_H3K27ac_H3K27s_4h.V2",
  "MTB_H3K27ac_H3K27s_Igg1.V2", "MTB_H3K27ac_H3K27s_Igg2.V2","MTB_H3K27ac_H3K27s_Igg5.V2" , "MTB_H3K27ac_H3K27s_Igg10.V2" , 
  "MTB_H3K27ac_H3K27ac_Pre+D.D","MTB_H3K27ac_H3K27ac_Post+D.D"   )
  to.plot$anno<-factor(to.plot$anno,levels = order)
  to.plot<-to.plot[!is.na(anno)]
  to.plot
 pdf("/Users/yang/data/scMTR_method.dev/plot/H3K27s_CUTTag.FRiP.pdf",width = 20,height=8)
 ggplot(to.plot[peaktype%in% c("H3K27me3:ChIP-seq:ref", "H3K27ac:ChIP-seq:ref")],aes(peaktype,FRiP*100,fill=anno))+
   geom_col(position = "dodge")+
   labs(x="Peaksets", y="FRiP (%)")+
   geom_text(data = to.plot[peaktype%in% c("H3K27me3:ChIP-seq:ref", "H3K27ac:ChIP-seq:ref")],aes(label=FRiP*100,x=peaktype, y=FRiP*100 ),position = position_dodge(0.9),color="black",size=4)+
   #facet_wrap(.~hist,scales = "free_y")+
   facet_grid(hist~.,scales = "free_y")+
   theme_bw()+
   theme(
     axis.text.x = element_text(angle=90,hjust=1,vjust = 1)
   )+
   plot.theme
 dev.off()
 
 
