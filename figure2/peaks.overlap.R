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
############################## peaks num, overlap ##############################
################################################################################

#de.novo.macs
peaks.macs<-list.files("/Users/yang/data/scMTR_method.dev/figure3.data",pattern="peaks.rds",full.names = T)
peaks.macs<-lapply(peaks.macs,function(i) readRDS(i) %>% 
                     as.data.table() %>% .[,.(chr=seqnames,start, end)]%>%
                     .[,anno:=basename(i) %>% sub(".peaks.rds","",.)]) %>% rbindlist()
#seacr
peaks.seacr <- list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/SEACR", pattern="tsv.normalized.bedgraph.stringent.be",full.name=T)
peaks.seacr<-lapply(peaks.seacr,function(i) fread(i,select = c(1:3)) %>% setnames(c("chr","start","end")) %>%
                      .[,anno:=basename(i) %>% sub(".fragments.tsv.normalized.bedgraph.stringent.bed","",.)]) %>% rbindlist()


#used

peaks.seacr <- list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/SEACR", pattern=".vs.IgG.stringent.bed",full.name=T)
peaks.seacr1<-lapply(peaks.seacr[c(1,3,4)],function(i) fread(i,select = c(1:3)) %>% setnames(c("chr","start","end")) %>%
                       .[,anno:=basename(i) %>% sub(".fragments.tsv.normalized.bedgraph.vs.IgG.stringent.bed","",.)]) %>% rbindlist()

peaks.seacr <- list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/SEACR", pattern=".relaxed.bed",full.name=T)
peaks.seacr2<-lapply(peaks.seacr[c(2,5)],function(i) fread(i,select = c(1:3)) %>% setnames(c("chr","start","end")) %>%
                       .[,anno:=basename(i) %>% sub(".fragments.tsv.normalized.bedgraph.vs.IgG.relaxed.bed","",.)]) %>% rbindlist()

peaks.seacr<-rbind(peaks.seacr1, peaks.seacr2)


#
#individual ct

peaks.ct <- list.files("/Users/yang/data/scMTR_method.dev/bulk.cuttag/fragments/SEACR", pattern="0.01",full.name=T)
peaks.ct<-lapply(peaks.ct[c(4,5,6,7,8)],function(i) fread(i,select = c(1:3)) %>% setnames(c("chr","start","end")) %>%
                   .[,anno:=basename(i) %>% sub(".fragments.tsv_seacr_top0.01.norm.peaks.stringent.bed","",.)]) %>% rbindlist()

#chip-seq
peaks.chip <- list.files("/Users/yang/data/scMTR_method.dev/ref.peaks/ref.peaks", pattern="bed.gz",full.name=T)
peaks.chip<-lapply(peaks.chip[c(3,4,5,6,7)],function(i) fread(i,select = c(1:3)) %>% setnames(c("chr","start","end")) %>%
                     .[,anno:=basename(i) %>% str_split(.,"[.]",simplify = T) %>% .[,1] ]) %>% rbindlist()

peaks.chip[,assay:="ChIP-seq"]
peaks.ct[,assay:="CUT&Tag"]
peaks.seacr[,assay:="SEACR"]

peaks.seacr[,anno:= sub("Def.endo.D0.","",anno)
            peaks.ct[,anno:=str_split(anno,"[.]",simplify = T) %>% .[,1] ]
            
            foo<-rbind(peaks.chip,peaks.ct,peaks.seacr) #peaks.macs
            foo$anno %>% unique
            foo$anno <-factor(foo$anno,levels = c( "H3K27me3","H3K27ac" ,"H3K4me1" , "H3K4me3", "H3K36me3" ,  "IgG"))
            foo$assay %>% unique
            foo$assay <-factor(foo$assay,levels = c("MACS2_narrow" , "MACS2_broad"  ,"SEACR","CUT&Tag",  "ChIP-seq"   ))
            
            #foo<-setkey(foo,anno,assay,chr,start,end)
            foo[,hist:=anno]
            foo[,anno:=paste0(anno,"_",assay)]
            foo<-setkey(foo,assay,hist,anno,chr,start,end)
            foo[,peakid:=paste0("peak_",.I),]
            foo.tmp<-copy(foo)
            foo.tmp<-split(foo.tmp,by=c("anno"))
            to.label<-foo[,.N,.(assay,hist)]
            pdf("/Users/yang/data/scMTR_method.dev/figure3.data/number.of.peaks.pdf",width=8,height = 5)
            ggplot(foo[,.N,.(assay,hist)],aes(hist,N,fill=hist)) +
              geom_col()+
              facet_grid(.~assay)+
              theme_bw()+
              scale_fill_manual(values=c(H3K27me3 =H3K27me3.color,H3K27ac=H3K27ac.color, H3K4me1= H3K4me1.color ,H3K4me3=H3K4me3.color ,H3K36me3=H3K36me3.color, IgG =IgG.color)) +
              labs(x="",y="Number of peaks")+
              geom_text(data = to.label, aes(label = N, y = N*1.1),size=4)+
              theme(
                #    strip.background = element_blank(),
                strip.text = element_text(size=16),
                #  strip.background = element_blank(),
                axis.text.x = element_text(size=16, angle = 45, vjust = 1,hjust=1),
                axis.text.y = element_text(size=16),
                axis.title = element_text(size=16),
                legend.justification = "center",
                legend.title = element_blank(),
                legend.text = element_blank(),
                # panel.border = element_blank(),
                # panel.grid.major = element_blank(),
                #  panel.grid = element_blank(),
                panel.background = element_blank()
              ) 
            dev.off()
            
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
            dev.off()
            
            pdf(paste0(in.file.dir,"/peaksets.overlap1.pdf"),width = 8,height=8)
            p1<-ggplot(data = tmp, aes(anno1, anno2, fill =  ratio.2))+
              geom_tile(color = "white")+
              scale_fill_viridis(discrete = F, begin = 0.1, end = 0.9, option = "magma", limit = c(0,100), space = "Lab", alpha = 0.8,name="Percentage") +
              theme_minimal()+
              theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                               size = 10, hjust = 1))+
              theme(axis.text.y = element_text(angle = 0, vjust = 0.5,
                                               size = 10, hjust = 1))+
              geom_text(aes(label = ratio.2), color = "white", size = 5) +
              theme(
                panel.grid = element_blank(),
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
            
            pdf(paste0(in.file.dir,"/peaksets.overlap2.pdf"),width = 8,height=8)
            
            p2<-ggplot(data = tmp, aes(anno2, anno1, fill =  ratio.1))+
              geom_tile(color = "white")+
              scale_fill_viridis(discrete = F, begin = 0.1, end = 0.9, option = "magma", limit = c(0,100), space = "Lab", alpha = 0.8,name="Percentage") +
              theme_minimal()+
              theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                               size = 10, hjust = 1))+
              theme(axis.text.y = element_text(angle = 0, vjust = 0.5,
                                               size = 10, hjust = 1))+
              geom_text(aes(label = ratio.1), color = "white", size = 5) +
              theme(
                panel.grid = element_blank(),
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
            
            
            pdf(paste0(in.file.dir,"/peaksets.overlap3.pdf"),width = 8,height=8)
            
            p2<-ggplot(data = tmp, aes(anno1, anno2, fill =  ratio.1))+
              geom_tile(color = "white")+
              scale_fill_viridis(discrete = F, begin = 0.1, end = 0.9, option = "magma", limit = c(0,100), space = "Lab", alpha = 0.8,name="Percentage") +
              theme_minimal()+
              theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                               size = 10, hjust = 1))+
              theme(axis.text.y = element_text(angle = 0, vjust = 0.5,
                                               size = 10, hjust = 1))+
              geom_text(aes(label = ratio.1), color = "white", size = 5) +
              theme(
                panel.grid = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.x = element_line(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.ticks = element_blank()#,
              ) +
              coord_fixed()
            print(p2)
            dev.off()
            
            saveRDS(signac,paste0(in.file.dir,"/seurat.rna.hist_bin5k.rds"))
            
            
            
            