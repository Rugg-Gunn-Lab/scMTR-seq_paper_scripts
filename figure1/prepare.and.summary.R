library(data.table)
library(purrr)
library(ggplot2)
library(viridis)
library(stringr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(patchwork)
set.seed(1234)
io<-list()

#input files:
in.file.dir <- "~/data/Sample_5742_200522YW_BMT2/Lane_8208_200522YW_BMT2/bulk"  
in.file.dir <- "~/data/scMTR_method.dev/Sample_5742_200522YW_BMT2/Lane_8208_200522YW_BMT2/bulk.paired.mapped"  

file.pattern <-  "L001_hs_sorted_rmdup.reads.split.bc.tsv.gz" #L001_hs_reads_summary.xls


files<-list.files(paste0(in.file.dir,"/reads.summary"),pattern="xls",full.names = T)#,recursive = T)
fread(files[1])

foo<-lapply(files, function(i) fread(i,select = c(1,2,3,5)) %>%  ## deduplicated: umi_chr_pos
              setnames(c("cellid", "all_reads", "mapped_reads", "deduplicated_reads")) %>%
              .[,c("mapping_rate","duplication_rate"):=list(mapped_reads/all_reads, 1-deduplicated_reads/mapped_reads)] %>% 
              #  .[,sub.lib:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,3] ] %>%
              .[,lib.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,4] ] %>%
              .[,sub.lib:=basename(i) %>% str_split(.,'_V1_|_V2_|_L001',simplify = T)%>% .[,2] ] 
)%>%
  rbindlist()

foo[,mapped_reads:=mapped_reads/2]
foo[,deduplicated_reads:=deduplicated_reads/2]
foo[,c("mapping_rate","duplication_rate"):=list(mapped_reads/all_reads, 1-deduplicated_reads/mapped_reads)]


foo$sub.lib%>% unique()
foo<-foo[deduplicated_reads>400000] 


#fread("~/data/Sample_5723_200422YW_BMT/Lane_8191_200422YW_BMT/bulk/lane8191_GGCATG_CACTGCAT_MTB_M8_L001_hs_sorted_rmdup.reads.split.bc.tsv.gz")
foo[ substring(cellid,10,11) %in% c("01") ,ab.type:="H3K27ac"]
foo[substring(cellid,10,11) %in% c("02")  ,ab.type:="H3K27me3"]
foo[substring(cellid,10,11) %in% c("09")  ,ab.type:="H3K27ac"]
foo[substring(cellid,10,11) %in% c("10")  ,ab.type:="H3K27me3"]
foo[ substring(cellid,10,11) %in% c("11") ,ab.type:="H3K27ac"]
foo[substring(cellid,10,11) %in% c("12")  ,ab.type:="H3K27me3"]

foo[,anno:=ifelse(substring(sub.lib,6,6)=="s",paste0("MTB_",ab.type,"_",sub.lib),paste0("STB_",ab.type,"_",sub.lib))]

ggplot(foo,aes(ab.type,log10(deduplicated_reads)))+
  geom_col()+
  facet_grid(.~sub.lib)

dev.off()

fwrite(foo,paste0((in.file.dir),"/reads.summary/rna.dna.reads.summary.tsv"),sep = "\t")
foo<-fread(paste0((in.file.dir),"/reads.summary/rna.dna.reads.summary.tsv"))
to.plot<-foo %>% melt()

to.plot[,unique(variable)]
to.plot[variable %in% c("all_reads","mapped_reads","deduplicated_reads"),value:=(value)/1000000]
to.plot[variable %in% c("all_reads","mapped_reads","deduplicated_reads"),unit:="million"]
to.plot[!variable %in% c("all_reads","mapped_reads","deduplicated_reads"),value:=(value)*100]
to.plot[!variable %in% c("all_reads","mapped_reads","deduplicated_reads"),unit:="%"]

means <- aggregate(value ~sub.lib+variable+unit+ ab.type+anno+lib.type, to.plot, mean) %>% as.data.table() %>% .[,value:=round(value*100)/100]

pdf(paste0((in.file.dir),"/reads.summary/mapping.reads.states.pdf"),width = 24,height=8)
ggplot(to.plot,aes(ab.type,value))+
  #  geom_jitter(aes(color=ab.type),alpha=0.6,size=0.1)+
  geom_col(position = "dodge",aes(fill=ab.type))+
  #   stat_summary(fun.y=mean, geom="point", shape=19, size=2,alpha=0.8)+
  geom_text(data = means, aes(label = value,y = value *0.8))+
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.6) +
  
  facet_grid(variable+unit~lib.type+sub.lib,scale="free_y")+
  theme_bw() +
  theme(
    #    strip.background = element_blank(),
    strip.text = element_text(size=12),
    axis.text.x = element_text(size=12, angle = 45, vjust = 0.5),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=14),
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    #  panel.grid = element_blank(),
    panel.background = element_blank()
  )   
dev.off()

pdf(paste0((in.file.dir),"/reads.summary/mapping.reads.states.STB.pdf"),width = 6,height=10)
ggplot(to.plot[substring(anno,1,3)=="STB" &sub.lib %in% c("H3K27me3","H3K27ac")],aes(ab.type,value))+
  #  geom_jitter(aes(color=ab.type),alpha=0.6,size=0.1)+
  geom_col(position = "dodge",aes(fill=ab.type))+
  #   stat_summary(fun.y=mean, geom="point", shape=19, size=2,alpha=0.8)+
  geom_text(data = means[substring(anno,1,3)=="STB" &sub.lib %in% c("H3K27me3","H3K27ac")], aes(label = value,y = value *0.8))+
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.6) +
  
  facet_grid(variable+unit~lib.type,scale="free_y")+
  theme_bw() +
  theme(
    #    strip.background = element_blank(),
    strip.text = element_text(size=12),
    axis.text.x = element_text(size=12, angle = 45, vjust = 0.5),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=14),
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    #  panel.grid = element_blank(),
    panel.background = element_blank()
  )   
dev.off()

########################################################################################################################################
############# after use perl to extract cellid and umi, useful forbi g file.   #############
############# prepare for fragments file: remove low quality cells/backgrounds ############# 
########################################################################################################################################

samples="H9"
#for (samples in c("HNES1","H9","E85")){
  if (samples %in% c("H9", "HNES1")){
    genomes<- "hs"
  } else if (samples == "E85") {
    genomes<-"mm10"
  }
  
  #input files:
  # single lane
  #in.file.dir <- "~/data/Sample_5640_050122YW_V3/Lane_8068_050122YW_V3"  
  # merged lanes

#  file.pattern <-  paste0("_L001_",genomes,"_sorted_merged_rmdup.reads.split.bc.tsv.gz")
  file.pattern <-  paste0("_L001_",genomes,"_sorted_rmdup.reads.split.bc.tsv.gz")
  
  #output files:
  out.file.dir <- in.file.dir #"/bi/home/wangy/projects/ctr_seq/08.matrix"
  out.file.folder <- "fragments"
  
  if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
    dir.create(paste0(out.file.dir,"/",out.file.folder))
  }
 
  
  reads.summary<-fread(paste0(((in.file.dir)),"/reads.summary/rna.dna.reads.summary.tsv"),sep = "\t")

  files<-list.files(in.file.dir,pattern = file.pattern,full.names = T, recursive = T) #%>% .[c(2,6,10,14)]
  chr.size<-fread("/Users/yang/scripts/tools/grch38.chrom.size.txt") %>% setnames(c("chr","len"))
  options(scipen=999) 
  for (i in files){
    # i<-files[1]
   # i<-files[15]
    out.file.name.prefix <- basename(i) %>% sub(".reads.split.bc.tsv.gz","",.)
    sub.libs<-basename(i) %>% str_split(.,'_V1_|_V2_|_L001',simplify = T)%>% .[,2]
    lib.types<-basename(i) %>% str_split(.,'_',simplify = T)%>% .[,4] 
    foo<-fread(i, fill=T, sep="\t", nThread=8 )%>% setnames(c("chr","start","end","umi","cellid","bc4"))
 #   foo[,unique(chr)]
    foo$umi<-NULL
    foo <- unique(foo)
#    foo[,num:=.N,.(chr,start,end,bc4)] 
#    foo[num>1]
    
    if (genomes == "hs"){
      chrs<-paste0("chr",c(1:22,"X","Y","M"))
    } else if (samples == "mm10") {
      chrs<-paste0("chr",c(1:19,"X","Y","M"))
    }
    foo<-foo[chr %in% chrs] 
    

    foo[bc4%in% c(1,9,11) ,ab.type:="H3K27ac"]
    foo[bc4%in% c(2,10,12)  ,ab.type:="H3K27me3"]

    tmp.bc<-foo[,.N,bc4]
    tmp.bc[,ratio:=N/sum(N)]
    foo<-foo[ab.type%in% reads.summary[sub.lib ==sub.libs &lib.type==lib.types,ab.type] &bc4 %in% tmp.bc[ratio>0.1,bc4]]
      
    to.plot<-foo[,.N,.(chr,ab.type)]
    to.plot$chr<-factor(to.plot$chr,levels=chrs)
    to.plot$ab.type<-factor(to.plot$ab.type,levels =
                              c("RNA","IgG","H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","H2AK119ub","H3K27ac_M","RNAPII","RNAPIIser5p","NANOG","OCT4")     
    )
    

    plot.width <- 4*length(foo[,unique(ab.type)])
    
    pdf(paste0(out.file.dir,"/",out.file.folder,"/", out.file.name.prefix,".uni.map.fragments.chr.pdf"),width=8,height = plot.width)
    
    p1<-ggplot(to.plot,aes(chr,N/1000,group=ab.type,fill=chr))+
      geom_col(position = "dodge")+
      facet_grid(ab.type~.,scale="free_y")+
      scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
      theme_bw(base_size = 18) +
      labs(x="",y="# of Unique Mapped Fragments (x 1000)",title = paste0(samples))+
      ggpubr::rotate_x_text(angle = 45)
    print(p1)
    dev.off()
    
    foo$bc4<-NULL
    foo<- merge(foo,reads.summary[sub.lib ==sub.libs &lib.type==lib.types,.(ab.type,lib.type,anno)],by="ab.type")
 #   foo.tmp<-foo[,.N,.(chr,start,end,anno,lib.type)] %>% split(by=c("anno","lib.type"),keep.by=F)
# for bulk assay, for one example:
#    foo[chr=="chr4" &start=="90833803" &end=="90834000"]
#    ab.type  chr    start      end      umi   cellid lib.type                   anno
#    1:  H3K27ac chr4 90833803 90834000 CTTCTGGT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    2:  H3K27ac chr4 90833803 90834000 AATCCGAT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    3:  H3K27ac chr4 90833803 90834000 TAAATTGT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    4:  H3K27ac chr4 90833803 90834000 TAAGGTGT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    5:  H3K27ac chr4 90833803 90834000 TCAACTAT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    6:  H3K27ac chr4 90833803 90834000 TAGCCCTT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    7:  H3K27ac chr4 90833803 90834000 TTCGTTGT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    8:  H3K27ac chr4 90833803 90834000 TCACTTTT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    9:  H3K27ac chr4 90833803 90834000 GTTGGCTT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    10:  H3K27ac chr4 90833803 90834000 GTATAAAT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    11:  H3K27ac chr4 90833803 90834000 GATCACGT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    12:  H3K27ac chr4 90833803 90834000 GCGGGAGT 01:01:01       V2  MTB_H3K27ac_H3K27s_4h
#    13: H3K27me3 chr4 90833803 90834000 ATCAGTTT 01:01:01       V2 MTB_H3K27me3_H3K27s_4h
# most of the umis are not similar, which means they are real unique cuts     
# so we should use UMI to unique reads, also it has been done at the begining step unique
    foo.tmp<-foo[,.N,.(chr,start,end,anno,lib.type)] %>% split(by=c("anno","lib.type"),keep.by=F)
    walk2(foo.tmp, paste0(out.file.dir,"/",out.file.folder,"/", names(foo.tmp),".fragments.tsv"), fwrite,sep="\t",col.names=F)
   
    pdf(paste0(out.file.dir,"/",out.file.folder,"/", out.file.name.prefix,".fragment.length.pdf"),width=8,height = plot.width)
   to.plot<-foo[,end-start,anno]
   p1<- ggplot(to.plot,aes(V1))+
      geom_histogram(binwidth = 10)+
      facet_grid(anno~.)+
      theme_bw(base_size = 18) +
      labs(x="Length of fragments",y="# Fragments",title = paste0(samples))#+
    print(p1)
    dev.off()
    
    
  #  tmp<-merge(foo,chr.size,by="chr")
  #  tmp[,start:=ifelse(start<100,0,start-100)]
  #  tmp[,end:=ifelse(end+100>len,len,end+100)]
  #  tmp$len<-NULL
  #  tmp[,summary(end-start)]
  #  tmp<-tmp[,.N,.(chr,start,end,anno,lib.type)] %>% split(by=c("anno","lib.type"),keep.by=F)
  #  walk2(tmp, paste0(out.file.dir,"/",out.file.folder,"/", names(tmp),".ext200.bed"), fwrite,sep="\t",col.names=F)
    
  }
  
  
  

########################## bedfiles for overlapping ####################
  
  in.file.dir <- "~/data/Sample_5742_200522YW_BMT2/Lane_8208_200522YW_BMT2/bulk"  
  
  samples="H9"
  #for (samples in c("HNES1","H9","E85")){
  if (samples %in% c("H9", "HNES1")){
    genomes<- "hs"
  } else if (samples == "E85") {
    genomes<-"mm10"
  }
  
  #input files:
  file.pattern <-  paste0("_L001_",genomes,"_sorted_rmdup.reads.split.bc.tsv.gz")
  
  #output files:
  out.file.dir <- in.file.dir #"/bi/home/wangy/projects/ctr_seq/08.matrix"
  out.file.folder <- "bedfiles"
  
  if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
    dir.create(paste0(out.file.dir,"/",out.file.folder))
  }
  
  
  reads.summary<-fread(paste0(((in.file.dir)),"/reads.summary/rna.dna.reads.summary.tsv"),sep = "\t")
  chr.size<-fread("/Users/yang/scripts/tools/grch38.chrom.size.txt") %>% setnames(c("chr","len"))
  
  files<-list.files(in.file.dir,pattern = file.pattern,full.names = T,recursive = T)
    options(scipen=999) 
  for (i in files){
    # i<-files[3]
    out.file.name.prefix <- basename(i) %>% sub(".reads.split.bc.tsv.gz","",.)
    sub.libs<-basename(i) %>% str_split(.,'_V1_|_V2_|_L001',simplify = T)%>% .[,2]
    lib.types<-basename(i) %>% str_split(.,'_',simplify = T)%>% .[,4] 
    foo<-fread(i, fill=T, sep="\t", nThread=8 )%>% setnames(c("chr","start","end","umi","cellid","bc4"))
    foo[,unique(chr)]
    
    if (genomes == "hs"){
      chrs<-paste0("chr",c(1:22,"X","Y","M"))
    } else if (samples == "mm10") {
      chrs<-paste0("chr",c(1:19,"X","Y","M"))
    }
    foo<-foo[chr %in% chrs] 
    
    
    foo[bc4%in% c(1,9,11) ,ab.type:="H3K27ac"]
    foo[bc4%in% c(2,10,12)  ,ab.type:="H3K27me3"]
    
    
    foo<-foo[ab.type%in% reads.summary[sub.lib ==sub.libs &lib.type==lib.types,ab.type]]
    
   
    foo$bc4<-NULL
    foo<- merge(foo,reads.summary[sub.lib ==sub.libs &lib.type==lib.types,.(ab.type,lib.type,anno)],by="ab.type")
    foo.tmp<-foo[,.(chr,start,end,anno,lib.type)] %>% split(by=c("anno","lib.type"),keep.by=F)
    walk2(foo.tmp, paste0(out.file.dir,"/",out.file.folder,"/", names(foo.tmp),".bed"), fwrite,sep="\t",col.names=F)
   
  }
  

