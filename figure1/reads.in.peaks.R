
# calculate counts per peak to assess the off-target, specificity, among histone marks in different conditions 
library(data.table)
library(purrr)
library(ggplot2)
library(viridis)
library(stringr)
library(bedtoolsr)
library(ggpointdensity)
source("~/scripts/utility.R")
library(R.utils)
library(wesanderson)
#wes_palettes

#################################################################
############### loading ENCODE chipseq ref peaks ################
#################################################################
peaks.files<-list.files("/Users/wangy/data/ref/encode.h9.ren.bin.tracks.peaks",pattern=".bed.gz",full.names = T,recursive = T)


#################################################################
############ loading fragment files from different assays########
#################################################################
in.file.dir <-"/Users/wangy/data/scMTR_method.dev/figure1.data/fragments"

bed.files<-list.files(in.file.dir,pattern = ".tsv.gz",full.names = T)%>%unique()

#output files:
out.file.dir <- dirname(in.file.dir) 
out.file.folder <- "ov.peaks"

if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
  dir.create(paste0(out.file.dir,"/",out.file.folder))
}

for (peaks in 1:2){ #length(peaks.files)
#  peaks<-1
  peakdata<-fread(peaks.files[peaks],select = c(1:4)) %>% setnames(c("chr","start","end","peakid")) %>% setkey(chr,start,end)
  for (beds in 1:length(bed.files)) { #

    
    ## calculate reads counts in peaksets
    foo<-bedtoolsr::bt.coverage(peakdata,
                              fread(bed.files[beds],select = c(1:3)) %>% setnames(c("chr",'start','end')) %>% setkey(chr,start,end))%>% 
  as.data.table() %>% 
  setnames(c("chr","start","end","peakid","counts","coverage.len","ref.len","coverage.rate"))%>% 
  .[,ref.peak:=basename(dirname(peaks.files[peaks]))]%>% 
.[,hist:=basename(bed.files[beds]) %>% sub(".fragments.tsv.gz","",.)] 
    
    out.file.name.prefix<- paste0("Refpeak_",basename(dirname(peaks.files[peaks])),"_VS_HIST_",basename(bed.files[beds]) %>% sub(".fragments.tsv.gz","",.))
png(paste0(out.file.dir,"/",out.file.folder,"/", out.file.name.prefix,".densityplot.png"),width = 500,height = 400)    
p1<-ggplot(foo,aes(counts,coverage.rate))+
#  geom_point(alpha=0.3)+
  geom_pointdensity() +
  ggtitle(label = paste0("Refpeak_",basename(dirname(peaks.files[peaks])),subtitle=paste0("n=",nrow(foo)," peaks, vs ", basename(bed.files[beds]) %>% sub(".fragments.tsv.gz","",.))))+
  scale_color_viridis() +
  theme_bw()+
  plot.theme
print(p1)
dev.off()

fwrite(foo, paste0(out.file.dir,"/",out.file.folder,"/", out.file.name.prefix,".tsv"),sep="\t")

system(paste0 ("gzip ", paste0(out.file.dir,"/",out.file.folder,"/", out.file.name.prefix,".tsv"))) ## require library(R.utils)


  }
}

#### get sequencing depth for normalization

in.file.dir <- paste0(out.file.dir,"/",out.file.folder) #"/bi/home/wangy/projects/ctr_seq/08.matrix"
ov.files<-list.files(in.file.dir,pattern = ".tsv.gz") #%>% gsub(".ext200","",.) %>%unique()
#ov.files<-list.files(in.file.dir,pattern = ".ext200.bed.gz.tsv.gz") 

tmp<- ov.files %>% data.table() %>% setnames("filename")

tmp[,refpeak:=str_split(filename,"[.]",simplify = T)%>% .[,1] %>% sub("Refpeak_","",.)]
#tmp[,source.lab:=str_split(filename,"[.]",simplify = T)%>% .[,3] ]
tmp[,hist:=str_split(filename,"_",simplify = T)%>% .[,6]] # %>% sub(".D|.D.bed.tsv|.D.rmdup.bed.tsv","",.) ]
tmp[,type:=str_split(filename,"_",simplify = T)%>% .[,5] ]
tmp[,anno:=str_split(filename,"_HIST_|.tsv.gz",simplify = T)%>% .[,2] %>%  substring(.,nchar(paste0(type,"_",hist,"_"))+1,100)]

bed.files<-list.files(paste0(dirname(in.file.dir),"/fragments"),pattern = ".tsv.gz",full.names = T)
tmp2<-lapply(bed.files,function(i)fread(i) %>% .[,filename:=basename((i))]) %>% 
    rbindlist() %>% .[,.N,filename]
fwrite(tmp2, paste0(in.file.dir,"/sequencing.depth.tsv"),sep="\t")
#  tmp2<-fread( paste0(in.file.dir,"/sequencing.depth.tsv"))

tmp[,samplename:=str_split(filename,"HIST_|.tsv.gz",simplify = T)%>% .[,2]] 

tmp<-merge(tmp[,.( filename,  refpeak, hist, type, anno,samplename)],tmp2[,.(samplename=filename%>% sub(".fragments.tsv.gz","",.),seqdepth=N)],by="samplename")

tmp[,in.file.dir:=in.file.dir]
fwrite(tmp, paste0(in.file.dir,"/plot.sample.metadata.tsv"),sep="\t")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#################### plot same antibody from different assays ################## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

metadata<-fread(paste0(in.file.dir,"/plot.sample.metadata.tsv"))

#output files:
out.file.dir <-  dirname(in.file.dir)#"/bi/home/wangy/projects/ctr_seq/08.matrix"
out.file.folder <- "ov.peaks.plot/ab.1v1.sample.2"

if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
  dir.create(paste0(out.file.dir,"/",out.file.folder),recursive = T)
}


for (plot.ref in metadata[,unique(refpeak)]){
#  for (plot.ab in metadata[,unique(hist)][c(1,2,4)]) {
  #plot.ref<-"h3k27ac"
  #plot.ab<-"H3K27me3"
  plot.ab<-  plot.ref
#  for (sourcelab in metadata[,unique(source.lab)]) {
   # sourcelab<-"r"
plot.files<-metadata[refpeak ==plot.ref  &tolower(hist)==plot.ab]#& source.lab== sourcelab



plot.files$anno<-factor(plot.files$anno,levels =
                          c( "H3K27ac.V2" ,  "H3K27ac.V1" , "H3K27me3.V2","H3K27me3.V1" ))#  "H3K27ac_r.V2"  , "H3K27me3_r.V2" ,  
                          #  "H3K27s_r1.V2" , "H3K27s_r2.V2", "H3K27s_0_5h.V2","H3K27s_1h.V2", "H3K27s_2h.V2", "H3K27s_4h.V2", 
                           # "H3K27s_Igg1.V2", "H3K27s_Igg2.V2", "H3K27s_Igg5.V2",   "H3K27s_Igg10.V2","H3K27ac_Pre+D.D", 
                            #"H3K27ac_Post+D.D","H3K27me3_Pre+D.D", "H3K27me3_Post+D.D" ))
setkey(plot.files,anno)
pairwise <- combn(plot.files[!is.na(anno),anno] %>% as.character, 2) %>%
  as.data.frame()

for (i in 1:ncol(pairwise)){
#  i<-1
sample.x<-pairwise[1,i]
sample.y<-pairwise[2,i]
#sample.x<- "STB_25k"
#sample.y<- "MTB"
foo<-lapply(c(sample.x,sample.y),function(i) fread(plot.files[anno==i,paste0(in.file.dir,"/",filename)],select=c(1:8)) %>% .[,anno:= i] ) %>% 
  rbindlist() %>% 
  merge(.,plot.files,by="anno")

foo[,counts:=counts/seqdepth*1000000]
to.plot<-dcast(foo, refpeak+peakid~anno,value.var = "counts")

select_cols <- c("refpeak","peakid",sample.x,sample.y)
to.plot<-to.plot[, ..select_cols] %>% setnames(c("refpeak","peakid","plot.x","plot.y"))

cors<-cor(x=to.plot$plot.x,y=to.plot$plot.y) 
xmax<-max(to.plot$plot.x)
ymax<-max(to.plot$plot.y)

out.file.name.prefix<-paste0("RefPeak.",plot.ref, "_PlotAb.", plot.ab,"_",sample.x,".vs.",sample.y)

#png(paste0(out.file.dir,"/",out.file.folder,"/", out.file.name.prefix,".png"),width = 800,height = 600,res = 150)    
pdf(paste0(out.file.dir,"/",out.file.folder,"/", out.file.name.prefix,".pdf"),width = 6,height = 4)    

model <- lm(plot.y ~ plot.x, data = to.plot) 
cod<-round( summary(model)$adj.r.squared,digits = 2)  #R squared

p1<-ggplot(to.plot,aes(plot.x,plot.y,color=refpeak))+
  geom_point(alpha=0.8,size=2) +
  scale_color_manual(breaks = c("h3k27me3","h3k27ac"),
                     values=c("#E6332A", "#365FAA"))+  #"#E69F00", 
  labs(x=paste0(sample.x," (RPM)"),y=paste0(sample.y," (RPM)"), title = paste0("Ref.Peak: ",plot.ref, ",", "Plot.Ab: ", plot.ab),subtitle = paste0(sample.x," vs ",sample.y,", N=",nrow(to.plot),"\n Cor=",round(cors,digits = 2),", adj.R2=",cod))+
  #  xlim(c(0,200))+
  # ylim(c(0,200))+
  #geom_ln(slope = 1)+
  theme_classic()+
  coord_fixed(ratio =xmax/ymax)+
  plot.theme
print(p1)
dev.off()

}
}


dev.off()


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#################### plot different antibody from same assay ################### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

metadata[,anno.type:=anno %>% sub("H3K27ac|H3K27me3","",.)]

for (samples in metadata[,unique(anno.type)]){
#  samples<- "H3K27s_0_5h.V2"
  # samples<-".V1"
 # for (sourcelab in metadata[,unique(source.lab)]) {
    # sourcelab<-"r"
    plot.files<-metadata[anno.type ==samples ]#& source.lab== sourcelab ]
    
  pairwise <- combn(c("H3K27me3","H3K27ac"), 2) %>%
    as.data.frame()
    
    for (i in 1:ncol(pairwise)){
      #i=1
      ab.y<-pairwise[1,i]
      ab.x<-pairwise[2,i]
      #<- "STB_25k"
      #sample.y<- "MTB"
      plot.files<-metadata[anno.type ==samples  & refpeak %in% tolower(c(ab.x,ab.y)) & hist %in% c(ab.x,ab.y)]
      
      foo<-lapply(1:nrow(plot.files),function(i) fread(plot.files[i,paste0(in.file.dir,"/",filename)],select=c(1:8)) %>% .[,filename:= plot.files[i,filename]] ) %>% 
        rbindlist() %>% 
        merge(.,plot.files,by="filename")
      
      foo[,counts:=counts/seqdepth*1000000]
      to.plot<-dcast(foo, refpeak+peakid~hist,value.var = "counts")
      
      select_cols <- c("refpeak","peakid",ab.x,ab.y)
      to.plot<-to.plot[, ..select_cols] %>% setnames(c("refpeak","peakid","plot.x","plot.y"))
      
      cors<-cor(x=to.plot$plot.x,y=to.plot$plot.y) 
      xmax<-max(to.plot$plot.x)
      ymax<-max(to.plot$plot.y)
      
      out.file.name.prefix<-paste0("RefPeak.",ab.x,".vs.",ab.y, "_PlotAb.", ab.x,".vs.",ab.y,"_",samples)
      #"_",sourcelab
      if (!file.exists(paste0(out.file.dir,"/",dirname(out.file.folder),"/ab.2v2"))){
        dir.create(paste0(out.file.dir,"/",dirname(out.file.folder),"/ab.2v2"))
      }

      to.label<-to.plot[,.N,refpeak]%>% .[,ypos:=.I]
      
      
      model <- lm(plot.y ~ plot.x, data = to.plot) 
    #  cod<-round( summary(model)$adj.r.squared,digits = 2)  #R squared
      cod<-format(summary(model)$adj.r.squared, scientific = T,digits = 2)
      
      png(paste0(out.file.dir,"/",dirname(out.file.folder),"/ab.2v2/", out.file.name.prefix,".png"),width = 800,height = 600,res = 150)    
     # pdf(paste0(out.file.dir,"/",dirname(out.file.folder),"/ab.2v2/", out.file.name.prefix,".pdf"),width = 6,height = 4)    
      p1<-ggplot(to.plot,aes(plot.x,plot.y))+
        geom_point(aes(color=refpeak),alpha=0.6,size=2) +
        scale_color_manual(breaks = c("h3k27me3","h3k27ac"),
                           values=c("#E6332A", "#365FAA"))+  #"#E69F00", 
        #   geom_pointdensity() +
        labs(x=paste0(ab.x," (RPM)"),y=paste0(ab.y," (RPM)"), title = paste0(samples),subtitle = paste0("Ab: ",ab.x,".vs.",ab.y,",\n Cor=",round(cors,digits = 2),", adj.R2=",cod))+
        
        #  xlim(c(0,200))+
        # ylim(c(0,200))+
  #    geom_text(data = to.label,aes(label=paste0(refpeak,":",N),x= xmax*0.6 , y= ymax*(1-ypos*0.1)),color="black",size=4)+
  #      geom_abline(slope = 1,linetype="dashed",alpha=0.8)+
        theme_classic()+
        coord_fixed(ratio =xmax/ymax)+
        plot.theme
      print(p1)
      dev.off()
      
      
      }}






