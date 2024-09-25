
library(alluvial)
library(ggalluvial)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)

################################################################################
###################### alluvial plot of chromatin states #######################
################################################################################

options(scipen = 9) # avoid 100000->1e+5
mm10<-fread("/Users/wangy/data/ref/mm10.chrom.sizes.txt")%>%  setnames(c("chr", "pos"))
chrs <- paste0("chr", c(1:19, "X", "Y"))

bin<-1000 #5kb
#bin<-500 # 

windows <- mm10[chr %in% chrs, .(start = seq(1, max(pos), by = bin)), chr] %>%
  .[, end := start + bin-2] %>%
  setkey(chr, start, end)
# as chromhmm bins are connected eg bin1 0-9400, bin2 9400-9800, so end subtract extra 1 to avoid false overlap
windows[,id:=paste0(chr,":",start,"-",end+1)]
setkey(windows,chr,start,end)


inputdir<-"/Users/wangy/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments/chromhmm_1000/states_15"

files<-list.files(inputdir,pattern = "txt",full.names = T)
fread(files[1])
to.plot<-as.data.frame(fread(files[1])) %>% tibble::column_to_rownames("State (Emission order)")
pdf(paste0(inputdir,"/emission.stats.pdf"),width = 6,height = 10)
pheatmap::pheatmap(to.plot[c(5,6,14,11,12,8,9,15,1,4,2,3,7,10,13),],
                   color = colorRampPalette(c("white","#00007F"
                                              ))(50) ,
                   cluster_rows  = F,
                   main = "Emission stats",
                  display_numbers = T, 
                  number_color="white",
                   cellwidth = 20,
                   cellheight = 20,
                   border_color = "black")
dev.off()

grep("^RPL|^RPS",hs.geneanno$symbol)
files[grep("overlap.txt",files)]

foo<-lapply(files[grep("overlap.txt",files)], function(f) fread(f,select = c(1:2))%>% 
              setnames(c("states","genome"))%>% 
              .[, anno:=basename(f)%>% str_split(.,"_overlap.txt",simplify = T)%>% .[,1]]
            
)%>% rbindlist()

foo[,mean(genome),states]

dcast(foo,states~anno,value.var = "genome") %>% as.data.frame() %>% tibble::column_to_rownames("states")
fread(files[5])

pdf(paste0(inputdir,"/genome.percentage.pdf"),width = 6,height = 10)
pheatmap::pheatmap(foo[!states=="Base",mean(genome),states] %>% 
                     as.data.frame() %>% 
                     tibble::column_to_rownames("states"),
                   color = colorRampPalette(c("white","red"
                   ))(100) ,
                   cluster_rows  = F,
                   cluster_cols = F,
                   main = "Genome %",
                   display_numbers = T,                   cellwidth = 20,
                   cellheight = 20,
                   border_color = "black")
dev.off()

######################## keep all posterior probability ########################

#
inputdir<-"/Users/wangy/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments/chromhmm_downsample"
files<-list.files(inputdir,pattern = "segments.bed",full.names = T,recursive = T)

for (sampl in  c("EPI","PE","TE")) {#"D0",
  #sampl<-"EPI",
  for (downsample in c("downsample1","downsample2","downsample3")){
  files<-list.files(paste0(inputdir,"/",downsample),pattern = "segments.bed",full.names = T)
  foo<-lapply(files[grep(sampl,files)], function(f) fread(f,select = c(1:4),skip = 1)%>% 
                setnames(c("chr","start","end","states"))%>% 
                .[, anno:=basename(f)%>% str_split(.,"_",simplify = T)%>% .[,1]]
              
  )%>% rbindlist()
  
  foo<-foverlaps(foo,windows,nomatch = 0L)
  
  fwrite(foo[,.(id,states)],paste0(inputdir,"/",downsample,"/",sampl,".kept.all.posterior.probability.tsv"),sep = "\t")
}
}

################################ change color ##################################

files<-list.files(inputdir,pattern = "segments.bed.gz",full.names = T,recursive = T)

foo<-lapply(files[], function(f) fread(f,skip=1)%>% 
              .[, anno:=basename(f)%>% str_split(.,"_",simplify = T)%>% .[,1]] %>%
              .[, downsample:=basename(dirname(f))]
            
)%>% rbindlist()

foo[,V5:=0]
foo[,V6:="."]
foo[,V7:=V2]
foo[,V8:=V3]

foo[V4 %in% c("E11","E12"), V9:="#EE766F"]
foo[V4 %in% c("E8","E9"), V9:="#C49B05"]
foo[V4 %in% c("E5","E6"), V9:="#59B031"] 

foo[V4 %in% c("E14"), V9:="#36B28F"]
foo[V4 %in% c("E15"), V9:="#11B6EB"]
foo[V4 %in% c("E1"), V9:="#9188C0"] #Bivalent enhancer
foo[V4 %in% c("E4"), V9:="grey60"] 
foo[V4 %in% c( "E2","E3","E7","E10","E13"), V9:="white"]

library(purrr)
foo.tmp<-foo%>% split(by=c("anno","downsample"),keep.by = F)
walk2(foo.tmp , paste0(inputdir,"/", names(foo.tmp),".dense.recolored.bed"), fwrite,sep="\t",col.names=F)
#manual add: 
# track name="D0_15" description="D0_15 (Emission ordered)" visibility=1 itemRgb="On"


############################ add  chromatin states #############################

inputdir<-"/Users/wangy/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments/chromhmm_downsample"

files<-list.files(inputdir,pattern = ".kept.all.posterior.probability.tsv",full.names = T,recursive = T)
foo<-lapply(files, function(f) fread(f)%>% 
              .[, anno:=basename(f)%>% str_split(.,".kept.all.posterior.probability.tsv",simplify = T)%>% .[,1]]%>%
              .[, downsample:=basename(dirname(f))]
            
)%>% rbindlist()

foo[states %in% c("E11","E12"), state.type:="Active.enhancer"]
foo[states %in% c( "E2","E3","E7","E10","E13"), state.type:="Unmarked"]
foo[states %in% c("E8","E9"), state.type:="Primed.enhancer"]
foo[states %in% c("E5","E6"), state.type:="Active.Promoter"] 
foo[states %in% c("E4"), state.type:="Mixed"] 
foo[states %in% c("E14"), state.type:="Active.Transcription"]
foo[states %in% c("E1"), state.type:="Heterochromatin"]
foo[states %in% c("E15"), state.type:="Repressed.Polycomb"] #Bivalent enhancer

##########################plot dynamics of all states ##########################

to.plot<-dcast(foo,id~anno+downsample,value.var = "state.type")

for (i in to.plot[,unique(EPI_downsample1)]) {
  # i<-"Heterochromatin.Repeats"
  kept.id<-foo[state.type==i,unique(id)]
  sub.to.plot.tmp<-to.plot[id %in% kept.id]
  
 # sub.to.plot.tmp<-to.plot[D0 %in% c(i,j) |D1 %in% c(i,j)|D2 %in% c(i,j)|D3 %in% c(i,j)]
  sub.to.plot.tmp<-as.data.frame(sub.to.plot.tmp[,.(freq=.N),.( EPI_downsample1,PE_downsample1,PE_downsample2,PE_downsample3,TE_downsample1,TE_downsample2,TE_downsample3)] )
  pdf(paste0(inputdir,"/alluvial.plot.states",i,".pdf"),width = 8,height = 4)
#  pdf(paste0(inputdir,"/alluvial.plot.states.promoters.pdf"),width = 8,height = 4)
  p1<-ggplot(sub.to.plot.tmp,
             aes(y=freq,
                 axis1= EPI_downsample1,
                 axis2 = PE_downsample1, 
                 axis3 = PE_downsample2,
                 axis4 = PE_downsample3,
                 axis5 = TE_downsample1,
                 axis6 = TE_downsample2,
                 axis7 = TE_downsample3
             ) )+
    geom_alluvium(aes(fill =as.factor (EPI_downsample1)), width = 1/2) +
    geom_stratum(aes(fill = as.factor (EPI_downsample1)), width = 1/2) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
    theme_void()
  print(p1)
  dev.off()
  fwrite(sub.to.plot.tmp,paste0(inputdir,"/alluvial.states",i,".tsv"),sep="\t")
}




