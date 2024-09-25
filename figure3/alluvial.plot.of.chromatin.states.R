
################################################################################
###################### alluvial plot of chromatin states #######################
################################################################################

library(alluvial)
library(ggalluvial)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)

options(scipen = 9) # avoid 100000->1e+5
hsGRCh38<-fread("/Users/wangy/data/ref/hg38.chrom.sizes.txt")%>%  setnames(c("chr", "pos"))
chrs <- paste0("chr", c(1:22, "X", "Y"))

bin<-200 #5kb
#bin<-500 # 

windows <- hsGRCh38[chr %in% chrs, .(start = seq(1, max(pos), by = bin)), chr] %>%
  .[, end := start + bin-2] %>%
  setkey(chr, start, end)

# as chromhmm bins are connected eg bin1 0-9400, bin2 9400-9800, so end subtract extra 1 to avoid false overlap
windows[,id:=paste0(chr,":",start,"-",end+1)]
setkey(windows,chr,start,end)

inputdir<-"/Users/wangy/data/scMTR_method.dev/figure4.data/chromhmm_adj/states_15"

files<-list.files(inputdir,pattern = "txt",full.names = T)
fread(files[13])
pdf(paste0(inputdir,"/emission.stats.pdf"),width = 6,height = 10)
pheatmap::pheatmap(as.data.frame(fread(files[13]))%>% tibble::column_to_rownames("State (Emission order)"),
                   color = colorRampPalette(c("white",  "grey80","#00007F"
                                              ))(100) ,
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

for (sampl in  c("D0","D1","D2","D3")) {#"D0",
  #sampl<-"D0",
  files<-list.files(inputdir,pattern = "dense.bed",full.names = T)
  foo<-lapply(files[grep(sampl,files)], function(f) fread(f,select = c(1:4),skip = 1)%>% 
                setnames(c("chr","start","end","states"))%>% 
                .[, anno:=basename(f)%>% str_split(.,"_dense.bed",simplify = T)%>% .[,1]]
              
  )%>% rbindlist()
  foo<-foverlaps(foo,windows,nomatch = 0L)
  fwrite(foo[,.(id,states)],paste0(inputdir,"/",sampl,".kept.all.posterior.probability.tsv"),sep = "\t")
}


files<-list.files(inputdir,pattern = ".kept.all.posterior.probability.tsv.gz",full.names = T)
foo<-lapply(files, function(f) fread(f)%>% 
              .[, anno:=basename(f)%>% str_split(.,".kept.all.posterior.probability.tsv.gz",simplify = T)%>% .[,1]]
            
)%>% rbindlist()

################################ change color ##################################
files<-list.files(inputdir,pattern = "dense.bed",full.names = T)
foo<-lapply(files, function(f) fread(f,skip=1)%>% 
              .[, anno:=basename(f)%>% str_split(.,".dense.bed",simplify = T)%>% .[,1]]
            
)%>% rbindlist()
foo[V4 %in% c(6,7), V9:="#EE766F"]
foo[V4 %in% c(9,10), V9:="#C49B05"]
foo[V4 %in% c(12,13), V9:="#59B031"]
foo[V4 %in% c(5), V9:="#36B28F"]
foo[V4 %in% c(14), V9:="#11B6EB"]
foo[V4 %in% c(1,3), V9:="#9188C0"]
foo[V4 %in% c(2,4,8,11,15), V9:="#FFFFFF"]
library(purrr)
foo.tmp<-foo%>% split(by=c("anno"),keep.by = F)
walk2(foo.tmp , paste0(inputdir,"/", names(foo.tmp),".dense.recolored.bed"), fwrite,sep="\t",col.names=F)

foo[states %in% c(1,3), state.type:="Transcription"]
foo[states %in% c( 2,4,8,11,15), state.type:="Unmarked"]
foo[states %in% c(5), state.type:="Primed enhancer"]
foo[states %in% c(6,7), state.type:="Active enhancer"]
foo[states %in% c(9,10), state.type:="Active promoter"]
foo[states %in% c(12,13), state.type:="Bivalent promoter"]
foo[states %in% c(14), state.type:="Repressed Polycomb"] #Bivalent enhancer

to.plot<-dcast(foo,id~anno,value.var = "state.type")

##########################plot dynamics of all states ##########################

sub.to.plot<-to.plot[!(D0=="Unmarked" & D1=="Unmarked" & D2=="Unmarked" &D3=="Unmarked")]
sub.to.plot.tmp<-as.data.frame(sub.to.plot[,.(freq=.N),.( D0,D1,D2,D3)])# %>%.[freq>=nrow(sub.to.plot.tmp)/1000]) #)#
sub.to.plot.tmp$freq %>% .[.>=1000] %>% hist()
pdf(paste0(inputdir,"/alluvial.plot.rm.allunmarked.pdf"),width = 8,height = 4)
p1<-ggplot(sub.to.plot.tmp[sub.to.plot.tmp$freq>=100,],
           aes(y=freq,
               axis1=D0,
               axis2 = D1, 
               axis3 = D2,
               axis4=D3) )+
  geom_alluvium(aes(fill = (D3)), width = 1/2) +
  geom_stratum(aes(fill =  (D3)), width = 1/2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  #geom_label(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  scale_x_discrete(limits = c("D0", "D1","D2","D3"),
                   expand = c(0.15, 0.05)) +
  #  scale_fill_manual(values = c("#e84f35", "#8abf57", "#f5ae94","#b1b1b1", "#8386c0")) +
  theme_void()
print(p1)
dev.off()


fwrite(sub.to.plot.tmp,paste0(inputdir,"/alluvial.rm.allunmarked.tsv"),sep="\t")


######################## plot dynamics of promoter #############################

sub.to.plot<-to.plot[D0%in% c("Active enhancer","Bivalent promoter") | D1%in% c("Active promoter","Bivalent promoter")  | D2%in% c("Active promoter","Bivalent promoter")  |D3%in% c("Active promoter","Bivalent promoter") ]
sub.to.plot.tmp<-as.data.frame(sub.to.plot[,.(freq=.N),.( D0,D1,D2,D3)])# %>%.[freq>=nrow(sub.to.plot.tmp)/1000]) #)#
sub.to.plot.tmp$freq %>% .[.>=1000] %>% hist()
pdf(paste0(inputdir,"/alluvial.plot.promoter.pdf"),width = 8,height = 4)
p1<-ggplot(sub.to.plot.tmp[sub.to.plot.tmp$freq>=100,],
           aes(y=freq,
               axis1=D0,
               axis2 = D1, 
               axis3 = D2,
               axis4=D3) )+
  geom_alluvium(aes(fill = (D3)), width = 1/2) +
  geom_stratum(aes(fill =  (D3)), width = 1/2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  #geom_label(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  scale_x_discrete(limits = c("D0", "D1","D2","D3"),
                   expand = c(0.15, 0.05)) +
  #  scale_fill_manual(values = c("#e84f35", "#8abf57", "#f5ae94","#b1b1b1", "#8386c0")) +
  theme_void()
print(p1)
dev.off()

######################## plot dynamics of promoter #############################

sub.to.plot<-to.plot[D0%in% c("Active promoter","Primed enhancer") | D1%in% c("Active promoter","Primed enhancer")  | D2%in% c("Active promoter","Primed enhancer")  |D3%in% c("Active promoter","Primed enhancer") ]
sub.to.plot.tmp<-as.data.frame(sub.to.plot[,.(freq=.N),.( D0,D1,D2,D3)])# %>%.[freq>=nrow(sub.to.plot.tmp)/1000]) #)#
sub.to.plot.tmp$freq %>% .[.>=1000] %>% hist()
pdf(paste0(inputdir,"/alluvial.plot.enhancer.pdf"),width = 8,height = 4)
p1<-ggplot(sub.to.plot.tmp[sub.to.plot.tmp$freq>=100,],
           aes(y=freq,
               axis1=D0,
               axis2 = D1, 
               axis3 = D2,
               axis4=D3) )+
  geom_alluvium(aes(fill = (D3)), width = 1/2) +
  geom_stratum(aes(fill =  (D3)), width = 1/2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  #geom_label(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  scale_x_discrete(limits = c("D0", "D1","D2","D3"),
                   expand = c(0.15, 0.05)) +
  #  scale_fill_manual(values = c("#e84f35", "#8abf57", "#f5ae94","#b1b1b1", "#8386c0")) +
  theme_void()
print(p1)
dev.off()



