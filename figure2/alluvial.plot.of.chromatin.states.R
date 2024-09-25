
library(alluvial)
library(ggalluvial)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)

################################################################################
####################### alluvial plot of chromatin states ###################### 
################################################################################

options(scipen = 9) # avoid 100000->1e+5
hsGRCh38<-fread("/Users/yang/data/ref/hg38.chrom.sizes.txt")%>%  setnames(c("chr", "pos"))
chrs <- paste0("chr", c(1:22, "X", "Y"))

bin<-200 #5kb
#bin<-500 # 

windows <- hsGRCh38[chr %in% chrs, .(start = seq(1, max(pos), by = bin)), chr] %>%
  .[, end := start + bin-2] %>%
  setkey(chr, start, end)

# as chromhmm bins are connected eg bin1 0-9400, bin2 9400-9800, so end subtract extra 1 to avoid false overlap
windows[,id:=paste0(chr,":",start,"-",end+1)]
setkey(windows,chr,start,end)

inputdir<-"/Users/yang/data/scMTR_method.dev/figure3.data/chromhmm_adj/states_15"

test<-fread("/Users/yang/data/scMTR_method.dev/figure3.data/chromhmm_adj/states_15/POSTERIOR/H9_MT_sc_15_chr14_posterior.txt")
windows[chr=="chr14"]

test[,.(start =(.I * 200 +1-200))]

inputdir<-"/Users/yang/data/scMTR_method.dev/figure3.data/chromhmm_adj/states_15"

files<-list.files(inputdir,pattern = "txt",full.names = T)
fread(files[1])
pdf(paste0(inputdir,"/emission.stats.pdf"),width = 6,height = 10)
pheatmap::pheatmap(as.data.frame(fread(files[1]))%>% tibble::column_to_rownames("State (Emission order)"),
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

# remove posterior probablity less than 0.75 bins 


for (sampl in  c("MT_Bulk","CT")) {
  #sampl<-"MT_sc",

files<-list.files(paste0(inputdir,"/POSTERIOR"),pattern = "txt",full.names = T)
foo<-lapply(files[grep(sampl,files)],function(f)fread(f) %>%
              .[,chr:=basename(f) %>% str_split (.,"_",simplify = T)  %>% .[,ifelse(sampl=="CT",4,5)]] %>% 
              .[,start :=(.I * 200 +1-200)] %>% 
              .[,end :=start +200 -2] %>%
              .[,id:=paste0(chr,":",start,"-",end+1)]
) %>% 
  rbindlist()

sel.bin<-melt(foo[,-(c("chr","start","end"))])
kept.bin<-sel.bin[value>=0.75,unique(id)]

files<-list.files(inputdir,pattern = "dense.bed",full.names = T)
foo<-lapply(files[grep(sampl,files)], function(f) fread(f,select = c(1:4),skip = 1)%>% 
         setnames(c("chr","start","end","states"))%>% 
         .[, anno:=basename(f)%>% str_split(.,"_dense.bed",simplify = T)%>% .[,1]]
       
       )%>% rbindlist()

foo<-foverlaps(foo,windows,nomatch = 0L)

test<-foo[!id %in% kept.bin]
pdf(paste0(inputdir,"/",sampl,".less.75.posterior.probability.bar.pdf"),width=5,height=3)
ggplot(test,aes(states))+ geom_bar() +scale_y_log10()+theme_bw()
dev.off()

foo<-foo[id %in% kept.bin]

fwrite(foo[,.(id,states)],paste0(inputdir,"/",sampl,".kept.75.posterior.probability.tsv"),sep = "\t")
}

files<-list.files(inputdir,pattern = ".tsv.gz",full.names = T)
foo<-lapply(files, function(f) fread(f)%>% 
              .[, anno:=basename(f)%>% str_split(.,".kept.75.posterior.probability.tsv.gz",simplify = T)%>% .[,1]]
            
)%>% rbindlist()

foo[states %in% c(1), state.type:="Repressed Polycomb"]
foo[states %in% c(2), state.type:="Repressed Polycomb"] #Bivalent enhancer
foo[states %in% c(3,4), state.type:="Bivalent promoter"]
foo[states %in% c( 10,12,15), state.type:="Unmarked"]
foo[states %in% c(5:8), state.type:="Active promoter"]
foo[states %in% c(9,11,13), state.type:="Active enhancer"]
foo[states %in% c(14), state.type:="Transcription"]
foo[states==14]
#to.plot<-dcast(foo,chr+start+end+id~anno,value.var = "state.type")
to.plot<-dcast(foo,id~anno,value.var = "state.type")
sub.to.plot<-to.plot[!(CT=="Unmarked" |MT_Bulk=="Unmarked" |MT_sc=="Unmarked")]
sub.to.plot<-to.plot[!(CT=="Unmarked")]

sub.to.plot[is.na(CT)]
sub.to.plot[is.na(MT_Bulk)]
sub.to.plot<-sub.to.plot[!(is.na(MT_sc)|is.na(CT)|is.na(MT_Bulk))]
sub.to.plot<-sub.to.plot[!(is.na(MT_sc)|is.na(CT))]

#sub.to.plot<-sub.to.plot[,.N,.(CT,MT_Bulk,MT_sc)]


foo<-fread(files[1],select = c(1:4),skip = 1 )%>% setnames(c("chr","start","end","states"))

is_alluvia_form((sub.to.plot), axes = 5:7, silent = TRUE)
is_alluvia_form((sub.to.plot), axes = 1:4, silent = TRUE)

sub.to.plot[,.(freq=.N),.(  CT     ,    MT_Bulk      ,     MT_sc)]
pdf(paste0(inputdir,"/alluvial.plot.states.pdf"),width = 8,height = 4)
ggplot(as.data.frame(sub.to.plot[,.(freq=.N),.( CT, MT_Bulk,     MT_sc)]),
       aes(y=freq,
           axis1=CT,
           axis2 = MT_Bulk, 
           axis3 = MT_sc) )+
  geom_alluvium(aes(fill = CT), width = 1/5) +
  geom_stratum(width=1/10) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), width = 1/5) +
  scale_x_discrete(limits = c("CT", "MT_Bulk","MT_sc"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = c("#e84f35", "#8abf57", "#f5ae94","#b1b1b1", "#8386c0")) +
 theme_void()
dev.off()





