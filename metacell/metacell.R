### load metacell metadata info

metacell<-fread("/Users/wangy/Desktop/metacell/cell2metacell_assignment.txt.gz") %>% setnames(c("cellid","metacell"))
# merge metacell into fragments file, for seurat object generation
### load fragments

in.file.dir <- "/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments"  

out.file.dir <- (in.file.dir) #"/bi/home/wangy/projects/ctr_seq/08.matrix"
out.file.folder <- "fragments"

if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
  dir.create(paste0(out.file.dir,"/",out.file.folder))
}

files<-list.files(in.file.dir,pattern = ".fragments.tsv.gz",full.names = T)

read.summary<-fread("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/rna.dna.kept.summary.tsv")
#read.summary <-fread("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/clean.all.reads.summary.Def.endo1.tsv.gz")
read.summary<-read.summary[lib.type=="RNA"& deduplicated_reads>=2000]

#"H3K27ac",
for (fi in c("H3K27me3","H3K36me3","H3K4me3","H3K4me1","rbIgG")){
  # fi= "H3K27ac" #H3K27ac|
  fi=files[grep(fi,files)]
  foo<-lapply(fi, function(i) fread(i, fill=T, sep="\t", nThread=8 )%>% 
                setnames(c("chr","start","end","cellid","num")) %>%  ## deduplicated: umi_chr_pos
                .[,sample.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,1] ] %>% #4
                .[,ab.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,3] %>% str_split(.,'[.]',simplify = T)%>% .[,1] ]%>%
                .[,sub.lib:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,2] ] #5
  )%>%
    rbindlist()
  
  foo<-merge(foo,read.summary[,.(cellid=new.cellid,sample.type)],by=c("cellid"))
 # foo[,.N,cellid]
  foo<-merge(foo,metacell,by="cellid")
  #1:   Def.endo.D2 43:36:06:S8 1972853
  #2:   Def.endo.D1 43:36:06:S8    3849
  #3:   Def.endo.D3 43:36:06:S8    3359
  #test<-foo[,.N,.(sample.type.y  ,  metacell)]
  # test<-foo[chr=="chr6"& start>142299881 & end <142304021]
  
  #foo[,metacell:=sub(":",".",metacell)]
  
  foo.tmp<-foo[,.N,.(chr,start,end,metacell,ab.type )] %>% setkey(chr,start,end) %>% split(by=c( "ab.type"),keep.by = F)
  
  walk2(foo.tmp, paste0(out.file.dir,"/",out.file.folder,"/", names(foo.tmp),".metacell.fragments.tsv"), fwrite,sep="\t",col.names=F)
  
}


