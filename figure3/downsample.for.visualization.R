library(data.table)
library(purrr)
library(stringr)
library(Seurat)

################################################################################
################## downsample for genome browser visualization ##################
################################################################################
read.summary<-fread("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/rna.dna.kept.summary.tsv")
sample1<-read.summary[lib.type=="RNA",sample(new.cellid,1000),sample.type]
sample2<-read.summary[lib.type=="RNA",sample(new.cellid,1000),sample.type]
sample3<-read.summary[lib.type=="RNA",sample(new.cellid,1000),sample.type]

sample4<-read.summary[lib.type=="RNA",sample(new.cellid,2000),sample.type]
sample5<-read.summary[lib.type=="RNA",sample(new.cellid,5000),sample.type]



sample1<-read.summary[new.cellid%in% sample1$V1]
sample2<-read.summary[new.cellid%in% sample2$V1]
sample3<-read.summary[new.cellid%in% sample3$V1]

sample4<-read.summary[new.cellid%in% sample4$V1]
sample5<-read.summary[new.cellid%in% sample5$V1]

fwrite(sample1, "/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/downsample.sample1.tsv",sep="\t")
fwrite(sample2, "/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/downsample.sample2.tsv",sep="\t")
fwrite(sample3, "/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/downsample.sample3.tsv",sep="\t")

fwrite(sample4, "/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/downsample.sample4.tsv",sep="\t")
fwrite(sample5, "/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/downsample.sample5.tsv",sep="\t")


sample1[lib.type=="RNA",.N,sub.lib]

sub.samples<-rbind(sample4[,sample:="sample4"],sample5[,sample:="sample5"])

files<-list.files("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/RNA",pattern = "rmdup.reads.split.bc.tsv.gz",full.names = T)

for (fi in files[c(1:8)]){
  #  fi=files[7]
  foo<-lapply(fi, function(i) fread(i, fill=T, sep="\t", nThread=8 )%>% 
                setnames(c("chr","start","end","umi","cellid","bc4")) %>%  ## deduplicated: umi_chr_pos
                .[,sample.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,5] ] %>% #4
                .[,lib.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,5] ]%>% #3
                .[,mapping.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,9] ]%>%
                .[,sub.lib:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,6] ] #5
  )%>%
    rbindlist()
  
  sublib<-basename(fi) %>% str_split(.,'_',simplify = T)%>% .[,6]

  foo<-foo[chr %in% paste0("chr",c(1:22,"X","Y","M"))]
  foo[,new.cellid:=paste0(cellid,":",sub.lib) ]
  foo[,cellid:=paste0(cellid,":0",bc4) ]

  foo<-merge(foo,sub.samples[lib.type=="RNA",.(cellid,new.cellid,sample,sample.type,ab.type)],by=c("new.cellid","cellid"))
  
  foo.tmp<-foo[,.N,.(chr,start,end,umi,new.cellid,cellid, sample,sample.type=sample.type.y, ab.type )]#%>% split(by=c( "ab.type"),keep.by = F)
  
  fwrite(foo.tmp,paste0("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/RNA/downsample.",sublib,".RNA.fragments.tsv"), sep="\t",col.names=T)
  #walk2(foo.tmp, paste0(out.file.dir,"/",out.file.folder,"/",sample.name,"_", names(foo.tmp),".fragments.tsv"), fwrite,sep="\t",col.names=F)
}

# merge RNA and split by downsample 
files<-list.files("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/RNA",pattern = "downsample",full.names = T)
fread(files[1])

  foo<-lapply(files, function(i) fread(i, fill=T, sep="\t", nThread=8 )
  )%>%
    rbindlist()
  
  foo.tmp<-foo[,.N,.(chr,start,end,umi,new.cellid,cellid, sample,sample.type, ab.type )]%>% split(by=c( "sample",'sample.type',"ab.type"),keep.by = F)
  walk2(foo.tmp, paste0("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/RNA/Down", names(foo.tmp),".fragments.tsv"), fwrite,sep="\t",col.names=F)
  

  ### DNA
  files<-list.files("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments",pattern = "HIST_",full.names = T)
  
  fread(files[1])
  foo<-lapply(files, function(i) fread(i, fill=T, sep="\t", nThread=8 )%>% 
                setnames(c("chr","start","end","new.cellid","N")) %>%  ## deduplicated: umi_chr_pos
                .[,ab.type:=basename(i) %>% str_split(.,'_|[.]',simplify = T)%>% .[,3] ]%>% #3
                .[,sub.lib:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,2] ] #5
  )%>%
    rbindlist()
  setkey(foo,chr,start,end)
  
  foo<-merge(foo,sub.samples[lib.type=="DNA",.(new.cellid,sample,sample.type,ab.type)],by=c("new.cellid","ab.type"))
  foo.tmp<-foo[,.N,.(chr,start,end,new.cellid, sample,sample.type, ab.type )]%>% split(by=c( "sample",'sample.type',"ab.type"),keep.by = F)
  walk2(foo.tmp, paste0("/Users/wangy/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/RNA/Down", names(foo.tmp),".fragments.tsv"), fwrite,sep="\t",col.names=F)
  
  rm(foo.tmp)
    