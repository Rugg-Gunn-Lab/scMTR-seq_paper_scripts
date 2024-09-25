library(data.table)
library(purrr)
library(ggplot2)
library(viridis)
library(stringr)
#library(tidyverse)
library(Seurat)

################################################################################  
########################## clean DNA fragments #################################  
################################################################################  

in.file.dir <- "/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments"  
genomes<- "mm10"
file.pattern <-  paste0(genomes,"_sorted_rmdup.reads.split.bc.uniq.tsv.gz")

#output files:
out.file.dir <- dirname(in.file.dir) #"/bi/home/wangy/projects/ctr_seq/08.matrix"
out.file.folder <- "fragments"

if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
  dir.create(paste0(out.file.dir,"/",out.file.folder))
}

files<-list.files(in.file.dir,pattern = file.pattern,full.names = T)
for (fi in files[1:4]){
#  fi=files[1]
foo<-lapply(fi, function(i) fread(i, fill=T, sep="\t", nThread=8 )%>% 
              setnames(c("chr","start","end","umi","cellid","bc4")) %>%  ## deduplicated: umi_chr_pos
              .[,sample.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,1] ] %>% #4
              .[,lib.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,2] ]%>% #3
              .[,mapping.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,4] ]%>%
              .[,sub.lib:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,3] ] #5
)%>%
  rbindlist()

foo<-foo[chr %in% paste0("chr",c(1:19,"X","Y"))]

foo$umi<-NULL
foo<-unique(foo)

foo[bc4 %in% c("1"), ab.type:="H3K4me1"]
foo[bc4 %in% c("2"), ab.type:="H3K4me3"]
foo[bc4 %in% c("3"), ab.type:="H3K27ac"]
foo[bc4 %in% c("4"), ab.type:="H3K27me3"]
foo[bc4 %in% c("5"), ab.type:="H3K9me3"]

foo[bc4 %in% c("6"), ab.type:="H3K36me3"]
foo[bc4 %in% c("7"), ab.type:="IgG1"]
foo[bc4 %in% c("8"), ab.type:="IgG2"]
foo<-foo[!is.na(ab.type)]

foo[,new.cellid:=paste0(cellid,":",sub.lib) ]
# rm no ab sample
tmp<-foo[!is.na(ab.type),.N,.(cellid,sub.lib,sample.type,ab.type)]%>% setorder(N)
tmp[,new.cellid:=paste0(cellid,":",sub.lib) ]
tmp[,all.reads.per.cell:=sum(N),new.cellid]
tmp<-dcast(tmp,new.cellid+all.reads.per.cell~ab.type,value.var = "N")
cutoff=500
tmp[,cutoff:=cutoff]
fwrite(tmp,paste0(in.file.dir,"/",sample.name,".cell.meta.info.allcell.tsv"),sep="\t" )

foo<-foo[new.cellid %in% tmp[all.reads.per.cell>=cutoff,new.cellid]]
foo.tmp<-foo[,.N,.(chr,start,end,new.cellid ,  ab.type )]%>% split(by=c( "ab.type"),keep.by = F)

walk2(foo.tmp, paste0(out.file.dir,"/",out.file.folder,"/",sample.name,"_", names(foo.tmp),".fragments.tsv"), fwrite,sep="\t",col.names=F)
}

################################################################################  
################### Annotate lineages from matched RNA data ####################  
################################################################################  


scRNA2<-readRDS( "/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/scRNA.data.seurat.processed.rds")
DimPlot(scRNA2, reduction = "umap", label = TRUE)
exp.celltype<-scRNA2@meta.data%>% as.data.frame() %>% tibble::rownames_to_column("cellid") %>% as.data.table
exp.celltype[,.N,RNA_snn_res.0.1]
exp.celltype[RNA_snn_res.0.1 %in% c(0),celltype:="TE"]
exp.celltype[RNA_snn_res.0.1==1,celltype:="PE"]
exp.celltype[RNA_snn_res.0.1==2,celltype:="EPI"]
exp.celltype[RNA_snn_res.0.1==3,celltype:="unknown"]
exp.celltype[,new.cellid:=paste0(substring(cellid,4,12),substring(cellid,1,2))]

read.summary<-fread( "/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/reads.summary/rna.dna.kept.summary.tsv.gz")
read.summary[mapping.type=="mm10"]

read.summary<-exp.celltype[new.cellid %in% read.summary[mapping.type=="mm10",new.cellid],.(new.cellid,celltype,sample_bio.id )]
read.summary[,.N,.(celltype,sample_bio.id )] %>% dcast(.,celltype~sample_bio.id ,value.var = "N")

in.file.dir <- "/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments"  

out.file.dir <- (in.file.dir) #"/bi/home/wangy/projects/ctr_seq/08.matrix"
out.file.folder <- "fragments"

if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
  dir.create(paste0(out.file.dir,"/",out.file.folder))
}

files<-list.files(in.file.dir,pattern = ".fragments.tsv",full.names = T)
options(scipen=999) 
for (fi in c("H3K27ac","H3K27me3","H3K36me3","H3K4me3","H3K4me1","H3K9me3")){ #,"IgG1","IgG2"
  # fi= "H3K27ac" #H3K27ac|
  fi=files[grep(fi,files[c(1:24)])]
  foo<-lapply(fi, function(i) fread(i, fill=T, sep="\t", nThread=8 )%>% 
                setnames(c("chr","start","end","cellid","num")) %>%  ## deduplicated: umi_chr_pos
                .[,sample.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,1] ] %>% #4
                .[,ab.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,3] %>% str_split(.,'[.]',simplify = T)%>% .[,1] ]%>%
                .[,sub.lib:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,2] ] #5
  )%>%
    rbindlist()

  foo<-merge(foo,read.summary[,.(cellid=new.cellid,celltype)],by=c("cellid"))
  ids<- foo[,paste0(cellid,"_", ab.type,".",celltype)] %>% unique()
  length(ids)
  foo[,cellid:=paste0(cellid,"_", ab.type,".",celltype)]
  setkey(foo,chr,start,end)
  foo.merged<- lapply(ids, function(i) bedtoolsr::bt.merge(foo[cellid == i,.(chr,start,end,cellid)])%>% 
                        as.data.table() %>% .[,cellid:=i] ) %>% rbindlist() %>% setnames(c("chr","start","end","cellid")) %>% 
    .[,ab.type:= str_split(cellid,'_|[.]',simplify = T)%>% .[,2] ] %>% 
    .[,celltype:= str_split(cellid,'[.]',simplify = T)%>% .[,2] ] %>% 
    .[,cellid:= str_split(cellid,'_',simplify = T)%>% .[,1] ] 
  foo.merged[,start:=start-1000]
  foo.merged[,end:=end+1000]
  foo.tmp<-foo.merged[,.N,.(chr,start,end,cellid , celltype,  ab.type )]%>% split(by=c("celltype", "ab.type"),keep.by = F)
  
#  walk2(foo.tmp, paste0(out.file.dir,"/",out.file.folder,"/", names(foo.tmp),".fragments.tsv"), fwrite,sep="\t",col.names=F)
  walk2(foo.tmp, paste0(out.file.dir,"/",out.file.folder,"/", names(foo.tmp),".ext2000.fragments.tsv"), fwrite,sep="\t",col.names=F)
  
}


############### merge different lineages 
in.file.dir <- "/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments"  
#in.file.dir <- "/Users/yang/data/scMTR_method.dev/Sample_5995_231117YW_MTR_E45/paired/fragments"  

out.file.dir <- (in.file.dir) #"/bi/home/wangy/projects/ctr_seq/08.matrix"
out.file.folder <- "fragments"

if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
  dir.create(paste0(out.file.dir,"/",out.file.folder))
}

files<-list.files(in.file.dir,pattern = ".fragments.tsv.gz",full.names = T)

#read.summary<-fread("/Users/yang/data/scMTR_method.dev/Sample_5995_231117YW_MTR_E45/paired/reads.summary/rna.dna.reads.kept.summary.tsv")
#read.summary <-fread("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/reads.summary/clean.all.reads.summary.Def.endo1.tsv.gz")
options(scipen=999) 

files<-files[c(1:24)]

for (fi in c("H3K27ac","H3K27me3","H3K36me3","H3K4me3","H3K4me1","H3K9me3")){
  # fi= "H3K4me1" #H3K27ac|
  foo<-lapply(files[grep(fi,files)], function (i) fread(i, fill=T, sep="\t", nThread=8 )%>% 
                setnames(c("chr","start","end","cellid","num")) %>%  ## deduplicated: umi_chr_pos
                .[,sample.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,1] ] %>% #4
                .[,ab.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,2]  ] )%>% 
    rbindlist() %>% setkey(chr,start,end)
  fwrite(foo, paste0(out.file.dir,"/",out.file.folder,"/merged.", fi,".fragments.tsv"),sep="\t",col.names = F)
}



################################################################################  
####################### Downsample for visualization  ##########################  
################################################################################  

################################### DNA data ###################################    

###down sampling TE to 947 cells
in.file.dir <- "/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments"  

out.file.dir <- (in.file.dir) #"/bi/home/wangy/projects/ctr_seq/08.matrix"
out.file.folder <- "downsample"

if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
  dir.create(paste0(out.file.dir,"/",out.file.folder))
}

files<-list.files(in.file.dir,pattern = ".fragments.tsv.gz",full.names = T)
options(scipen=999) 
fi= "H3K27ac" #H3K27ac|
fi=files[grep(fi,files)]
foo<-lapply(fi[c(1,3,4)], function(i) fread(i, fill=T, sep="\t", nThread=8 )%>% 
              setnames(c("chr","start","end","cellid","num")) %>%  ## deduplicated: umi_chr_pos
              .[,sample.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,1] ] %>% #4
              .[,ab.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,2]  ]
            
)%>%
  rbindlist()

foo[,.N,.(cellid,sample.type)] %>% .[,.N,sample.type]
id1<-foo[sample.type=="TE",unique(cellid)] %>% sample(.,947)
id2<-foo[sample.type=="TE",unique(cellid)] %>% sample(.,947)
id3<-foo[sample.type=="TE",unique(cellid)] %>% sample(.,947)

files<-list.files("/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments/downsample",pattern="TE",full.names = T)
id1<-fread(files[grep("downsample1",files)][1]) %>% .[,unique(V4)]
id2<-fread(files[grep("downsample2",files)][1]) %>% .[,unique(V4)]
id3<-fread(files[grep("downsample3",files)][1]) %>% .[,unique(V4)]

dir.create(paste0(out.file.dir,"/",out.file.folder,"/downsample1"))
dir.create(paste0(out.file.dir,"/",out.file.folder,"/downsample2"))
dir.create(paste0(out.file.dir,"/",out.file.folder,"/downsample3"))
files<-list.files(in.file.dir,pattern = ".fragments.tsv.gz",full.names = T)
files<-files[grep("TE",files)]

for (fi in c("H3K27ac","H3K27me3","H3K36me3","H3K4me3","H3K4me1","H3K9me3","IgG1","IgG2")){
  # fi= "H3K4me1" #H3K27ac|
  i=files[grep(fi,files)] 
  foo<-fread(i, fill=T, sep="\t", nThread=8 )%>% 
                setnames(c("chr","start","end","cellid","num")) %>%  ## deduplicated: umi_chr_pos
                .[,sample.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,1] ] %>% #4
                .[,ab.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,2]  ]
  setkey(foo,chr,start,end)
  foo[cellid%in% id1,]
  
  fwrite(foo[cellid%in% id1,c(1:5)], paste0(out.file.dir,"/",out.file.folder,"/downsample1","/TE.", fi,".downsample.fragments.tsv"),sep="\t",col.names = F)
  fwrite(foo[cellid%in% id2,c(1:5)], paste0(out.file.dir,"/",out.file.folder,"/downsample2","/TE.", fi,".downsample.fragments.tsv"),sep="\t",col.names = F)
   fwrite(foo[cellid%in% id3,c(1:5)], paste0(out.file.dir,"/",out.file.folder,"/downsample3","/TE.", fi,".downsample.fragments.tsv"),sep="\t",col.names = F)
 
}

###############down sampling PE to 947 cells
in.file.dir <- "/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments"  
#in.file.dir <- "/Users/yang/data/scMTR_method.dev/Sample_5995_231117YW_MTR_E45/paired/fragments"  

out.file.dir <- (in.file.dir) #"/bi/home/wangy/projects/ctr_seq/08.matrix"
out.file.folder <- "downsample"

if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
  dir.create(paste0(out.file.dir,"/",out.file.folder))
}

files<-list.files(in.file.dir,pattern = ".fragments.tsv.gz",full.names = T)
options(scipen=999) 
fi= "H3K27ac" #H3K27ac|
fi=files[grep(fi,files)]
foo<-lapply(fi[c(3)], function(i) fread(i, fill=T, sep="\t", nThread=8 )%>% 
              setnames(c("chr","start","end","cellid","num")) %>%  ## deduplicated: umi_chr_pos
              .[,sample.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,1] ] %>% #4
              .[,ab.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,2]  ]
            
)%>%
  rbindlist()

foo[,.N,.(cellid,sample.type)] %>% .[,.N,sample.type]
id1<-foo[sample.type=="PE",unique(cellid)] %>% sample(.,947)
id2<-foo[sample.type=="PE",unique(cellid)] %>% sample(.,947)
id3<-foo[sample.type=="PE",unique(cellid)] %>% sample(.,947)

files<-list.files("/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments/downsample",pattern="PE",full.names = T)
id1<-fread(files[grep("downsample1",files)][1]) %>% .[,unique(V4)]
id2<-fread(files[grep("downsample2",files)][1]) %>% .[,unique(V4)]
id3<-fread(files[grep("downsample3",files)][1]) %>% .[,unique(V4)]

files<-list.files(in.file.dir,pattern = ".fragments.tsv.gz",full.names = T)
files<-files[grep("PE",files)]

for (fi in c("H3K27ac","H3K27me3","H3K36me3","H3K4me3","H3K4me1","H3K9me3","IgG1","IgG2")){
  # fi= "H3K4me1" #H3K27ac|
  i=files[grep(fi,files)] 
  foo<-fread(i, fill=T, sep="\t", nThread=8 )%>% 
    setnames(c("chr","start","end","cellid","num")) %>%  ## deduplicated: umi_chr_pos
    .[,sample.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,1] ] %>% #4
    .[,ab.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,2]  ]
  setkey(foo,chr,start,end)
  
  foo[cellid%in% id1,]
  fwrite(foo[cellid%in% id1,c(1:5)], paste0(out.file.dir,"/",out.file.folder,"/downsample1","/PE.", fi,".downsample.fragments.tsv"),sep="\t",col.names = F)
  fwrite(foo[cellid%in% id2,c(1:5)], paste0(out.file.dir,"/",out.file.folder,"/downsample2","/PE.", fi,".downsample.fragments.tsv"),sep="\t",col.names = F)
  fwrite(foo[cellid%in% id3,c(1:5)], paste0(out.file.dir,"/",out.file.folder,"/downsample3","/PE.", fi,".downsample.fragments.tsv"),sep="\t",col.names = F)
 
}
############## rename epi 
in.file.dir <- "/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments"  
#in.file.dir <- "/Users/yang/data/scMTR_method.dev/Sample_5995_231117YW_MTR_E45/paired/fragments"  

out.file.dir <- (in.file.dir) #"/bi/home/wangy/projects/ctr_seq/08.matrix"
out.file.folder <- "downsample"

if (!file.exists(paste0(out.file.dir,"/",out.file.folder))){
  dir.create(paste0(out.file.dir,"/",out.file.folder))
}

files<-list.files(in.file.dir,pattern = ".fragments.tsv.gz",full.names = T)
files<-files[grep("EPI",files)]

for (fi in c("H3K27ac","H3K27me3","H3K36me3","H3K4me3","H3K4me1","H3K9me3","IgG1","IgG2")){
  # fi= "H3K4me1" #H3K27ac|
  i=files[grep(fi,files)] 
  foo<-fread(i, fill=T, sep="\t", nThread=8 )%>% 
    setnames(c("chr","start","end","cellid","num")) %>%  ## deduplicated: umi_chr_pos
    .[,sample.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,1] ] %>% #4
    .[,ab.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,2]  ]
  setkey(foo,chr,start,end)
  fwrite(foo, paste0(out.file.dir,"/",out.file.folder,"/downsample1","/EPI.", fi,".downsample.fragments.tsv"),sep="\t",col.names = F)
  fwrite(foo, paste0(out.file.dir,"/",out.file.folder,"/downsample2","/EPI.", fi,".downsample.fragments.tsv"),sep="\t",col.names = F)
  fwrite(foo, paste0(out.file.dir,"/",out.file.folder,"/downsample3","/EPI.", fi,".downsample.fragments.tsv"),sep="\t",col.names = F)
  
}


################################### RNA data ###################################    
in.file.dir <- "/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments"  
files<-list.files(in.file.dir,pattern="_mm10_sorted_rmdup.reads.split.bc.tsv.gz",full.names = T)
fread(files[1])
foo<-lapply(files, function(i) fread(i, fill=T, sep="\t", nThread=8 )%>% 
              setnames(c("chr","start","end","umi","cellid","num")) %>%  ## deduplicated: umi_chr_pos
              .[,sub.lib:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,3] ]
            
)%>%
  rbindlist()
foo<-foo[chr%in% paste0("chr",c(1:19,"X","Y"))]

foo[,cellid:=paste0(cellid,":",sub.lib)]
# RT only keep 1,2 for E45
foo<-foo[num%in% c(1,2)]

for (cells in c("EPI","PE","TE")){
  files<-list.files("/Users/yang/data/scMTR_method.dev/240501YW_MTR_E45_rep/fragments/fragments/fragments/downsample",pattern=cells,full.names = T,recursive = T)
  id1<-fread(files[grep("downsample1",files)][1]) %>% .[,unique(V4)]
  id2<-fread(files[grep("downsample2",files)][1]) %>% .[,unique(V4)]
  id3<-fread(files[grep("downsample3",files)][1]) %>% .[,unique(V4)]
  
  fwrite(foo[cellid%in% id1,c(1:5)], paste0(out.file.dir,"/",out.file.folder,"/downsample1/",cells,".downsample.RNA.fragments.tsv"),sep="\t",col.names = F)
  fwrite(foo[cellid%in% id2,c(1:5)], paste0(out.file.dir,"/",out.file.folder,"/downsample2/",cells,".downsample.RNA.fragments.tsv"),sep="\t",col.names = F)
  fwrite(foo[cellid%in% id3,c(1:5)], paste0(out.file.dir,"/",out.file.folder,"/downsample3/",cells,".downsample.RNA.fragments.tsv"),sep="\t",col.names = F)
  
}




