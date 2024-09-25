library(RColorBrewer)
library(ggplot2)
library(purrr)
library(data.table)
library(stringr)

## down samples and perform Correlations
## 3500, 2000,1500, 1000,900,800,700,600,500,400,300, 200, 100,50 cells


################################################################################
#####################  downsample different number of cells ####################
################################################################################

in.file.dir <- "/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments"

file.pattern <- "sorted.bed.gz.tbi"
fragment.list<-list.files(in.file.dir,pattern=file.pattern,full.names = T) %>% sub(".tbi","",.)

scRNA<-readRDS("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/RNA/scRNA.data.seurat.rds")
scRNA<-subset(scRNA,subset= sample_bio.id=="Def.endo.D0")
scRNA$CB<-paste0(substring(colnames(scRNA),4,12),substring(colnames(scRNA),1,2)) 
rownames(scRNA@meta.data)<-scRNA$CB
RNAcounts<-scRNA@assays$RNA@counts
colnames(RNAcounts) <-paste0(substring(colnames(RNAcounts),4,12),substring(colnames(RNAcounts),1,2)) 
#7479
frag.path<-fragment.list[1]
kept.cells<-unique(fread(frag.path)%>% .$V4)
#7484
kept.cells<-kept.cells[kept.cells %in%colnames(RNAcounts)]
#7479

rm(scRNA)
rm(RNAcounts)

select.cells<-c()
for (i in 1:5) {
  for (sample.size in c ( 3500, 2000,1500, 1000,900,800,700,600,500,400,300, 200, 100,50)) {
  select.cells<-rbind(select.cells, data.table(cellid=sample(kept.cells, size = sample.size),
                                               rep=paste0("rep",i), 
                                               sample.size=paste0(sample.size,"cells") ) )
  }
}

out.file<-"/Users/yang/data/scMTR_method.dev/figure3.data/downsamples"
dir.create(out.file)
fwrite(select.cells,paste0(out.file,"/downsample.cells.metadata.tsv"),sep="\t")


select.cells <-fread(paste0(out.file,"/downsample.cells.metadata.tsv"))

fread(fragment.list[1])

foo<-lapply(fragment.list,function(i) fread(i) %>% 
         setnames(c("chr","start","end","cellid","num")) %>% 
         .[,ab.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,4] ]) %>% rbindlist()

for (reps in paste0("rep",1:5)) {
  for (sample.siz in unique(select.cells$sample.size)) {

      foo.tmp<-foo[cellid %in% select.cells[rep==reps & sample.size== sample.siz, cellid] , .(chr,start,end,cellid,num,ab.type)] %>% split(by="ab.type",drop=T)
      walk2(foo.tmp,  paste0(out.file,"/",names(foo.tmp),"_",sample.siz,"_", reps,".fragments.tsv"),fwrite, sep="\t",col.names=F )
      system(paste0("gzip ",paste0(out.file,"/",names(foo.tmp),"_",sample.siz,"_", reps,".fragments.tsv")) )
  }
  }

################################################################################
####################  load samples and perform correlations  ###################
################################################################################

out.file<-"/Users/yang/data/scMTR_method.dev/figure3.data/downsamples/bin_cors"
dir.create(out.file)
files1<-list.files("/Users/yang/data/scMTR_method.dev/bulk.cuttag/fragments/bin_counts",pattern = ".gz", recursive = T,full.names = T)
files2<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/bin_counts",pattern = ".gz", recursive = T,full.names = T)
files3<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/bulk/fragments/bin_counts",pattern = ".gz", recursive = T,full.names = T)

for (reps in c(paste0("rep",1:5))) {
files4<-list.files("/Users/yang/data/scMTR_method.dev/figure3.data/downsamples/bin_counts",pattern = reps, recursive = T,full.names = T)


for (binsize in c ( "5000.tsv", "20000.tsv")) {
  #binsize<-"5000.tsv"
  foo1<-lapply(files1[grep(binsize,files1)], function(i) fread(i,select = c(1,2,3)) %>%  ## deduplicated: umi_chr_pos
                 setnames(c("chrom", "bin", "hist")) %>%
                 .[,sample.type:= "CUT&Tag"] %>% #4
                 .[,ab.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,1] ]%>% #3
                 .[,bin.type:= basename(dirname(i))]
  )%>%
    rbindlist()
  
  foo2<-lapply(files2[grep(binsize,files1)], function(i) fread(i,select = c(1,2,3)) %>%  ## deduplicated: umi_chr_pos
                 setnames(c("chrom", "bin", "hist")) %>%
                 .[,sample.type:= "scMTR-seq"] %>% #4
                 .[,ab.type:=basename(i) %>% str_split(.,'[.]',simplify = T)%>% .[,4] ]%>% #3
                 .[,bin.type:= basename(dirname(i))]
  )%>%
    rbindlist()
  
  foo3<-lapply(files3[grep(binsize,files3)], function(i) fread(i,select = c(1,2,3)) %>%  ## deduplicated: umi_chr_pos
                 setnames(c("chrom", "bin", "hist")) %>%
                 .[,sample.type:= "BMT-seq"] %>% #4
                 .[,ab.type:=basename(i) %>% str_split(.,'[.]|_',simplify = T)%>% .[,2] ]%>% #3
                 .[,bin.type:= basename(dirname(i))]
  )%>%
    rbindlist()
  
  foo4<-lapply(files4[grep(binsize,files4)], function(i) fread(i,select = c(1,2,3)) %>%  ## deduplicated: umi_chr_pos
                 setnames(c("chrom", "bin", "hist")) %>%
                 .[,sample.type:= reps] %>% #4
                 .[,ab.type:=basename(i) %>% str_split(.,'.fragments.tsv',simplify = T)%>% .[,1] ]%>% #3
                 .[,bin.type:= basename(dirname(i))]
  )%>%
    rbindlist()
  
  foo<-rbind(foo1,foo2,foo3,foo4)
  #foo[!ab.type =="RNA"]
  foo<-foo[!grep("IgG",ab.type)]
  
  #foo<-foo[!grep("IgG",ab.type)]

  foo[,repd:=.N,.(chrom,bin)]
  to.plot<-foo[,.N,repd]
  
  pdf(paste0(out.file, "/bin.",binsize,"_",reps,"_rep.num.pdf" ) ,width=6,height=4)
  p1<- ggplot(to.plot,aes(repd,N))+
    geom_col()+
    labs(x="# of reproducible bins",y="# of bins")+
    theme_bw()
  print(p1)
  dev.off()
  
  rep.num=5
  used.met="pearson"
  test<-dcast(foo[repd>=rep.num ], chrom+bin~sample.type+ab.type,value.var = "hist")
      
  used.obs<- "pairwise.complete.obs"# "complete.obs"
  data_normalized <- test[,-(c("chrom","bin"))]
  M = cor( data_normalized %>% as.data.frame() , use = used.obs,method = used.met)
  fwrite(melt(M), paste0(out.file, "/bin.",binsize,"_",used.obs,"_cor.",used.met,"_cov.",rep.num,"_",reps,".tsv" ),sep="\t" )
      
  pdf(paste0(out.file, "/bin.",binsize,"_",used.obs,"_cor.",used.met,"_cov.",rep.num,"_",reps,".heatmap.pdf" ) ,width=28,height=28)
  p2<-pheatmap::pheatmap(M ,
                        display_numbers = T,
                        clustering_method = "ward.D2",
                        cellwidth = 20,
                        cellheight = 20,
                        color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100) )
           print(p2)
      dev.off()
      
      
      }
    
}

################################################################################
#####  plot downsample corrleation with bulk/aggreated full sc datasets  #######
################################################################################
files<-list.files(out.file, pattern = ".tsv", full.names = T )

to.plot<-lapply(files[-grep("pdf",files)],function(i)fread(i) %>% 
         setnames(c("sample1","sample2","cors")) %>% 
         .[,binsize:=basename(i) %>% str_split("[.]",simplify = T) %>% .[,2] ]%>% 
         .[,reps:= basename(i) %>% str_split("[.]|_",simplify = T) %>% .[,11] ]
         ) %>% rbindlist()

to.plot[,unique(sample1)]
to.plot[,hist1:=sample1 %>% str_split("_",simplify = T) %>% .[,2]]
to.plot[,hist2:=sample2 %>% str_split("_",simplify = T) %>% .[,2]]

to.plot[,type1:=sample1 %>% str_split("_",simplify = T) %>% .[,1]]
to.plot[,type2:=sample2%>% str_split("_",simplify = T) %>% .[,1]]

to.plot[,cell1:=sample1 %>% str_split("_",simplify = T) %>% .[,3]]
to.plot[,cell2:=sample2 %>% str_split("_",simplify = T) %>% .[,3]]
to.plot[,unique(cell1)]

tmp.to.plot<-to.plot[hist1==hist2 & type1 %in% "CUT&Tag" & !type2%in% c("CUT&Tag","BMT-seq")]
tmp.to.plot[reps==type2,type2:=cell2]


tmp.to.plot[,unique(type2)]
tmp.to.plot[type2=="scMTR-seq",type2:="7479"]
tmp.to.plot[,type2:=sub("cells","",type2)]

tmp.to.plot$type2<-as.numeric(tmp.to.plot$type2)
pdf(paste0(out.file, "/downsamples.vs.CUT&Tag.cor.pointplot.pdf" ) ,width=8,height=8)
ggplot(tmp.to.plot,aes(type2, cors)) +
  geom_point(aes(color=binsize),size=0.2)+
  facet_grid(hist1~.)+
  labs(title = "downsamples vs CUT&Tag")+
  theme_bw()
   dev.off()
   
 tmp.to.plot<-to.plot[hist1==hist2 & type1 %in% "scMTR-seq" & !type2%in% c("CUT&Tag","BMT-seq","scMTR-seq")]
 tmp.to.plot[reps==type2,type2:=cell2]
 
 
 tmp.to.plot[,unique(type2)]
 tmp.to.plot[type2=="scMTR-seq",type2:="7479"]
 tmp.to.plot[,type2:=sub("cells","",type2)]
 
 tmp.to.plot$type2<-as.numeric(tmp.to.plot$type2)
 pdf(paste0(out.file, "/downsamples.vs.scMTR-seq.cor.pointplot.pdf" ) ,width=8,height=8)
 ggplot(tmp.to.plot,aes(type2, cors)) +
   geom_point(aes(color=binsize),size=0.2)+
   facet_grid(hist1~.)+
   labs(title = "downsamples vs scMTR-seq")+
   theme_bw()
   dev.off()
   
