
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(purrr)
library(data.table)
library(stringr)
################################################################################# 
######## to generate bin counts, use bash script to generate bin counts #########
################################################################################# 
out.file<-"/Users/yang/data/scMTR_method.dev/figure3.data/bin_cors"

files1<-list.files("/Users/yang/data/scMTR_method.dev/bulk.cuttag/fragments/bin_counts",pattern = ".gz", recursive = T,full.names = T)
files2<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/paired/fragments/fragments/bin_counts",pattern = ".gz", recursive = T,full.names = T)
files3<-list.files("/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/bulk/fragments/bin_counts",pattern = ".gz", recursive = T,full.names = T)
files1[grep("1000.tsv",files1)]
files2[grep("1000.tsv",files2)]
files3[grep("1000.tsv",files3)]

################################################################################# 
######## calculate correlation with different bin size, rep num of bins #########
#################################################################################
for (binsize in c ("200.tsv", "500.tsv","1000.tsv", "5000.tsv", "10000.tsv" ,"20000.tsv", "50000.tsv","100000.tsv", "200000.tsv")) {
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
  foo<-rbind(foo1,foo2,foo3)
  #foo[!ab.type =="RNA"]
  foo<-foo[!grep("IgG",ab.type)]
  foo<-foo[!grep("IgG",ab.type)]
  foo[,repd:=.N,.(chrom,bin)]
  to.plot<-foo[,.N,repd]
  
  pdf(paste0(out.file, "/bin.",binsize,"_rep.num.pdf" ) ,width=6,height=4)
 p1<- ggplot(to.plot,aes(repd,N))+
    geom_col()+
    labs(x="# of reproducible bins",y="# of bins")+
    theme_bw()
 print(p1)
  dev.off()
  
 for (rep.num in c(3:15) ) {
 #   rep.num=5
    for (used.met in c("pearson", "spearman")){
    #used.met="pearson"
    test<-dcast(foo[repd>=rep.num ], chrom+bin~sample.type+ab.type,value.var = "hist")
    
    used.obs<- "pairwise.complete.obs"# "complete.obs"
    # data_normalized <- scale(test[,-(c("chrom","bin"))])
    data_normalized <- test[,-(c("chrom","bin"))]
    M = cor( data_normalized %>% as.data.frame() , use = used.obs,method = used.met)
     pdf(paste0(out.file, "/bin.",binsize,"_",used.obs,"_cor.",used.met,"_cov.",rep.num,".heatmap.pdf" ) ,width=8,height=8)
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
  }

