library(data.table)
library(ggplot2)
files<-list.files("/Users/yang/data/Sample_5769_SLX-22093_220722YW_BMT3",pattern = "mtx",full.names = T)
tmp<-merge(fread(files[1],select = c(1:2)) %>% setnames(c("gene", "RT1")), 
      fread(files[2],select = c(1:2)) %>% setnames(c("gene", "RT2")),by=c("gene"),all=T) %>% 
  merge(fread(files[3],select = c(1:2)) %>% setnames(c("gene", "RT3")), by=c("gene"),all=T)%>% 
  merge(fread(files[4],select = c(1:2)) %>% setnames(c("gene", "RT4")), by=c("gene"),all=T) %>% melt()


tmp[is.na(value),value:=0]
tmp[,rpm:=value/sum(value)*1000000,variable]


#BiocManager::install("PerformanceAnalytics")
library("PerformanceAnalytics")

tmp[,log2rpm:=log2(rpm+1)]
chart.Correlation(dcast(tmp,gene~variable,value.var = "rpm") %>% as.data.frame()%>% tibble::column_to_rownames("gene"), histogram=TRUE, pch=19)
chart.Correlation(dcast(tmp,gene~variable,value.var = "log2rpm") %>% as.data.frame()%>% tibble::column_to_rownames("gene"), histogram=TRUE, pch=19)

tmp[,cv:=sd(log2rpm)/mean(log2rpm),gene]
tmp[,means:=mean(log2rpm),gene]
tmp[variable=="RT1"]
tmp[substring(gene,1,3)=="MT-"]
tmp[,sum(value),substring(gene,1,3)=="MT-"]
tmp[log2rpm>15]
tmp[gene%in%c("CD24","PHC1","TCF4","USP44","DUSP6","ZIC2")] %>% dcast(.,gene~variable,value.var = "log2rpm")
tmp[gene%in%c("POU5F1","SOX2","NANOG","DPPA5","KLF4","DPPA3","TFCP2L1","DEPTOR","CBFA2T2")] %>% dcast(.,gene~variable,value.var = "log2rpm")

tmp[rpm>0 &rpm <=1,.N,variable]
to.plot<-rbind(tmp[rpm>1 &rpm <=10,.N,variable] %>% .[,group:="1<rpm<=10"] %>% melt(),
tmp[rpm>10 &rpm <=100,.N,variable]%>% .[,group:="10<rpm<=100"]%>% melt(),
tmp[rpm>100 &rpm <=1000,.N,variable]%>% .[,group:="100<rpm<=1000"]%>% melt(),
tmp[rpm>1000 ,.N,variable]%>% .[,group:="rpm>1000"]%>% melt())# %>% dcast(.,group~variable,value.var = "value")

pdf("/Users/yang/data/Sample_5769_SLX-22093_220722YW_BMT3/number.of.genes.RNA.bulk.pdf",width=4,height=6)
ggplot(to.plot,aes(x=variable,y=value))+
  geom_bar(aes(fill=variable),position = "dodge",stat = "summary")+
  facet_grid(group~.,scales = "free_y")+
  labs(y="Number of genes",x="")+
  theme_bw()+ theme(
    #    strip.background = element_blank(),
    strip.text = element_text(size=12),
    axis.text.x = element_text(size=12),
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


tmp[cv>1 &means>1]
tmp[gene=="UNC45B"]
tmp[gene=="SOX2"]
tmp[gene=="POU5F1"]
tmp[gene=="NANOG"]
tmp[gene=="GAPDH"]
tmp[gene=="CDX2"]
tmp[gene=="GATA3"]
tmp[gene=="TEAD1"]
tmp[gene=="TEAD2"]
tmp[gene=="TEAD3"]
tmp[gene=="TEAD4"]

fread(files[4],select = c(1:2)) %>% setnames(c("gene", "RT4"))

fread(files[2])
fread(files[3])
