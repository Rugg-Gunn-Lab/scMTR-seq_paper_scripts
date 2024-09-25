
library(data.table)
library(purrr)
library(ggplot2)
library(viridis)
library(stringr)

## 1. as for different sample we have different RT, so we can separate rna cellid based on sub:BC1:BC2:BC3:RT1 or sub:BC1:BC2:BC3:RT2, even sub:BC1:BC2:BC3 is the same;
## 2. to match DNA data, we can only use sub:BC1:BC2:BC3 as cell id;
## 3. we use this RT information to estimate the collision rate of cellids with sub:BC1:BC2:BC3

################################################################################
############################ collision plot ####################################
################################################################################
in.file.dir <- "/Users/yang/data/scMTR_method.dev/Sample_5927_SLX-23063_051623YW_MTR_Def.endo.DNA/RNA"  

files<-list.files(paste0(in.file.dir ,"/reads.summary"),recursive = T,pattern="xls",full.names = T)#recursive = T,

foo<-lapply(files, function(i) fread(i,select = c(1,2,3,5)) %>%  ## deduplicated: umi_chr_pos
              setnames(c("cellid", "all_reads", "mapped_reads", "deduplicated_reads")) %>%
              .[,c("mapping_rate","duplication_rate"):=list(mapped_reads/all_reads, 1-deduplicated_reads/mapped_reads)] %>% 
              #  .[,sub.lib:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,3] ] %>%
              .[,sample.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,4] ] %>% #4
              .[,lib.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,5] ]%>% #3
              .[,mapping.type:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,9] ]%>%
              .[,sub.lib:=basename(i) %>% str_split(.,'_',simplify = T)%>% .[,6] ] #5
)%>%
  rbindlist()

foo<-foo[deduplicated_reads>0]

foo[lib.type=="R",lib.type:="RNA"]

foo[lib.type=="RNA"& substring(cellid,10,11) %in% c("01"),sample.type:="Def.endo.D0"]
foo[lib.type=="RNA"& substring(cellid,10,11) %in% c("02") ,sample.type:="Def.endo.D1"]
foo[lib.type=="RNA"& substring(cellid,10,11) %in% c("03"),sample.type:="Def.endo.D2"]
foo[lib.type=="RNA"& substring(cellid,10,11) %in% c("04"),sample.type:="Def.endo.D3"]

foo[,new.cellid:=paste0(substring(cellid,1,8),":",sub.lib)]

kept.cells<-foo[deduplicated_reads>=2000,unique(new.cellid)]

foo<-foo[new.cellid %in% kept.cells& !sample.type=="HIST"]

to.plot<-dcast(foo,new.cellid~sample.type,value.var ="deduplicated_reads" )
to.plot[Def.endo.D0<2000 &Def.endo.D1<2000,type.D0.D1:="low"]
to.plot[(Def.endo.D0>=2000 |Def.endo.D1>=2000)& Def.endo.D0/(Def.endo.D0+Def.endo.D1) > 0.8,type.D0.D1:="Specific1"]
to.plot[(Def.endo.D0>=2000 |Def.endo.D1>=2000)& Def.endo.D1/(Def.endo.D0+Def.endo.D1) > 0.8,type.D0.D1:="Specific2"]
to.plot[(Def.endo.D0>=2000 |Def.endo.D1>=2000)& is.na(type.D0.D1),type.D0.D1:="collision"]

to.label<-to.plot[,.N,type.D0.D1]
dev.off()
out.file<-"/Users/yang/data/scMTR_method.dev/figure4.data"

pdf( paste0(out.file,"/collision.plot.rna.type.D0.D1.pdf"),width=5,height=4)
ggplot(to.plot,aes(Def.endo.D0 ,Def.endo.D1))+
  geom_point(size=0.1,aes(color=type.D0.D1))+
  scale_color_manual(values =c(low="white",Specific1="#E6332A",Specific2="#365FAA",collision="grey60") )+
  labs(title  = paste0("Xs=",to.label[type.D0.D1=="Specific1",N], ", Ys=",to.label[type.D0.D1=="Specific2",N], ", collision cells=",to.label[type.D0.D1=="collision",N] ))+
  theme_bw()
dev.off()

to.plot[Def.endo.D0<2000 &Def.endo.D2<2000,type.D0.D2:="low"]
to.plot[(Def.endo.D0>=2000 |Def.endo.D2>=2000)& Def.endo.D0/(Def.endo.D0+Def.endo.D2) > 0.8,type.D0.D2:="Specific1"]
to.plot[(Def.endo.D0>=2000 |Def.endo.D2>=2000)& Def.endo.D2/(Def.endo.D0+Def.endo.D2) > 0.8,type.D0.D2:="Specific2"]
to.plot[(Def.endo.D0>=2000 |Def.endo.D2>=2000)& is.na(type.D0.D2),type.D0.D2:="collision"]

to.label<-to.plot[,.N,type.D0.D2]

pdf( paste0(out.file,"/collision.plot.rna.type.D0.D2.pdf"),width=5,height=4)

ggplot(to.plot,aes(Def.endo.D0 ,Def.endo.D2))+
  geom_point(size=0.1,aes(color=type.D0.D2))+
  scale_color_manual(values =c(low="white",Specific1="#E6332A",Specific2="#365FAA",collision="grey60") )+
  labs(title  = paste0("Xs=",to.label[type.D0.D2=="Specific1",N], ", Ys=",to.label[type.D0.D2=="Specific2",N], ", collision cells=",to.label[type.D0.D2=="collision",N] ))+
  theme_bw()
dev.off()

to.plot[Def.endo.D0<2000 &Def.endo.D3<2000,type.D0.D3:="low"]
to.plot[(Def.endo.D0>=2000 |Def.endo.D3>=2000)& Def.endo.D0/(Def.endo.D0+Def.endo.D3) > 0.8,type.D0.D3:="Specific1"]
to.plot[(Def.endo.D0>=2000 |Def.endo.D3>=2000)& Def.endo.D3/(Def.endo.D0+Def.endo.D3) > 0.8,type.D0.D3:="Specific2"]
to.plot[(Def.endo.D0>=2000 |Def.endo.D3>=2000)& is.na(type.D0.D3),type.D0.D3:="collision"]

to.label<-to.plot[,.N,type.D0.D3]
pdf( paste0(out.file,"/collision.plot.rna.type.D0.D3.pdf"),width=5,height=4)
ggplot(to.plot,aes(Def.endo.D0 ,Def.endo.D3))+
  geom_point(size=0.1,aes(color=type.D0.D3))+
  scale_color_manual(values =c(low="white",Specific1="#E6332A",Specific2="#365FAA",collision="grey60") )+
  labs(title  = paste0("Xs=",to.label[type.D0.D3=="Specific1",N], ", Ys=",to.label[type.D0.D3=="Specific2",N], ", collision cells=",to.label[type.D0.D3=="collision",N] ))+
  theme_bw()
dev.off()


to.plot[Def.endo.D1<2000 &Def.endo.D3<2000,type.D1.D3:="low"]
to.plot[(Def.endo.D1>=2000 |Def.endo.D3>=2000)& Def.endo.D1/(Def.endo.D1+Def.endo.D3) > 0.8,type.D1.D3:="Specific1"]
to.plot[(Def.endo.D1>=2000 |Def.endo.D3>=2000)& Def.endo.D3/(Def.endo.D1+Def.endo.D3) > 0.8,type.D1.D3:="Specific2"]
to.plot[(Def.endo.D1>=2000 |Def.endo.D3>=2000)& is.na(type.D1.D3),type.D1.D3:="collision"]

to.label<-to.plot[,.N,type.D1.D3]
pdf( paste0(out.file,"/collision.plot.rna.type.D1.D3.pdf"),width=5,height=4)

ggplot(to.plot,aes(Def.endo.D1 ,Def.endo.D3))+
  geom_point(size=0.1,aes(color=type.D1.D3))+
  scale_color_manual(values =c(low="white",Specific1="#E6332A",Specific2="#365FAA",collision="grey60") )+
  labs(title  = paste0("Xs=",to.label[type.D1.D3=="Specific1",N], ", Ys=",to.label[type.D1.D3=="Specific2",N], ", collision cells=",to.label[type.D1.D3=="collision",N] ))+
  theme_bw()
dev.off()

to.plot[Def.endo.D1<2000 &Def.endo.D2<2000,type.D1.D2:="low"]
to.plot[(Def.endo.D1>=2000 |Def.endo.D2>=2000)& Def.endo.D1/(Def.endo.D1+Def.endo.D2) > 0.8,type.D1.D2:="Specific1"]
to.plot[(Def.endo.D1>=2000 |Def.endo.D2>=2000)& Def.endo.D2/(Def.endo.D1+Def.endo.D2) > 0.8,type.D1.D2:="Specific2"]
to.plot[(Def.endo.D1>=2000 |Def.endo.D2>=2000)& is.na(type.D1.D2),type.D1.D2:="collision"]

to.label<-to.plot[,.N,type.D1.D2]
pdf( paste0(out.file,"/collision.plot.rna.type.D1.D2.pdf"),width=5,height=4)

ggplot(to.plot,aes(Def.endo.D1 ,Def.endo.D2))+
  geom_point(size=0.1,aes(color=type.D1.D2))+
  scale_color_manual(values =c(low="white",Specific1="#E6332A",Specific2="#365FAA",collision="grey60") )+
  labs(title  = paste0("Xs=",to.label[type.D1.D2=="Specific1",N], ", Ys=",to.label[type.D1.D2=="Specific2",N], ", collision cells=",to.label[type.D1.D2=="collision",N] ))+
  theme_bw()

dev.off()

to.plot[Def.endo.D2<2000 &Def.endo.D3<2000,type.D2.D3:="low"]
to.plot[(Def.endo.D2>=2000 |Def.endo.D3>=2000)& Def.endo.D2/(Def.endo.D2+Def.endo.D3) > 0.8,type.D2.D3:="Specific1"]
to.plot[(Def.endo.D2>=2000 |Def.endo.D3>=2000)& Def.endo.D3/(Def.endo.D2+Def.endo.D3) > 0.8,type.D2.D3:="Specific2"]
to.plot[(Def.endo.D2>=2000 |Def.endo.D3>=2000)& is.na(type.D2.D3),type.D2.D3:="collision"]


to.label<-to.plot[,.N,type.D2.D3]
pdf( paste0(out.file,"/collision.plot.rna.type.D2.D3.pdf"),width=5,height=4)

ggplot(to.plot,aes(Def.endo.D2 ,Def.endo.D3))+
  geom_point(size=0.1,aes(color=type.D2.D3))+
  scale_color_manual(values =c(low="white",Specific1="#E6332A",Specific2="#365FAA",collision="grey60") )+
  labs(title  = paste0("Xs=",to.label[type.D2.D3=="Specific1",N], ", Ys=",to.label[type.D2.D3=="Specific2",N], ", collision cells=",to.label[type.D2.D3=="collision",N] ))+
  theme_bw()
dev.off()
