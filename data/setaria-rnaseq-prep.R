library(tm)
library(ggplot2)
library(gridExtra)
library(plyr)
library(shiny)
library(reshape2)
library(reshape)
library(stringr)
library(matrixStats)


setwd("/Users/mgehan/Documents/github/setaria-expression-shiny/data/")

experiment<-read.table(file='experiment.information.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)
saveRDS(experiment,"experiment.information.rds")

leaf<-read.table(file='all_genes_lists/Sv_leaf_annotation.csv',sep=',',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE, quote="", fill=TRUE, na.strings=c(""," ","NA"))
inflor<-read.table(file='all_genes_lists/Sviridis_inflorescence_development_annotation.csv',sep=',',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE, quote="", fill=TRUE, na.strings=c(""," ","NA"))


leaf.sum<-leaf[, lapply(.SD, sum), by=c()]

leaf.sub<-subset(leaf,select=-c(TranscriptID,maize_classic.x, classic_description.x,
                                WikiGene.description,InterPro.description,
                                Best.hit.arabi.name,arabi.symbol,arabi.defline,
                                Best.hit.rice.name,rice.symbol,rice.defline,
                                Maizev4.gene,Zmays_284.gene,Bradi.gene,
                                Sobic.gene,Seita.gene,Pfam,Panther,
                                KOG,KEGG.ec,KO,GO))

leaf.sub.sum<-ddply(leaf.sub,.(GeneID),summarise,S1_1=sum(S1_1),
                    S1_2=sum(S1_2),S1_3=sum(S1_3),S1_4=sum(S1_4),S2_1=sum(S2_1),
                    S2_2=sum(S2_2),S2_3=sum(S2_3),S2_4=sum(S2_4),S3_1=sum(S3_1),
                    S3_2=sum(S3_2),S3_3=sum(S3_3),S3_4=sum(S3_4),S4_1=sum(S4_1),
                    S4_2=sum(S4_2),S4_3=sum(S4_3),S4_4=sum(S4_4))

leaf.sub.sum$means1<-NA
leaf.sub.sum$means2<-NA
leaf.sub.sum$means3<-NA
leaf.sub.sum$means4<-NA

leaf.sub.sum$means1<-rowMeans(subset(leaf.sub.sum, select = c(S1_1, S1_2, S1_3,S1_4)), na.rm = TRUE)
leaf.sub.sum$means2<-rowMeans(subset(leaf.sub.sum, select = c(S2_1, S2_2, S2_3,S2_4)), na.rm = TRUE)
leaf.sub.sum$means3<-rowMeans(subset(leaf.sub.sum, select = c(S3_1, S3_2, S3_3,S3_4)), na.rm = TRUE)
leaf.sub.sum$means4<-rowMeans(subset(leaf.sub.sum, select = c(S4_1, S4_2, S4_3,S4_4)), na.rm = TRUE)

leaf.sub.sum$sd1<-NA
leaf.sub.sum$sd2<-NA
leaf.sub.sum$sd3<-NA
leaf.sub.sum$sd4<-NA

leaf.sub.sum$sd1<-apply(subset(leaf.sub.sum, select = c(S1_1, S1_2, S1_3,S1_4)),1, sd)
leaf.sub.sum$sd2<-apply(subset(leaf.sub.sum, select = c(S2_1, S2_2, S2_3,S2_4)),1, sd)
leaf.sub.sum$sd3<-apply(subset(leaf.sub.sum, select = c(S3_1, S3_2, S3_3,S3_4)),1, sd)
leaf.sub.sum$sd4<-apply(subset(leaf.sub.sum, select = c(S4_1, S4_2, S4_3,S4_4)),1, sd)

leaf.anno<-subset(leaf,select=c(GeneID,maize_classic.x, classic_description.x,
                                WikiGene.description,InterPro.description,
                                Best.hit.arabi.name,arabi.symbol,arabi.defline,
                                Best.hit.rice.name,rice.defline,
                                Maizev4.gene,Zmays_284.gene,Bradi.gene,
                                Sobic.gene,Seita.gene,Pfam,Panther,
                                KOG,KEGG.ec,KO,GO))
leaf.anno.unique<-unique(leaf.anno)
colnames(leaf.anno.unique)<-c("Gene_ID","maize_classic.x", "classic_description.x",
                              "WikiGene.description","InterPro.description",
                              "Best.hit.arabi.name","arabi.symbol","arabi.define",
                              "Best.hit.rice.name","rice.define",
                              "Maizev4.gene","Zmays_284.gene","Bradi.gene",
                              "Sobic.gene","Seita.gene","Pfam","Panther",
                              "KOG","KEGG.ec","KO","GO")

inflor.anno<-subset(inflor,select=c(Gene_ID,maize_classic.x, classic_description.x,
                                  WikiGene.description,InterPro.description,
                                  Best.hit.arabi.name,arabi.symbol,arabi.define,
                                  Best.hit.rice.name,rice.define,
                                  Maizev4.gene,Zmays_284.gene,Bradi.gene,
                                  Sobic.gene,Seita.gene,Pfam,Panther,
                                  KOG,KEGG.ec,KO,GO))
inflor.anno.unique<-unique(inflor.anno)

all.anno<-merge(leaf.anno.unique,inflor.anno.unique,by=c("Gene_ID","maize_classic.x", "classic_description.x",
                                                         "WikiGene.description","InterPro.description",
                                                         "Best.hit.arabi.name","arabi.symbol","arabi.define",
                                                         "Best.hit.rice.name","rice.define",
                                                         "Maizev4.gene","Zmays_284.gene","Bradi.gene",
                                                         "Sobic.gene","Seita.gene","Pfam","Panther",
                                                         "KOG","KEGG.ec","KO","GO"), all=TRUE)

inflor.sub<-subset(inflor,select=-c(maize_classic.x, classic_description.x,
                                   WikiGene.description,InterPro.description,
                                   Best.hit.arabi.name,arabi.symbol,arabi.define,
                                   Best.hit.rice.name,rice.symbol,rice.define,
                                   Maizev4.gene,Zmays_284.gene,Bradi.gene,
                                   Sobic.gene,Seita.gene,Pfam,Panther,
                                   KOG,KEGG.ec,KO,GO,Short_name,Short_description))

leaf.sub.sum1<-subset(leaf.sub.sum,select=c(GeneID,means1,means2,means3,means4,sd1,sd2,sd3,sd4))
datamerge<-merge(leaf.sub.sum1,inflor.sub,by.x="GeneID",by.y="Gene_ID",all=TRUE)
datamerge.empty<-subset(datamerge,select=-c(GeneID))
datamerge.nona<-datamerge[rowSums(is.na(datamerge.empty)) != ncol(datamerge.empty), ]

datamerge.nona$percents1<- NA
datamerge.nona$percents2<- NA
datamerge.nona$percents3<- NA
datamerge.nona$percents4<- NA
datamerge.nona$percent10d<- NA
datamerge.nona$percent12d<- NA
datamerge.nona$percent14d<- NA
datamerge.nona$percent15d<- NA
datamerge.nona$percent16d<- NA
datamerge.nona$percent18d<- NA

datamerge.nona$percents1<- (datamerge.nona$sd1/datamerge.nona$means1)*100
datamerge.nona$percents2<- (datamerge.nona$sd2/datamerge.nona$means2)*100
datamerge.nona$percents3<- (datamerge.nona$sd3/datamerge.nona$means3)*100
datamerge.nona$percents4<- (datamerge.nona$sd4/datamerge.nona$means4)*100
datamerge.nona$percent10d<- (datamerge.nona$X10DAS_SD/datamerge.nona$X10DAS_mean)*100
datamerge.nona$percent12d<- (datamerge.nona$X12DAS_SD/datamerge.nona$X12DAS_mean)*100
datamerge.nona$percent14d<- (datamerge.nona$X14DAS_SD/datamerge.nona$X14DAS_mean)*100
datamerge.nona$percent15d<- (datamerge.nona$X15DAS_SD/datamerge.nona$X15DAS_mean)*100
datamerge.nona$percent16d<- (datamerge.nona$X16DAS_SD/datamerge.nona$X16DAS_mean)*100
datamerge.nona$percent18d<- (datamerge.nona$X18DAS_SD/datamerge.nona$X18DAS_mean)*100

all.anno1<-ddply(all.anno, .(Gene_ID), head, n = 1)
data.anno<-merge(datamerge.nona,all.anno1,by.x="GeneID",by.y="Gene_ID",all.x=TRUE)

saveRDS(data.anno,"setaria-all-data.rds")

plot<-readRDS("../setaria.plot.rds")
ggplot(plot, aes(x=leaf.segment, y=mean, group=factor(GeneID), colour=factor(GeneID))) +
  geom_line()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  geom_point()+
  theme_bw()+
  labs(title="Leaf Segment Expression", y="TPM", x="Leaf Segment")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right")
