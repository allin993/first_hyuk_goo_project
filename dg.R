library(DESeq2)
BiocManager::install("airway")
library(airway)
data(gse)# airway가 제공하는 데이터 이름 
gse# rna seq airway smooth muscle control/ dexamethasone 4 samples
assay(gse)
head(assay(gse),3)
gse
colData(gse)
rowRanges(gse)
gse$cell<-gse$donor
gse$dex<-gse$condition
colData(gse)


library(DESeq2)
dds<-DESeqDataSet(gse, design= ~cell +dex)
dds
dds<-DESeq(dds)# filtering 필요없는 유전자 제거 (8sample 0 or 1 , )
assay(dds)[,1]

rowSums(assay(dds))==0

assay(dds)[rowSums(assay(dds))>=1,]
dds<-dds[rowSums(assay(dds))>1,]
#variance stablizing transformation
#발현량이 많을 수록 variance가 크다. variance를 동등하게 만들어주는 함수가 필요
vsd<-vst(dds) #clustereing을 하기 위해서 사용 
plotPCA(vsd) #batch  effect 고려하기 위해 cell로 보정. +앞에 보정 사항 적어 놓음.


dds<-DESeq(dds)
res<-results(dds)
res

#basemean: 평균 발현량, log2FC: 대조군 대비 실험군의 변화, p-value: 얼마나 유의

res.05<-results(dds,alpha=0.05)
resLFC1<-results(dds, lfcThreshold = 1)
summary(resLFC1)
res[res$padj<=0.1,]

#volcano plot
df<-data.frame(subset(res,!is.na(padj)))
df$sig<-NA
df$sig<-replace(df$sig, df$padj<0.01& df$log2FoldChange>0,"up")

# 행의 평균값으로 heatmap! 

install.packages("pheatmap")
topVarGenes<-head(order(rowVars(assay(vsd)),decreasing = T),20)
mat<-assay(vsd)[topVarGenes,]
mat<-mat-rowMeans(mat)
anno<-as.data.frmaee(colData(vsd))









