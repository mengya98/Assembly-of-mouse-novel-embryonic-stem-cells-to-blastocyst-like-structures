setwd('G:/Mouse_embryo/mblastoid_plus_blastocyst')
library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(harmony)

E_A_E4.75ICM<-read.table('G:/Mouse_embryo/data/E_A_E4.75ICM.xls',row.names = 1,header = T)
E_B_E4.75ICM<-read.table('G:/Mouse_embryo/data/E_B_E4.75ICM.xls',row.names = 1,header = T)
E_B_E4.75TE<-read.table('G:/Mouse_embryo/data/E_B_E4.75TE.xls',row.names = 1,header = T)
E_C_E4.75ICM<-read.table('G:/Mouse_embryo/data/E_C_E4.75ICM.xls',row.names = 1,header = T)
E_D_E4.5<-read.table('G:/Mouse_embryo/data/E_D_E4.5.xls',row.names = 1,header = T)
EM14<-read.table('G:/Mouse_embryo/data/EM_RNA_14.umi_counts.xls',row.names = 1,header = T)
EM15<-read.table('G:/Mouse_embryo/data/EM_RNA_15.umi_counts.xls',row.names = 1,header = T)
EM16<-read.table('G:/Mouse_embryo/data/EM_RNA_16.umi_counts.xls',row.names = 1,header = T)
EM17<-read.table('G:/Mouse_embryo/data/EM_RNA_17.umi_counts.xls',row.names = 1,header = T)

EM1<-read.table('G:/Mouse_embryo/data/EM_RNA_1.umi_counts.xls',row.names = 1,header = T)
EM2<-read.table('G:/Mouse_embryo/data/EM_RNA_2.umi_counts.xls',row.names = 1,header = T)
EM3<-read.table('G:/Mouse_embryo/data/EM_RNA_3.umi_counts.xls',row.names = 1,header = T)
EM4<-read.table('G:/Mouse_embryo/data/EM_RNA_4.umi_counts.xls',row.names = 1,header = T)
EM5<-read.table('G:/Mouse_embryo/data/EM_RNA_5.umi_counts.xls',row.names = 1,header = T)
EM6<-read.table('G:/Mouse_embryo/data/EM_RNA_6.umi_counts.xls',row.names = 1,header = T)
EM7<-read.table('G:/Mouse_embryo/data/EM_RNA_7.umi_counts.xls',row.names = 1,header = T)
EM9<-read.table('G:/Mouse_embryo/data/EM_RNA_9.umi_counts.xls',row.names = 1,header = T)
EM10<-read.table('G:/Mouse_embryo/data/EM_RNA_10.umi_counts.xls',row.names = 1,header = T)
EM12<-read.table('G:/Mouse_embryo/data/EM_RNA_12.umi_counts.xls',row.names = 1,header = T)
EM13<-read.table('G:/Mouse_embryo/data/EM_RNA_13.umi_counts.xls',row.names = 1,header = T)

em<- CreateSeuratObject(counts = cbind(EM14,EM15,EM16,EM17,E_A_E4.75ICM,E_B_E4.75ICM,E_B_E4.75TE,E_C_E4.75ICM,E_D_E4.5,
                                       EM1,EM2,EM3,EM4,EM5,EM6,EM7,EM9,EM10,EM12,EM13), project = "em", min.cells = 5) #%>%
em@meta.data$Sample <- c(rep("EM14", ncol(EM14)),rep("EM15", ncol(EM15)),rep("EM16", ncol(EM16)),rep("EM17", ncol(EM17)),
                         rep('E_A_E4.75ICM',ncol(E_A_E4.75ICM)),rep('E_B_E4.75ICM',ncol(E_B_E4.75ICM)),rep('E_B_E4.75TE',ncol(E_B_E4.75TE)),rep('E_C_E4.75ICM',ncol(E_C_E4.75ICM)),
                         rep('E_D_E4.5',ncol(E_D_E4.5)),
                         rep("EM1", ncol(EM1)), rep("EM2", ncol(EM2)),rep("EM3", ncol(EM3)),rep("EM4", ncol(EM4)),
                         rep("EM5", ncol(EM5)),rep("EM6", ncol(EM6)),rep("EM7", ncol(EM7)),rep("EM9", ncol(EM9)),
                         rep("EM10", ncol(EM10)),rep("EM12", ncol(EM12)),rep("EM13", ncol(EM13)))
em@meta.data$Origin <- c(rep("Normal", ncol(EM14)),rep("Normal", ncol(EM15)),rep("Normal", ncol(EM16)),rep("Normal", ncol(EM17)),
                         rep('Normal',ncol(E_A_E4.75ICM)),rep('Normal',ncol(E_B_E4.75ICM)),rep('Normal',ncol(E_B_E4.75TE)),rep('Normal',ncol(E_C_E4.75ICM)),
                         rep('Normal',ncol(E_D_E4.5)),
                         rep("Positive", ncol(EM1)), rep("Positive", ncol(EM2)),rep("Positive", ncol(EM3)),rep("Positive", ncol(EM4)),
                         rep("Positive", ncol(EM5)),rep("Positive", ncol(EM6)),rep("Positive", ncol(EM7)),rep("Positive", ncol(EM9)),
                         rep("Positive", ncol(EM10)),rep("Positive", ncol(EM12)),rep("Positive", ncol(EM13)))

em[["percent.mito"]] <- PercentageFeatureSet(object = em, pattern = "^mt.")

pdf("em.FeatureScatter.pdf",width = 24,height = 14)
VlnPlot(object = em, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 2)
plot1 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "percent.mito") 
plot2 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
CombinePlots(plots = list(plot1,plot2))
dev.off()

em <- subset(x = em, subset = nFeature_RNA > 2500 & nCount_RNA < 750000 &percent.mito < 5 )#nCount_RNA < 600000
pdf("em.afterFeatureScatter.pdf",width = 24,height = 14)
VlnPlot(object = em, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 2)
plot1 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "percent.mito") 
plot2 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
CombinePlots(plots = list(plot1,plot2))
dev.off()

em <- NormalizeData(object = em, normalization.method = "LogNormalize", scale.factor = 1e6)
em <- FindVariableFeatures(object = em,selection.method = 'vst', nfeatures = 2000)
em.genes <- rownames(x = em)
em <- ScaleData(object = em, features = em.genes)
em <- RunPCA(object = em, features = VariableFeatures(object = em))
pdf("em.sample(去batch前).pdf",width = 20,height = 8)
#options(repr.plot.height =5, repr.plot.width = 15)
plot1 <- DimPlot(object = em, reduction = "pca", pt.size = 1, group.by = "Sample")
plot2 <- VlnPlot(object = em, features = "PC_1", group.by = "Sample", pt.size = 1)
CombinePlots(plots = list(plot1, plot2))
dev.off()

em <- JackStraw(object = em, num.replicate = 1000,dims = 30) 
em <- ScoreJackStraw(object = em, dims = 1:30)
pdf("em.pc.pdf",width = 10,height = 6)
JackStrawPlot(object = em, dims = 1:30)
ElbowPlot(object = em,ndims = 30)
dev.off()


#harmony
em <- em %>%
  RunHarmony("Sample", plot_convergence = TRUE)
pdf("em.sample(去batch后,harmony).pdf",width = 20,height = 8)
plot1 <- DimPlot(object = em, reduction = "harmony", pt.size = 1, group.by = "Sample")
plot2 <- VlnPlot(object = em, features = "harmony_1", group.by = "Sample", pt.size = 1)
CombinePlots(plots = list(plot1, plot2))
dev.off()


em <- FindNeighbors(object = em, reduction="harmony",dims = 1:11)
em <- FindClusters(object = em, resolution = 0.3)
em <- RunUMAP(object = em,reduction="harmony",dims = 1:11)

pdf("./em.sample.pdf",width = 5.5,height = 4)
DimPlot(em, reduction = "umap", group.by = "Sample", pt.size = 0.5)
dev.off()
pdf("./em.origin.pdf",width = 4.9,height = 4)
DimPlot(em, reduction = "umap", group.by = "Origin", pt.size = 0.5)
dev.off()

###文章图###
em$cell_type<-''
em$cell_type[em$seurat_clusters%in%c(5)]='bl-like_TE'
em$cell_type[em$seurat_clusters%in%c(1)]='bl-like_ICM'
em$cell_type[em$seurat_clusters%in%c(0,4)]='bl-like_PE'

em_normal<-read.table('G:/Mouse_embryo/blastocyst/em.normal.meta.xls',header=T)
meta<-data.frame(em@meta.data)
meta_em_normal<-data.frame(em@meta.data)
meta_em_normal$identity<-rownames(meta_em_normal)
em_normal<-em_normal[em_normal$identity%in%c(meta_em_normal$identity),]
em_normal<-em_normal[,c(1,9)]
meta_em_normal<-merge(meta_em_normal,em_normal,by='identity',all=T)
meta_em_normal$cell_type.y[meta_em_normal$cell_type.y%in%c('Morula')]<-'TE'
meta_em_normal[760,10]<-''
meta_em_normal$cell_type.y[meta_em_normal$cell_type.y%in%c('ICM')]<-'Control_ICM'
meta_em_normal$cell_type.y[meta_em_normal$cell_type.y%in%c('PE')]<-'Control_PE'
meta_em_normal$cell_type.y[meta_em_normal$cell_type.y%in%c('TE')]<-'Control_TE'

for (i in 1:nrow(meta_em_normal)) {
  if (meta_em_normal[i,10]=='') {
    meta_em_normal[i,12]<-meta_em_normal[i,11]
  }
  if (meta_em_normal[i,10]!='') {
    meta_em_normal[i,12]<-meta_em_normal[i,10]
  }
}
colnames(meta_em_normal)[12]<-'cell_type'
rownames(meta_em_normal)<-meta_em_normal$identity
meta_em_normal<-meta_em_normal[,c(2:9,12)]
meta_em_normal<-meta_em_normal[colnames(em),]
em@meta.data<-meta_em_normal

cols1<-c("mediumpurple1","darkorchid1",'lightblue',"gold1","darkorange","darkorange4")#
em$cell_type<-factor(em$cell_type,levels=c('Control_ICM','Control_PE','Control_TE','bl-like_ICM','bl-like_PE','bl-like_TE'))
pdf('./em_celltype.pdf',width=5.1,height=4)
DimPlot(em,reduction = "umap", pt.size = 0.5, group.by = 'cell_type',label = F,cols = cols1)+
  theme_classic()+
  labs(x = "UMAP1", y = "UMAP2",title = '')
dev.off()

te_markers_intersect2<-read.table('G:/Mouse_embryo/mblastoid/TE_markers_intersect.xls',header=T)

#VLNPLOT#
all_markers<-c("Pou5f1","Gata4","Cdx2","Nanog","Sox17","Eomes","Sox2",'Pdgfra',"Elf5",'Klf2','Foxa2','Gata2')
pdf("./em.StackedVlnPlot.pdf", width = 6, height = 10)
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = 'cell_type',... )  + 
    xlab("") + ylab(feature) + ggtitle("") + theme_classic()+
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_line(color ='black'),
          axis.text.y = element_text(size = rel(1),color ='black'), 
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5,color ='black'),
          plot.margin = plot.margin , title = element_blank())
  return(p)
}
## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(color ='black',angle = 45,hjust = 1,vjust = 1), axis.ticks.x = element_line(color ='black'))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 3)
  return(p)
}
StackedVlnPlot(em, features = c(all_markers,te_markers_intersect2$gene),pt.size=0,cols=cols1)
dev.off()

#featureplot#
for (i in 1:length(all_markers)) {
  pdf(paste0('./featureplot/',all_markers[i],'.pdf'),width=4.5,height=4)
  p<-FeaturePlot(em, reduction = "umap", all_markers[i],pt.size = 0.5)+
    theme_classic()+
    labs(x = "UMAP1", y = "UMAP2",title = all_markers[i])+
    scale_colour_distiller(palette='Reds',direction = 1)+
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}

gene_features <-  unique(c("Pou5f1","Nanog","Sox2",'Klf2','Zfp42','Esrrb','Tdgf1','Gdf3','Utf1','Fgf4','Tbx3',
                           "Gata4","Gata6","Sox17",'Pdgfra','Col4a1','Dab2','Foxa2','Grb10',
                           "Cdx2","Eomes","Elf5",'Gata2','Gata3','Tead3','Krt8','Tfap2c','Krt18','Tacstd2','Tead4','Dppa1',
                           te_markers_intersect2$gene))
pdf("./heatmap_featuregenes.pdf",width = 7, height = 6)
DoHeatmap(object = em, features = gene_features,group.bar=T,angle=30,size=2.5,group.by = 'cell_type',
          group.colors =cols1)+
  theme(panel.border = element_blank())
dev.off()
