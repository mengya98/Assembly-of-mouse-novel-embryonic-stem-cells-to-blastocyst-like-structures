setwd("G:/Mouse_embryo/blastocyst")
library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(patchwork)

E_A_E4.75ICM<-read.table('G:/Mouse_embryo/data/E_A_E4.75ICM.xls',row.names = 1,header = T)
E_B_E4.75ICM<-read.table('G:/Mouse_embryo/data/E_B_E4.75ICM.xls',row.names = 1,header = T)
E_B_E4.75TE<-read.table('G:/Mouse_embryo/data/E_B_E4.75TE.xls',row.names = 1,header = T)
E_C_E4.75ICM<-read.table('G:/Mouse_embryo/data/E_C_E4.75ICM.xls',row.names = 1,header = T)
E_C_latemorula<-read.table('G:/Mouse_embryo/data/E_C_latemorula.xls',row.names = 1,header = T)
E_D_E4.5<-read.table('G:/Mouse_embryo/data/E_D_E4.5.xls',row.names = 1,header = T)
E_E_earlymorula<-read.table('G:/Mouse_embryo/data/E_E_earlymorula.xls',row.names = 1,header = T)

EM14<-read.table('G:/Mouse_embryo/data/EM_RNA_14.umi_counts.xls',row.names = 1,header = T)
EM15<-read.table('G:/Mouse_embryo/data/EM_RNA_15.umi_counts.xls',row.names = 1,header = T)
EM16<-read.table('G:/Mouse_embryo/data/EM_RNA_16.umi_counts.xls',row.names = 1,header = T)
EM17<-read.table('G:/Mouse_embryo/data/EM_RNA_17.umi_counts.xls',row.names = 1,header = T)

em<- CreateSeuratObject(counts = cbind(EM14,EM15,EM16,EM17,E_A_E4.75ICM,E_B_E4.75ICM,E_B_E4.75TE,E_C_E4.75ICM,E_C_latemorula,E_D_E4.5,E_E_earlymorula), project = "em", min.cells = 5) #%>%

em@meta.data$batch <- c(rep("EM14", ncol(EM14)),rep("EM15", ncol(EM15)),rep("EM16", ncol(EM16)),rep("EM17", ncol(EM17)),
                        rep('E_A_E4.75ICM',ncol(E_A_E4.75ICM)),rep('E_B_E4.75ICM',ncol(E_B_E4.75ICM)),rep('E_B_E4.75TE',ncol(E_B_E4.75TE)),rep('E_C_E4.75ICM',ncol(E_C_E4.75ICM)),
                        rep('E_C_latemorula',ncol(E_C_latemorula)),rep('E_D_E4.5',ncol(E_D_E4.5)),rep('E_E_earlymorula',ncol(E_E_earlymorula)))

em[["percent.mito"]] <- PercentageFeatureSet(object = em, pattern = "^mt.")

pdf("em.FeatureScatter.pdf",width = 24,height = 14)
VlnPlot(object = em, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 2)
plot1 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "percent.mito") 
plot2 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
CombinePlots(plots = list(plot1,plot2))
dev.off()

em <- subset(x = em, subset = nFeature_RNA > 2500 & nCount_RNA < 750000 &percent.mito < 5 )#nCount_RNA < 600000
em
em <- NormalizeData(object = em, normalization.method = "LogNormalize", scale.factor = 1e6)
em <- FindVariableFeatures(object = em,selection.method = 'vst', nfeatures = 2000)
em.genes <- rownames(x = em)
em <- ScaleData(object = em, features = em.genes)
em <- RunPCA(object = em, features = VariableFeatures(object = em))

pdf("em.batch(去batch前).pdf",width = 20,height = 8)
#options(repr.plot.height =5, repr.plot.width = 15)
plot1 <- DimPlot(object = em, reduction = "pca", pt.size = 1, group.by = "batch")
plot2 <- VlnPlot(object = em, features = "PC_1", group.by = "batch", pt.size = 1)
CombinePlots(plots = list(plot1, plot2))
dev.off()

em <- JackStraw(object = em, num.replicate = 1000,dims = 20) 
em <- ScoreJackStraw(object = em, dims = 1:20)
pdf("em.pc.pdf",width = 10,height = 6)
JackStrawPlot(object = em, dims = 1:20)
ElbowPlot(object = em,ndims = 20)
dev.off()

#harmony
em <- em %>%
  RunHarmony("batch", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(em, 'harmony')
pdf("em.batch(去batch后,harmony).pdf",width = 20,height = 8)
plot1 <- DimPlot(object = em, reduction = "harmony", pt.size = 1, group.by = "batch")
plot2 <- VlnPlot(object = em, features = "harmony_1", group.by = "batch", pt.size = 1)
CombinePlots(plots = list(plot1, plot2))
dev.off()

em <- FindNeighbors(object = em, dims = 1:12)
em <- FindClusters(object = em, resolution = 0.3)
em <- RunTSNE(object = em,dims = 1:12)

em <- FindNeighbors(object = em, reduction="harmony",dims = 1:12)
em <- FindClusters(object = em, resolution = 0.3)
em <- RunUMAP(object = em,reduction="harmony",dims = 1:12)

em$cell_type<-''
em$cell_type[em$seurat_clusters%in%c(0,2,3)]='TE'
em$cell_type[em$seurat_clusters%in%c(4)]='ICM'
em$cell_type[em$seurat_clusters%in%c(5)]='PE'
em$cell_type[em$seurat_clusters%in%c(1,6)]='Morula'

em$cell_type<-factor(em$cell_type,levels=c('ICM','PE','TE','Morula'))
pdf('em_celltype.pdf',width=4.8,height=4)
DimPlot(em,reduction = "umap", pt.size = 1, group.by = 'cell_type',label = F)+
  theme_classic()+
  labs(x = "UMAP1", y = "UMAP2",title = '')
dev.off()

all_markers<-c("Pou5f1","Gata4","Cdx2","Nanog","Sox17","Eomes","Sox2",'Pdgfra',"Elf5",'Klf2','Foxa2','Gata2')
#vlnplot#
pdf("em.StackedVlnPlot.pdf", width = 6, height = 4)
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = 'cell_type',... )  + 
    xlab("") + ylab(feature) + ggtitle("") + theme_classic()+#coord_flip()+
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
    theme(axis.text.x=element_text(color ='black'), axis.ticks.x = element_line(color ='black'))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 3)
  return(p)
}
StackedVlnPlot(em, features = all_markers,pt.size=0)
dev.off()

for (i in 1:length(all_markers)) {
  pdf(paste0('./featureplot/',all_markers[i],'.pdf'),width=4.5,height=4)
  p<-FeaturePlot(em, reduction = "umap", all_markers[i],pt.size = 1)+
    theme_classic()+
    #theme(axis.text= element_blank(), axis.ticks=element_blank())+
    labs(x = "UMAP1", y = "UMAP2",title = all_markers[i])+
    scale_colour_distiller(palette='Reds',direction = 1)+
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}

#heatmap#
te_markers<-read.table('G:/Mouse_embryo/TE Marker.txt')
gene_features <-  unique(c("Pou5f1","Nanog","Sox2",'Klf2','Zfp42','Esrrb','Tdgf1','Gdf3','Utf1','Fgf4','Tbx3',
                           "Gata4","Gata6","Sox17",'Pdgfra','Col4a1','Dab2','Foxa2','Grb10',
                           "Cdx2","Eomes","Elf5",'Gata2','Gata3','Tead3','Krt8','Tfap2c','Krt18','Tacstd2','Tead4','Dppa1',
                           te_markers$V1))
pdf("heatmap_featuregenes.pdf",width = 6, height = 10)
DoHeatmap(object = em, features = gene_features,group.bar=T,angle=0,size=3,group.by = 'cell_type')+
   theme(panel.border = element_blank())
dev.off()
