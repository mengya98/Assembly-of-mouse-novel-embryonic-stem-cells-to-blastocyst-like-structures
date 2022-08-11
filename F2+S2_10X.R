setwd('G:/Mouse_embryo/10X/')
library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)

em.data<-Read10X(data.dir = "./rawdata/raw_feature_bc_matrix/")
em<-CreateSeuratObject(counts = em.data,Project="em",min.cells = 3,min.features = 200)
em

em[["percent.mito"]] <- PercentageFeatureSet(object = em, pattern = "^mt-")
pdf("em.FeatureScatter.pdf",width = 20,height = 10)
VlnPlot(object = em, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0.5)
plot1 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "percent.mito") 
plot2 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
CombinePlots(plots = list(plot1,plot2))
dev.off()

em <- subset(x = em, subset = nFeature_RNA > 2500 & nCount_RNA < 40000 & percent.mito < 5 )#
em
pdf("em.FeatureScatter.afterQC.pdf",width = 20,height = 10)
VlnPlot(object = em, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0.5)
plot1 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "percent.mito") 
plot2 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
CombinePlots(plots = list(plot1,plot2))
dev.off()

em <- NormalizeData(object = em, normalization.method = "LogNormalize", scale.factor = 1e6)
em <- FindVariableFeatures(object = em,selection.method = 'vst', nfeatures = 2000)
em.genes <- rownames(x = em)
em <- ScaleData(object = em, features = em.genes)
em <- RunPCA(object = em, features = VariableFeatures(object = em))

em <- JackStraw(object = em, num.replicate = 1000,dims = 30) 
em <- ScoreJackStraw(object = em, dims = 1:30)
pdf("em.pc.pdf",width = 10,height = 6)
JackStrawPlot(object = em, dims = 1:30)
ElbowPlot(object = em,ndims = 30)
dev.off()

em <- FindNeighbors(object = em,dims = 1:10)
em <- FindClusters(object = em, resolution = 0.4)
em <- RunUMAP(object = em,dims = 1:10)


###celltype umap###
em$cell_type<-''
em$cell_type[em$seurat_clusters%in%c(4,5,6)]='ICM'
em$cell_type[em$seurat_clusters%in%c(7,8)]='Intermediate'
em$cell_type[em$seurat_clusters%in%c(0,1,2,3)]='PE'

em$cell_type<-factor(em$cell_type,levels=c('ICM','PE','Intermediate'))
pdf('./em_celltype.pdf',width=5.1,height=4)
DimPlot(em,reduction = "umap", pt.size = 0.5, group.by = 'cell_type',label = F)+
  theme_classic()+
  labs(x = "UMAP1", y = "UMAP2",title = '')#+
dev.off()

te_markers_intersect<-read.table('G:/Mouse_embryo/te_markers_intersect.xls',header=T)
gene_features <-  unique(c("Pou5f1","Nanog","Sox2",'Klf2','Zfp42','Esrrb','Tdgf1','Gdf3','Utf1','Fgf4','Tbx3',
                           "Gata4","Gata6","Sox17",'Pdgfra','Col4a1','Dab2','Foxa2','Grb10',
                           "Cdx2","Eomes",'Gata3','Tead3',te_markers_intersect$gene))
pdf("./heatmap_featuregenes.pdf",width = 6, height = 6)
#pdf("heatmap_normal.pdf",width = 6, height = 12)
DoHeatmap(object = em, features = gene_features,group.bar=T,angle=0,size=3,group.by = 'cell_type',hjust = 0.5)+
  #scale_color_manual(name='Cluster',values=cols1,limits=c("0,Neutrophil-1","1,Neutrophil-2","2,Neutrophil-3","3,Neutrophil-4","4,B cell","5,Monocyte-1","6,Monocyte-2","7,Erythroidal","8,DC","9,Neutrophil-5","10,B cell","11,Lymphocyte(NK&T)","12,Basophil","13,B cell"))+
  theme(panel.border = element_blank())
dev.off()

all_markers<-c("Pou5f1","Gata4","Cdx2","Nanog","Sox17","Eomes","Sox2",'Pdgfra',"Elf5",'Klf2','Foxa2','Gata2','Id2','Lgals1','S100a11','Tmsb10')
#vlnplot#
pdf("./em.StackedVlnPlot.pdf", width = 6, height = 4)
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
  p<-FeaturePlot(em, reduction = "umap", all_markers[i],pt.size = 0.5)+
    theme_classic()+
    #theme(axis.text= element_blank(), axis.ticks=element_blank())+
    labs(x = "UMAP1", y = "UMAP2",title = all_markers[i])+
    scale_colour_distiller(palette='Reds',direction = 1)+
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}

##DEG##
Idents(em)<-em$cell_type
em.markers <- FindAllMarkers(em, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
em.markers<-em.markers[em.markers$p_val_adj <= 0.05 & em.markers$pct.1 >= 0.5 ,]
em.markers<-em.markers[order(em.markers$avg_log2FC,em.markers$cluster,decreasing = T),]
em.markers<-em.markers[order(em.markers$cluster,decreasing = T),]
write.table(em.markers,file = "./em.markers.xls", row.names = F, col.names = T, quote = F, sep = "\t")

