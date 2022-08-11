setwd('G:/Mouse_embryo/mblastoid')
library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(harmony)


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

em<- CreateSeuratObject(counts = cbind(EM1,EM2,EM3,EM4,EM5,EM6,EM7,EM9,EM10,EM12,EM13), project = "em", min.cells = 5) #%>%
em@meta.data$sample <- c(rep("EM1", ncol(EM1)), rep("EM2", ncol(EM2)),rep("EM3", ncol(EM3)),rep("EM4", ncol(EM4)),
                         rep("EM5", ncol(EM5)),rep("EM6", ncol(EM6)),rep("EM7", ncol(EM7)),rep("EM9", ncol(EM9)),
                         rep("EM10", ncol(EM10)),rep("EM12", ncol(EM12)),rep("EM13", ncol(EM13)))

em[["percent.mito"]] <- PercentageFeatureSet(object = em, pattern = "^mt.")

pdf("em.FeatureScatter.pdf",width = 24,height = 14)
VlnPlot(object = em, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 2)
plot1 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "percent.mito") 
plot2 <- FeatureScatter(object = em, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
CombinePlots(plots = list(plot1,plot2))
dev.off()

em <- subset(x = em, subset = nFeature_RNA > 2500 & nCount_RNA < 750000 &percent.mito < 5 )#nCount_RNA < 600000
em <- NormalizeData(object = em, normalization.method = "LogNormalize", scale.factor = 1e6)
em <- FindVariableFeatures(object = em,selection.method = 'vst', nfeatures = 2000)
em.genes <- rownames(x = em)
em <- ScaleData(object = em, features = em.genes)
em <- RunPCA(object = em, features = VariableFeatures(object = em))
pdf("em.sample(去batch前).pdf",width = 20,height = 8)
#options(repr.plot.height =5, repr.plot.width = 15)
plot1 <- DimPlot(object = em, reduction = "pca", pt.size = 1, group.by = "sample")
plot2 <- VlnPlot(object = em, features = "PC_1", group.by = "sample", pt.size = 1)
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
  RunHarmony("sample", plot_convergence = TRUE)
pdf("em.sample(去batch后,harmony).pdf",width = 20,height = 8)
plot1 <- DimPlot(object = em, reduction = "harmony", pt.size = 1, group.by = "sample")
plot2 <- VlnPlot(object = em, features = "harmony_1", group.by = "sample", pt.size = 1)
CombinePlots(plots = list(plot1, plot2))
dev.off()

em <- FindNeighbors(object = em, reduction="harmony",dims = 1:9)
em <- FindClusters(object = em, resolution = 1.2)
em <- RunUMAP(object = em,reduction="harmony",dims = 1:9)

pdf("./em.sample.pdf",width =4.7,height = 4)
DimPlot(em, reduction = "umap", group.by = "sample", pt.size = 1)
dev.off()

all_markers<-c("Pou5f1","Gata4","Cdx2","Nanog","Sox17","Eomes","Sox2",'Pdgfra','Klf2','Foxa2','Tead3','Tead2','Msx1','Pik3r3')
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

###文章图###
em$cell_type<-''
em$cell_type[em$seurat_clusters%in%c(2,4,5,8,9)]='ICM'
em$cell_type[em$seurat_clusters%in%c(6)]='TE'
em$cell_type[em$seurat_clusters%in%c(0,1,3,7,10)]='PE'

em$cell_type<-factor(em$cell_type,levels=c('ICM','PE','TE'))
pdf('./em_celltype.pdf',width=4.7,height=4)
DimPlot(em,reduction = "umap", pt.size = 1, group.by = 'cell_type',label = F)+
  theme_classic()+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  #theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_text(color = 'black'))
dev.off()

#heatmap#
te_markers<-read.table('G:/Mouse_embryo/TE Marker.txt')
te_markers2<-read.table('G:/Mouse_embryo/blastocyst/em.markers.xls',header=T)
te_markers2<-te_markers2[te_markers2$cluster=='TE',]
te_markers<-unique(c(te_markers$V1,te_markers2$gene))
Idents(em)<-em$cell_type
em.markers <- FindAllMarkers(em, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)##test.use = "roc",

te_markers_intersect<-em.markers[em.markers$gene %in% te_markers,]
te_markers_intersect<-te_markers_intersect[te_markers_intersect$cluster=='TE',]
te_markers_intersect<-te_markers_intersect[te_markers_intersect$p_val_adj <= 0.05 & te_markers_intersect$pct.1 >= 0.5 & te_markers_intersect$avg_log2FC >= 0.5,]
te_markers_intersect<-te_markers_intersect[order(te_markers_intersect$avg_log2FC,decreasing = T),]
write.table(te_markers_intersect,file = "./TE_markers_intersect.xls", row.names = F, col.names = T, quote = F, sep = "\t")
gene_features<-unique(c("Pou5f1","Nanog","Sox2",'Klf2','Zfp42','Esrrb','Tdgf1','Gdf3','Utf1','Fgf4','Tbx3',
                           "Gata4","Gata6","Sox17",'Pdgfra','Col4a1','Dab2','Foxa2','Grb10',
                           "Cdx2","Eomes","Elf5",'Gata2','Gata3','Tead3','Krt8','Tfap2c','Krt18','Tacstd2','Tead4','Dppa1',
                           te_markers_intersect$gene))
pdf("./heatmap_featuregenes.pdf",width = 6, height = 15)
DoHeatmap(object = em, features = gene_features,group.bar=T,angle=0,size=3,group.by = 'cell_type',hjust = 0.5)+
    theme(panel.border = element_blank())
dev.off()


pdf("./em.StackedVlnPlot.pdf", width = 7, height = 40)
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
StackedVlnPlot(em, features = c(all_markers,te_markers_intersect$gene),pt.size=0)
dev.off()

#count the number/proportion of cells expressing Nanog/Sox17/Cdx2/Sox2/Pdgfra/Tead3#
umap<-data.frame(em@reductions$umap@cell.embeddings)
umap$Nanog<-em@assays$RNA@data['Nanog',]
umap$Sox17<-em@assays$RNA@data['Sox17',]
umap$Cdx2<-em@assays$RNA@data['Cdx2',]
umap$Sox2<-em@assays$RNA@data['Sox2',]
umap$Pdgfra<-em@assays$RNA@data['Pdgfra',]
umap$Tead3<-em@assays$RNA@data['Tead3',]
meta<-data.frame(em@meta.data)
umap$celltype<-meta$cell_type


em.markers1 <- em.markers[em.markers$pct.1>=0.25 & em.markers$p_val_adj<=0.05,]
em.markers1<-em.markers1[order(em.markers1$avg_log2FC,em.markers1$cluster,decreasing = T),]
em.markers1<-em.markers1[order(em.markers1$cluster,decreasing = F),]
write.table(em.markers1,file = "./em.markers.xls", row.names = F, col.names = T, quote = F, sep = "\t")
