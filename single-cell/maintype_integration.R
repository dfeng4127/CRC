library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer) 

setwd('./CRC/Single_cellSeq/A&T/AT_integ/')
# preset
# for all A and T sample 
data <- Read10X(data.dir = "./filtered_feature_bc_matrix")
data <- CreateSeuratObject(counts = data, project = "N1", min.cells = 3, min.features = 200)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)  

# integration
integ <- list(N1=N1,N2=N2,N3=N3,N4=N4,N5=N5,N6=N6,N7=N7,N8=N8,
              T1=T1,T2=T2,T3=T3,T4=T4,T5=T5,T6=T6,T7=T7,T8=T8)


#  standard normalization and variable feature selection
integ <- lapply(X = integ, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE,scale.factor = 100000)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = integ)
integ <- lapply(X = integ, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# reference can be changed 
anchors <- FindIntegrationAnchors(object.list = integ, reference = c(2, 7, 10, 15), reduction = "rpca", dims = 1:20) 
ATmerge <- IntegrateData(anchorset = anchors, dims = 1:50)

ATmerge <- ScaleData(ATmerge, verbose = FALSE)
ATmerge <- RunPCA(ATmerge, verbose = FALSE)
ATmerge <- FindNeighbors(ATmerge, dims = 1:10)
ATmerge <- FindClusters(ATmerge, resolution = 0.5)
ATmerge <- RunUMAP(ATmerge, dims = 1:15)

# add information
Idents(ATmerge)<-'orig.ident'
ATmerge$orig.ident <- factor(ATmerge$orig.ident,levels = c('N1','T1','N2','T2','N3','T3','N4','T4',
                                                           'N5','T5','N6','T6','N7','T7','N8','T8'))
new.cluster.ids <- c('NKM','NKM','NKM','NKM','NKM','NKM','KM','KM',
                     'NKM','NKM','NKM','NKM','KM','KM','KM','KM')
names(new.cluster.ids) <- levels(ATmerge)
ATmerge <- RenameIdents(ATmerge, new.cluster.ids)
ATmerge$kras<-Idents(ATmerge)

Idents(ATmerge) <-'orig.ident'
new.cluster.ids <- c('IV','IV','I','I','III','III','I','I',
                     'II','II','II','II','III','III','I','I')
names(new.cluster.ids) <- levels(ATmerge)
ATmerge <- RenameIdents(ATmerge, new.cluster.ids)
ATmerge$grade<- Idents(ATmerge)
ATmerge$grade <- factor(ATmerge$grade,levels = c('I','II','III','IV'))

Idents(ATmerge)<-'orig.ident'
new.cluster.ids <- c('A','T','A','T','A','T','A','T',
                     'A','T','A','T','A','T','A','T')
names(new.cluster.ids) <- levels(ATmerge)
ATmerge <- RenameIdents(ATmerge, new.cluster.ids)
ATmerge$stim<-Idents(ATmerge)

Idents(ATmerge)<-'orig.ident'
new.cluster.ids <- c('N1','N1','N2','N2','N5','N5','N9','N9',
                     'N10','N10','N11','N11','N12','N12','N13','N13')
names(new.cluster.ids) <- levels(ATmerge)
ATmerge <- RenameIdents(ATmerge, new.cluster.ids)
ATmerge$patient<-Idents(ATmerge)

# Find marker
Idents(ATmerge)<-'seurat_clusters'
AT_markers <- FindAllMarkers(ATmerge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# annotation
Idents(ATmerge)<-'seurat_clusters'
DimPlot(object = ATmerge, group.by='seurat_clusters',pt.size = 0.2,label = T)

mark <- c("EPCAM", "THY1", "COL3A1", "PECAM1",
          "CD3D", "CD3E", "CD3G","KLRD1", "PRF1", "GNLY",'NKG7',
          "CD79A", "CD79B", "KIT",'CPA3',"ITGAX", "CD68", "CD14",'IGHA1','IGHG2','S100A8')
DotPlot(merge_sub, features = mark,cols = c('#FFFFFF','#FF0033') ) + RotatedAxis()+theme_test()+
  theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

new.cluster.ids <- c('BC','TC','NKT','TC','Plasma','NKT',
                     'Plasma','Epithelial','Myeloid','Fibro','Epithelial',
                     'Neutrophil','Myeloid','Endo','BC','Mast') 
names(new.cluster.ids) <- levels(ATmerge)
ATmerge <- RenameIdents(ATmerge, new.cluster.ids)
ATmerge$celltype <- Idents(ATmerge)

# proportion
cell.prop<-as.data.frame(prop.table(table(Idents(ATmerge), ATmerge$stim)))
colnames(cell.prop)<-c("Celltype","Stim","Proportion")
ggplot(data=cell.prop,aes(x=Stim,y=Proportion,fill= Celltype))+
  geom_bar(stat="identity",position="fill")+ggtitle("")+theme_bw()
DimPlot(object = ATmerge,raster = T,split.by='stim',pt.size = 2,group.by = 'celltype',label = F,cols =pal_npg('nrc')(10)) +
  RotatedAxis()+theme_test() 
