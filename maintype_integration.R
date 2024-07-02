library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer) 

# READ  rds if necessary
setwd('/Users/user/Documents/上海十院/实验数据/CRC/Single_cellSeq/A&T/AT_integ/')
A3 <- readRDS("./A3.rds")
A2 <- readRDS("./A2.rds")
A5 <- readRDS("./A5.rds")
A9 <- readRDS("./A9.rds")
A10 <- readRDS("./A10.rds")
A11 <- readRDS("./A11.rds")
A12 <- readRDS("./A12.rds")
A13 <- readRDS("./A13.rds")
T3 <- readRDS("./T3.rds")
T2 <- readRDS("./T2.rds")
T5 <- readRDS("./T5.rds")
T9 <- readRDS("./T9.rds")
T10 <- readRDS("./T10.rds")
T11 <- readRDS("./T11.rds")
T12 <- readRDS("./T12.rds")
T13 <- readRDS("./T13.rds")

integ <- list(A2=A2,A3=A3,A5=A5,A9=A9,A10=A10,A11=A11,A12=A12,A13=A13,
                T2=T2,T3=T3,T5=T5,T9=T9,T10=T10,T11=T11,T12=T12,T13=T13)

# joint.bcs <- intersect(colnames(A2), colnames(A3))

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
ATmerge$orig.ident <- factor(ATmerge$orig.ident,levels = c('A2','T2','A3','T3','A5','T5','A9','T9',
                                                                 'A10','T10','A11','T11','A12','T12','A13','T13'))
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
new.cluster.ids <- c('N2','N2','N3','N3','N5','N5','N9','N9',
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
          "CD79A", "CD79B", "KIT",'CPA3',"ITGAX", "CD68", "CD14",
          'IGHA1','IGHG2','S100A8','S100A9')
# plot
DotPlot(merge_sub, features = mark,cols = c('#FFFFFF','#FF0033') ) + RotatedAxis()+theme_test()+
  theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

new.cluster.ids <- c('BC','TC','NKT','TC','Plasma','NKT',
                     'Plasma','Epithelial','Myeloid','Fibro','Epithelial',
                     'Myeloid','Myeloid','unkonwn1','Endo','BC',
                     'unkonwn2','Mast','unkonwn3') 
names(new.cluster.ids) <- levels(ATmerge)
ATmerge <- RenameIdents(ATmerge, new.cluster.ids)
ATmerge$celltype <- Idents(ATmerge)
DimPlot(object = ATmerge, group.by='celltype',pt.size = 0.2,label = T)


# subset cluster exclude BT
Idents(ATmerge)<-'celltype'
merge_sub <- subset(ATmerge, idents = c('BC','TC','Plasma','NKT','Myeloid',
                                                'Epithelial','Fibro','Endo','Mast'))
merge_sub$celltype <- Idents(merge_sub)
merge_sub <- RunUMAP(merge_sub, dims = 1:10)

DimPlot(merge_sub, reduction = "umap", label = TRUE, pt.size = 0.1) + 
  theme_test() + theme(legend.position = 'bottom')
DimPlot(merge_sub, reduction = "umap", label = TRUE, pt.size = 0.1,split.by = 'stim') + 
  theme_test() + theme(legend.position = 'bottom')
DimPlot(merge_sub, reduction = "umap", label = TRUE, pt.size = 0.1) 

DimPlot(object = merge_sub, split.by='grade',pt.size = 0.2)+ 
  theme_test() + theme(legend.position = 'bottom')
DimPlot(object = merge_sub, group.by='grade',pt.size = 0.2)+ 
  theme_test() + theme(legend.position = 'bottom')

DimPlot(object = merge_sub, group.by='celltype',pt.size = 0.2) 

DimPlot(object = merge_sub, group.by='seurat_clusters',pt.size = 0.2,label = T)
pdf('AT_celltype.pdf',width = 8,height = 4)
DotPlot(merge_sub, features = mark ) + RotatedAxis()+theme_test()+
  theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_colour_gradientn (colours = c ("#0099CC", "white", "#CC0033"))
dev.off()

mer <- subset(x = merge_sub, downsample = 1000)

FeaturePlot(mer, features = c("EPCAM",  "COL3A1", "PECAM1","CD3D", 
                                    'NKG7',"CD79A","KIT",'CPA3', "CD68",'S100A8'),ncol = 5,
            keep.scale='all',cols =  c("lightgrey", "#FF9900") ,label.size = 0.1 ,pt.size=0.1)
FeaturePlot(merge_sub, features = c("CD3D",'NKG7',
                                    "CD68", 'IGHA1'),ncol = 2,
            keep.scale='all',cols =  c("lightgrey", "#FF9900") ,label.size = 3)
FeaturePlot(mer, features = c("EPCAM",  "COL3A1",'NKG7',"CD79A"),
            keep.scale='all',cols =  c("lightgrey", "#FF9900") ,label.size = 0.1,ncol = 2,
            pt.size=0.1) 
DimPlot(object = merge_sub, group.by='grade',pt.size = 0.2,label = F) + theme_test()
DimPlot(object = mtx, group.by='celltype',label = T,  label.size = 3)
FeaturePlot(merge_sub, features = c('IL34','SORT1'),raster = T,pt.size = 2)

VlnPlot(merge_sub, features = c('GRN'),group.by = 'stim',pt.size = 0.1)
FeaturePlot(Fibro, features = c('VEGFR1R2'))


Idents(merge_sub) <- 'celltype'
VlnPlot(object = mer, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,  pt.size = 0)

# change color
library(ggsci)

# down_mer <- subset(x = merge_sub, downsample = 4000)
DimPlot(object = merge_sub,group.by = 'celltype',label = F,cols =c(palette = 'Set1'),raster = T) 
DimPlot(object = merge_sub,split.by='stim',raster = T,pt.size = 2,group.by = 'celltype',label = F,cols =pal_aaas('default')(10))+
  RotatedAxis()+theme_test()
pdf('AT_celltype.pdf',width = 8,height = 4)
DimPlot(object = merge_sub,raster = T,split.by='stim',pt.size = 2,group.by = 'celltype',label = F,cols =pal_npg('nrc')(10)) +
  RotatedAxis()+theme_test() 
dev.off()
# saveRDS(merge_sub, file = "/Users/user/Documents/上海十院/实验数据/CRC/Single_cellSeq/merge_sub.rds")

# test
pdf('ctsl_stim.pdf',width = 8,height = 4)
FeaturePlot(merge_sub, features = c('IGFBP3','TMEM219'),
            keep.scale='all' ,raster = T,pt.size = 1.5)
dev.off()
VlnPlot(object = merge_sub, features = c("IGFBP3"), group.by = 'celltype', pt.size = 0.01)

pdf('c3_cfd.pdf',width = 8,height = 4)
FeaturePlot(merge_sub, features = c('IGFBP3'),
            keep.scale='all',label.size = 0.1,ncol = 2,raster = T,pt.size = 2) 
dev.off()


# fibro
Fibro <- readRDS("~/Documents/上海十院/实验数据/CRC/Single_cellSeq/rds_celltype/Fibro.rds")
DimPlot(Fibro, reduction = "umap",label = T)
VlnPlot(object = Fibro, features = c('COL3A1','FAP'),pt.size = 0.01)
pdf(file = './c3.pdf' ,height =8, width = 12)
FeaturePlot(Fibro, features = c('PTGDS','FBLN1','CFD','GSN'),pt.size =0.5)
dev.off()
DimPlot(object = Fibro, split.by ='sampling_site',pt.size = 0.2,label = T)
VlnPlot(object = Fibro, features = c('FAP'),pt.size = 0.01)
saveRDS(Fibro, file = "/Users/user/Documents/上海十院/实验数据/CRC/Single_cellSeq/database/E-MTAB-8410/Fibro.rds")

mark <- c( 'FAP','TWIST1',"MYH11",'ACTA2',"TAGLN",'IL6',"WNT2B","RSPO3",
           'WNT5B','THBS4','ADAMDEC1')
mark1 <- c('WNT2B','WNT4','BMP4','BMP5','WNT5A','WNT5B','FAP')
mark2 <- c('PI16','COL15A1','BMP4','BMP5','WNT5A','WNT5B','FAP')
pdf(file = './dot.pdf' ,height = 4, width = 7)
DotPlot(Fibro, features = mark,cols = c('#FFFFFF','#FF0033') ) + RotatedAxis()+theme_test()+
  theme(legend.position = 'top',axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))
dev.off()
Idents(Fibro) <-'sampling_site'
FeaturePlot(mtx, features = c('PDGFRA'),pt.size = 0.05)
Idents(Fibro) <-'seurat_clusters'
Allmarkers <- FindAllMarkers(Fibro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- Allmarkers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
# annotation
new.cluster.ids <- c('crypt fibroblasts','CAF1','crypt fibroblasts','myofibroblasts1','THBS4+ fibroblasts','ADAMDEC1+ fibroblasts',
                     'CAF2','iCAF','myofibroblasts1','myofibroblasts2','villus fibroblasts') 
names(new.cluster.ids) <- levels(Fibro)
Fibro <- RenameIdents(Fibro, new.cluster.ids)
Fibro$celltype <- Idents(Fibro)
pdf(file = './do.pdf' ,height = 5, width = 5)
DimPlot(object = Fibro, group.by='celltype',pt.size = 0.2,label = T)+ 
  theme_test() + theme(legend.position = 'bottom')
dev.off()
pdf(file = './do.pdf' ,height = 5, width = 5)
DimPlot(Fibro, reduction = "umap",label = T)+ 
  theme_test() + theme(legend.position = 'bottom')
dev.off()
DimPlot(Fibro, reduction = "umap",split.by = 'sampling_site', label = T)+ 
  theme_test() + theme(legend.position = 'bottom')
FeaturePlot(Fibro, features = c('SULF1'),pt.size = 0.05,min.cutoff = 0)

# percent analysis

Idents(merge_sub)<-'celltype'
merge_no_epi <- subset(merge_sub, idents = c('BC','TC','NKT','Plasma','Myeloid',
                                                'Fibro','Mast','Endo'))
table(Idents(merge_sub))
prop.table(table(Idents(merge_sub)))
prop.table(table(Idents(merge_no_epi),merge_no_epi$orig.ident))
cell.prop<-as.data.frame(prop.table(table(Idents(merge_sub), merge_sub$orig.ident)))
colnames(cell.prop)<-c ("Cluster","ID","Proportion") 

pdf('id_pro.pdf',width = 5,height = 4)
ggplot(data=cell.prop,aes(x=factor(ID,levels = c('A2','A3','A5','A9','A10','A11','A12','A13',
                                                 'T2','T3','T5','T9','T10','T11','T12','T13')) ,y=Proportion,fill=Cluster))+
  geom_bar(stat="identity",position="fill")+scale_fill_nejm()+ggtitle("")+theme_bw()+xlab('Patient')
dev.off()

# AT CELLNUMBER
Idents(merge_sub)<-'orig.ident'
a <- table(Idents(merge_sub))
a <- data.frame(a)
b <- c('P2','P2','P3','P3','P5','P5','P9','P9','P10','P10','P11','P11','P12','P12','P13','P13')
d <- c('A','T')
c <- cbind(a,b)
c<- cbind(c,d)
library(ggplot2)
library(plyr)
ggplot(c,aes(x=b,y=Freq,fill=d))+geom_bar(stat = 'identity')+xlab(' ')+ylab('Cells Number')+
  scale_x_discrete(limits=c('P11','P9','P3','P2','P5','P10','P12','P13'))+
  scale_fill_manual(values = c('#99CC00','#CC0033'))+theme_bw()+
  theme(panel.grid = element_blank()) +geom_text(aes(label=Freq))


# grade
table(Idents(merge_no_epi))
prop.table(table(Idents(merge_no_epi)))
prop.table(table(Idents(merge_no_epi),merge_no_epi$grade))
cell.prop<-as.data.frame(prop.table(table(Idents(merge_no_epi), merge_no_epi$grade)))
colnames(cell.prop)<-c("Celltype","Grade","Proportion")

pdf('grade_pro.pdf',width = 3.5,height = 4)

ggplot(data=cell.prop,aes(x=Grade,y=Proportion,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+scale_fill_aaas()+ggtitle("")+theme_bw()
dev.off()

# AT
Idents(merge_sub)<-'celltype'
table(Idents(merge_sub))
prop.table(table(Idents(merge_sub)))
prop.table(table(Idents(merge_sub),merge_sub$stim))
cell.prop<-as.data.frame(prop.table(table(Idents(merge_sub), merge_sub$stim)))
colnames(cell.prop)<-c("Celltype","Stim","Proportion")

pdf('at_pro.pdf',width = 3,height = 4)
ggplot(data=cell.prop,aes(x=Stim,y=Proportion,fill= Celltype))+
  geom_bar(stat="identity",position="fill")+ggtitle("")+theme_bw()
dev.off()

table(Idents(merge_no_epi))
prop.table(table(Idents(merge_no_epi)))
prop.table(table(Idents(merge_no_epi),merge_no_epi$stim))
cell.prop<-as.data.frame(prop.table(table(Idents(merge_no_epi), merge_no_epi$stim)))
colnames(cell.prop)<-c("Celltype","Stim","Proportion")
ggplot(data=cell.prop,aes(x=Stim,y=Proportion,fill= Celltype))+
  geom_bar(stat="identity",position="fill")+ggtitle("")+theme_bw()

# kras
table(Idents(merge_no_epi))
prop.table(table(Idents(merge_no_epi)))
prop.table(table(Idents(merge_no_epi),merge_no_epi$kras))
cell.prop<-as.data.frame(prop.table(table(Idents(merge_no_epi), merge_no_epi$kras)))
colnames(cell.prop)<-c("Celltype","KRAS","Proportion")
ggplot(data=cell.prop,aes(x=KRAS,y=Proportion,fill= Celltype))+
  geom_bar(stat="identity",position="fill")+ggtitle("")+theme_bw()


# filtrate save
Idents(merge_sub)<-'stim'
a <- merge_sub@active.ident
Idents(merge_sub)<-'orig.ident'
b <- merge_sub@active.ident
Idents(merge_sub)<-'celltype'
c <- merge_sub@active.ident
Idents(merge_sub)<-'grade'
d <- merge_sub@active.ident
Idents(merge_sub)<-'kras'
e <- merge_sub@active.ident

a <-as.matrix(a)
b <-as.matrix(b)
c <-as.matrix(c)
d <-as.matrix(d)
e <-as.matrix(e)

# cbind all 
# ????
f<- cbind(a,b)
f<- cbind(f,c)
f<- cbind(f,d)
f<- cbind(f,e)

write.csv(f,file = 'C:/Users/df/Desktop/info.csv')

## Add genetic info
Idents(merge_sub)<-'orig.ident'

merge_sub$orig.ident <- factor(merge_sub$orig.ident,levels = c('A2','T2','A3','T3','A5','T5','A9','T9',
                                                 'A10','T10','A11','T11','A12','T12','A13','T13'))
new.cluster.ids <- c('NKM','NKM','NKM','NKM','NKM','NKM','KM','KM',
                     'NKM','NKM','NKM','NKM','KM','KM','KM','KM')
names(new.cluster.ids) <- levels(merge_sub)
merge_sub <- RenameIdents(merge_sub, new.cluster.ids)
merge_sub$kras<-Idents(merge_sub)
saveRDS(merge_sub, file = "/Volumes/E ???/SHSY/Single_cellSeq/A&T/AT_integ/merge_sub.rds")
DimPlot(object = merge_sub, group.by='kras',pt.size = 0.2,label = T)

a<- DotPlot(merge_sub,features = 'C3',cols = c("green","red"),group.by = 'celltype')+coord_flip()+theme_bw()+
  theme(panel.grid.major.x = element_blank(),axis.text.x = element_text(size = 12,angle = 45,hjust = 1) )

b <- a$data %>% select('pct.exp','avg.exp','id')
b$tumor <-  'crc'
write_csv(b,file = '/Users/user/Documents/上海十院/实验数据/CRC/Single_cellSeq/database/c3/crc_c3.csv')

pdf(file = './c3v.pdf' ,height = 4, width = 6)
VlnPlot(object = merge_sub, features = c('ITGB2','ITGAM','C3AR1','C5AR1'),group.by = 'celltype',pt.size = 0.01)
dev.off()

library(CellChat)
Idents(merge_sub)<-'celltype'
pdf(file = './c3v.pdf' ,height = 4, width = 6)
StackedVlnPlot(merge_sub, c('ITGB2','ITGAM','ITGAX','C3AR1','C5AR1'),
               angle.x = 45,pt.size = 0.0001)

library(CellChat)
Idents(merge_sub)<-'celltype'
cy <- cx[c(21:31),]
pdf(file = './c3v.pdf' ,height = 4, width = 6)
StackedVlnPlot(merge_sub, cy$id,
               angle.x = 45)
dev.off()


VlnPlot(merge_sub, features = c('CCL20'),group.by = 'celltype',pt.size = 0,split.by = 'stim')
VlnPlot(merge_sub, features = c('CXCL8'),group.by = 'celltype',pt.size = 0,split.by = 'stim')
FeaturePlot(merge_sub, features = c('CCL20','CXCL8'),split.by = 'stim')

# annotation
Idents(merge_sub)<-'seurat_clusters'
DimPlot(object = merge_sub, group.by='seurat_clusters',pt.size = 0.2,label = T)

mark <- c("EPCAM", "THY1", "COL3A1", "PECAM1",
          "CD3D", "CD3E", "CD3G","KLRD1", "PRF1", "GNLY",'NKG7',
          "CD79A", "CD79B", "KIT",'CPA3',"ITGAX", "CD68", "CD14",'IGHA1','IGHG2','S100A8')
DotPlot(merge_sub, features = mark,cols = c('#FFFFFF','#FF0033') ) + RotatedAxis()+theme_test()+
  theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

new.cluster.ids <- c('BC','TC','NKT','TC','Plasma','NKT',
                     'Plasma','Epithelial','Myeloid','Fibro','Epithelial',
                     'Neutrophil','Myeloid','Endo','BC','Mast') 
names(new.cluster.ids) <- levels(merge_sub)
merge_sub <- RenameIdents(merge_sub, new.cluster.ids)
merge_sub$celltype <- Idents(merge_sub)
DimPlot(object = merge_sub, group.by='celltype',pt.size = 0.2,label = T)
saveRDS(merge_sub, file = "~/Documents/上海十院/实验数据/CRC/Single_cellSeq/merge_sub.rds")

# plot cell number
df <- matrix(runif(n= 32 , min= 1 , max= 20 ), nrow= 2 )
ha <- HeatmapAnnotation(
  cell = anno_barplot(
    c(5775,8685,4743,7534,1958,9247,1838,1963,
      4762,5882,4708,8784,4378,11418,4401,2577), 
    bar_width = 1, 
    gp = gpar(col = "white", fill = "#FFE200"), 
    border = FALSE,
    height = unit(2, "cm")
  )
  
)
pdf('./heat0.pdf',width = 6,height = 4)
Heatmap(
  df, 
  name = "mat", 
  top_annotation = ha,
  cluster_columns = F
)
dev.off()

Idents(merge_sub) <- 'celltype'
pdf('./vln1.pdf',width = 7,height = 4)
VlnPlot(object = merge_sub, features = c('CFD'),split.by = 'stim',pt.size = -0.01,ncol = 1,
        cols = c('#99CC00','#CC0033'))
dev.off()

# epi
Epithelial <- readRDS("~/Documents/上海十院/实验数据/CRC/Single_cellSeq/rds_celltype/Epithelial.rds")
Idents(Epithelial) <- 'stim'
Allmarkers <- FindAllMarkers(Epithelial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- Allmarkers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
DimPlot(object = Epithelial, split.by ='stim',pt.size = 0.2,label = T)
ggsave('./myc.pdf',width = 8,height = 4)
FeaturePlot(Epithelial, features = c('PLAC8'),split.by = 'stim',pt.size = 0.5,min.cutoff = 0)

Idents(merge_sub)<-'celltype'
sub_crc <- subset(merge_sub,downsample=3000)
crc_matrix <-as.matrix(sub_crc@assays[["RNA"]]@data)
list=list(ca=inter121311109532)
gsva_matrix <- gsva(crc_matrix, list, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
sub_crc$ca <- data.frame(t(gsva_matrix)) 

VlnPlot(object = ST2, features = 'HALLMARK-MYOGENESIS',group.by = 'interface',pt.size = -0.01,ncol = 1)


data2 <- as.data.frame(t(gsva_matrix)) 
Epithelial@meta.data<- cbind(Epithelial@meta.data,type=data2)

# all marker
# Find marker
Idents(Epithelial) <- 'stim'

markers <- FindAllMarkers(Epithelial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top25 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 25, order_by = avg_log2FC)
DoHeatmap(Epithelial, features = top25$gene) 
ggsave('~/Desktop/at.pdf',width = 5,height = 7)


library(CellChat)
genefilt <- top25 %>% filter(cluster=='Endo') %>% select(gene) 
genefilt_mtx <- top25_mtx %>% filter(cluster=='Endo') %>% select(gene) # top25_mtx from mtx data
StackedVlnPlot(merge_sub,features = genefilt$gene) 
StackedVlnPlot(mtx,features = genefilt_mtx$gene)

intersect(genefilt$gene,genefilt_mtx$gene)

# filter gene by specific in vlnplot and featureplot in mergesub and mtx 
BC <- union(genefilt$gene[c(1:10,12:20,23)],genefilt_mtx$gene[c(1:10)])
TC <- union(genefilt$gene[c(1:22)],genefilt_mtx$gene[c(1,4,23,25)])
NKT <- union(genefilt$gene[c(1:22)],genefilt_mtx$gene[c(1:8,13,14,16,19,21,22,24)])
Plasma <- union(genefilt$gene[c(1,3:9,11:19,21)],genefilt_mtx$gene[c(3,6:8,10:11,14:16,19,22,24)])
Epithelial <- union(genefilt$gene[c(1:9,11:12,14:15,17:20)],genefilt_mtx$gene[c(1:10)])
Myeloid <- union(genefilt$gene[c(1:25)],genefilt_mtx$gene[c(2)])
Fibro <- union(genefilt$gene[c(1:16,18:21,25)],genefilt_mtx$gene[c(1:2,4:10)])
Neutrophil <- genefilt$gene[c(1:9,11:17,19:25)]
Endo <- union(genefilt$gene[c(1:5,7:11,13:24)],genefilt_mtx$gene[c(1:10)])
Mast <- union(genefilt$gene[c(1:15,17:20)],genefilt_mtx$gene[c(2:12)])
# merge for save
geneset <- list(B_cell=BC,T_cell=TC,NKT=NKT,Epithelial=Epithelial,Plasma=Plasma,Fibroblast=Fibro,
                Endothelial=Endo,Macrophage=Myeloid,Neutrophil=Neutrophil,Mast=Mast)

alltype_gene <- c(BC,TC,NKT,Plasma,Epithelial,Myeloid,Fibro,Neutrophil,Endo,Mast)# 262 genes
save(top25,top25_mtx,geneset,alltype_gene,file = './celltype_geneset.RData')

FeaturePlot(merge_sub, features = genefilt$gene[c(18:25)],raster = T,pt.size = 3)
FeaturePlot(merge_sub, features = genefilt_mtx$gene[c(7:10)],raster = T,pt.size = 3)

StackedVlnPlot(merge_sub,features = Endo,pt.size = 0.1)
FeaturePlot(merge_sub, features = 'C3AR1',raster = T,split.by = 'stim',pt.size = 3)


# macro
ccl <- c('CCL8','CXCL10','CXCL8','CXCL2','IL10RA',
         'CXCL3','TNFSF10','TNFSF13B','CCL2','LIF','CXCL11',
         'CSF1','IL12RB2','CXCL1','CCL3L3','TNF')
FeaturePlot(macro, features = 'CCR2',split.by = 'stim', pt.size = 0.5,min.cutoff = 0)
VlnPlot(macro, features = c('C3AR1'),group.by  = 'stim',pt.size = 0.05)
FeaturePlot(merge_sub, features = c('CTSB','CTSL'),split.by = 'stim',raster = T,pt.size = 2,min.cutoff = 0)
FeaturePlot(macro, features = c('CD1C','C1S'),split.by = 'stim',raster = T,pt.size = 3,min.cutoff = 0)

DefaultAssay(Epithelial) <- 'RNA'
FeaturePlot(Epithelial, features = c('IGLC3','LYPD8'),split.by  = 'stim',pt.size = 0.05)




