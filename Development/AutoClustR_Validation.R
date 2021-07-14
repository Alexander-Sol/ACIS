#AutoClustR Validation

library(BiocManager)
BiocManager::install("GEOquery")
library(GEOquery)
install.packages("aricode")
library(aricode)

#Load in the mouse brain dataset
#This doesn't seem to be annotated
test <- read.table("C:/Users/asolivai/Desktop/R_Files/AutoClustR/MouseBrainCells/GSE60361_C1-3005-Expression.txt/GSE60361_C1-3005-Expression.txt")
test.head <- test[1, ]
colnames(test) <- test.head
rownames(test) <- make.names(test[ ,1], unique = T)
nchar(test.head[2])
mb.annotations <- substring(test.head, 12, 14)
mb.letters <- substring(test.head[2:3006], 12, 12)
unique(mb.letters)

series.matrix <- getGEO(filename = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/MouseBrainCells/GSE60361_series_matrix.txt.gz")
miniml <- getGEO(filename = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/MouseBrainCells/GSE60361_family.soft.gz")


#Can't figure out annotation, going to go with GSE72056, melanoma tumors
melanoma <- getGEO(GEO = "GSE72056")

mel.test <- read.table("C:/Users/asolivai/Desktop/R_Files/AutoClustR/Melanoma/GSE72056_melanoma_single_cell_revised_v2.txt/GSE72056_melanoma_single_cell_revised_v2.txt",
                       header = T)
row.names(mel.test) <- make.names(mel.test[ ,1 ], unique = T)
mel.test <- mel.test[ ,-1]
tumor.no <- as.numeric(mel.test[1, ])
malignancy <- as.numeric(mel.test[2, ])
annotations <- as.numeric(mel.test[3, ])
mel.test <- mel.test[-c(1,2,3),]

table(annotations)
mel.seurat <- CreateSeuratObject(mel.test, project = "MelanomaTest")
VlnPlot(mel.seurat, features = "mt.per")
mel.seurat$tumor <- tumor.no
mel.seurat$annotations <- annotations
mel.seurat$malignancy <- malignancy
mel.seurat <- SCTransform(mel.seurat)
mel.seurat <- AutoClustR(mel.seurat)
DimPlot(mel.seurat)
class(mel.seurat@active.ident)
MARI(mel.seurat@active.ident, mel.seurat$annotations)

VlnPlot(mel.seurat, "RPL27")

ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

ClusterPurity(mel.seurat@active.ident, mel.seurat$annotations)



#Selected the wrong number of pcs without prefiltering. I'm gonna prefilter based on RPL27 exp and see what happens

mel.seurat <- housekeeping.filter(mel.seurat, housekeeping.gene = "RPL27")
VlnPlot(mel.seurat, "nCount_RNA")
mean(mel.seurat$nCount_RNA)
sd(mel.seurat$nCount_RNA)
subset(mel.seurat)
mel.seurat <- subset(mel.seurat, hk > mean(mel.seurat$hk) - 2.5*sd(mel.seurat$hk) &
                   hk < mean(mel.seurat$hk) + 2.5*sd(mel.seurat$hk))
mel.seurat <- SCTransform(mel.seurat)
mel.seurat <- Predict_nPCs(mel.seurat, period = 5)
mel.seurat <- AutoClustR(mel.seurat, period = 5)
DimPlot(mel.seurat)
NMI(mel.seurat@active.ident, mel.seurat$annotations)

#using a period of 5 for pc estimation matches with a visual elbow selection, but it significantly reduces ARI and NMI of the clustering solution



#Load in Pancreas dataset
panc.1 <- read.csv("C:/Users/asolivai/Desktop/R_Files/AutoClustR/PancreaticIslet1/GSM2230757_human1_umifm_counts.csv/GSM2230757_human1_umifm_counts.csv")
panc.1 <- t(panc.1)
colnames(panc.1) <- panc.1[2, ]
panc.1.annotation <- panc.1[3,]
table(as.character(panc.1.annotation))
panc.1 <- panc.1[-c(1,2,3), ]
panc.1.seurat <- CreateSeuratObject(panc.1, project = "HumanPancreas2")
rm(panc.1)

panc.1.seurat <- 

VlnPlot(panc.1.seurat, "nCount_RNA")
plot(panc.1.seurat@assays$RNA["RPL27A"], panc.1.seurat$nCount_RNA)
mean(log(p1s$nCount_RNA))
sd(log(p1s$nCount_RNA))
density(p1s$logcounts) %>% plot()
p1s <- panc.1.seurat
p1s$logcounts <- log(p1s$nCount_RNA)
p1s <- subset(p1s, logcounts > mean(p1s$logcounts) - 2*sd(p1s$logcounts) &
              logcounts < mean(p1s$logcounts) + 2*sd(p1s$logcounts))
p1s <- SCTransform(p1s)
p1s$annotations <- panc.1.annotation
p1s$na <- !is.na(p1s$annotations)
p1s <- subset(p1s, na)
p1s <- AutoClustR(p1s)
ari.algo.1 <- ARI(p1s@active.ident, p1s$annotations)
nmi.algo.1 <- NMI(p1s@active.ident, p1s$annotations)
sum(is.na(p1s$annotations))
#Running again but using SLM in Findclusters
p1s.louv.3 <- AutoClustR(p1s)
ari.algo.2 <- ARI(p1s.louv.2@active.ident, p1s$annotations)
ari.algo1.2 <- ARI(p1s.louv.3@active.ident, p1s$annotations)
nmi.algo.2 <- NMI(p1s.slm@active.ident, p1s$annotations)
table(p1s@active.ident)
table(p1s$annotations)
b <- NA
is.na(b)

#Load in Pancreas Dataset #2
panc.2 <- read.csv("C:/Users/asolivai/Desktop/R_Files/AutoClustR/PancreaticIslet2/GSM2230758_human2_umifm_counts.csv.gz")
panc.2 <- t(panc.2)
colnames(panc.2) <- panc.2[2, ]
panc.2.annotation <- panc.2[3,]
table(as.character(panc.2.annotation))
panc.2 <- panc.2[-c(1,2,3), ]
panc.2.seurat <- CreateSeuratObject(panc.2, project = "HumanPancreas2")

rm(panc.2)



#Taken from old code
pbmc3 <- Read10X("C:/Users/asolivai/Desktop/R_Files/AutoClustR/pbmc3k_filtered_gene_bc_matrices.tar/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
p3 <- CreateSeuratObject(pbmc3, project = "PBMC_3000")
occ <- Read10X("C:/Users/asolivai/Desktop/R_Files/AutoClustR/Menon2019_GSE137537/GSE137537/", gene.column = 1)
retina <- CreateSeuratObject(occ, project = "Menon2019")
retina$annotation <- occular.ann$Labels

#For IKAP, need to normalize, compute pct mito, then pass to IKAP.
filterPlot <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  #obj <- subset(obj, percent.mt < 15 & nFeature_RNA > 200 & 
  #               nCount_RNA < 36500 & nCount_RNA > 1000)
  return(obj)
}

retina <- filterPlot(retina)
retina <- SCTransform(retina)


VlnPlot(retina, features = "percent.mt")

occular.ann <- read.table("C:/Users/asolivai/Desktop/R_Files/AutoClustR/Menon2019_GSE137537/GSE137537/GSE137537_sample_annotations.tsv", header = TRUE)


head(occular.ann$Labels)
table(occular.ann$Labels)

write.table(occular.ann$Barcode, file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/Menon2019_GSE137537/GSE137537/barcodes.tsv",
            sep = '\t', col.names = FALSE, quote = FALSE, row.names = FALSE)

#Testing AutoClustR

day.60.sub <- subset(day60.seurat, cells = sample(Cells(day60.seurat), 
                                                  replace = F,
                                                  size = 2000))
d60.sub.auto <- AutoClustR(day.60.sub)







