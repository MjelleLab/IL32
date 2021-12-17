R-Codes


################################ Differentially expressed genes INA-6 ###############################
require(clusterProfiler)
require(org.Hs.eg.db)
require(edgeR)
require(ggplot2)



mat.cellines <- read.table("Matrix_CellLines.csv",row.names=1) #Count table of WT and KO

INA <- DGEList(mat.cellines)
keep <- rowSums(INA$counts>1) >= dim(INA)[2]/2
INA <- INA[keep,]
method <- "TMM" # c("TMM","RLE","upperquartile","none")
y <- calcNormFactors(INA, method=method)
groups <- c("WT","WT","KO1","KO1","KO2","KO2")
groups.df <- as.data.frame(groups)
colnames(groups.df) <- "INA6"
des <- model.matrix(~0 + groups)  
v <- voom(y,plot = T)
fit <- lmFit(v, design=des)
fit <- eBayes(fit)
#fit <- treat(fit)
contrasts <- makeContrasts(KO1KO2_vs_WT=(0.5*(groupsKO1+groupsKO2))-groupsWT,
                           levels=des) 

fit2 <- contrasts.fit(fit, contrasts=contrasts)
fit2 <- eBayes(fit2)

colSums(decideTests(fit2)!=0)
colSums(decideTests(fit2, method="glob")!=0)

KO1KO2_vs_WT <- topTable(fit2, coef="KO1KO2_vs_WT", sort.by="P", n=Inf)
KO1KO2_vs_WT$ENSEMBL <- rownames(KO1KO2_vs_WT)

KO1KO2_vs_WT.names <- bitr(rownames(KO1KO2_vs_WT), fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"),OrgDb = org.Hs.eg.db)
KO1KO2_vs_WT.names.merge <- merge(KO1KO2_vs_WT,KO1KO2_vs_WT.names,by="ENSEMBL")

KO1KO2_vs_WT.names.merge$treshold <- ifelse((KO1KO2_vs_WT.names.merge$logFC < 0), "Down","Up")

require(ggrepel)
pdf(file="Vulcano_INA.pdf",width = 2.9,height = 1.9,useDingbats = F, compress = T)
ggplot(KO1KO2_vs_WT.names.merge, aes(logFC ,-log10(adj.P.Val),colour=treshold)) + geom_point(alpha=0.6,size=0.8) +
  theme_bw(base_size=8) +  labs(x="logFC (KO vs WT)", y="Adjusted p-value (-log10)") +
  geom_text_repel(data=subset(KO1KO2_vs_WT.names.merge,-log10(adj.P.Val) > 8.6 & logFC < 0 | -log10(adj.P.Val)>8.87 & logFC >0),segment.size = 0.2,
                  aes(logFC ,label=SYMBOL  ),size=1.2,force=0, hjust=1) +
  scale_colour_manual(values = c("Up"= "#e6550d", "Down"= "#3182bd"),name="") +
  geom_hline(yintercept=1.3, colour="#bdbdbd") 

dev.off()

require(ggbio)
pca.INA <- prcomp(as.data.frame(t(cpm(INA$counts,log=T))))
pdf(file="INA_PCA.pdf",width = 2.7,height = 2.2,useDingbats = F,compress = T)
autoplot(pca.INA, data = groups.df, colour = 'INA6',frame = FALSE,alpha=0.75) + theme_bw(base_size=8) +scale_color_manual(values =c("#e6550d","#e6550d","#4292c6","#4292c6","#636363","#636363"))
dev.off()


######################### END ###########################################


###################### #Differentially expressed genes in MMRF CoMMpass ###############


## read count matrix containing baseline samples: 

df.baseline.orderInClinics <- read.table(file="MMRF_CoMMpass_IA13a_E74GTF_Salmon_Gene_Counts.txt",header=T)

dataComp.dge <- DGEList(df.baseline.orderInClinics) 
keep <- rowSums(edgeR::cpm(dataComp.dge)>1) >= dim(dataComp.dge)[2]/8
dataComp.dge <- dataComp.dge[keep,]
dataComp.dge.cpm <- cpm(dataComp.dge$counts,log=T)
IL32Factor <- as.numeric(dataComp.dge.cpm[grep("ENSG00000008517",rownames(dataComp.dge.cpm)),])
IL32Factor <- IL32Factor
quantile(IL32Factor,seq(0,1,0.1))
plotDensities(as.numeric(as.character(IL32Factor)))
IL32Factor <- as.numeric(as.character(IL32Factor))
IL32Factor[as.numeric(IL32Factor) > 1.5197169] <- "High"
IL32Factor[as.numeric(IL32Factor) <= 1.5197169] <- "Low"
IL32Factor <- factor(IL32Factor)
table(IL32Factor)

dataComp.dge$samples$lib.size <- colSums(dataComp.dge$counts)
method <- "TMM" # c("TMM","RLE","upperquartile","none")
dataComp.dge <- calcNormFactors(dataComp.dge, method=method)
cc <- ifelse(IL32Factor =="High", "pink", "green3")
des <- model.matrix(~IL32Factor)
colnames(des) <- c("Intercept", "cond")
v <- voom(dataComp.dge,plot = F,design=des)
fit <- lmFit(v, des)
fit <- eBayes(fit)
contrasts <- makeContrasts(cond=cond, levels=des) 
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
colSums(decideTests(fit2) != 0)
High_vs_Low_IL32 <- topTable(fit2,coef="cond",sort.by="P", adjust.method="BH", n=Inf)



######################## END ####################################






####################### Survival MMRF CoMMpass ##########################



fitSurv.os <- survfit(Surv(Time.mmrf,Status.mmrf) ~ IL32Factor, 
                      data = as.data.frame((df.baseline.orderInClinics.surv.t)))


coxph(Surv(Time.mmrf,Status.mmrf) ~ IL32Factor, 
      data = as.data.frame((df.baseline.orderInClinics.surv.t)))


# df.baseline.orderInClinics.surv.t is the data frame containing the survival data
# Time.mmrf is the survival time
# Status.mmrf is the status of the patient
# IL32Factor is a vector defining which patients have high (top Q10) and low (low Q90) IL32-expression. The vector is similar to the IL32Factor object above. 


require(survminer)
ggsurvplot(fitSurv.os ,title="", xlab="Time (Yrs)", ylab="Overall survival prbability",
                             pval.size = 3,
                             pval.coord = c(1000,1),
                             size=0.4,
                             legend = "right",
                             censor.size=2,
                             break.time.by = 365,
                             pval ="8.9e-5",#"p=0.41",
                             palette = c("#4292c6","#e6550d"), 
                             risk.table = T,
                             xscale=365.25,
                             xlim=c(0,5*365),
                             legend.title = "Expression",legend.labs = c("Low (Q90)","High (Q10)"))
######################## END ####################################



###################### Gene ontology MMRF Compass ##########################
require(org.Hs.eg.db)
require(tidyverse)
require(clusterProfiler)


High_vs_Low_IL32.sign <- High_vs_Low_IL32[High_vs_Low_IL32$adj.P.Val<0.05,]
### High_vs_Low_IL32 is the abovementioned limma-voom results from  the differential expression analysis of High vs Low IL32
High_vs_Low_IL32.sign.up <- High_vs_Low_IL32.sign[High_vs_Low_IL32.sign$logFC>0,]
High_vs_Low_IL32.sign.up <- High_vs_Low_IL32.sign.up[High_vs_Low_IL32.sign.up$AveExpr>0,]
High_vs_Low_IL32.sign.down <- High_vs_Low_IL32.sign[High_vs_Low_IL32.sign$logFC< 0,]
High_vs_Low_IL32.sign.down <- High_vs_Low_IL32.sign.down[High_vs_Low_IL32.sign.down$AveExpr>0,]


High_vs_Low_IL32.uni <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(High_vs_Low_IL32),columns=c("ENTREZID"), keytype="ENSEMBL") 
High_vs_Low_IL32.sign.up.annotate <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(High_vs_Low_IL32.sign.up),columns=c("ENTREZID"), keytype="ENSEMBL") 
High_vs_Low_IL32.sign.up.annotate.enrich <- enrichGO(gene= High_vs_Low_IL32.sign.up.annotate$ENTREZID,universe= High_vs_Low_IL32.uni$ENTREZID,OrgDb         = org.Hs.eg.db,
                                              ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable      = TRUE)

High_vs_Low_IL32.sign.up.annotate.enrich.simplify <- clusterProfiler::simplify(High_vs_Low_IL32.sign.up.annotate.enrich) 
High_vs_Low_IL32.sign.up.annotate.enrich.df <- as.data.frame(High_vs_Low_IL32.sign.up.annotate.enrich.simplify)
High_vs_Low_IL32.sign.up.annotate.enrich.df$log10P <- -log10(High_vs_Low_IL32.sign.up.annotate.enrich.df$p.adjust)
High_vs_Low_IL32.sign.up.annotate.enrich.df.order <- High_vs_Low_IL32.sign.up.annotate.enrich.df[order(High_vs_Low_IL32.sign.up.annotate.enrich.df$log10P,decreasing = T),]
High_vs_Low_IL32.sign.up.annotate.enrich.df.order.top <- head(High_vs_Low_IL32.sign.up.annotate.enrich.df.order,20)

ggplot(High_vs_Low_IL32.sign.up.annotate.enrich.df.order.top, aes(x = reorder(Description,log10P), y = log10P)) + 
  geom_bar(stat="identity",fill="#005824") +
  theme_bw(base_size = 8) +
  theme(legend.key.size = unit(0.4, "cm")) +
  coord_flip() + 
  labs(y=expression(-Log['10'](Adj.~p-value)), x="")


High_vs_Low_IL32.sign.down.annotate <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(High_vs_Low_IL32.sign.down),columns=c("ENTREZID"), keytype="ENSEMBL") 
High_vs_Low_IL32.sign.down.annotate.enrich <- enrichGO(gene= High_vs_Low_IL32.sign.down.annotate$ENTREZID,universe= High_vs_Low_IL32.uni$ENTREZID,OrgDb         = org.Hs.eg.db,
                                                ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.5,qvalueCutoff  = 0.5,readable      = TRUE)

High_vs_Low_IL32.sign.down.annotate.enrich.simplify <- clusterProfiler::simplify(High_vs_Low_IL32.sign.down.annotate.enrich) 
High_vs_Low_IL32.sign.down.annotate.enrich.df <- as.data.frame(High_vs_Low_IL32.sign.down.annotate.enrich.simplify)
High_vs_Low_IL32.sign.down.annotate.enrich.df$log10P <- -log10(High_vs_Low_IL32.sign.down.annotate.enrich.df$p.adjust)
High_vs_Low_IL32.sign.down.annotate.enrich.df.order <- High_vs_Low_IL32.sign.down.annotate.enrich.df[order(High_vs_Low_IL32.sign.down.annotate.enrich.df$log10P,decreasing = T),]

ggplot(High_vs_Low_IL32.sign.down.annotate.enrich.df.order.top, aes(x = reorder(Description,log10P), y = log10P)) + 
  geom_bar(stat="identity",fill="#005824") +
  theme_bw(base_size = 8) +
  theme(legend.key.size = unit(0.4, "cm")) +
  coord_flip() + 
  labs(y=expression(-Log['10'](Adj.~p-value)), x="")

######################## END ####################################



###################### Single cell data analysis GSE106218 #######################
library(Seurat)
require(ggplot2)


all.counts <-  read.table(file="GSE106218_GEO_processed_MM_raw_TPM_matrix.txt", header = T,sep="\t")

#Count matrix available for download:   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106218


all.counts.sce <- all.counts
all.counts.sce <- all.counts.sce[grep("ERCC",rownames(all.counts.sce),invert = T),]
all.counts.sce.names <- sapply(strsplit(colnames(all.counts.sce),"_"), `[`, 1)
sc.myeloma <- CreateSeuratObject(counts = all.counts.sce)
sc.myeloma <- NormalizeData(object = sc.myeloma)
sc.myeloma <- FindVariableFeatures(object = sc.myeloma)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc.myeloma), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sc.myeloma)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2
all.genes <- rownames(sc.myeloma)
sc.myeloma <- ScaleData(sc.myeloma, features = all.genes)
sc.myeloma <- RunPCA(sc.myeloma, features = VariableFeatures(object = sc.myeloma))


sc.myeloma <- FindNeighbors(sc.myeloma, dims = 1:10)
sc.myeloma <- FindClusters(sc.myeloma, resolution = 0.5)
head(Idents(sc.myeloma), 5)
sc.myeloma <- RunUMAP(sc.myeloma, dims = 1:10)

DimPlot(sc.myeloma, reduction = "umap", group.by = "orig.ident",pt.size = 0.01)  +
  scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#CCCC00","#b15928"),name = "Sample") + 
  theme_classic(base_size = 8) + theme(legend.key.size = unit(.2, "cm")) + labs(title="")

DimPlot(sc.myeloma, reduction = "umap",pt.size = 0.01) + 
  scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"),name="Cluster") +
  theme_classic(base_size = 8) +theme(legend.key.size = unit(.2, "cm"))

FeaturePlot(pbmc, features = c("IL32"),pt.size = 0.01) + labs(color = "IL32\nexpression",title = "") + 
  theme_classic(base_size = 8) + theme(legend.key.size = unit(.2, "cm"))


########## Gene ontology single cell ############

pbmc[18942,]@meta.data$nCount_RNA > 265.707
pbmc@meta.data$IL32_Group <- ifelse(test = pbmc[18942,]@meta.data$nCount_RNA > 265.707,yes="High",no = "low")
pbmc@meta.data$IL32_Expressed <- ifelse(test = pbmc[18942,]@meta.data$nCount_RNA > 0,yes="High",no = "low") #265.707

Idents(pbmc_groups) <- pbmc_groups@meta.data$IL32_Expressed
pbmc_groups@meta.data$IL32_Group_Patients <- paste(pbmc_groups@meta.data$orig.ident,pbmc_groups@meta.data$IL32_Group,sep="_")
Idents(pbmc_groups) <- pbmc_groups@meta.data$IL32_Group_Patients  

High_Low_IL32.MM17 <- FindMarkers(pbmc_groups, ident.1 = "MM17_High", ident.2 = "MM17_low", min.pct = 0)
head(High_Low_IL32.MM17)
High_Low_IL32.MM02EM <- FindMarkers(pbmc_groups, ident.1 = "MM02EM_High", ident.2 = "MM02EM_low", min.pct = 0)
head(High_Low_IL32.MM02EM)
High_Low_IL32.MM33 <- FindMarkers(pbmc_groups,clusters="MM33", ident.1 = "MM33_High", ident.2 = "MM33_low", min.pct = 0)
head(High_Low_IL32.MM33)


High_Low_IL32.MM17.up <- High_Low_IL32.MM17[High_Low_IL32.MM17$avg_log2FC>0.5,]
High_Low_IL32.MM02EM.up <- High_Low_IL32.MM02EM[High_Low_IL32.MM02EM$avg_log2FC>0.50,]
High_Low_IL32.MM33.up <- High_Low_IL32.MM33[High_Low_IL32.MM33$avg_log2FC>0.5,]

High_Low_IL32.rbind <- rbind(High_Low_IL32.MM17.up,High_Low_IL32.MM02EM.up,High_Low_IL32.MM33.up)
High_Low_IL32.rbind.rmDUP <- High_Low_IL32.rbind[!duplicated(High_Low_IL32.rbind$Gene),]

High_Low_IL32.up.annotate <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(High_Low_IL32.rbind.rmDUP),columns=c("ENTREZID"), keytype="SYMBOL") 
High_Low_IL32.up.universe <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(pbmc),columns=c("ENTREZID"), keytype="SYMBOL") 

High_Low_IL32.up.annotate.up.enrich <- enrichGO(gene= High_Low_IL32.up.annotate$ENTREZID,universe= High_Low_IL32.up.universe$ENTREZID,OrgDb         = org.Hs.eg.db,
                                                ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable      = TRUE)

dotplot(High_Low_IL32.up.annotate.up.enrich)
High_Low_IL32.up.annotate.up.enrich.df <- as.data.frame(High_Low_IL32.up.annotate.up.enrich)
High_Low_IL32.up.annotate.up.enrich.df$log10P <- -log10(High_Low_IL32.up.annotate.up.enrich.df$p.adjust)
High_Low_IL32.up.annotate.up.enrich.df.order <- High_Low_IL32.up.annotate.up.enrich.df[order(High_Low_IL32.up.annotate.up.enrich.df$log10P,decreasing = T),]
High_Low_IL32.up.annotate.up.enrich.df.order.top <- head(High_Low_IL32.up.annotate.up.enrich.df.order,20)
High_Low_IL32.up.annotate.up.enrich.df.order.top

require(tidyverse)

####### Gene ontology plot #######
ggplot(High_Low_IL32.up.annotate.up.enrich.df.order.top, aes(x = reorder(Description,log10P), y = log10P)) + 
  geom_bar(stat="identity",fill="#005824") +
  theme_bw(base_size = 8) +
  theme(legend.key.size = unit(0.4, "cm")) +
  coord_flip() + 
  labs(y=expression(-Log['10'](Adj.~p-value)), x="")

############################################# END #########################################################

