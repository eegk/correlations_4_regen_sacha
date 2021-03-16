#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
libs<-c("corrplot","ROTS","Hmisc","tcR","tidyr","limma","NMF","tsne","Rtsne" ,"UpSetR",
        "AnnotationDbi","RColorBrewer","papmap","sva","GO.db","fitdistrplus","ff","plyranges",
        "annotables","Rsamtools","GenomicFeatures","ggarrange","pheatmap","rms",
        "dplyr","DESeq2","ggplot2","ggrepel","edgeR","doParallel","ggradar","colorspace",
        "reshape2","rmarkdown","org.Hs.eg.db","treemapify",
        "variancePartition","ggpubr","factoextra","Mfuzz", "universalmotif","ggbio","GenomicRanges",
        "scales","tibble","RColorBrewer","tidyr","ggplot2","reshape2","circlize",
        "colorspace","Vennerable","enrichR","cowplot","data.table",
        "ROCR","mlbench","caretEnsemble","gbm","rpart","UpSetR","SCENIC","AUCell","RcisTarget","plyr",
        "tidyverse","hrbrthemes","fmsb","colormap","viridis","survminer","survival","ggalluvial")
lapply(libs, require, character.only = TRUE) ; rm(libs)
library(tidyverse)
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

### Initial analysis; Re-done downstream.
#############################################################################################################################################
setwd("/Users/gonzae34/Documents/projects_gnjatic/correlation_matrix_regeneron/")
my_correlation_matrix <- readxl::read_excel("/Users/gonzae34/Documents/projects_gnjatic/correlation_matrix_regeneron/For Edgar.xlsx",sheet = 1)
my_correlation_matrix <- as.data.frame(my_correlation_matrix)
my_correlation_matrix <- t(my_correlation_matrix) 

colnames(my_correlation_matrix) <- my_correlation_matrix[1,]
my_correlation_matrix <- my_correlation_matrix[-1,]

my_correlation_matrix <- as.data.frame(my_correlation_matrix)
dim(my_correlation_matrix)

my_correlation_matrix[is.na(my_correlation_matrix)] <- "NA"
table(is.na(my_correlation_matrix))

my_correlation_matrix <- my_correlation_matrix[,-2]
test <- melt(my_correlation_matrix,id="x")
test$value <- as.numeric(test$value)

table(is.na(test$value)) ### 2067 1533
#############################################################################################################################################

### Correlations between markers
#############################################################################################################################################
# P1_aSMA has no values whatsoever, Thus is removed from analysis.
Markers <- levels(as.factor(as.character(test$x)))
Markersv2 <- Markers[-1]

### Loop to correlate samples between markers
rm(my_storage) ; my_storage <- list() ; count=1

for (xxi in 1:length(Markersv2) ){
  
  for( xxy in 1:length(Markersv2)){
    
    my_test <- cor.test ( test$value[test$x %in% Markersv2[xxi]], test$value[test$x %in% Markersv2[xxy]], method="spearman")
      
    my_storage[[count]] <- data.frame( p.value = my_test$p.value , S = my_test$statistic , Rho = my_test$estimate ,
    
                                                                          A =  Markersv2[xxi] , 
                                       
                                                                          B =  Markersv2[xxy] )
    
    count=count+1  } }
### Store results into a data frame class object
my_storage <- do.call(rbind,my_storage)
my_storage$adj.p.value <- p.adjust(my_storage$p.value)
my_storage$nLogFDR<- -log10(my_storage$adj.p.value)
my_storage$nLogFDR[my_storage$nLogFDR> 5] <- 5

### For a figure the most interesting is to sort the data into a matrix.
### organize data for hierarchical clustering of correlation values
out <- pivot_wider(data=my_storage[,c("A","B","Rho")], names_from=A, values_from=Rho) # %>% unnest() 
out <- as.data.frame(out)
rownames(out) <- out$B
colnames(out)
out <- pheatmap(as.data.frame(out[,-1]), angle_col = 90)
dev.off()
dev.off()
### Copy table and organize based on hierarchical clustering
testXYY <- my_storage
testXYY$A <- factor(testXYY$A, levels=out$tree_row$labels[out$tree_row$order])
testXYY$B <- factor(testXYY$B, levels=out$tree_col$labels[out$tree_col$order])
### Select only significant values for plot
ix <- which(testXYY$nLogFDR < 1.3)
###
pdf(file="/Users/gonzae34/Downloads/figure_for_sacha_regeneron.pdf",width = 8, height = 6)
ggplot(testXYY[,], aes(x=A, y=B, color= Rho, size=nLogFDR )) + 
  geom_point(shape=16, colour = "black", aes(size = max(nLogFDR))) +
  geom_point(shape=16, colour = "white", aes(size = 0.8*max(nLogFDR))) +
  geom_point(shape=16, aes(size = 0.81*nLogFDR)) + 
  geom_text(aes(label=round(Rho,1)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='Rho',limits=c(-1,1)) + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_classic2() + rotate_x_text(angle=45) + 
  labs(x ='', y='', title='') 
dev.off()
###

#############################################################################################################################################

out_rho <- pivot_wider(data=my_storage[,c("A","B","Rho")], names_from=A, values_from=Rho)
out_rho <- as.data.frame(out_rho)
write.csv(file="out_rho.csv",out_rho)


out_p.value <- pivot_wider(data=my_storage[,c("A","B","p.value")], names_from=A, values_from=p.value) 
out_p.value <- as.data.frame(out_p.value)
write.csv(file="out_p.value.csv",out_p.value)

out_FDR <- pivot_wider(data=my_storage[,c("A","B","adj.p.value")], names_from=A, values_from=adj.p.value) 
out_FDR <- as.data.frame(out_FDR)
write.csv(file="out_fdr.csv",out_FDR)

###

Markers
### ex0
cor.test ( test$value[which(test$x %in% "P2_Trem2")], test$value[which(test$x %in% "P1_CD3")], method="spearman")
cor.test ( test$value[test$x %in% "P2_Trem2"], test$value[test$x %in% "P1_CD3"], method="pearson")
summary(lm(test$value[which(test$x %in% "P2_Trem2")] ~ test$value[which(test$x %in% "P1_CD3")]))

ggplot() + aes(x=test$value[which(test$x %in% "P2_Trem2")], y=test$value[which(test$x %in% "P1_CD3")]) + geom_point() + theme_bw() +
  labs(x="P2_Trem2",y="P1_CD3") + geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 2) + 
  stat_smooth( method = 'lm', color = 'purple')
### ex1
cor.test ( test$value[which(test$x %in% "P2_Trem2")], test$value[which(test$x %in% "P1_CD68")], method="spearman")
cor.test ( test$value[test$x %in% "P2_Trem2"], test$value[test$x %in% "P1_CD68"], method="pearson")
summary(lm(test$value[which(test$x %in% "P2_Trem2")] ~ test$value[which(test$x %in% "P1_CD68")]))

ggplot() + aes(x=test$value[which(test$x %in% "P2_Trem2")], y=test$value[which(test$x %in% "P1_CD68")]) + geom_point() + theme_bw() +
  labs(x="P2_Trem2",y="P1_CD68") + geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 2) + 
  stat_smooth( method = 'lm', color = 'purple') 
### ex2
cor.test ( test$value[which(test$x %in% "P2_Trem2")], test$value[which(test$x %in% "P2_CD1c")], method="spearman")
cor.test ( test$value[test$x %in% "P2_Trem2"], test$value[test$x %in% "P2_CD1c"], method="pearson")
summary(lm(test$value[which(test$x %in% "P2_Trem2")] ~ test$value[which(test$x %in% "P2_CD1c")]))

ggplot() + aes(x=test$value[which(test$x %in% "P2_Trem2")], y=test$value[which(test$x %in% "P2_CD1c")]) + geom_point() + theme_bw() +
  labs(x="P2_Trem2",y="P2_CD1c") + geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 2) + 
  stat_smooth( method = 'lm', color = 'purple') 
### ex3
cor.test ( test$value[which(test$x %in% "P2_DC-LAMP")], test$value[which(test$x %in% "P1_CD68")], method="spearman")
cor.test ( test$value[(test$x %in% "P2_DC-LAMP")], test$value[(test$x %in% "P1_CD68")], method="spearman")
cor.test ( test$value[test$x %in% "P2_DC-LAMP"], test$value[test$x %in% "P1_CD68"], method="pearson")
summary(lm(test$value[which(test$x %in% "P2_DC-LAMP")] ~ test$value[which(test$x %in% "P1_CD68")]))

ggplot() + aes(x=test$value[which(test$x %in% "P2_DC-LAMP")], y=test$value[which(test$x %in% "P1_CD68")]) + geom_point() + theme_bw() +
  labs(x="P2_DC-LAMP",y="P1_CD68") + geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 2) + 
  stat_smooth( method = 'lm', color = 'purple')

### ex4
cor.test ( test$value[which(test$x %in% "P1_NKp46")], test$value[which(test$x %in% "P2_CD3")], method="spearman")
cor.test ( test$value[(test$x %in% "P1_NKp46")], test$value[(test$x %in% "P2_CD3")], method="spearman")
cor.test ( test$value[test$x %in% "P1_NKp46"], test$value[test$x %in% "P2_CD3"], method="pearson")
summary(lm(test$value[which(test$x %in% "P1_NKp46")] ~ test$value[which(test$x %in% "P2_CD3")]))

ggplot() + aes(x=test$value[which(test$x %in% "P1_NKp46")], y=test$value[which(test$x %in% "P2_CD3")]) + geom_point() + theme_bw() +
  labs(x="P1_NKp46",y="P2_CD3") + geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 2) + 
  stat_smooth( method = 'lm', color = 'purple')

cor( test$value[which(test$x %in% "P1_NKp46")], test$value[which(test$x %in% "P2_CD3")], method="spearman")

#DC-LAMP / CD68 and TREM2 / CD1c are anticorrelated at -0.8 with Prism’s calculations, while you are finding only -0.5 and -0.3 respectively.
#Conversely, between P1_NKp46 and P2_CD3, you have -0.4 with very strong p value, while I have -0.1 with a p value of 0.6…
sort(test$value[which(test$x %in% "P1_CD3")])
ix <- which(!(test$value %in% "NA"))
testx <- test[ix,]
ix <- which(!is.na(testx$value ))
testx <- test[ix,]
cor.test( test$value[which(test$x %in% "P1_NKp46")], test$value[which(test$x %in% "P2_CD3")], method="spearman")
cor.test ( sample(testx$value[which(testx$x %in% "P1_NKp46")],119), sample(testx$value[which(testx$x %in% "P2_CD3")],119), method="spearman")
cor.test ( testx$value[which(testx$x %in% "P1_NKp46")], sample(testx$value[which(testx$x %in% "P2_CD3")],119), method="spearman")
cor( test$value[which(test$x %in% "P1_NKp46")], test$value[which(test$x %in% "P2_CD3")], method="spearman",use="pairwise.complete.obs")
test$value[which(test$x %in% "P1_NKp46" & (test$x %in% "NA"))]

########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

### re-analysis with more details and updated data.
########################################################################################################################################################

### load data
my_correlation_matrix <- readxl::read_excel("/Users/gonzae34/Documents/projects_gnjatic/correlation_matrix_regeneron/For Edgar_sites_b.xlsx",sheet = 1)

### organize data
my_correlation_matrix <- as.data.frame(my_correlation_matrix)
my_correlation_matrix <- t(my_correlation_matrix) 

colnames(my_correlation_matrix) <- my_correlation_matrix[1,]

my_metadata <- my_correlation_matrix[1:6,]
my_metadata <- my_metadata[,-c(1,2)]
my_metadata <- my_metadata[-c(6),]
my_metadata <- as.data.frame(my_metadata)
my_metadata <- t(my_metadata)
my_metadata <- as.data.frame(my_metadata)
colnames(my_metadata)<-c("PatientID","code1","code2","code3","Class")

my_correlation_matrix <- my_correlation_matrix[-c(1:6),]

my_correlation_matrix <- as.data.frame(my_correlation_matrix)
dim(my_correlation_matrix)

my_correlation_matrix[is.na(my_correlation_matrix)] <- "NA"
table(is.na(my_correlation_matrix))

my_correlation_matrix <- as.data.frame(my_correlation_matrix)
my_correlation_matrix$a <- paste(my_correlation_matrix$a,my_correlation_matrix$h,sep="___")
my_correlation_matrix <- my_correlation_matrix[,-2]

test <- melt(my_correlation_matrix,id="a")
test$value <- as.numeric(test$value)

table(is.na(test$value)) ### 2067 1533

colnames(test)[1]<-"x"
test$markers_roi <- test$x
test$x <- tidyr::separate(data.frame( test$x ), 1, sep="___", c("a","b"))$a

### Prepare values and correlate
Markers <- levels(as.factor(as.character(test$x)))
Markersv2 <- Markers[-1]

### loop for correlations
rm(my_storage) ; my_storage <- list() ; count=1

for (xxi in 1:length(Markersv2) ){
  for( xxy in 1:length(Markersv2)){
    my_test <- cor.test ( test$value[test$x %in% Markersv2[xxi]], test$value[test$x %in% Markersv2[xxy]], method="spearman")
    my_storage[[count]] <- data.frame( p.value = my_test$p.value , S = my_test$statistic , Rho = my_test$estimate ,
                                       A =  Markersv2[xxi] , 
                                       B =  Markersv2[xxy] )
    count=count+1  } }

### Store results
my_storage <- do.call(rbind,my_storage)
my_storage$adj.p.value <- p.adjust(my_storage$p.value)
my_storage$nLogFDR<- -log10(my_storage$adj.p.value)
my_storage$nLogFDR[my_storage$nLogFDR> 5] <- 5

### Prepare dendrograms
out_mat <- pivot_wider(data=my_storage[,c("A","B","Rho")], names_from=A, values_from=Rho) # %>% unnest() 
out_mat <- as.data.frame(out_mat)
rownames(out_mat) <- out_mat$B
colnames(out_mat)

### produce the figure
head(my_metadata)
out <- pheatmap(as.data.frame(out_mat[,-1]),angle_col = 90)
dev.off()
dev.off()

### extract clusters
marker_clustering <- data.frame( cluster = cutree(out$tree_row, k=3), name = names(cutree(out$tree_row, k=3)))
### store clusters
write.csv(file="marker_clustering.csv",marker_clustering)

### store figure
pdf(file="markers_clustering_dendrogram.pdf", width = 5, height = 6)
plot(out$tree_row)
abline(h=1.9, col="red", lty=2, lwd=2)
dev.off()
###

library(ggdendro)

### alternative plot
### https://stackoverflow.com/questions/43794870/plotting-a-clustered-heatmap-with-dendrograms-using-rs-plotly
### dx <- dendro_data(as.dendrogram(hclust(dist(t(as.matrix(scale(out_mat[,-1]))))))) ; rm(dx)

dendro_data_markers <- dendro_data(out$tree_row)
# helper function for creating dendograms
ggdend <- function(df) {
  ggplot() + geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
             labs(x = "", y = "") + theme_minimal() +
             theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) }

dx <- ggdend(dendro_data_markers$segments) + coord_flip() + scale_y_reverse()

###
testXYY <- my_storage
testXYY$A <- factor(testXYY$A, levels=out$tree_row$labels[out$tree_row$order])
testXYY$B <- factor(testXYY$B, levels=out$tree_col$labels[out$tree_col$order])
###

###
ix <- which(testXYY$nLogFDR < 1.3)
###
pdf(file="figure1_Spearman_rho_Markers_across_ROIs.pdf",width = 11, height = 5)

ggarrange(
dx + theme(plot.margin=unit(c(t=1,r=-5,b=1,l=1),"pt")) ,

ggplot(testXYY[,], aes(x=A, y=B, color= Rho, size=nLogFDR )) + 
  geom_point(shape=16, colour = "black", aes(size = max(nLogFDR))) +
  geom_point(shape=16, colour = "white", aes(size = 0.8*max(nLogFDR))) +
  geom_point(shape=16, aes(size = 0.81*nLogFDR)) + 
  geom_text(aes(label=round(Rho,1)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='Rho',limits=c(-1,1)) + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_classic2() + rotate_x_text(angle=45) + 
  labs(x ='', y='', title='') + 
  theme(plot.margin=unit(c(t=1,r=1,b=1,l=-10),"pt"))
, ncol=2, nrow=1, align = 'hv' )

dev.off()
########################################################################################################################################################


### Correlation between ROIs
########################################################################################################################################################
Markers <- levels(as.factor(as.character(test$markers_roi)))
Markersv2 <- Markers
Markersv2 <- Markersv2[-grep("aSMA",Markersv2)]

rm(my_storage_mkr) ; my_storage_mkr <- list() ; count=1

for (xxi in 1:length(Markersv2) ){
  for( xxy in 1:length(Markersv2)){
    my_test <- cor.test ( test$value[test$markers_roi %in% Markersv2[xxi]], test$value[test$markers_roi %in% Markersv2[xxy]], method="spearman")
    my_storage_mkr[[count]] <- data.frame( p.value = my_test$p.value , S = my_test$statistic , Rho = my_test$estimate ,
                                       A =  Markersv2[xxi] , 
                                       B =  Markersv2[xxy] )
    
    count=count+1  } }

###
my_storage_mkr <- do.call(rbind,my_storage_mkr)
my_storage_mkr$adj.p.value <- p.adjust(my_storage_mkr$p.value)
my_storage_mkr$nLogFDR<- -log10(my_storage_mkr$adj.p.value)
my_storage_mkr$nLogFDR[my_storage_mkr$nLogFDR> 5] <- 5
###
out_mat <- pivot_wider(data=my_storage_mkr[,c("A","B","Rho")], names_from=A, values_from=Rho) # %>% unnest() 
out_mat <- as.data.frame(out_mat)
rownames(out_mat) <- out_mat$B
colnames(out_mat)
###
write.csv(file="marker_correlation_per_roi.csv",out_mat)
###
pdf(file="figure2_Spearman_rho_Markers_across_ROIs.pdf",width = 26,height = 24)
out <- pheatmap(as.data.frame(out_mat[,-1]),angle_col = 90)
dev.off()
dev.off()
###
#marker_clustering_per_roi <- data.frame( cluster = as.factor(cutree(out$tree_row, k=6)) )
#out <- pheatmap(as.data.frame(out_mat[,-1]),angle_col = 90, annotation_col =  marker_clustering_per_roi)
###
marker_clustering_per_roi <- data.frame( cluster = cutree(out$tree_row, k=6), name = names(cutree(out$tree_row, k=6)))
###
write.csv(file="marker_clustering_per_roi.csv",marker_clustering_per_roi)
###
pdf(file="markers_clustering_per_roi.pdf",width = 26,height = 8)
plot(out$tree_row)
#abline(h=1.9, col="red", lty=2, lwd=2)
dev.off()
###
#out$tree_row$order
#out$tree_row$order

### too much for ggplot to handle something comprehensive
#testXYY <- my_storage_mkr
#testXYY$A <- factor(testXYY$A, levels=out$tree_row$labels[out$tree_row$order])
#testXYY$B <- factor(testXYY$B, levels=out$tree_col$labels[out$tree_col$order])
###
#testXYY$nLogFDR<- -log10(testXYY$adj.p.value)
#testXYY$nLogFDR[testXYY$nLogFDR> 5] <- 5
###
#ix <- which(testXYY$nLogFDR > 1.3)
###
#pdf(file="figure_for_sacha_regeneron_per_roi_v2.pdf",width = 24, height = 22)
#ggplot(testXYY[ix,], aes(x=A, y=B, color= Rho, size=nLogFDR )) + 
#  geom_point(shape=16, colour = "black", aes(size = max(nLogFDR))) +
#  geom_point(shape=16, colour = "white", aes(size = 0.8*max(nLogFDR))) +
#  geom_point(shape=16, aes(size = 0.81*nLogFDR)) + 
  #geom_point(data=testXYY[-ix,],shape=16, aes(size = 0.81*nLogFDR)) + 
  #geom_point(data=testXYY[ix,],shape=4, aes(size = 0.81*nLogFDR)) + 
#  geom_text(aes(label=round(Rho,1)),size=2,color="black") +
#  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='Rho',limits=c(-1,1)) + #,breaks = c(-2.8, 0, 1, 5)) +
#  theme_classic2() + rotate_x_text(angle=45) + 
#  labs(x ='', y='', title='') 
#dev.off()

###
########################################################################################################################################################
###

### Organize data for multiple testing
########################################################################################################################################################
dim(test) ### flat version of my_correlation_matrix
dim(my_correlation_matrix) # table from excel with data per ROI

### inspect
### my_correlation_matrix[1:5,1:5]

### set rownames based on marker x ROI
rownames(my_correlation_matrix) <- my_correlation_matrix[,1]
### remove column with names
my_correlation_matrix <- my_correlation_matrix[,-1]


### Convert the matrix to a numeric matrix
data <- my_correlation_matrix %>% mutate_at(c(1:20), as.numeric)

### inspect
dim(data)
table(is.na(data)) ### This looks correct 2067, 1533

###
########################################################################################################################################################
###

### identify clusters with clara algorithm using 50 sample boostrap
clara.res <- clara(t(data), 6, samples = 50, pamLike = TRUE, correct.d=TRUE)
### summary
print(clara.res)
### extract cluster assignment
dd <- cbind(t(data), cluster = clara.res$cluster)
head(dd, n = 4)
dd <- as.data.frame(dd)
table(dd$cluster) ###

### data with missing value
pdf(file="raw_data_including_missing_data.pdf",width = 10,height = 22)
pheatmap(data[,order(clara.res$cluster)], cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 90)
dev.off()

### Clean space
rm(out,out_mat,dx,dd,clara.res,testXYY,my_test)

### copy original
data_with_NAs <- data

### replace to ceros
data[is.na(data)] <- 0
### remove aSMA
data <- data[-grep("aSMA",rownames(data)),]


### PCA analysis
########################################################################################################################################################

### PCA
mydata <- prcomp( t(data) , scale=TRUE, center = TRUE )

### SCREE
Scree.plot <- fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw()
eigen <- get_eigenvalue(mydata)
pca1 <- format(round(eigen$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(eigen$variance.percent[2], 2), nsmall = 2) ;Scree.plot ; rm(eigen)

my_metadata$PatientID
prcomp_obj_mssm_pca <- as.data.frame(mydata$x)
prcomp_obj_mssm_pca$PatientID <- rownames(prcomp_obj_mssm_pca)

prcomp_obj_mssm_pca <- merge(prcomp_obj_mssm_pca,my_metadata,by="PatientID")
###
pdf(file="pca_CEROS_filled.pdf",width = 8,height = 6)
ggplot(data=prcomp_obj_mssm_pca, aes(x=PC1, y=PC2, color = Class)) + 
  geom_point() + 
  geom_label_repel(aes(label=PatientID),force = 10, max.iter = 5000,show.legend = FALSE) +
  theme_classic2() + 
  labs(x=paste("PC1:",pca1,"%",sep=""),title="NA=Cero | 20 Markers", y=paste("PC1:", pca2,"%", sep=""), color="Fluid / Origin") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(4.0)),
        plot.title=element_text(size=rel(4.8)) ) + 
  theme(legend.position="bottom") + 
  guides(color=guide_legend(ncol=2), shape = guide_legend(override.aes = list(size = rel(4.0)))) 
dev.off()
###


### PCA
#data <- data[-grep("aSMA",rownames(data)),]
mydata <- prcomp( t(na.omit(data_with_NAs)) , scale=TRUE, center = TRUE )

### SCREE
Scree.plot <- fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw()
eigen <- get_eigenvalue(mydata)
pca1 <- format(round(eigen$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(eigen$variance.percent[2], 2), nsmall = 2) ;Scree.plot ; rm(eigen)

my_metadata$PatientID
prcomp_obj_mssm_pca <- as.data.frame(mydata$x)
prcomp_obj_mssm_pca$PatientID <- rownames(prcomp_obj_mssm_pca)

prcomp_obj_mssm_pca <- merge(prcomp_obj_mssm_pca,my_metadata,by="PatientID")
###
pdf(file="pca_without_filled.pdf",width = 8,height = 6)
ggplot(data=prcomp_obj_mssm_pca, aes(x=PC1, y=PC2, color = Class)) + 
  geom_point() + 
  geom_label_repel(aes(label=PatientID),force = 10, max.iter = 5000,show.legend = FALSE) +
  theme_classic2() + 
  labs(x=paste("PC1:",pca1,"%",sep=""),title="NA.Omit | 6 Markers", y=paste("PC1:", pca2,"%", sep=""), color="Fluid / Origin") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(4.0)),
        plot.title=element_text(size=rel(4.8)) ) + 
  theme(legend.position="bottom") + 
  guides(color=guide_legend(ncol=2), shape = guide_legend(override.aes = list(size = rel(4.0)))) 
dev.off()
###
########################################################################################################################################################


### linear model
########################################################################################################################################################
###
my_metadata$Class <- gsub("-","",my_metadata$Class)
###
identical(my_metadata$PatientID, colnames(data))
###

### Design matrix for the model
design <- model.matrix( ~ 0 + Class, data = my_metadata )
### names
colnames(design) <- gsub("Class","",colnames(design))
### contrast for comparisons
contr.matrix <- makeContrasts("C vs LN" = LiverC - LiverLN, 
                              "C vs N" = LiverC - LiverN, 
                              "C vs T" = LiverC - LiverT, 
                              "LN vs N" = LiverLN - LiverN, 
                              "LN vs T" = LiverLN - LiverT, 
                              "N vs T" = LiverN - LiverT, 
                              levels = colnames(design) )

vfit <- lmFit(data, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

### Store data
lmfreq_results <- list()
for ( i in 1:ncol(summary(decideTests(efit))) ) { lmfreq_results[[i]] <- topTable(efit, coef=i,n=Inf,adjust.method="fdr") 
lmfreq_results[[i]]$Comparison <- colnames(summary(decideTests(efit)))[i]
lmfreq_results[[i]]$Marker_ROI <- rownames( topTable(efit, coef=i,n=Inf,adjust.method="fdr") ) }
regeneron_dm_res <- do.call(rbind,lmfreq_results)
rownames(regeneron_dm_res) <- NULL

regeneron_dm_res$nLogFDR <- -log10(regeneron_dm_res$adj.P.Val)
table(regeneron_dm_res$adj.P.Val<0.05)


### filtering
regeneron_dm_res$Marker <- tidyr::separate(data.frame( regeneron_dm_res$Marker_ROI ), 1, sep="___", c("a","b"))$a

for_figure <- regeneron_dm_res[which(regeneron_dm_res$adj.P.Val<0.05),]
for_figure <- for_figure[,c(1,7,8,9,10)]
for_figure <- for_figure[order(for_figure$Marker_ROI),]

for_figure$condition <- paste(for_figure$Comparison,for_figure$Marker_ROI,sep="___")
for_figure <- for_figure[,c(1,4,6)]

### this means results per marker using Marker instead of Marker ROI
#for_figure2<- cbind( aggregate(for_figure$logFC~for_figure$condition, FUN=mean),
#                     aggregate(for_figure$nLogFDR~for_figure$condition, FUN=mean) )
#for_figure2 <- for_figure2[,-3]
#colnames(for_figure2) <- c("MX","logFC","nLogFDR")

for_figure2 <- for_figure
colnames(for_figure2) <- c("logFC","nLogFDR","MX")

for_figure2$Comparison <- tidyr::separate(data.frame( for_figure2$MX ), 1, sep="___", c("a","b"))$a
#for_figure2$Marker <- tidyr::separate(data.frame( for_figure2$MX ), 1, sep="___", c("a","b"))$b
for_figure2$Marker <- paste(tidyr::separate(data.frame( for_figure2$MX ), 1, sep="___", c("a","b"))$b,tidyr::separate(data.frame( for_figure2$MX ), 1, sep="___", c("a","b","c"))$c,sep="___")

###
out <- pivot_wider(data=for_figure2[,c("Marker","Comparison","logFC")], names_from=Marker, values_from=logFC) # %>% unnest() 
out <- as.data.frame(out)
rownames(out) <- out$Comparison
colnames(out)
out[is.na(out)]<-0
out <- pheatmap(as.data.frame(out[,-1]))
dev.off()
dev.off()
###
testXYY <- for_figure2
testXYY$Comparison <- factor(testXYY$Comparison, levels=out$tree_row$labels[out$tree_row$order])
testXYY$Marker<- factor(testXYY$Marker, levels=out$tree_col$labels[out$tree_col$order])

pdf(file="Markers_between_groups_moderate_t.pdf",width = 15,height = 4)
ggplot(testXYY, aes(x=Marker, y=Comparison, color=logFC, size=nLogFDR )) + 
  geom_point(shape=16, colour = "black", aes(size = max(nLogFDR))) +
  geom_point(shape=16, colour = "white", aes(size = 0.8*max(nLogFDR))) +
  geom_point(shape=16, aes(size = 0.81*nLogFDR)) + 
  #geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') + #,,limits=c(-1,1), breaks = c(-2.8, 0, 1, 5)) +
  theme_classic2() + rotate_x_text(angle=45) + 
  labs(x ='', y='', title='') 
dev.off()
dev.off()


#### melt data for figures violin plots
########################################################################################################################################################
data_melted <- data
data_melted$Markers <- rownames(data_melted)
data_melted <- (melt(data_melted))

colnames(data_melted)[2] <- "PatientID"
data_melted <- merge(data_melted,my_metadata,by="PatientID")
data_melted$Marker <- tidyr::separate(data.frame( data_melted$Markers ), 1, sep="___", c("a","b"))$a

head(data_melted)
table(data_melted$Class)
head(regeneron_dm_res)
table(regeneron_dm_res$Comparison)

### Comparisons
my_comparisons <- list( c("LiverC", "LiverLN"),
                        c("LiverC", "LiverN"),
                        c("LiverC", "LiverT"),
                        c("LiverN","LiverLN"),
                        c("LiverN","LiverT"),
                        c("LiverT","LiverLN"))
### Marker
ix <- which(data_melted$Marker %in% "P1_CD68" & data_melted$value > 0)

pdf(file="mixed_linear_model_outlines_4_P1_CD68_ROIs.pdf", width = 18, height = 5)
ggboxplot(data=data_melted[ix,], x="Class" , y="value",color='Class') +
    facet_wrap(~Markers,ncol=10, scales = 'free') +
    theme_bw() +  #geom_quasirandom() + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    #scale_fill_manual(values = c("firebrick","steelblue")) +
    #scale_color_manual(values = c("firebrick","firebrick","steelblue","steelblue","steelblue")) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE ) +
    theme(legend.position = "none") + 
    labs(x ='', y='Value', title='',color="",fill="") + theme(text=element_text(size=rel(5.5)),strip.text=element_text(size=rel(1.5),face="bold"))
dev.off()
    
library(ggbeeswarm)

pdf(file="mixed_linear_model_outlines_4_P1_CD68_Marker.pdf", width = 5, height=6)    
ggboxplot(data=data_melted[ix,], x="Class" , y="value",color='Class') +
      theme_bw() + geom_quasirandom(size=1,alpha=0.5) + 
      theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
      stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE ) +
      theme(legend.position = "none") + 
      labs(x ='', y='P1_CD68 Value', title='', color="", fill="") + theme(text=element_text(size=rel(5.5)),strip.text=element_text(size=rel(1.5),face="bold"))
dev.off()


###
########################################################################################################################################################
###

### Correlate full sample profiles
########################################################################################################################################################
Markers <- levels(as.factor(as.character(test$variable)))
Markersv2 <- Markers

rm(my_storage_pts) ; my_storage_pts <- list() ; count=1

for (xxi in 1:length(Markersv2) ){
  for( xxy in 1:length(Markersv2)){
    my_test <- cor.test ( test$value[test$variable %in% Markersv2[xxi]], test$value[test$variable %in% Markersv2[xxy]], method="spearman")
    my_storage_pts[[count]] <- data.frame( p.value = my_test$p.value , S = my_test$statistic , Rho = my_test$estimate ,
                                           A =  Markersv2[xxi] , 
                                           B =  Markersv2[xxy] )
    
    count=count+1  } }
###
my_storage_pts <- do.call(rbind,my_storage_pts)
my_storage_pts$adj.p.value <- p.adjust(my_storage_pts$p.value)
my_storage_pts$nLogFDR<- -log10(my_storage_pts$adj.p.value)
my_storage_pts$nLogFDR[my_storage_pts$nLogFDR> 5] <- 5
###
out_mat_pt <- pivot_wider(data=my_storage_pts[,c("A","B","Rho")], names_from=A, values_from=Rho) # %>% unnest() 
out_mat_pt <- as.data.frame(out_mat_pt)
rownames(out_mat_pt) <- out_mat_pt$B
colnames(out_mat_pt)
###
#write.csv(file="correlation_per_pts.csv",out_mat_pt)
###
#pdf(file="figure_for_sacha_regeneron_per_pts_v1.pdf",width = 6, height = 6)
out <- pheatmap(as.data.frame(out_mat_pt[,-1]),angle_col = 90)
dev.off()
dev.off()
###
ann_col <- data.frame(Class=my_metadata$Class)
rownames(ann_col) <- my_metadata$PatientID
###
pdf(file="figure_for_sacha_regeneron_per_pts_v3.pdf",width = 7, height = 5)
out <- pheatmap(as.data.frame(out_mat_pt[,-1]),angle_col = 90,annotation_col = ann_col,annotation_row = ann_col)
dev.off()
dev.off()

### too ugly
#testXYY <- my_storage_pts
#testXYY$A <- factor(testXYY$A, levels=out$tree_row$labels[out$tree_row$order])
#testXYY$B <- factor(testXYY$B, levels=out$tree_col$labels[out$tree_col$order])
###
#testXYY$nLogFDR<- -log10(testXYY$adj.p.value)
#testXYY$nLogFDR[testXYY$nLogFDR> 5] <- 5
###
#ix <- which(testXYY$nLogFDR < 1.3)
###
#pdf(file="figure_for_sacha_regeneron_per_pts_v2.pdf",width = 8, height = 6)
#ggplot(testXYY[,], aes(x=A, y=B, color= Rho, size=nLogFDR )) + 
#  geom_point(shape=16, colour = "black", aes(size = max(nLogFDR))) +
#  geom_point(shape=16, colour = "white", aes(size = 0.8*max(nLogFDR))) +
#  geom_point(shape=16, aes(size = 0.81*nLogFDR)) + 
#  geom_text(aes(label=round(Rho,1)),size=2,color="black") +
#  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='Rho',limits=c(-1,1)) + #,breaks = c(-2.8, 0, 1, 5)) +
#  theme_classic2() + rotate_x_text(angle=45) + 
#  labs(x ='', y='', title='') 
#dev.off()
###

###
########################################################################################################################################################
###

### my_correlation_matrix[1:5,1:5]
### data[1:5,1:5]
### test[1:5,]

test$ROI <- tidyr::separate(data.frame( test$markers_roi ), 1, sep="___", c("a","b"))$b
test$Marker <- tidyr::separate(data.frame( test$markers_roi ), 1, sep="___", c("a","b"))$a

colnames(test)[2]<- "PatientID"

my_extended_data_long <- merge(test, my_metadata, by="PatientID")

my_comparisons <- list( c("LiverC", "LiverLN"),c("LiverC", "LiverN"),
                        c("LiverC", "LiverT"), c("LiverLN", "LiverN"), 
                        c("LiverLN", "LiverT"),c("LiverN", "LiverT") )

library(ggbeeswarm)
pdf(file="figure_for_sacha_regeneron_boxplots_marker_scale_free.pdf",width = 16, height = 8)
ggboxplot(data=my_extended_data_long, x="Class",y="value") + theme_minimal() +
  facet_wrap(~Marker,ncol=9,scales = 'free') +
  geom_quasirandom(size=1)+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) 
dev.off()

pdf(file="figure_for_sacha_regeneron_boxplots_marker_scale_fixed.pdf",width = 16, height = 8)
ggboxplot(data=my_extended_data_long, x="Class",y="value") + theme_minimal() +
  facet_wrap(~Marker,ncol=9) +
  geom_quasirandom(size=1)+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) 
dev.off()

my_extended_data_long$ROI <- tidyr::separate(data.frame( my_extended_data_long$markers_roi ), 1, sep="___", c("a","b"))$b

pdf(file="figure_for_sacha_regeneron_boxplots_marker_scale_free_ROI.pdf",width = 16, height = 70)
ggboxplot(data=my_extended_data_long, x="Class",y="value") + theme_bw() +
  facet_wrap(~Marker+ROI,ncol=9,scales = 'free') +
  geom_quasirandom(size=1)+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1))  
dev.off()

### Selected correlations
########################################################################################################################################################
###

head(test)







#"P2_CD66b" & "P2_CLEC9A"
########################################################################################################################################################
ix <- which(test$Marker %in% c("P2_CD66b","P2_CLEC9A"))

out <- as.data.frame(pivot_wider(data=test[ix,c(2,3,6)], names_from=Marker, values_from=value)  %>% unnest())
cor.test(out$P2_CD66b,out$P2_CLEC9A,method = "spear")


pdf(file="anticorrelation_between_markers_case1.pdf", width = 6, height = 5)
ggplot(data=out, aes(x=P2_CD66b, y=P2_CLEC9A)) + 
  geom_point(size=1) + ylim(0,max(out$P2_CLEC9A)) +
  stat_density2d(alpha=0.5) + 
  geom_rug(alpha=0.01,color="black")  +
  geom_abline(intercept = 0, slope = 1, color='firebrick', linetype = 2,size=1) + 
  stat_smooth(method = 'lm', color = 'steelblue') + 
  theme_bw() +
  labs(x="P2_CD66b", y="P2_CLEC9A", title="") +
  annotate(geom="text", x=4.2, y=2, label="Spearman Rho: -0.47\np.value:4.214e-08",color="steelblue",size=6) +
  theme(text=element_text(size=21))
dev.off()


#"P2_CD66b" & "P1_FOXP3"
########################################################################################################################################################
ix <- which(test$Marker %in% c("P2_CD66b","P1_FOXP3"))

out <- as.data.frame(pivot_wider(data=test[ix,c(2,3,6)], names_from=Marker, values_from=value)  %>% unnest())
cor.test(out$P2_CD66b,out$P1_FOXP3,method = "spear")


pdf(file="anticorrelation_between_markers_case2.pdf", width = 6, height = 5)
ggplot(data=out, aes(x=P2_CD66b, y=P1_FOXP3)) + 
  geom_point(size=1) + ylim(0,max(out$P1_FOXP3)) +
  stat_density2d(alpha=0.5) + 
  geom_rug(alpha=0.01,color="black")  +
  geom_abline(intercept = 0, slope = 1, color='firebrick', linetype = 2,size=1) + 
  stat_smooth(method = 'lm', color = 'steelblue') + 
  theme_bw() +
  labs(x="P2_CD66b", y="P1_FOXP3", title="") +
  annotate(geom="text", x=3, y=6, label="Spearman Rho: -0.45\np.value:2.313e-08",color="steelblue",size=6) +
  theme(text=element_text(size=21))
dev.off()

#"P1_CD68" & "P2_DC-LAMP"
########################################################################################################################################################
ix <- which(test$Marker %in% c("P1_CD68","P2_DC-LAMP"))

out <- as.data.frame(pivot_wider(data=test[ix,c(2,3,6)], names_from=Marker, values_from=value)  %>% unnest())
cor.test(out$P1_CD68,out$`P2_DC-LAMP`,method = "spear")

pdf(file="anticorrelation_between_markers_case3.pdf", width = 6, height = 5)
ggplot(data=out, aes(x=P1_CD68, y=`P2_DC-LAMP`)) + 
  geom_point(size=1) + ylim(0,max(out$`P2_DC-LAMP`)) +
  stat_density2d(alpha=0.5) + 
  geom_rug(alpha=0.01,color="black")  +
  geom_abline(intercept = 0, slope = 1, color='firebrick', linetype = 2,size=1) + 
  stat_smooth(method = 'lm', color = 'steelblue') + 
  theme_bw() +
  labs(x="P1_CD68", y="P2_DC-LAMP", title="") +
  annotate(geom="text", x=20, y=5, label="Spearman Rho: -0.51\np.value:1.845e-09",color="steelblue",size=6) +
  theme(text=element_text(size=21))
dev.off()


#"P1_NKp46" & "P2_CD3"
########################################################################################################################################################
ix <- which(test$Marker %in% c("P1_NKp46","P2_CD3"))

out <- as.data.frame(pivot_wider(data=test[ix,c(2,3,6)], names_from=Marker, values_from=value)  %>% unnest())
cor.test(out$P1_NKp46,out$P2_CD3,method = "spear")

pdf(file="anticorrelation_between_markers_case4.pdf", width = 6, height = 5)
ggplot(data=out, aes(x=P1_NKp46, y=P2_CD3)) + 
  geom_point(size=1) + ylim(0,max(out$P2_CD3)) +
  stat_density2d(alpha=0.5) + 
  geom_rug(alpha=0.01,color="black")  +
  geom_abline(intercept = 0, slope = 1, color='firebrick', linetype = 2,size=1) + 
  stat_smooth(method = 'lm', color = 'steelblue') + 
  theme_bw() +
  labs(x="P1_NKp46", y="P2_CD3", title="") +
  annotate(geom="text", x=3, y=20, label="Spearman Rho: -0.41\np.value:4.538e-03",color="steelblue",size=6) +
  theme(text=element_text(size=21))
dev.off()


########################################################################################################################################################
save.image(file="analysis.regeneron.histology.RData")
########################################################################################################################################################



