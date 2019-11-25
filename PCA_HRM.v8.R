##################################################
### PCA and Clustering
##################################################

### Please cite:
# A High Resolution Melting Analysis (HRM)-based Genotyping Toolkit for Chilling Requirement in Peach (Prunus persica)
# (Will be publish in 2019)

### update
## v.8	2019.11.25


### 1 make absolute fluorescence into percentage
#Set up: (1)Data import (2)Create a data frame to store normalized data
input_file = readline('Enter the file name: ') ###Don't forget .csv!
exported.data.file = read.csv(input_file,header = T)### import the raw melt curve

### 2 select melting region
max_lim = readline('Enter the upper limit of the melt region: ')
min_lim = readline('Enter the lower limit of the melt region: ')
exported.data.file <- subset(exported.data.file, Temperature>min_lim & Temperature<max_lim)  # Excluding the pre- and post-melt region
normfluor =exported.data.file #to store the normalized data

### 3 make absolute fluorescence into percentage
for (i in 2:length(exported.data.file[1,]))
{sample.fluorescence.norm <-(exported.data.file[,i])/max(exported.data.file[,i])
normfluor[,i]=sample.fluorescence.norm}
normfluor

### 4 make the median curve
medcurve <- apply(normfluor, 1, median)

### 5 use median curve to normalize percentage melt curves
for(i in 2:length(normfluor[1,])){
  normfluor[,i]=normfluor[,i]-medcurve}
normfluor

### 6 Principal component analysis (PCA)
PCA.analysis.file <-prcomp(normfluor[,2:28], scale = FALSE) 
plot(PCA.analysis.file,type="line",  main="PCs of HRM") # Pareto plot

### 7 Clustering
library(mclust)
cluster.data <- Mclust(PCA.analysis.file$rotation[,1:3], G = 3)##CHANGE HOW MANY PCs TO BE SELECTED HERE!!!!!!!!!!!
df.PCA.MM1 <- as.data.frame(PCA.analysis.file$rotation[,1:3])##CHANGE HOW MANY PCs TO BE SELECTED HERE!!!!!!!!!!!

## report the variant call
cluster.data$classification

## preview plots
#PC1 vs PC2
plot(df.PCA.MM1[,1:2], bg=cluster.data$classification,pch=21,xlab='PC1', ylab='PC2')
identify(df.PCA.MM1[,1:2], labels = rownames(df.PCA.MM1))
#PC2 vs PC3, ignore this if only 2 PCs are selected
plot(df.PCA.MM1[,2:3], bg=cluster.data$classification,pch=21,xlab='PC2', ylab='PC3')
identify(df.PCA.MM1[,2:3], labels = rownames(df.PCA.MM1))
#PC1 vs PC3, ignore this if only 2 PCs are selected
plot(df.PCA.MM1[,c(1,3)], bg=cluster.data$classification, pch=21, xlab='PC1', ylab='PC3')
identify(df.PCA.MM1[,c(1,3)], labels = rownames(df.PCA.MM1))
#3D plot
library(plot3D)
points3D(x=df.PCA.MM1$PC1,y=df.PCA.MM1$PC2,z=df.PCA.MM1$PC3, 
         colvar=cluster.data$classification, 
         xlab = 'PC1', ylab = 'PC2', zlab = 'PC3', pch = 16,
         colkey = F, col=c('black', 'red', 'green')
         )

identify(df.PCA.MM1[,1:3], labels = rownames(df.PCA.MM1))

### 8 Data export
#PC1 vs PC2
pdf(paste(input_file, "PC1 vs PC2.pdf", sep=" ")) 
plot(df.PCA.MM1[,1:2], bg=cluster.data$classification,pch=21,xlab='PC1', ylab='PC2')
dev.off()
#PC2 vs PC3, ignore this if only 2 PCs are selected
pdf(paste(input_file, "PC2 vs PC3.pdf", sep=" ")) # Enter the file name
plot(df.PCA.MM1[,2:3], bg=cluster.data$classification,pch=21,xlab='PC2', ylab='PC3')
dev.off()
#PC1 vs PC3, ignore this if only 2 PCs are selected
pdf(paste(input_file, "PC1 vs PC3.pdf", sep=" ")) # Enter the file name
plot(df.PCA.MM1[,c(1,3)], bg=cluster.data$classification, pch=21, xlab='PC1', ylab='PC3')
dev.off()
#3D plot
pdf(paste(input_file, "3D.pdf", sep=" ")) # Enter the file name
library(plot3D)
points3D(x=df.PCA.MM1$PC1,y=df.PCA.MM1$PC2,z=df.PCA.MM1$PC3, 
         colvar=cluster.data$classification, 
         xlab = 'PC1', ylab = 'PC2', zlab = 'PC3',pch=16,
         colkey = F, col=c('black', 'red', 'green'))
dev.off()

#Output the variant calls
write.csv(cluster.data$classification, file = paste(input_file, "variant call.csv", sep=" "))

