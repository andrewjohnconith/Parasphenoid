#install.packages("/Users/home/Dropbox/Amherst/Courses/SOURCE/Undergraduate/PostDoc/Week3/geomorph_3.0.4.tar.gz", repos = NULL, type="source")
library(geomorph); library(abind)

#Set WD
setwd("/Users/aconith/Library/CloudStorage/OneDrive-DePaulUniversity/DePaul/Research/Research/ParaSWriting/Github/HybridGMM")

#Read in data and TPS file
PSInds<-readland.tps("HybridPSsemi1_20230717.tps", specID = "imageID", readcurves = T)

#Remove additional LMs
RMvec<-c(7,9,31:50,111:130)
PSInds1<-PSInds[-RMvec,,]

#Procrustes
Curves<-read.csv("PSCurveslide4.csv") #3 semi lm sets
Y.gpa<-gpagen(A = PSInds1, curves = Curves, ProcD = T)
#Y.gpa<-gpagen(A = PSInds[c(-7,-9),,], curves = Curves) #Removes lateral vomer LMs

psPCA<-gm.prcomp(Y.gpa$coords)
plot(psPCA, pch=19, cex=0.8)
text(psPCA$x, pos = 4, label = dimnames(Y.gpa$coords)[[3]], cex=0.4)
picknplot.shape(plot(psPCA),method = "points", mag = 2)

write.csv(psPCA$x[,1:5], "PsPCScores_20240701.csv")

####Allometry####
gdf <- geomorph.data.frame(Y.gpa, species = dimnames(Y.gpa$coords)[[3]])

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova <- procD.lm(coords~Csize, data = gdf, iter = 999, RRPP=TRUE) 
summary(HybAnova) 

shape.resid <- arrayspecs(HybAnova$residuals, p=dim(Y.gpa$coords)[1], k=dim(Y.gpa$coords)[2]) # size-adjusted residuals

#allometry-free shapes
PSShape <- shape.resid + array(Y.gpa$consensus, dim(shape.resid))

#Allo Free PCA
PSPCA_allofree<-gm.prcomp(A = PSShape)

par(pty='s')
plot(PSPCA_allofree, axis1 = 1, axis2 = 2, pch=19)
text(PSPCA_allofree$x, pos = 4, label = dimnames(PSShape)[[3]])
picknplot.shape(plot(PSPCA_allofree))

open3d()
#PC1
plotRefToTarget(M1 = PSShape[,,"LFxTRC5137ps"], M2 = PSShape[,,"LFxTRC5085ps"],
                method = "vector", mag = 2)

#PC2
plotRefToTarget(M1 = PSShape[,,"LFxTRC5078ps"], M2 = PSShape[,,"LFxTRC5110ps"],
                method = "vector", mag = 2)


#writeland.tps(A = PSShape, file = "PSnoallo_20230414.tps")

#Save PC Scores
#setwd("/Users/aconith/Dropbox/Amherst/Post-Doc/Muscles/Parasphenoid/GMM")
#write.csv(cbind.data.frame(Y.gpa$Csize, psPCA$x[,1:3], PSPCA_allofree$x[,1:3]), "PsPCScores_20240606.csv")


####PS wireframe####
TPStoMLogikaPS <- function(A){
  #A - An array (p x k x n) containing landmark coordinates for a set of specimens 
  Indviduals<-length(dimnames(A)[[3]])
  Landmarks<-dim(A)[1]
  Dimensions<-dim(A)[2]
  Name<-dimnames(A)[[3]]
  Names<-NULL; for(i in 1:length(Name)){Names[i]<-paste0("'", Name[i])}
  
  MyList<-list(paste("[individuals]"), Indviduals, paste("\n[landmarks]"), Landmarks, paste("\n[dimensions]"), Dimensions, paste("\n[names]"), Names, paste("\n[rawpoints]"))
  lapply(MyList, write, "MorphologikaFile.txt", append=TRUE, ncolumns=1)
  
  for(i in 1:length(Names)){
    #loop through Y.gpa$coords
    #if dimensions=3, then ncolumns=3
    #if dimensions=2, then ncolumns=2
    NewList<-list(paste0(Names[i]), t(A[,,i]), paste0("\n"))
    lapply(NewList, write, "MorphologikaFile.txt", append=TRUE, ncolumns=3)
  }
}

TPStoMLogikaPS(A = Y.gpa$coords)
