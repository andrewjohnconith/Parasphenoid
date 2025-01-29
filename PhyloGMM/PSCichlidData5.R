####F5 Cichlid Data####
install.packages("geomorph")
#Load packages
library(geomorph); library(plyr); library(abind); library(geiger); library(mvMORPH)
setwd("/Users/aconith/Library/CloudStorage/OneDrive-DePaulUniversity/DePaul/Research/Research/ParaSWriting/Github/PhyloGMM/")

#Read in data
AllInds<-readland.tps("PS_20240625.tps", specID = "imageID", readcurves = T)
AllData<-read.csv("Reduced_TropheopsCatalogueEd4.csv")
curveslide<-read.csv("PSCurveslide4.csv")

#Perform Procrustes
Y.gpa<-gpagen(AllInds, curves = curveslide, ProcD = T)

gdf <- geomorph.data.frame(Y.gpa, species = dimnames(Y.gpa$coords)[[3]])

MyDepthVec<-as.factor(AllData$Depth)
MyColVec<-c("black","gray")

#Perform PCA
par(pty='s')
CichlidPCA<-gm.prcomp(gdf$coords)

plot(CichlidPCA, pch=19, axis1 = 1, axis2 = 2, col=MyColVec[unclass(MyDepthVec)])
text(CichlidPCA$x[,c(1,2)], pos = 2, label = rownames(CichlidPCA$x), cex=0.5)
picknplot.shape(plot(CichlidPCA))

write.csv(cbind(Y.gpa$Csize,CichlidPCA$x[,1:5]),"CichlidPCA_AllInds.csv")

#Check for shape differences between habitats
Demap<-as.factor(AllData$Depth)
names(Demap)<-rownames(AllData$ID)
Fit<-procD.lm(f1 = coords~Demap, data = gdf)
summary(Fit)


####Species Means####
NaturalData<-read.csv("Reduced_TropheopsCatalogueEd4.csv")
#Coordinates#
NaturalData$PhyloName<-as.factor(NaturalData$PhyloName)
#Create mean GPA from all specimens
CichlidMeanShapes<-vector("list", length(levels(NaturalData$PhyloName)))

for(i in 1:length(levels(NaturalData$PhyloName))){
  
  CichlidNames <- levels(NaturalData$PhyloName)[i]
  NamesRows <- which(NaturalData$PhyloName == CichlidNames)
  
  CatTPS <- Y.gpa$coords[,,paste0(as.character(NaturalData[NamesRows,"SpecimenID"]))]
  
  #If you have a single entry need to convert to array (R thinks it is a matrix)
  if(class(CatTPS)[1] != "array"){
    CatTPS<-array(CatTPS, dim=c(88,3,1))
  }
  
  Cichlid.gpa<-gpagen(CatTPS, curves = curveslide, ProcD = T)
  CichlidMeanShapes[[i]] <- mshape(Cichlid.gpa$coords)
  
}

#Convert from list to array
CichlidMeanShapes <-array(unlist(CichlidMeanShapes), dim=c(88,3,length(levels(NaturalData$PhyloName))))

#Add species names
CichlidDimName<-vector(mode="list", 3)
CichlidDimName[[3]]<-levels(NaturalData$PhyloName)
dimnames(CichlidMeanShapes)<-CichlidDimName

CichlidMeanShapesY<-gpagen(CichlidMeanShapes, curves = curveslide, ProcD = T)

morphmean <- aggregate(cbind(as.factor(Depth),Csize)~PhyloName,data=NaturalData,FUN="mean",na.rm=TRUE,na.action=NULL)
CichlidMeanShapesY$Csize<-morphmean$Csize

rownames(morphmean)<-morphmean$PhyloName

#write.csv(morphmean, "morphmean2.csv")

####PCA####
PS_CichlidMeanShapes<-CichlidMeanShapesY$coords
TaxaIDData<-read.csv("morphmean2.csv", row.names = 1)

CichlidPCA2<-gm.prcomp(PS_CichlidMeanShapes)
MyDepthVec<-as.factor(TaxaIDData$Depth)
MyColVec<-c("black","gray")

par(pty='s')
plot(CichlidPCA2, col=MyColVec[unclass(MyDepthVec)], pch=19, axis1 = 2, axis2 = 3)
text(CichlidPCA2$x[,c(1,2)], pos = 4, label = rownames(TaxaIDData), cex = 0.5)
picknplot.shape(plot(CichlidPCA2, col=MyColVec[unclass(MyDepthVec)], pch=19, axis1 = 2, axis2 = 3), method = "points", mag = 1)

write.csv(CichlidPCA2$x[,1:5],"CichlidPCA_MeanInds.csv")


####PCM####
Tree<-read.tree("Tree-CompMeth-Ultra_20181101_bb.tre")

Pruning<-treedata(Tree, TaxaIDData)
Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

FData$PC1<-as.numeric(as.character(FData$PC1))
FData$PC2<-as.numeric(as.character(FData$PC2))

RM_Cichlid<-dimnames(PS_CichlidMeanShapes)[[3]][dimnames(PS_CichlidMeanShapes)[[3]]%in%rownames(FData)]
PS_CichlidMeanShapes_Tr<-PS_CichlidMeanShapes[,,RM_Cichlid]


PS_CichlidMeanShapes_Tr<-PS_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(PS_CichlidMeanShapes_Tr)[[3]]))]

gdf <- geomorph.data.frame(coords = PS_CichlidMeanShapes_Tr, species = dimnames(PS_CichlidMeanShapes_Tr)[[3]], Csize = as.numeric(FData$Csize))

##Checking for allometry##
  #No signal for allometry
HybAnova <- procD.pgls(coords~Csize, phy = FTree, data = gdf) 
summary(HybAnova) 

##
#Perform PCA
MyDepthVec<-as.factor(FData$Depth)
MyColVec<-c("black","gray")

par(pty='s')
plot(x = FData$PC1, y = FData$PC2, col=MyColVec[unclass(MyDepthVec)], pch=19)
text(x = FData$PC1, y = FData$PC2, pos = 4, label = rownames(FData), cex = 0.5)

####PCMs####
Demap<-as.factor(FData$Depth)
names(Demap)<-rownames(FData)

DePC1<-FData$PC1
names(DePC1)<-rownames(FData)

Fit1<-procD.pgls(f1 = DePC1~Demap, phy = FTree)
summary(Fit1)

##
DePC2<-FData$PC2
names(DePC2)<-rownames(FData)

Fit2<-procD.pgls(f1 = DePC2~Demap, phy = FTree)
summary(Fit2)
#
DePC3<-as.numeric(FData$PC3)
names(DePC3)<-rownames(FData)

Fit3<-procD.pgls(f1 = DePC3~Demap, phy = FTree)
summary(Fit3)


###Full shape
Fit1<-procD.pgls(f1 = PS_CichlidMeanShapes_Tr~Demap, phy = FTree)
summary(Fit1)



####Transformed Trees####
####PCMs####
Demap<-as.factor(FData$Depth)
names(Demap)<-rownames(FData)

DePC1<-FData$PC1
names(DePC1)<-rownames(FData)

Fit1<-procD.pgls(f1 = DePC1~Demap, phy = FTree0_5)
summary(Fit1)

##
DePC2<-FData$PC2
names(DePC2)<-rownames(FData)

Fit2<-procD.pgls(f1 = DePC2~Demap, phy = FTree0_5)
summary(Fit2)
#
DePC3<-as.numeric(FData$PC3)
names(DePC3)<-rownames(FData)

Fit3<-procD.pgls(f1 = DePC3~Demap, phy = FTree0_5)
summary(Fit3)


###Full shape
Fit1<-procD.pgls(f1 = PS_CichlidMeanShapes_Tr~Demap, phy = FTree0_5)
summary(Fit1)


####Transformed Trees####
####PCMs####
Demap<-as.factor(FData$Depth)
names(Demap)<-rownames(FData)

DePC1<-FData$PC1
names(DePC1)<-rownames(FData)

Fit1<-procD.pgls(f1 = DePC1~Demap, phy = FTree1_5)
summary(Fit1)

##
DePC2<-FData$PC2
names(DePC2)<-rownames(FData)

Fit2<-procD.pgls(f1 = DePC2~Demap, phy = FTree1_5)
summary(Fit2)
#
DePC3<-as.numeric(FData$PC3)
names(DePC3)<-rownames(FData)

Fit3<-procD.pgls(f1 = DePC3~Demap, phy = FTree1_5)
summary(Fit3)


###Full shape
Fit1<-procD.pgls(f1 = PS_CichlidMeanShapes_Tr~Demap, phy = FTree1_5)
summary(Fit1)

##Evo Models
source("CichlidDepth_SOURCE.R")

Files<-MakeFiles()


#Your tree and data is stored as 'Tree' and 'Data' respectively.
Tree<-Files[[1]]
Data<-Files[[2]]
str(Data)

#For cichlids as listed as factors
Data[,1]<-as.numeric(as.character(Data[,1]))
Data[,4]<-as.numeric(as.character(Data[,4]))
Data[,5]<-as.numeric(as.character(Data[,5]))
Data[,6]<-as.numeric(as.character(Data[,6]))
Data[,7]<-as.numeric(as.character(Data[,7]))
Data[,8]<-as.numeric(as.character(Data[,8]))
str(Data)

#Geiger Evolutionary Models#
#Single variable
MyModels<-GeigerModels(Tree, Data, "PC1")
#View results
MyModels

#OUCH Evolutionary Models#
#Single Variable
MyOUCHModels<-OuchModels(Tree, Data, "PC1", "Depth")
#View results
MyOUCHModels


#OUwie Evolutionary Models#
#Generate trees with stocastically mapped characters.
MySimTrees<-SimmapTrees(Tree, Data, "Depth", 1000)

#Perform the OU analysis with your recently generated trees.
#Run the analysis over multiple cores to speed it up.
MyOUwieModels<-OUwieModelsMC(MySimTrees, Data, "PC1", "Depth", "OUM")
#View results
MyOUwieModels

##Depth##
FinalOutput<-CollateOUwieNormDiet(MyOUwieModels)
FinalOutput
apply(FinalOutput, MARGIN = 2, FUN = median)
apply(FinalOutput, MARGIN = 2, FUN = range)
apply(FinalOutput, MARGIN = 2, FUN = quantile, probs=c(0.025,0.975))


#Descriptions of the models you can run in place of "OUM".#
#model=OU1  a single peak Ornstein-Uhlenbeck model across the entire tree.
#model=OUM  a multi-peak Ornstein-Uhlenbeck model with different optima (theta) for each regime.
#model=OUMV  a multi-peak Ornstein-Uhlenbeck model with different optima (theta) and different Brownian rate parameter (sigma2) for each regime.
#model=OUMA  a multi-peak Ornstein-Uhlenbeck model with different optima (theta) and different strength of selection parameter (alpha) for each regime.
#OUMA didn't seem to run to completion.
#model=OUMVA  a multi-peak Ornstein-Uhlenbeck model with different optima (theta), different Brownian rate parameter (sigma2), and different strength of selection parameter (alpha) for each regime.
#OUMVA rarely has enough information to work.

akaike.weights(c(-100.31395, -97.68538, -100.16865, -115.8604))


AncestralMap(Tree, Data, "PC1")

####mvMORPH####
#Simmap
library(RRphylo)
FTree1_5<-rescaleRR(FTree,RR=NULL,height=NULL,trend=NULL,delta=NULL,kappa=NULL,lambda=1.5)
FTree0_5<-rescaleRR(FTree,RR=NULL,height=NULL,trend=NULL,delta=NULL,kappa=NULL,lambda=0.5)

DepthVec<-as.factor(Data$Depth)
names(DepthVec)<-rownames(Data)
SimMapTreesSYM_Two<-make.simmap(tree = Tree, x = DepthVec, model = "SYM", nsim=1000)
SimMapTreesSYM_Two1_5<-make.simmap(tree = FTree1_5, x = DepthVec, model = "SYM", nsim=1000)
SimMapTreesSYM_Two0_5<-make.simmap(tree = FTree0_5, x = DepthVec, model = "SYM", nsim=1000)

##Save
save(SimMapTreesSYM_Two, file = "SimMapTreesSYM_Two.rda")

###Trait data setup###
CData <- Data[,4:6]

##single models##
#BM
BM1SampledOutput <- mvBM(Tree, CData, model="BM1", scale.height=TRUE)

BM1SampledOutput1_5 <- mvBM(FTree1_5, CData, model="BM1", scale.height=TRUE)
BM1SampledOutput0_5 <- mvBM(FTree0_5, CData, model="BM1", scale.height=TRUE)

#Single Peak OU#
OU1SampledOutput <- mvOU(Tree, CData, model="OU1", scale.height=TRUE)

OU1SampledOutput1_5 <- mvOU(FTree1_5, CData, model="OU1", scale.height=TRUE)
OU1SampledOutput0_5 <- mvOU(FTree0_5, CData, model="OU1", scale.height=TRUE)

#Single Peak OU# ?
#EBSampledOutput <- mvEB(Tree, CData, scale.height=TRUE, method="pic")

##Habitat multi-peak models##
#Run multivariate OU models
OUMSampledOutputTwo <- mclapply(SimMapTreesSYM_Two, function(x) {
  mvOU(x, CData, model="OUM", scale.height=TRUE)
})
#save(OUMSampledOutputTwo, file="OUMSampledOutputTwo.rda")


OUMSampledOutputTwo1_5 <- mclapply(SimMapTreesSYM_Two1_5, function(x) {
  mvOU(x, CData, model="OUM", scale.height=TRUE)
})
#save(OUMSampledOutputTwo1_5, file="OUMSampledOutputTwo1_5.rda")

OUMSampledOutputTwo0_5 <- mclapply(SimMapTreesSYM_Two0_5, function(x) {
  mvOU(x, CData, model="OUM", scale.height=TRUE)
})
#save(OUMSampledOutputTwo0_5, file="OUMSampledOutputTwo0_5.rda")


SimMapTreesSYM_Two1_5


save(OUMSampledOutputTwo, file="OUMSampledOutputTwo.rda")
rm(OUMSampledOutputTwo)


#Extract model output data
BM1SampledOutput[c(1,3)]

BM1SampledOutput1_5[c(1,3)]
BM1SampledOutput0_5[c(1,3)]

OU1SampledOutput[c(1,3)]

OU1SampledOutput1_5[c(1,3)]
OU1SampledOutput0_5[c(1,3)]

OUMSampledOutputTwoRBIND <- do.call("rbind.data.frame", lapply(OUMSampledOutputTwo, '[', c(1,3)))
colMeans(OUMSampledOutputTwoRBIND)

OUMSampledOutputTwo1_5RBIND <- do.call("rbind.data.frame", lapply(OUMSampledOutputTwo1_5, '[', c(1,3)))
colMeans(OUMSampledOutputTwo1_5RBIND)

OUMSampledOutputTwo0_5RBIND <- do.call("rbind.data.frame", lapply(OUMSampledOutputTwo0_5, '[', c(1,3)))
colMeans(OUMSampledOutputTwo0_5RBIND)
OUMSampledOutputTwo0_5RBIND

#CI
l.model.loglik <- lm(LogLik ~ 1, OUMSampledOutputTwoRBIND)
l.model.aicc <- lm(AICc ~ 1, OUMSampledOutputTwoRBIND)

# Calculate the confidence interval
confint(l.model.loglik, level=0.95)
confint(l.model.aicc, level=0.95)

#CI 1.5
l.model.loglik <- lm(LogLik ~ 1, OUMSampledOutputTwo1_5RBIND)
l.model.aicc <- lm(AICc ~ 1, OUMSampledOutputTwo1_5RBIND)

# Calculate the confidence interval
confint(l.model.loglik, level=0.95)
confint(l.model.aicc, level=0.95)

#CI 0.5
l.model.loglik <- lm(LogLik ~ 1, OUMSampledOutputTwo0_5RBIND)
l.model.aicc <- lm(AICc ~ 1, OUMSampledOutputTwo0_5RBIND)

# Calculate the confidence interval
confint(l.model.loglik, level=0.95)
confint(l.model.aicc, level=0.95)


scor.wts.ps <-akaike.weights(c(-372.6258, -366.0162, -375.1166))
scor.wts.ps1.5 <-akaike.weights(c(-353.2936,-370.5366,-374.9716))
scor.wts.ps0.5 <-akaike.weights(c(-378.3578,-366.1009,-375.2124))

