#rm(list=ls())

require(geiger); require(qpcR); require(phytools); require(ouch); require(ape); require(mvMORPH)
require(OUwie); require(picante); require(parallel); require(surface); require(auteur)

MakeFiles<-function(){
    
checkWD=function(){
  cat(getwd(),"\n")
  if(readline("Is this the correct working directory? [y/n]")=="n")
  {
    cat("Exiting script. Please change working directory and run again. \n")
    break
  }
}

importTree=function(){ #interactive
  print(list.files())
  as.numeric(readline("What index value is your tree? \n"))->treeIndex
  read.tree(list.files()[treeIndex])
}

importData=function(){ #interactive
  print(list.files())
  as.numeric(readline("What index value is your dataset? \n"))->dataIndex
  as.matrix(read.table(list.files()[dataIndex], header=T, row.names=1, na.strings="NA", sep=","))
}

	checkWD()
    TreeFile<-importTree()
    DataFile<-importData()

Pruning<-treedata(TreeFile, DataFile)

CentTreePrune<-Pruning$phy
CentTreePrune<-ladderize(CentTreePrune)
CentDataPrune<-Pruning$data

CentDataPrune<- CentDataPrune[match(CentTreePrune$tip.label, rownames(CentDataPrune)),]
CentDataPrune<-as.data.frame(CentDataPrune)

CentTreePrune<-ladderize(CentTreePrune,T)
plot(CentTreePrune, cex=0.5)
axisPhylo()

Files<-list(CentTreePrune, CentDataPrune)

}

Standardize <- function(Data, Trait, Std) {

 NewTrait <- Data[,which(colnames(Data)==Trait)]/Data[,which(colnames(Data)==Std)]
	
 Data[,which(colnames(Data)==Trait)] <- NewTrait

  Data

}

GeigerModels<-function(Tree, Data, Trait){
	
Skull<-Data[,Trait]
names(Skull)<-rownames(Data)

CentTreePrune<-Tree

brown_skull<-fitContinuous(CentTreePrune, Skull, model="BM")
eb_skull<-fitContinuous(CentTreePrune, Skull, model="EB")
ou_skull<-fitContinuous(CentTreePrune, Skull, model="OU")

model.skull<-matrix(,3,4,dimnames = list(c("Brownian Motion", "Early Burst", "Ornstein-Uhlenbeck"),c("log likelihood", "AICc", "Delta AICc", "AICc Weights")))

model.skull[,1]<-c(brown_skull$opt$lnL, eb_skull$opt$lnL, ou_skull$opt$lnL)
model.skull[,2]<-c(brown_skull$opt$aicc, eb_skull$opt$aicc, ou_skull$opt$aicc)

aic.all.skull<-as.matrix(model.skull[,2])

scor.wts.skull<-akaike.weights(aic.all.skull)

model.skull[,3]<-scor.wts.skull$deltaAIC
model.skull[,4]<-scor.wts.skull$weights
model.skull
}


OuchModels<-function(Tree, Data, Trait, Regime){

BatTreePrune<-Tree
BatDataPrune<-Data

BatTreePruneOU<-ape2ouch(BatTreePrune)

DFBatTreePruneOU<-as(BatTreePruneOU,"data.frame")

BatDataPrune$labels<-row.names(BatDataPrune)

BatDataPruneOU<-merge(DFBatTreePruneOU, BatDataPrune, by="labels", all=T)

row.names(BatDataPruneOU)<-BatDataPruneOU$nodes

RUN_DFBatTreePruneOU <-ouchtree(nodes= BatDataPruneOU$nodes, ancestors= BatDataPruneOU$ancestors, times= BatDataPruneOU$times, labels= BatDataPruneOU$labels)

#Brownian#
brown_posize<-brown(data= BatDataPruneOU[Trait], tree= RUN_DFBatTreePruneOU)
Brown_Pos<-summary(brown_posize)

#Single peak OU#
BatDataPruneOU$ou1<-as.factor("single")
ou1_PC1<-hansen(data= BatDataPruneOU[Trait], tree=RUN_DFBatTreePruneOU, regimes= BatDataPruneOU["ou1"], sqrt.alpha=1, sigma=1, fit=TRUE)
OU1_Sing<-summary(ou1_PC1)

#plot(ou1_PC1)


#Multi-peak OU coded by diet#
#Ancestors derived from ancestral state reconstruction.
select_anc<-ace(Data[,Regime], Tree, type="discrete")


#Needs to be automated... Done#
anc_node<-select_anc$lik.anc

anc_nodes<-NULL
for(i in 1:dim(anc_node)[1]){

anc_nodes[i]<-order(anc_node[i,],decreasing=T)[1]

}

myanctree<-BatTreePrune
myanctree$node.label<-anc_nodes

myancoutree<-ape2ouch(myanctree)
myancdataframe<-as(myancoutree, "data.frame")

Tips<-length(Tree$tip.label)
TheAncestors<-Tree$Nnode

diet_ancou<-c(as.vector(myancdataframe$labels[1:TheAncestors]), BatDataPruneOU[,Regime][TheAncestors+1:Tips])

BatDataPruneOU$dietanc<-as.factor(diet_ancou)

ouanc_diet<-hansen(data= BatDataPruneOU[Trait],tree=RUN_DFBatTreePruneOU, regimes= BatDataPruneOU["dietanc"], sqrt.alpha=1, sigma=1, fit=TRUE)
OU_Anc<-summary(ouanc_diet)

#plot(ouanc_diet, cex=0.6)


#Contrasting OU models#
OU.skull<-matrix(,3,4,dimnames = list(c("Brownian", "Single Peak", "Multi-Peak"),c("log likelihood", "AICc", "Delta AICc", "AICc Weights")))

OU.skull[,1]<-c(Brown_Pos$loglik, OU1_Sing$loglik, OU_Anc$loglik)
OU.skull[,2]<-c(Brown_Pos$aic.c, OU1_Sing$aic.c, OU_Anc$aic.c)

aic.OU.skull<-as.matrix(OU.skull[,2])

OU.wts.skull<-akaike.weights(aic.OU.skull)

OU.skull[,3]<-OU.wts.skull$deltaAIC
OU.skull[,4]<-OU.wts.skull$weights
OU.skull
}

SimmapTrees<-function(Tree, Data, Regime, TreeSim=5){
X<-Data[,Regime]
names(X)<-rownames(Data)

mtrees<-make.simmap(tree=Tree, x=X, nsim=TreeSim)
mtrees
}


PlotSimTreeDiet<-function(SimTrees, FontSize=0.5){
	cols<-setNames(c("blue","red"),c(1, 2))
	
	par(mfrow=c(2,2))
		for(i in 1:4){
			plotSimmap(SimTrees[[i]], cols, fsize=FontSize, pts=F)
		}
	par(mfrow=c(1,1))
}


PlotSimTreeHabitat<-function(SimTrees, FontSize=0.5){
	cols<-setNames(c("blue", "red", "yellow", "darkgreen"),c(1, 2, 3, 4))
	
	par(mfrow=c(2,2))
		for(i in 1:4){
			plotSimmap(SimTrees[[i]], cols, fsize=FontSize, pts=F)
		}
	par(mfrow=c(1,1))
}


OUwieModels<-function(SimTree, Data, Trait, Regime, Model){
Data$Species<-rownames(Data)
BatData<-Data[,c("Species", Regime, Trait)]

write.table(BatData, "TempFile.txt", row.names=F)
TempData<-read.table("TempFile.txt", header=T)

OUMVSampledOutput = lapply(SimTree, function(x) {
	OUwie(x, TempData, model = Model, simmap.tree = TRUE)
 })
}

OUwieModelsMC<-function(SimTree, Data, Trait, Regime, Model){
Data$Species<-rownames(Data)
BatData<-Data[,c("Species", Regime, Trait)]

write.table(BatData, "TempFile.txt", row.names=F)
TempData<-read.table("TempFile.txt", header=T)

OUMVSampledOutput = mclapply(SimTree, function(x) {
	OUwie(x, TempData, model = Model, simmap.tree = TRUE)
 })
}


CollateOUwieNormDiet<-function(OUOutput){

OUwieResults<-matrix(NA,length(OUOutput),8)

for(i in 1:length(OUOutput)){
	
	OUwieResults[i,1] <-OUOutput[[i]]$loglik
	OUwieResults[i,2] <-OUOutput[[i]]$AICc
	OUwieResults[i,3] <-OUOutput[[i]]$solution[1,1]
	OUwieResults[i,4] <-OUOutput[[i]]$solution[1,2]
	OUwieResults[i,5] <-OUOutput[[i]]$solution[2,1]
	OUwieResults[i,6] <-OUOutput[[i]]$solution[2,2]
	OUwieResults[i,7] <-OUOutput[[i]]$theta[1]
	OUwieResults[i,8] <-OUOutput[[i]]$theta[2]


}

colnames(OUwieResults)<-c("log likelihood", "AICc", "Alpha Peak 1", "Alpha Peak 2", "Rate Peak 1", "Rate Peak 2", "Theta Peak 1", "Theta Peak 2")
write.csv(OUwieResults, "OUwieResultsDiet.csv")
OUwieResults


}


CollateOUwieEigDiet<-function(OUOutput){
Less<-NULL
Temps<-list()
for(i in 1:length(OUOutput)){
Temps[[i]]<-OUOutput[[i]]$eigval

Less[i]<-sum(Temps[[i]]<0)
}

Remove<-which(Less!=0)

MyOUwieModelFinal<-OUOutput[-Remove]


OUwieResults<-matrix(NA,length(MyOUwieModelFinal),8)

for(i in 1:length(MyOUwieModelFinal)){
	
	OUwieResults[i,1] <-MyOUwieModelFinal[[i]]$loglik
	OUwieResults[i,2] <-MyOUwieModelFinal[[i]]$AICc
	OUwieResults[i,3] <-MyOUwieModelFinal[[i]]$solution[1,1]
	OUwieResults[i,4] <-MyOUwieModelFinal[[i]]$solution[1,2]
	OUwieResults[i,5] <-MyOUwieModelFinal[[i]]$solution[2,1]
	OUwieResults[i,6] <-MyOUwieModelFinal[[i]]$solution[2,2]
	OUwieResults[i,7] <-MyOUwieModelFinal[[i]]$theta[1]
	OUwieResults[i,8] <-MyOUwieModelFinal[[i]]$theta[2]



	}

colnames(OUwieResults)<-c("log likelihood", "AICc", "Alpha Peak 1", "Alpha Peak 2", "Rate Peak 1", "Rate Peak 2", "Theta Peak 1", "Theta Peak 2")
write.csv(OUwieResults, "OUwieResultsDiet.csv")
OUwieResults


}



CollateOUwieNormHabitat<-function(OUOutput){

OUwieResults<-matrix(NA,length(OUOutput),14)

for(i in 1:length(OUOutput)){
	
	OUwieResults[i,1] <-OUOutput[[i]]$loglik
	OUwieResults[i,2] <-OUOutput[[i]]$AICc
	OUwieResults[i,3] <-OUOutput[[i]]$solution[1,1]
	OUwieResults[i,4] <-OUOutput[[i]]$solution[1,2]
	OUwieResults[i,5] <-OUOutput[[i]]$solution[1,3]
	OUwieResults[i,6] <-OUOutput[[i]]$solution[1,4]
	OUwieResults[i,7] <-OUOutput[[i]]$solution[2,1]
	OUwieResults[i,8] <-OUOutput[[i]]$solution[2,2]
	OUwieResults[i,9] <-OUOutput[[i]]$solution[2,3]
	OUwieResults[i,10] <-OUOutput[[i]]$solution[2,4]
	OUwieResults[i,11] <-OUOutput[[i]]$theta[1]
	OUwieResults[i,12] <-OUOutput[[i]]$theta[2]
	OUwieResults[i,13] <-OUOutput[[i]]$theta[3]
	OUwieResults[i,14] <-OUOutput[[i]]$theta[4]



	}

colnames(OUwieResults)<-c("log likelihood", "AICc", "Alpha Peak 1", "Alpha Peak 2", "Alpha Peak 3", "Alpha Peak 4", "Rate Peak 1", "Rate Peak 2", "Rate Peak 3", "Rate Peak 4", "Theta Peak 1", "Theta Peak 2", "Theta Peak 3", "Theta Peak 4")
write.csv(OUwieResults, "OUwieResultsHabitat.csv")
OUwieResults

}

CollateOUwieEigHabitat<-function(OUOutput){
Less<-NULL
Temps<-list()
for(i in 1:length(OUOutput)){
Temps[[i]]<-OUOutput[[i]]$eigval

Less[i]<-sum(Temps[[i]]<0)
}

Remove<-which(Less!=0)

MyOUwieModelFinal<-OUOutput[-Remove]


OUwieResults<-matrix(NA,length(MyOUwieModelFinal),14)

for(i in 1:length(MyOUwieModelFinal)){
	
	OUwieResults[i,1] <-MyOUwieModelFinal[[i]]$loglik
	OUwieResults[i,2] <-MyOUwieModelFinal[[i]]$AICc
	OUwieResults[i,3] <-MyOUwieModelFinal[[i]]$solution[1,1]
	OUwieResults[i,4] <-MyOUwieModelFinal[[i]]$solution[1,2]
	OUwieResults[i,5] <-MyOUwieModelFinal[[i]]$solution[1,3]
	OUwieResults[i,6] <-MyOUwieModelFinal[[i]]$solution[1,4]
	OUwieResults[i,7] <-MyOUwieModelFinal[[i]]$solution[2,1]
	OUwieResults[i,8] <-MyOUwieModelFinal[[i]]$solution[2,2]
	OUwieResults[i,9] <-MyOUwieModelFinal[[i]]$solution[2,3]
	OUwieResults[i,10] <-MyOUwieModelFinal[[i]]$solution[2,4]
	OUwieResults[i,11] <-MyOUwieModelFinal[[i]]$theta[1]
	OUwieResults[i,12] <-MyOUwieModelFinal[[i]]$theta[2]
	OUwieResults[i,13] <-MyOUwieModelFinal[[i]]$theta[3]
	OUwieResults[i,14] <-MyOUwieModelFinal[[i]]$theta[4]


	}

colnames(OUwieResults)<-c("log likelihood", "AICc", "Alpha Peak 1", "Alpha Peak 2", "Alpha Peak 3", "Alpha Peak 4", "Rate Peak 1", "Rate Peak 2", "Rate Peak 3", "Rate Peak 4", "Theta Peak 1", "Theta Peak 2", "Theta Peak 3", "Theta Peak 4")
write.csv(OUwieResults, "OUwieResultsHabitat.csv")
OUwieResults


}



AncestralMap<-function(Tree, Data, Trait){
	AncRecon <-Data[,Trait]
	names(AncRecon)<-rownames(Data)

	contMap(Tree, AncRecon)
}


SURFACE<-function(Tree, Data, Trait){
if(class(Trait) == "numeric"){
names(Trait)<-rownames(Data)
Trees<-Tree
Traits<-Trait

Traits<-as.data.frame(Traits)
Trees<-nameNodes(Trees)

olist<-convertTreeData(Trees, Traits)
otree<-olist[[1]]; odata<-olist[[2]]

fwd<-surfaceForward(otree, odata, aic_threshold = 0, exclude = 0, verbose = FALSE, plotaic = FALSE)
k<-length(fwd)
fsum<-surfaceSummary(fwd)  

bwd<-surfaceBackward(otree, odata, starting_model = fwd[[k]], aic_threshold = 0, only_best = FALSE, verbose = FALSE, plotaic = FALSE)
bsum<-surfaceSummary(bwd)
kk<-length(bwd)

surfaceTreePlot(Trees, bwd[[kk]], labelshifts = T, cex=0.5)

SurfaceOutput<-list(fsum, bsum)	
SurfaceOutput	
}

else{
	
Trees<-Tree
Traits<-Trait

Traits<-as.data.frame(Traits)
Trees<-nameNodes(Trees)

olist<-convertTreeData(Trees, Traits)
otree<-olist[[1]]; odata<-olist[[2]]

fwd<-surfaceForward(otree, odata, aic_threshold = 0, exclude = 0, verbose = FALSE, plotaic = FALSE)
k<-length(fwd)
fsum<-surfaceSummary(fwd)  

bwd<-surfaceBackward(otree, odata, starting_model = fwd[[k]], aic_threshold = 0, only_best = FALSE, verbose = FALSE, plotaic = FALSE)
bsum<-surfaceSummary(bwd)
kk<-length(bwd)

surfaceTreePlot(Trees, bwd[[kk]], labelshifts = T, cex=0.5)

SurfaceOutput<-list(fsum, bsum)	
SurfaceOutput		
	
}

}

#Disparity#
Disparity<-function(Tree, Data, Trait, Range=c(0,0.8), Sims=10000){

if(class(Trait) == "numeric"){
names(Trait)<-rownames(Data)
Trees<-Tree
Traits<-Trait

Traits<-as.data.frame(Traits)
Trees<-nameNodes(Trees)

DamDisp<-dtt(Trees, Traits, mdi.range=Range, nsim=Sims)

#dtt - Average disparity for each node.
#MDI - MDI Statistic is the area between the DTT for the data and the median of the simulations.

	Disp<-list(DamDisp$dtt, DamDisp$MDI)
	
	}
	
else{
Trees<-Tree
Traits<-Trait

Traits<-as.data.frame(Traits)
Trees<-nameNodes(Trees)

DamDisp<-dtt(Trees, Traits, mdi.range=Range, nsim=Sims)

#dtt - Average disparity for each node.
#MDI - MDI Statistic is the area between the DTT for the data and the median of the simulations.

	Disp<-list(DamDisp$dtt, DamDisp$MDI)

	}
}


PMorphospace<-function(Tree, Data, Regime, Trait1, Trait2){
	
	attach(Data)
	
	X<-cbind(Trait1,Trait2)
	rownames(X)<-rownames(Data)
	
	phylomorphospace(Tree, X)
	
	colour.code <- c("blue", "red", "yellow", "darkgreen", "purple", "grey", "cornflowerblue", "aquamarine")
	
	points(X[ ,1], X[ ,2], col=colour.code[unclass(Regime)], pch=20, cex=1.2)
	
	detach(Data)
	
}


Rate.Shifts <- function(Tree, Data, Trait, Ngens=10000) {

 AuteurData <- Data[,Trait]
 names(AuteurData) <- rownames(Data)

 rjmcmc.bm(Tree, AuteurData, summary=TRUE, reml=TRUE, ngen=Ngens, sample.freq=500, simplestart=T, fileBase=Trait)

 shifts.plot(phy=Tree, base.dir=paste("BM", Trait, "parameters", sep="."), burnin=0.25, legend=TRUE, edge.width=2)

}

CompareRates.multTrait <- function(phy, x, TraitCov=T, ms.err=NULL, ms.cov=NULL){
  #Compares LLik of R-matrix vs. LLik of R-matrix with constrained diagonal
  
  #TraitCov = TRUE assumes covariation among traits (default)
  #ms.err allows the incorporation of within-species measurement error. Input is a matrix of species (rows) by within-species variation for each trait (columns).
  #ms.cov allows the incorporation of within-species covariation between traits. Input is a matrix of species (rows) by within-species covariation for each pair of traits (columns). These must be provided in a specific order, beginning with covariation between trait 1 and the rest, then trait 2 and the rest, etc. For instance, for 4 traits, the columns are: cov_12, cov_13, cov_14, cov_23, cov_24 cov_34.
  
  #Some calculations adapted from 'evol.vcv' in phytools (Revell, 2012)
  
  library(MASS)
  x<-as.matrix(x)
  N<-nrow(x)
  p<-ncol(x)
  C<-vcv.phylo(phy)
  C<-C[rownames(x),rownames(x)]
  if (is.matrix(ms.err)){    
    ms.err<-as.matrix(ms.err[rownames(x),])}
  if (is.matrix(ms.cov)){    
    ms.cov<-as.matrix(ms.cov[rownames(x),])}
  
  #Cholesky decomposition function for diagonal-constrained VCV
  build.chol<-function(b){
    c.mat<-matrix(0,nrow=p,ncol=p)
    c.mat[lower.tri(c.mat)] <- b[-1]  
    c.mat[p,p]<-exp(b[1])
    c.mat[1,1]<-sqrt(sum((c.mat[p,])^2))
    if(p>2){
      for (i in 2:(p-1)){
        c.mat[i,i]<-ifelse( (c.mat[1,1]^2-sum((c.mat[i,])^2) )>0,
                            sqrt(c.mat[1,1]^2-sum((c.mat[i,])^2)), 0)
      }}
    return(c.mat) 
  }
  
  #Fit Rate matrix for all traits: follows code of L. Revell (evol.vcv)
  a.obs<-colSums(solve(C))%*%x/sum(solve(C))   
  D<-matrix(0,N*p,p)
  for(i in 1:(N*p)) for(j in 1:p) if((j-1)*N<i&&i<=j*N) D[i,j]=1.0
  y<-as.matrix(as.vector(x))
  one<-matrix(1,N,1)
  R.obs<-t(x-one%*%a.obs)%*%solve(C)%*%(x-one%*%a.obs)/N
  if (TraitCov==F)    #for TraitCov = F
  { R.obs<-diag(diag(R.obs),p)  }
  #Calculate observed likelihood with or without measurement error
  LLik.obs<-ifelse(is.matrix(ms.err)==TRUE, 
                   -t(y-D%*%t(a.obs))%*%ginv((kronecker(R.obs,C)+ diag(as.vector(ms.err))))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                     determinant((kronecker(R.obs,C)+ diag(as.vector(ms.err))))$modulus[1]/2 , 
                   -t(y-D%*%t(a.obs))%*%ginv(kronecker(R.obs,C))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                     determinant(kronecker(R.obs,C))$modulus[1]/2
  ) 
  
  #Fit common rate for all traits; search over parameter space   
  sigma.mn<-mean(diag(R.obs))   #reasonable start value for diagonal
  
  #Within-species measurement error matrix
  if(is.matrix(ms.err)){m.e<-diag(as.vector(ms.err))}
  
  #Within-species measurement error and trait covariation matrix
  if (is.matrix(ms.err) && is.matrix(ms.cov)){	
    within.spp<-cbind(ms.err,ms.cov)
    rc.label<-NULL
    for (i in 1:p){ rc.label<-rbind(rc.label,c(i,i)) }
    for (i in 1:p){
      for (j in 2:p){ if (i!=j && i<j){rc.label<-rbind(rc.label,c(i,j))} }}
    m.e<-NULL
    for (i in 1:p){
      tmp<-NULL
      for (j in 1:p){
        for (k in 1:nrow(rc.label)){
          if(setequal(c(i,j),rc.label[k,])==T) {tmp<-cbind(tmp,diag(within.spp[,k]))}
        }
      }
      m.e<-rbind(m.e,tmp)
    }
  }
  
  #likelihood optimizer for no trait covariation
  lik.covF<-function(sigma){  
    R<-R.obs
    diag(R)<-sigma
    LLik<-ifelse(is.matrix(ms.err)==TRUE, 
                 -t(y-D%*%t(a.obs))%*%ginv((kronecker(R,C)+ m.e))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant((kronecker(R,C)+ m.e))$modulus[1]/2 , 
                 -t(y-D%*%t(a.obs))%*%ginv(kronecker(R,C))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant(kronecker(R,C))$modulus[1]/2
    ) 
    if (LLik == -Inf) { LLikk <- -1e+10  }
    return(-LLik)
  }
  
  #likelihood optimizer with trait covariation
  lik.covT<-function(sigma){  
    low.chol<-build.chol(sigma)
    R<-low.chol%*%t(low.chol)
    
    LLik<-ifelse(is.matrix(ms.err)==TRUE, 
                 -t(y-D%*%t(a.obs))%*%ginv((kronecker(R,C)+ m.e))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant((kronecker(R,C)+ m.e))$modulus[1]/2 , 
                 -t(y-D%*%t(a.obs))%*%ginv(kronecker(R,C))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant(kronecker(R,C))$modulus[1]/2
    ) 
    if (LLik == -Inf)  {LLikk <- -1e+10  }
    return(-LLik)
  }
  
  ##Optimize for no trait covariation
  if (TraitCov==F)    
  { model1<-optim(sigma.mn,fn=lik.covF,method="L-BFGS-B",lower=c(0.0))}
  ##Optimize with trait covariation
  R.offd<-rep(0,(p*(p-1)/2))
  if (TraitCov==T)  
  {model1<-optim(par=c(sigma.mn,R.offd),fn=lik.covT,method="L-BFGS-B")}
  
  #### Assemble R.constrained
  if (TraitCov==F){R.constr<-diag(model1$par,p)}
  if (TraitCov==T){  
    chol.mat<-build.chol(model1$par)
    R.constr<-chol.mat%*%t(chol.mat)}
  
  if(model1$convergence==0)
    message<-"Optimization has converged."
  else
    message<-"Optim may not have converged.  Consider changing start value or lower/upper limits."
  LRT<- (-2*((-model1$value-LLik.obs)))
  LRT.prob<-pchisq(LRT, (p-1),lower.tail=FALSE) #df = Nvar-1
  AIC.obs<- -2*LLik.obs+2*p+2*p #(2p twice: 1x for rates, 1x for anc. states)
  AIC.common<- -2*(-model1$value)+2+2*p #(2*1: for 1 rate 2p for anc. states)
  return(list(Robs=R.obs, Rconstrained=R.constr,Lobs=LLik.obs,Lconstrained=(-model1$value),LRTest=LRT,Prob=LRT.prob,
              AICc.obs=AIC.obs,AICc.constrained=AIC.common,optimmessage=message))   
 }

Rate.Correlation <- function(Tree, Data, Trait1, Trait2) {
 DataM<-as.matrix(Data)

  Rates<-CompareRates.multTrait(Tree, DataM[,c(which(colnames(DataM)==Trait1), which(colnames(DataM)==Trait2))])
  Rates

}


Trait.Correlation <- function(Tree, Data, Trait1, Trait2) {
T1 <- Data[,Trait1]
names(T1) <- rownames(Data)
T1.pic <- pic(T1, Tree)

T2 <- Data[,Trait2]
names(T2) <- rownames(Data)
T2.pic<-pic(T2, Tree)

lm.pic <- lm(T2.pic ~ T1.pic -1) 

plot(T1.pic, T2.pic, pch=19)
abline(lm.pic, col="red")

#summary(lm.pic)
Pvalue.pic <- anova(lm.pic)$'Pr(>F)'[1]
Pvalue.pic
}


MultiVarModels <- function(Tree, Data, Trait1, Trait2, NSim=10, Model="OUM", Test=c(BM, EB, OU)){
attach(Data)

names(Habitat) <- Tree$tip.label
SimTree <- make.simmap(Tree, Habitat, model="ER", nsim=NSim)

names(DietFive) <- Tree$tip.label
SimTree <- make.simmap(Tree, DietFive, model="ER", nsim=2)

names(DietEight) <- Tree$tip.label
SimTree <- make.simmap(Tree, DietEight, model="ER", nsim=2)


DataOU<-cbind(Data$PC1, Data$PC2)

rownames(DataOU)<-Tree$tip.label


OUSampledOutput = lapply(SimTree, function(x) {
mvBM(x, DataOU, model=Model)
})

OUSampledOutput = lapply(SimTree, function(x) {
mvBM(x, DataOU, model="BMM")
})


OUSampledOutput = lapply(SimTree, function(x) {
mvEB(x, DataOU)
})


OUSampledOutput = lapply(SimTree, function(x) {
mvOU(x, DataOU, model="OU1")
})

OUSampledOutput = lapply(SimTree, function(x) {
mvOU(x, DataOU, model="OUM", scale.height=TRUE)
})
}