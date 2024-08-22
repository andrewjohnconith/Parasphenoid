rm(list = ls(all = TRUE)) 
#load qtl library
library(qtl)

setwd("/Users/aconith/Dropbox/Amherst/Post-Doc/Map/ReAnalyze/MapPhenos/Analysis/Bones/Parasphenoid/MQM/")
MyMap<-read.csv("F5Map_noNAIndsLGd.csv")
MyData<-read.csv("PsPCScores_20240626.csv")

Combine<-merge(MyData,MyMap,by="HybridID",all=T)

write.csv(Combine,"F5Map_wData_noNAIndsLGd_Ps.csv")

#load data
setwd("/Users/aconith/Dropbox/Amherst/Post-Doc/Map/ReAnalyze/MapPhenos/Analysis/Bones/Parasphenoid/MQM/")

F5All <- read.cross(format="csv",file="F5Map_wData_noNAIndsLGd_Ps.csv",
                    na.strings="NA",genotypes=c("AA","AB","BB"),
                    alleles=c("A","B"),convertXdata=TRUE)


gt <- geno.table(F5All)
todrop3 <- rownames(gt[gt$P.value < 1e-300,])
F5Omit <- drop.markers(F5All, todrop3)


augdata1<-fill.geno(F5Omit)
#write.cross(cross = augdata1, format = "csv", filestem = "AugCrossF5_PS-202406b")

##Read in augmented data
augdata1 <- read.cross(format="csv",file="AugCrossF5_PS-202406.csv",
                    na.strings="NA",genotypes=c("AA","AB","BB"),
                    alleles=c("A","B"),convertXdata=TRUE)


####PS####
#If you get the singular matrix, re-run the fill.geno step
scan1<-mqmscan(augdata1,pheno.col=2,plot=T,model="dominance")  #Step 1: scan without cofactors
summary(scan1)

#cofactorsindex<-NULL
find.markerindex(augdata1,find.marker(augdata1,chr=7,pos=20))
#find.markerindex(augdata1,find.marker(augdata1,chr="10b",pos=15))
#find.markerindex(augdata1,find.marker(augdata1,chr="13b",pos=0))
find.markerindex(augdata1,find.marker(augdata1,chr="13a",pos=15))
augdata1<-fill.geno(F5Omit)

#cofactorsindex<-c(794, 520, 419, 162, 254, 40, 636, 92, 166, 63, 503, 590, 286, 582, 443, 640, 658, 57, 55, 752, 737, 25, 815, 461, 835, 719, 337, 734, 66, 216, 71, 235, 317, 365, 548, 765, 642, 627, 288, 672, 12, 85, 100, 156, 448, 455, 468, 516, 694, 824) #v1
#cofactorsindex<-c(242, 395, 675, 11, 802, 355, 618, 685, 53, 176, 329, 519, 131, 591, 200, 344, 56, 501, 407, 19, 253, 358, 412, 809, 704, 598, 715, 447, 589, 512, 436, 165, 728, 748, 463, 400, 606, 784, 451, 92, 609, 257, 36, 364, 456, 472, 671) #v2
cofactorsindex<-c(354, 445, 426, 482, 165, 75, 812, 261, 276, 604, 408, 123, 175, 734, 311, 227, 112, 181, 38, 339, 699, 315, 307, 521, 650, 716, 639, 99, 770, 448, 406, 167, 729, 418, 582, 59, 337, 630, 498, 319, 789, 374, 35, 71, 85, 170, 293, 333, 457, 472)#v3
homemadecofactors<-mqmsetcofactors(augdata1,cofactors=c(cofactorsindex)) #Step 3: designate cofactors for multivariate scan.

homescan1<-mqmscan(cross=augdata1,cofactors=homemadecofactors,pheno.col=2,cofactor.significance=0.02,verbose=T,plot=T,model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan1)

find.marker(augdata1,chr=7,pos=20)
effectplot(F5All, pheno.col=2, mname1="scaffold_21_2195347",var.flag=c("group"))

#find.marker(augdata1,chr="10b",pos=15)
#effectplot(F5All, pheno.col=2, mname1="scaffold_75_1746024",var.flag=c("group"))

find.marker(augdata1,chr=12,pos=30)
effectplot(F5All, pheno.col=2, mname1="scaffold_9_9277793",var.flag=c("group"))

find.marker(augdata1,chr="13a",pos=15) #scaffold_44_2763411
find.marker(augdata1,chr="13a",pos=20) #scaffold_44_2763411
find.marker(augdata1,chr="13a",pos=35) #scaffold_62_1648691
find.marker(augdata1,chr="13a",pos=40)
effectplot(augdata1, pheno.col=2, mname1="scaffold_44_2763411",var.flag=c("group"))

find.markerindex(augdata1,find.marker(augdata1,chr=7,pos=20))
bayesint(homescan1,7,qtl.index=256,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.

find.markerindex(augdata1,find.marker(augdata1,chr=12,pos=30))
bayesint(homescan1,12,qtl.index=446,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.

find.markerindex(augdata1,find.marker(augdata1,chr="13a",pos=20))
bayesint(homescan1,"13a",qtl.index=459,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.

#Epistasis
effectplot(F5All, pheno.col = 2, mname1="scaffold_26_2173503")
effectplot(F5All, pheno.col = 2, mname1="scaffold_44_2763411", mname2="scaffold_21_2195347")
effectplot(F5All, pheno.col = 2, mname1="scaffold_21_2195347", mname2="scaffold_44_2763411")

effectplot(F5All, pheno.col = 2, mname1="scaffold_44_3044568", mname2="scaffold_21_2195347")
effectplot(F5All, pheno.col = 2, mname1="scaffold_21_2195347", mname2="scaffold_44_3044568")

results <- mqmpermutation(augdata1, pheno.col=2, scanfunction=mqmscan, cofactors=homemadecofactors, n.perm=10, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(results)
summary(sig.results)
mqmplot.permutations(results)




####Auto#### 9, 10, 11
#augdata1<-fill.geno(F5Omit)
augdata1<-augdata2
augdata2<-augdata1
auto.F2Complete<-auto.F2
auto.F2True<-auto.F2

gt <- geno.table(F5All)
todrop3 <- rownames(gt[gt$P.value < 1e-300,])
F5Omit <- drop.markers(F5All, todrop3)
augdata1<-fill.geno(F5Omit)

auto.F2 <- mqmautocofactors(augdata1,50) #assigns 250 cofactors across genome accounting for marker density (10 cofactors/chromosome)
#auto.F2 #this will be different each time you run it
homescan2<-mqmscan(cross=augdata1,cofactors=auto.F2,pheno.col=2,cofactor.significance=0.02,verbose=T,plot=T,model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan2)

results <- mqmpermutation(augdata1, pheno.col=2, scanfunction=mqmscan, cofactors=auto.F2, n.perm=10, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(results)
summary(sig.results)
mqmplot.permutations(results)



find.marker(augdata1,chr=7,pos=20)
effectplot(F5All, pheno.col=2, mname1="scaffold_21_2195347",var.flag=c("group"))

find.marker(augdata1,chr="10b",pos=15)
effectplot(F5All, pheno.col=2, mname1="scaffold_75_1746024",var.flag=c("group"))

find.marker(augdata1,chr="13a",pos=40)
effectplot(F5All, pheno.col=2, mname1="scaffold_38_2529240",var.flag=c("group"))
#scaffold_38_1591063	scaffold_38_608518	scaffold_38_2111260	scaffold_38_2348122	scaffold_38_2529240

effectplot(F5All, pheno.col = 2, mname1="scaffold_44_2763411")
effectplot(F5All, pheno.col = 2, mname1="scaffold_38_2348122", mname2="scaffold_21_2195347")
effectplot(F5All, pheno.col = 2, mname1="scaffold_21_2195347", mname2="scaffold_38_2348122")
effectplot(F5All, pheno.col = 2, mname1="scaffold_21_2195347", mname2="scaffold_75_1746024")
effectplot(F5All, pheno.col = 2, mname1="scaffold_75_1746024", mname2="scaffold_21_2195347")

results <- mqmpermutation(augdata1, pheno.col=8, scanfunction=mqmscan, cofactors=auto.F2, n.perm=100, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(results)
summary(sig.results)



find.markerindex(augdata1,find.marker(augdata1,chr="13a",pos=20))
bayesint(homescan2,"13a",qtl.index=459,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.
find.marker(augdata1,chr=15,pos=35)
find.marker(augdata1,chr=18,pos=20)
effectplot(F5All, pheno.col=4, mname1="scaffold_23_1502215",var.flag=c("group"))
effectplot(F5All, pheno.col=4, mname1="scaffold_6_9719709",var.flag=c("group"))

find.markerindex(augdata1,find.marker(augdata1,chr=12,pos=35))
bayesint(homescan2,12,qtl.index=427,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.
find.marker(augdata1,chr="13a",pos=20)
effectplot(F5All, pheno.col=4, mname1="scaffold_14_6883973",var.flag=c("group"))
effectplot(F5All, pheno.col=5, mname1="scaffold_14_6883973",var.flag=c("group"))


####Epistasis####
EpiFull <- calc.genoprob(F5All, step=2.5, err=0.001)
#Depth
outEpiFullDepth <- scantwo(EpiFull, pheno.col=10, verbose=FALSE)
outEpiFullDepthRedCOMP <- scantwo(EpiFull, pheno.col=7, chr=c(5,7,11,12,16,17,22), verbose=FALSE)
#save(outHeartEpiFullDepthP, file="outHeartEpiFullDepthP.rda")
#outHeartEpiFullRatio; #outHeartEpiFullSize; #outHeartEpiFullDepth

plot(outEpiFullDepthRedCOMP, lower="cond-int")

outEpiFullDepthP <- scantwo(EpiFull, n.perm=100, pheno.col=5, chr=c(1,5,6,8,9,"10a",11,15,16,17,18,19,21), verbose=T)
outEpiFullDepthPFull <- scantwo(EpiFull, n.perm=100, pheno.col=5, verbose=T)
summary(outEpiFullDepthPFull)

summary(outEpiFullDepthRedCOMP, thresholds = c(9.67,7.84,6.30,6.85,3.69)+1) #5%

summary(outEpiFullDepth, thresholds = c(9.67,7.84,6.30,6.85,3.69)) #5%
summary(outEpiFullDepth, thresholds = c(9.32,6.98,5.99,6.71,3.42)) #10%

mar<-find.marker(F5All,chr=c(7,17),pos=c(20,20))
#mar<-find.marker(F5All,chr=c(8,11),pos=c(22.5,10))
#mar<-find.marker(F5All,chr=c(15,18),pos=c(2.5,30))


effectplot(F5All, pheno.col = 10, mname1=mar[2])
geno.crosstab(cross = F5All, mname1 = mar[1], mname2 = mar[2])
