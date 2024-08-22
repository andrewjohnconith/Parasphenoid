library(qtl); library(plyr); library(plotrix)
setwd("/Users/aconith/Dropbox/Amherst/Post-Doc/Map/ReAnalyze/MapPhenos/Analysis/Bones/Parasphenoid/FineMap/LG13/")
PS<-read.csv("PsPCScores_20240626.csv")
Map<-read.csv("LG13_Stringent.csv")
PSMap<-merge(PS, Map, by="HybridID", all=T)

write.csv(PSMap, "LG13Map_PSpc.csv")

LG13 <- read.cross(format="csv",file="LG13Map_PSpcEd.csv",
                           na.strings="NA",genotypes=c("AA","AB","BB"),
                           alleles=c("A","B"),convertXdata=TRUE)

#Check sample sizes
apply(LG13$geno$`13`$data, 2, table)
  #

#ScafNames<-names(LG3$geno$`3`$map)
ScafNames<-rownames(geno.table(LG13))[as.logical(apply(X = geno.table(LG13,13)[,3:5]>=10,MARGIN = 1,FUN = mean)==1)]
EffResult<-matrix(data = NA, nrow = length(ScafNames), ncol = 6)
         
for (i in 1:length(ScafNames)){
  effect<-effectplot(LG13, pheno.col=2, mname1=ScafNames[i], var.flag=c("group"), draw = F) #effect plot at qtl
  EffResult[i,]<-c(effect$Means, effect$SEs)
}


CombinedDataEffectE<-cbind.data.frame(ScafNames,as.numeric(sub("^[^_]*_", "", ScafNames)), EffResult,EffResult[,3]-EffResult[,1])
colnames(CombinedDataEffectE)<-c("Scaff", "Pos", "AA_mean", "AB_mean", "BB_mean", "AA_se", "AB_se", "BB_se", "AA.BB")

CombinedDataEffectEo<-CombinedDataEffectE[order(CombinedDataEffectE$Pos),]

write.csv(CombinedDataEffectEo,"LG13FineMapDataT_PS_noaug.csv", quote = F, row.names = F)
#write.csv(CombinedDataEffectEo,"LG19FineMapDataT_WLResidual_Stringent_noaug.csv", quote = F, row.names = F)


##Plot
setwd("/Users/aconith/Library/CloudStorage/OneDrive-DePaulUniversity/DePaul/Research/Research/ParaSWriting/Github/FineMap/Reduced/LG13")
CombinedDataEffectEo<-read.csv("LG13FineMapDataT_PS_noaugRED.csv")
FSTs<-read.csv("MyLiftOverScaffold13RED.csv")

MapPP<-CombinedDataEffectEo
MapP<-MapPP[complete.cases(MapPP$AA.BB) & complete.cases(MapPP$AA_se) & complete.cases(MapPP$BB_se),]

#AA-BB
par(pty='s', mar = c(5,5,2,5))
plot(MapP$Pos, MapP$AA.BB+MapP$AA_se, xlim=c(9996786,21283499), ylim = range(c(MapP$AA.BB,MapP$AA.BB+MapP$AA_se,MapP$AA.BB+MapP$BB_se,MapP$AA.BB-MapP$AA_se,MapP$AA.BB-MapP$BB_se), na.rm = T), type='n', lwd=1, lty=2, xlab = 'LG13', ylab = 'Avg. Pheno. Effect')
polygon(c(MapP$Pos,rev(MapP$Pos)),c(MapP$AA.BB+MapP$AA_se,rev(MapP$AA.BB-MapP$AA_se)), col = "purple", border = "purple", lwd=2)
polygon(c(MapP$Pos,rev(MapP$Pos)),c(MapP$AA.BB+MapP$BB_se,rev(MapP$AA.BB-MapP$BB_se)), col = "dark blue", border = "dark blue", lwd=2)
lines(MapP$Pos, MapP$AA.BB, col="light blue", lwd=1.5)
par(new = T)
plot(FSTs$startLG, FSTs$Fst, cex=0.5, ylim=c(0,1), xlim=c(9996786,21283499), col="blue",  ylab="", xlab="", axes=F, pch=19)


#AB-BB
par(pty='s', mar = c(5,5,2,5))
plot(MapP$Pos, MapP$AB.BB+MapP$AB_se, xlim=c(9996786,21283499), ylim = range(c(MapP$AB.BB,MapP$AB.BB+MapP$AB_se,MapP$AB.BB+MapP$BB_se,MapP$AB.BB-MapP$AB_se,MapP$AB.BB-MapP$BB_se), na.rm = T), type='n', lwd=1, lty=2, xlab = 'LG13', ylab = 'Avg. Pheno. Effect')
polygon(c(MapP$Pos,rev(MapP$Pos)),c(MapP$AB.BB+MapP$AB_se,rev(MapP$AB.BB-MapP$AB_se)), col = "purple", border = "purple", lwd=2)
polygon(c(MapP$Pos,rev(MapP$Pos)),c(MapP$AB.BB+MapP$BB_se,rev(MapP$AB.BB-MapP$BB_se)), col = "dark blue", border = "dark blue", lwd=2)
lines(MapP$Pos, MapP$AB.BB, col="light blue", lwd=1.5)
par(new = T)
plot(FSTs$startLG, FSTs$Fst, cex=0.5, ylim=c(0,1), xlim=c(9996786,21283499), col="blue",  ylab="", xlab="", axes=F, pch=19)

AC_13<-read.csv("/Users/aconith/Library/CloudStorage/OneDrive-DePaulUniversity/DePaul/Research/Research/PS2/Liftover/AstCalLG13.csv")
plot(AC_13$wStart,AC_13$Lab, cex=0.5)


abline(h=0)

abline(v=11992751, col="black") ##13a length
abline(v=27951429, col="black")

abline(v=12104957, col="red") #Bayes Int
abline(v=19425681, col="red")
abline(v=19425731, col="orange") #Bayes Peak

abline(v=26978820, col="purple") #calm2
abline(v=26983366, col="purple") #calm2

abline(v=19425681, col="blue") #fine map peak

abline(v=12763810, col="gray") #Adam12
abline(v=12851929, col="gray") #Adam12
