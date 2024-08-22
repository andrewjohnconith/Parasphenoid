####Fine Map Plots####
library(ggplot2); library(dplyr)
setwd("/Users/aconith/Library/CloudStorage/OneDrive-DePaulUniversity/DePaul/Research/Research/ParaSWriting/Github/FineMap/Reduced/LG7")

####LG7####
FullMap<-read.csv("LG7FineMapDataT_PS_noaugRED1.csv")

FineMap_Plot <- ggplot(FullMap, aes(x=Pos, y=AA.BB)) +
  scale_color_manual(values=c("#ca0020")) +
  scale_fill_manual(values=c("#ca0020"))

FineMap_PlotSE <- FineMap_Plot + theme_bw() + annotate("rect", xmin = 43819868, xmax = 47752984, ymin = -Inf, ymax = Inf, fill = "#404040", alpha = 0.2) +
  #geom_ribbon(data = FullMap, alpha=0.25, aes(ymin=AA.BB-AA_se, ymax=AA.BB+AA_se), col=NA) +
  #geom_ribbon(data = FullMap, alpha=0.25, aes(ymin=AA.BB-BB_se, ymax=AA.BB+BB_se), col=NA) +
  geom_line(size=0.9) +
  labs(x = "Marker position on LG7 (Mbs)", y = "Average Phenotypic Effect") +
  theme(legend.position = "none", aspect.ratio = 1, legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))
  
FineMap_PlotSE +
  geom_vline(xintercept = 46644689, col="purple") + #DYM
  geom_vline(xintercept = 48423569, col="dark green") #NOTCH1a

####Fst LG7 Plot####
setwd("/Users/aconith/Library/CloudStorage/OneDrive-DePaulUniversity/DePaul/Research/Research/ParaSWriting/Github/FineMap/Reduced/LG7")
#Read in data
FullMap<-read.csv("MyLiftOverScaffold7RED1.csv")

#Insert a user defined alpha value
#To adjust the opacity range play with the
#'scale_alpha_continuous' and the 0 value below
FullMap$FstAlpha<-FullMap$Fst
FullMap$FstAlpha[which(FullMap$FstAlpha<0.6)]<-0

Zrange <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
FullMap$FstAlpha[which(FullMap$FstAlpha>0.6)]<-Zrange(FullMap$FstAlpha[which(FullMap$FstAlpha>0.6)])

FineMap_Fst <- filter(FullMap, !is.na(Fst)) %>%
  ggplot(aes(x=startLG, y=Fst))

FineMap_FstPoint <- FineMap_Fst + theme_bw() + annotate("rect", xmin = 43819868, xmax = 47753039, ymin = -Inf, ymax = Inf, fill = "#404040", alpha = 0.2) +
  geom_point(size=2) +
  #geom_point(size=2, aes(alpha=FstAlpha)) +
  #scale_alpha_continuous(range = c(0.2,1)) +
  labs(x = "Marker position on LG7 (Mbs)", y = "Fst") +
  theme(legend.position = "none", aspect.ratio = 1, legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))

FineMap_FstPointGene <- FineMap_FstPoint + 
  geom_hline(yintercept = 0.9, col="dark green", linetype="longdash") + #Z=2 position
  geom_hline(yintercept = 0.6, col="dark green", linetype="dashed") #Z=1 position

FineMap_FstPointGene



####LG13####
setwd("/Users/aconith/Library/CloudStorage/OneDrive-DePaulUniversity/DePaul/Research/Research/ParaSWriting/Github/FineMap/Reduced/LG13")
FullMap<-read.csv("LG13FineMapDataT_PS_noaugRED.csv")

FineMap_Plot <- ggplot(FullMap, aes(x=Pos, y=AA.BB)) +
  scale_color_manual(values=c("#ca0020")) +
  scale_fill_manual(values=c("#ca0020"))

FineMap_PlotSE <- FineMap_Plot + theme_bw() + annotate("rect", xmin = 12104957, xmax = 19425681, ymin = -Inf, ymax = Inf, fill = "#404040", alpha = 0.2) +
  #geom_ribbon(data = FullMap, alpha=0.25, aes(ymin=AA.BB-AA_se, ymax=AA.BB+AA_se), col=NA) +
  #geom_ribbon(data = FullMap, alpha=0.25, aes(ymin=AA.BB-BB_se, ymax=AA.BB+BB_se), col=NA) +
  geom_line(size=0.9) +
  labs(x = "Marker position on LG13 (Mbs)", y = "Average Phenotypic Effect") +
  theme(legend.position = "none", aspect.ratio = 1, legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))

FineMap_PlotSE +
  geom_vline(xintercept = 12763810, col="purple") #ADAM12

####Fst LG13 Plot####
setwd("/Users/aconith/Library/CloudStorage/OneDrive-DePaulUniversity/DePaul/Research/Research/ParaSWriting/Github/FineMap/Reduced/LG13")
#Read in data
FullMap<-read.csv("MyLiftOverScaffold13RED.csv")

#Insert a user defined alpha value
#To adjust the opacity range play with the
#'scale_alpha_continuous' and the 0 value below
FullMap$FstAlpha<-FullMap$Fst
FullMap$FstAlpha[which(FullMap$FstAlpha<0.6)]<-0

Zrange <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
FullMap$FstAlpha[which(FullMap$FstAlpha>0.6)]<-Zrange(FullMap$FstAlpha[which(FullMap$FstAlpha>0.6)])

FineMap_Fst <- filter(FullMap, !is.na(Fst)) %>%
  ggplot(aes(x=startLG, y=Fst))

FineMap_FstPoint <- FineMap_Fst + theme_bw() + annotate("rect", xmin = 12104957, xmax = 19425681, ymin = -Inf, ymax = Inf, fill = "#404040", alpha = 0.2) +
  geom_point(size=2) +
  #geom_point(size=2, aes(alpha=FstAlpha)) +
  #scale_alpha_continuous(range = c(0.2,1)) +
  labs(x = "Marker position on LG13 (Mbs)", y = "Fst") +
  theme(legend.position = "none", aspect.ratio = 1, legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))

FineMap_FstPointGene <- FineMap_FstPoint + 
  geom_hline(yintercept = 0.9, col="dark green", linetype="longdash") + #Z=2 position
  geom_hline(yintercept = 0.6, col="dark green", linetype="dashed") #Z=1 position

FineMap_FstPointGene

