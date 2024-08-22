####Fine Map Plots####
library(ggplot2); library(dplyr)
setwd("/Users/aconith/Library/CloudStorage/OneDrive-DePaulUniversity/DePaul/Research/Research/ParaSWriting/Github/FineMap")

####LG7####
FullMap<-read.csv("LG7FineMapDataT_PS_noaug.csv")

FineMap_Plot <- ggplot(FullMap, aes(x=Pos, y=AA.BB)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values=c("#ca0020")) +
  scale_fill_manual(values=c("#ca0020")) #+
  facet_wrap(~ Trait, scales = "free", labeller=labeller(Trait = MyLabels), nrow = 1)

FineMap_PlotSE <- FineMap_Plot + theme_bw() + annotate("rect", xmin = 43819923, xmax = 47752984, ymin = -Inf, ymax = Inf, fill = "#404040", alpha = 0.2) +
  geom_ribbon(data = FullMap, alpha=0.25, aes(ymin=AA.BB-AA_se, ymax=AA.BB+AA_se), col=NA) +
  geom_ribbon(data = FullMap, alpha=0.25, aes(ymin=AA.BB-BB_se, ymax=AA.BB+BB_se), col=NA) +
  geom_line(size=0.9) +
  labs(x = "Marker position on LG7 (Mbs)", y = "Average Phenotypic Effect") +
  theme(legend.position = "none", aspect.ratio = 1, legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))
  
FineMap_PlotSE +
  geom_vline(xintercept = 48423569, col="purple") + #DYM
  geom_vline(xintercept = 46670572, col="dark green") #NOTCH1a
  
####LG13####
FullMap<-read.csv("LG13FineMapDataT_PS_noaug.csv")

FineMap_Plot <- ggplot(FullMap, aes(x=Pos, y=AA.BB)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values=c("#ca0020")) +
  scale_fill_manual(values=c("#ca0020"))

FineMap_PlotSE <- FineMap_Plot + theme_bw() + annotate("rect", xmin = 12104957, xmax = 19425681, ymin = -Inf, ymax = Inf, fill = "#404040", alpha = 0.2) +
  geom_ribbon(data = FullMap, alpha=0.25, aes(ymin=AA.BB-AA_se, ymax=AA.BB+AA_se), col=NA) +
  geom_ribbon(data = FullMap, alpha=0.25, aes(ymin=AA.BB-BB_se, ymax=AA.BB+BB_se), col=NA) +
  geom_line(size=0.9) +
  labs(x = "Marker position on LG13 (Mbs)", y = "Average Phenotypic Effect") +
  theme(legend.position = "none", aspect.ratio = 1, legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))

FineMap_PlotSE +
  geom_vline(xintercept = 12763810, col="purple") #ADAM12
