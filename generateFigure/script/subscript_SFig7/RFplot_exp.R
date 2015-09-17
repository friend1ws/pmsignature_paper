
getDFForRF <- function(G, sampleInd, typeVec, typeName) {
  
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  targetChr <- paste("chr", c(1:22, "X", "Y"), sep="")
  seqlen_chr <- GenomeInfoDb::seqlengths(genome)[targetChr];
  
  genomeStartInd <- cumsum(as.numeric(seqlen_chr)) - as.numeric(seqlen_chr)
  names(genomeStartInd) <- targetChr
  
  absPos <- genomeStartInd[G@mutationPosition[G@mutationPosition[,3] == sampleInd, 1]] + G@mutationPosition[G@mutationPosition[,3] == sampleInd, 2]
  mutType <- factor(sapply(G@mutationPosition[G@mutationPosition[,3] == sampleInd, 4], 
                    function(x) {if (x %in% typeVec) typeName else "other"}), levels = c(typeName, "other"))
  
  mutOrder <- order(absPos)
  absPos_sort <- absPos[mutOrder]
  mutType_sort <- mutType[mutOrder]
  
  dist1 <- abs(absPos_sort - c(0, absPos_sort[1:(length(absPos_sort)-1)]))
  dist2 <- abs(absPos_sort - c(absPos_sort[2:length(absPos_sort)], 0))
  dist <- apply(cbind(dist1, dist2), 1, min)
  
  return(data.frame(mutind = 1:length(absPos_sort), pos = absPos_sort, type = mutType_sort, log_dist = log10(dist)))
  
}




if (!file.exists("../result/Fig_revise2")) {
  dir.create("../result/Fig_revise2")
}

##########

G <- readMPFile("../../analysis/AlexandrovEtAl/result/MPFormat/Glioma-Low-Grade.mp.txt.gz", numBases = 3)
mutNum <- getMutNum(G)
sampleInd <- which(mutNum[,2] == max(mutNum[,2]))

mutinfo <- getDFForRF(G, sampleInd, c(4, 8, 12, 16), "Np[C>[A]]pT")

ggplot(mutinfo, aes(x = mutind, y = log_dist, colour = type)) + 
  geom_point() +
  xlab("ordered mutation number") +
  ylab("minimum log distance") +
  ggtitle(G@sampleList[sampleInd]) + 
  labs(fill="#sample", colour="#sample") +
  theme_bw() +
  theme(plot.title = element_text(size=rel(1.5)),
        axis.title = element_text(size = rel(1.2)), 
        axis.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(face="bold", size=rel(1.2)));


ggsave("../result/Fig_revise2/LGG_rf.eps", width=15, height=6, units="in");




####


G <- readMPFile("../../analysis/AlexandrovEtAl/result/MPFormat/Breast.mp.txt.gz", numBases = 3);

# mutNum <- getMutNum(G)
# sampleInd <- which(mutNum[,2] == max(mutNum[,2]))
sampleInd <- which(G@sampleList == "PD4103a")

mutinfo <- getDFForRF(G, 16, c(13:16, 29:32, 45:48), "Tp[C>[AGT]]pN")


ggplot(mutinfo, aes(x = mutind, y = log_dist, colour = type)) + 
  geom_point() +
  xlab("ordered mutation number") +
  ylab("minimum log distance") +
  ggtitle(G@sampleList[sampleInd]) + 
  labs(fill="#sample", colour="#sample") +
  theme_bw() +
  theme(plot.title = element_text(size=rel(1.5)),
        axis.title = element_text(size = rel(1.2)), 
        axis.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(face="bold", size=rel(1.2)));

ggsave("../result/Fig_revise2/Breast_rf.eps", width=15, height=6, units="in");



