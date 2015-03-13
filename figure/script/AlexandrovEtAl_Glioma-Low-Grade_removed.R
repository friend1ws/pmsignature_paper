library(pmsignature)

G <- readMPFile("../../AlexandrovEtAl/result/MPFormat/Glioma-Low-Grade.mp.txt.gz", numBases = 5, trDir = TRUE);

##########
# membership plot
load("../../AlexandrovEtAl/result/Param_ind5_dir/Glioma-Low-Grade.3.Rdata");
visMembership(G, resultForSave[[1]], ylog = TRUE, colourBrewer = "Set2", toSample = 100);
ggsave("../../supp/Glioma-Low-Grade.membership.eps", height = 4, width = 15, pointsize = 9);
##########


##########
# remove the sample with maximum number of mutation and check the signatures
mutNum <- getMutNum(G);
mutMaxInd <- which(mutNum[,2] == max(mutNum[,2]));
G@countData[3, G@countData[2,] == mutMaxInd] <- 0;
BG <- readBGFile(G);

Param2 <- getPMSignature(G, K = 2, BG = BG);
Param3 <- getPMSignature(G, K = 3, BG = BG);
Param4 <- getPMSignature(G, K = 4, BG = BG);
Param5 <- getPMSignature(G, K = 5, BG = BG);
Params <- c(Param2, Param3, Param4, Param5);

alignOrder <- alignmentSignature(Params);


tempMar <- par("mar");
par(mar = 0.4 * tempMar);

for (l1 in 1:length(Params)) {
  
  mat <- matrix(c(1:l1, rep(0, 4 - l1)), 1, 4, byrow = TRUE);
  layout(mat);
   
  for (l2 in 1:l1) {
    aInd <- which(alignOrder[[l1]] == l2);
    visPMSignature(Params[[l1]], aInd);
  }
  
  outputName <- paste("../../supp/Glioma-Low-Grade.removed_signature_K", (l1 + 1), ".eps", sep="");
  dev.copy2eps(file=outputName, height = 3, width = 15, pointsize = 9);
}





