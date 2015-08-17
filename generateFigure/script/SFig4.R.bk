# script for generating subfigures in the Supplementary Figure 4.
# data generated in "../../analysis/AlexandrovEtAl/result/" is necessary.



if (!file.exists("../result/SFig4")) {
dir.create("../result/SFig4")
}


library(pmsignature)

# dummy signature estimation for generating class
inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt", package="pmsignature");
G <- readMPFile(inputFile, numBases = 3, type = "full");
Param <- getPMSignature(G, K = 2);

# reading the signature from Alexandrov et al. Nature 2013
nature2013_sig_raw <- read.table("../data/AlexandrovEtAl_signatures.txt", header=T, sep="\t");
nature2013_sig <- t(nature2013_sig_raw[,4:30]);
colnames(nature2013_sig) <- nature2013_sig_raw[,3];

Param@signatureFeatureDistribution[1,1,] <- nature2013_sig[which(rownames(nature2013_sig) == "Signature.14"),];
Param@signatureFeatureDistribution[2,1,] <- nature2013_sig[which(rownames(nature2013_sig) == "Signature.14"),];

visPMSignature(Param, 1, isScale = TRUE)
ggsave("../result/SFig4/example96_Glioma-Low-Grade.eps", width=12, height=3, units="in");


##########

G <- readMPFile("../../analysis/AlexandrovEtAl/result/MPFormat/Glioma-Low-Grade.mp.txt.gz", numBases = 5, trDir = TRUE);


##########
# membership plot
load("../../analysis/AlexandrovEtAl/result/Param_ind5_dir/Glioma-Low-Grade.3.Rdata");
visMembership(G, resultForSave[[1]], ylog = FALSE, colourBrewer = "Set2", toSample = 100);
ggsave("../result/SFig4/Glioma-Low-Grade_membership.eps", height = 4, width = 12, pointsize = 9);

visMembership(G, resultForSave[[1]], ylog = TRUE, colourBrewer = "Set2", toSample = 100);
ggsave("../result/SFig4/Glioma-Low-Grade_membership_log.eps", height = 4, width = 12, pointsize = 9);
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
    visPMSignature(Params[[l1]], aInd, charSize = 1, isScale = TRUE);
  }
  
  outputName <- paste("../result/SFig4/Glioma-Low-Grade_removed_signature_K", (l1 + 1), ".eps", sep="");
  dev.copy2eps(file=outputName, height = 2.5, width = 15, pointsize = 9);
}


##########
# original signature

#' create the matrix showing the relationships between cancer type 
#' and the specified numer of mutation signatures
type2sigNum <- read.table("../data/AlexandrovEtAl_sigNum.txt", sep=",", header=FALSE);

K <- type2sigNum[type2sigNum[,1] == "Glioma-Low-Grade", 2];
inputRdata <- paste("../../analysis/AlexandrovEtAl/result/Param_ind5_dir/Glioma-Low-Grade.", K, ".Rdata", sep="");
  
load(inputRdata);

mat <- matrix(c(1:(K - 1), rep(0, K - 1)), 1, 4, byrow = TRUE);
layout(mat);

for (k in 1:(K - 1)) {
  visPMSignature(resultForSave[[1]], k, charSize = 1, isScale = TRUE);
}

outputName <- paste("../result/SFig4/Glioma-Low-Grade_original_signature_K", K, ".eps", sep = "");
dev.copy2eps(file=outputName, height = 2.5, width = 15, pointsize = 9);





