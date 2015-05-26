# script for generating subfigures for Figure 1.
# data generated in "../../analysis/AlexandrovEtAl/result/" is necessary.


if (!file.exists("../result/Fig1")) {
  dir.create("../result/Fig1")
}

#' generating the subfigure of Figure 1 (mutation signature example of the full model)
##########
inputFile <- "../../analysis/AlexandrovEtAl/result/MPFormat/Liver.mp.txt.gz";
G <- readMPFile(inputFile, numBases = 5, trDir = TRUE, type = "full");
BG_prob <- readBGFile(G);
Param <- getPMSignature(G, K = 3, BG = BG_prob);

if (sum(Param@signatureFeatureDistribution[1,1,(256 * 10 + 1):(256 * 11)]) > 0.2) {
  visPMSignature(Param, 1);
} else {
  visPMSignature(Param, 2);
}

ggsave("../result/Fig1/example3076_TCR.eps", width=12, height=5, units="in");
##########


#' generating the subfiture of Figure 1 (mutation signature example of the independent model)
##########
G <- readMPFile(inputFile, numBases = 5, trDir = TRUE, type = "independent");
BG_prob <- readBGFile(G);
Param <- getPMSignature(G, K = 3, BG = BG_prob);


if (Param@signatureFeatureDistribution[1, 1, 5] > 0.6) {
  visPMSignature(Param, 1, isScale = TRUE, alpha = 2); 
} else {
  visPMSignature(Param, 2, isScale = TRUE, alpha = 2);
}

outputName <- "../result/Fig1/example18_TCR.eps";
dev.copy2eps(file=outputName, height = 6.66, width = 10);
##########



