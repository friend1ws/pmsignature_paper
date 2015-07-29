# script for generating figures for revision
# data generated in "../../analysis/AlexandrovEtAl/result/" is necessary.


if (!file.exists("../result/Fig_revise1")) {
  dir.create("../result/Fig_revise1")
}


# dummy learning
inputFile <- "../../analysis/AlexandrovEtAl/result/MPFormat/Uterus.mp.txt.gz";
G <- readMPFile(inputFile, numBases = 3, trDir = FALSE, type = "full");
BG_prob <- readBGFile(G);
Param <- getPMSignature(G, K = 3, BG = BG_prob);

# substitute the signature obtained in Alexandrov et al., Nature, 2013
nature2013_sig_raw <- read.table("../data/AlexandrovEtAl_signatures.txt", header=T, sep="\t");
nature2013_sig <- t(nature2013_sig_raw[,4:30]);
colnames(nature2013_sig) <- nature2013_sig_raw[,3];
Param@signatureFeatureDistribution[1,1,] <- nature2013_sig["Signature.10",]

visPMSignature(Param, 1)

ggsave("../result/Fig_revise1/POLE_AlexandrovEtAl.eps", width=12, height=3, units="in");



inputFile <- "../../analysis/AlexandrovEtAl/result/MPFormat/Colorectum.mp.txt.gz";
G <- readMPFile(inputFile, numBases = 5, trDir = TRUE, type = "independent");
BG_prob <- readBGFile(G);
Param <- getPMSignature(G, K = 5, BG = BG_prob);


Sig <- Param@signatureFeatureDistribution
for (k in 1:4) {
  if (Sig[k,1,1] > 0.5 & Sig[k,3,4] > 0.5 & Sig[k,4,4] > 0.5) {
    POLE1 <- k
  }
  if (Sig[k,1,3] > 0.5 & Sig[k,3,4] > 0.5 & Sig[k,4,3] > 0.5) {
    POLE2 <- k
  }
}

visPMSignature(Param, POLE1)
ggsave("../result/Fig_revise1/POLE_pmsignature_1.eps", width=5, height=3, units="in")

visPMSignature(Param, POLE2)
ggsave("../result/Fig_revise1/POLE_pmsignature_2.eps", width=5, height=3, units="in")

