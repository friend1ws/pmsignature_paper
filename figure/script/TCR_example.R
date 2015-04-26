##########
# new figure 1(a)

inputFile <- "../../AlexandrovEtAl/result/MPFormat/Liver.mp.txt.gz";
G <- readMPFile(inputFile, numBases = 5, trDir = TRUE, type = "full");
BG_prob <- readBGFile(G);
Param <- getPMSignature(G, K = 3, BG = BG_prob);

if (sum(Param@signatureFeatureDistribution[1,1,(256 * 10 + 1):(256 * 11)]) > 0.2) {
  visPMSignature(Param, 1);
} else {
  visPMSignature(Param, 2);
}

ggsave("../../manuscript/example3076_TCR.eps", width=12, height=5, units="in");


# new figure 1(b)
G <- readMPFile(inputFile, numBases = 5, trDir = TRUE, type = "independent");
BG_prob <- readBGFile(G);
Param <- getPMSignature(G, K = 3, BG = BG_prob);


if (Param@signatureFeatureDistribution[1, 1, 5] > 0.6) {
  visPMSignature(Param, 1, isScale = TRUE, alpha = 2); 
} else {
  visPMSignature(Param, 2, isScale = TRUE, alpha = 2);
}

ggsave("../../manuscript/example18_TCR.eps", width=12, height=5, units="in");




