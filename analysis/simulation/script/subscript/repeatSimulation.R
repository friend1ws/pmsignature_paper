library(pmsignature)

mutationNum <- as.integer(commandArgs()[5]);
sampleNum <- as.integer(commandArgs()[6]);
alpha <- as.numeric(commandArgs()[7]);
gamma <- as.numeric(commandArgs()[8]);
outputFile <- commandArgs()[9];

fdim <- c(6, 4, 4, 4, 4);

Dists <- matrix(0, 100, 4);

for (i in 1:100) {
    G_sim <- makeSimData(type = "independent", numBases = 5, trDir = FALSE, K = 5, sampleNum = sampleNum, mutationNum = mutationNum, param_alpha = alpha, param_gamma = gamma, isBG = TRUE);
    bgprob <- readBGFile(G_sim[[1]]);
    Param_est <- getPMSignature(G_sim[[1]], K = 5, BG = bgprob);
    Dists[i, ] <- getCosDistance(Param_est@signatureFeatureDistribution, G_sim[[2]], Param_est@possibleFeatures);
}

write.table(Dists, file=outputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t");


