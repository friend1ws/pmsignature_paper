library(pmsignature)

mutationNum <- 25;
featureNum <- as.integer(commandArgs()[5]);
sampleNum <- as.integer(commandArgs()[6]);
alpha <- as.numeric(commandArgs()[7]);
gamma <- as.numeric(commandArgs()[8]);
outputFile <- commandArgs()[9];

fdim <- c(6, 4, 4, 4, 4, rep(2, featureNum));

Dists <- matrix(0, 100, 4);

for (i in 1:100) {
  G_sim <- makeSimData(type = "custom", K = 4, fdim = fdim, sampleNum = sampleNum, mutationNum = mutationNum, param_alpha = alpha, param_gamma = gamma);
  Param_est <- getPMSignature(G_sim[[1]], K = 4);
  Dists[i, ] <- getCosDistance(Param_est@signatureFeatureDistribution, G_sim[[2]], Param_est@possibleFeatures);
}


write.table(Dists, file=outputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
