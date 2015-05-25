library(pmsignature);
library(grid);


G_sim <- makeSimData(type = "independent", numBases = 5, trDir = FALSE, K = 5, sampleNum = 25, mutationNum = 100, isBG = TRUE);
brprob <- readBGFile(G_sim[[1]]);


value <- c();
sigNum <- c();
type <- c();

for (K in 2:8) {
  
  Param_est <- getPMSignature(G_sim[[1]], K = K, BG = brprob);
  Param_boot <- bootPMSignature(G_sim[[1]], Param0 = Param_est, bootNum = 100, BG = brprob);
  
  if (K == 5) {
    Param_est_5 <- Param_est;
  }
  
  sigNum <- c(sigNum, K, K, K);
  
  # likelihood
  value <- c(value, slot(Param_est, "loglikelihood"));
  type <- c(type, "log-likelihood");
  
  # bootstrap standard error
  bF <- Param_boot[[1]];
  bFs <- matrix(0, 100, K - 1);
  for(b in 1:100) {
    for(k in 1:(K - 1)) {
      bFs[b, k] <- sum(bF[b,k,,]);
    }
  }
  value <- c(value, sum(bFs) / ((K - 1) * 100));
  type <- c(type, "bootstrap error");
  
  # correlation
  Q <- slot(Param_est, "sampleSignatureDistribution");
  value <- c(value, max(cor(Q) - 2 * diag(rep(1, K))));
  type <- c(type, "max. corr.");
  
}


type <- factor(type, levels = c("log-likelihood", "bootstrap error", "max. corr."));
simulation_stat <- data.frame(value = value, type = type, sigNum = sigNum);


ggplot(simulation_stat, aes(x = sigNum, y = value)) +
  geom_line(stat = "identity") + 
  facet_wrap(~type, scales = "free") +
  xlab("#signature") +
  theme_bw() +
  theme(panel.margin = unit(1, "lines"),
        plot.title = element_text(size=rel(2)),
        axis.title = element_text(size = rel(1.2)), 
        axis.text.x = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(1.2), angle = 90, hjust = 0.5),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_line(colour="grey60", linetype="dashed"),
        strip.text = element_text(face="bold", size=rel(1.2)));

ggsave("../result/SText1/simulation_stat.eps", width=15, height=5, units = "in");




Param_true <- Param_est_5
Param_true@signatureFeatureDistribution <- G_sim[[2]];

#' export eps width=1000, height=400

removeList <- rep(1, 4);
relvect <- c();
par(mar=c(1, 0, 2, 0));
par(mfrow=c(2, 4));
for (i in 1:4) {
  visPMSignature(Param_true, i, charSize = 1)
  mtext(paste("true sig.", i),
        outer = FALSE,     # don't use outer margins
        side = 3,          # write at the top
        cex = 0.8,         # relative size of the character
        line = 0.5,        # write 0.5 lines outwards．
        col = "black");
  
  tempMin <- 1e10;
  F_true_i <- slot(Param_true, "signatureFeatureDistribution")[i,,];
  for (j in 1:4) {
    F_est_j <- slot(Param_est_5, "signatureFeatureDistribution")[j,,]
    dist_F_ij <- sum((F_true_i - F_est_j)^2);
    if (removeList[j] == 1 & dist_F_ij < tempMin) {
      tempJ <- j;
      tempMin <- dist_F_ij;
    }
  }
  relvect <- c(relvect, tempJ);
  removeList[tempJ] <- 0;
  
}

for (i in 1:length(relvect)) {
  visPMSignature(Param_est_5, relvect[i], charSize = 1)
  mtext(paste("estimated. sig.", i),
        outer = FALSE,     # don't use outer marginsく
        side = 3,          # write at the top
        cex = 0.8,         # relative size of the character
        line = 0.5,        #  write 0.5 lines outwards．．
        col = "black");
  
}


dev.copy2eps(file="../result/SText1/simulation_signature_example.eps", height = 5.5, width = 15);


