library(ggplot2);
library(dplyr);


mutNum <- c(10, 25, 50, 100, 250, 500, 1000);
sampleNum <- c(10, 25, 50, 100);
alpha <- c(0.5, 1, 2);
gamma <- c(0.5, 1, 2);

MutNum <- rep(0, 25200);
SampleNum <- rep(0, 25200);
Alpha <- rep(0, 25200);
Gamma <- rep(0, 25200);
Dist <- rep(0, 25200);

Ind <- 1:100;

for (i1 in 1:7) {
  for (i2 in 1:4) {
    for (i3 in 1:3) {
      for (i4 in 1:3) {
        a <- read.table(paste("../../../simulation/result/cosineDist/", mutNum[i1], "_", sampleNum[i2], "_", alpha[i3], "_", gamma[i4], ".txt", sep=""), sep="\t", header=FALSE);

        Dist[Ind] <- rowSums(a) / 4;
        MutNum[Ind] <- mutNum[i1];
        SampleNum[Ind] <- sampleNum[i2];
        Alpha[Ind] <- alpha[i3];
        Gamma[Ind] <- gamma[i4];
        
        Ind <- Ind + 100;
      }
    }
  }
}



simulation.result <- data.frame(mutNum=MutNum, sampleNum=SampleNum, alpha=Alpha, gamma=Gamma, distance=Dist);
summary_simulation.result <- simulation.result  %.% group_by(mutNum, sampleNum, alpha, gamma) %.% summarise(cosine_similality = mean(distance), seDist = sd(distance));
summary_simulation.result$seDist2 <- sapply(summary_simulation.result$cosine_similality + summary_simulation.result$seDist, min, 1) - summary_simulation.result$cosine_similality;


summary_simulation.result.1_1 <- summary_simulation.result %.% filter(alpha == 1 & gamma == 1);
summary_simulation.result.1_1$seDist2 <- sapply(summary_simulation.result.1_1$cosine_similality + summary_simulation.result.1_1$seDist, min, 1) - summary_simulation.result.1_1$cosine_similality;


pd <- position_dodge(0.2);
ggplot(summary_simulation.result, aes(x=factor(mutNum), y=cosine_similality, ymax=1, ymin = 0.3, fill=factor(sampleNum), colour=factor(sampleNum), group=factor(sampleNum))) +
  geom_line(position=pd) +
  geom_point(position=pd, size=2, shape=21) +
  geom_errorbar(aes(ymin=cosine_similality - seDist, ymax=cosine_similality + seDist2), width=0.01, position=pd) + 
  facet_grid(gamma ~ alpha, labeller = label_both) +
  xlab("#mutation") +
  ylab("mean cosine similarity") +
  labs(fill="#sample", colour="#sample") +
  theme_bw() +
  theme(plot.title = element_text(size=rel(2)),
        axis.title = element_text(size = rel(1.2)), 
        axis.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        strip.text= element_text(face="bold", size=rel(1.2)));

ggsave("../../result/simulation_result.eps", width=15, height=6);


G_sim <- makeSimData(type = "independent", numBases = 5, trDir = FALSE, K = 5, sampleNum = 25, mutationNum = 100, isBG = TRUE);
brprob <- readBGFile(G_sim[[1]]);

Ls <- rep(0, 7);
Bs <- rep(0, 7);
Cs <- rep(0, 7);

for (K in 2:8) {
  
  Param_est <- getPMSignature(G_sim[[1]], K = K, BG = brprob);
  Param_boot <- bootPMSignature(G_sim[[1]], Param0 = Param_est, bootNum = 100, BG = brprob);

  Ls[K - 1] <- slot(Param_est, "loglikelihood");
  bF <- Param_boot[[1]];
  bFs <- matrix(0, 100, K - 1);
  for(b in 1:100) {
    for(k in 1:(K - 1)) {
      bFs[b, k] <- sum(bF[b,k,,]);
    }
  }
  Bs[(K - 1)] <- sum(bFs) / ((K - 1) * 100);
  
  Q <- slot(Param_est, "sampleSignatureDistribution");
  
  Cs[(K - 1)] <- max(cor(Q) - 2 * diag(rep(1, K)));
}


#' export eps width=900, height=300
par(mfrow=c(1, 3));
plot(2:8, Ls, xlab="#signature", ylab="log-likelihood", type="l");
grid();
plot(2:8, Bs, xlab="#signature", ylab="bootstrap-error", type="l");
grid();
plot(2:8, Cs, xlab="#signature", ylab="maximum cor. efficient", type="l");
grid(); 

outputName <- "signature_summary.eps";
dev.copy2eps(file=outputName, height = 3, width = 9, pointsize = 12);
par(.pardefault);



Param_true <- Param_est
Param_true@signatureFeatureDistribution <- G_sim[[2]];

#' export eps width=1000, height=400

removeList <- rep(1, 4);
relvect <- c();
par(mar=c(1, 0, 2, 0));
par(mfrow=c(2, 4));
for (i in 1:4) {
  visPMSignature(Param_true, i)
  mtext(paste("true sig.", i),
        outer = FALSE,     # don't use outer margins
        side = 3,          # write at the top
        cex = 0.8,         # relative size of the character
        line = 0.5,        # write 0.5 lines outwards．
        col = "black");
  
  tempMin <- 1e10;
  F_true_i <- slot(Param_true, "signatureFeatureDistribution")[i,,];
  for (j in 1:4) {
    F_est_j <- slot(Param_est, "signatureFeatureDistribution")[j,,]
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
  visPMSignature(Param_est, relvect[i])
  mtext(paste("estimated. sig.", i),
        outer = FALSE,     # don't use outer marginsく
        side = 3,          # write at the top
        cex = 0.8,         # relative size of the character
        line = 0.5,        #  write 0.5 lines outwards．．
        col = "black");

}


outputName <- "signature_summary.eps";
dev.copy2eps(file=outputName, height = 3, width = 9, pointsize = 12);
par(.pardefault);
