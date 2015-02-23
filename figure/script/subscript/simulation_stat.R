
value <- c();
sigNum <- c();
type <- c();


for(K in 2:6) {
  fileName <- paste("../../../UTUC/result/Param_ind/", K, ".Rdata", sep="");
  load(fileName);
  
  sigNum <- c(sigNum, K, K, K);
  
  # likelihood
  value <- c(value, resultForSave[[1]]@loglikelihood);
  type <- c(type, "log-likelihood");
  
  # bootstrap standard error
  bF <- resultForSave[[2]][[1]];
  bFs <- matrix(0, 100, K - 1);
  for(b in 1:100) {
      for(k in 1:(K - 1)) {
        bFs[b, k] <- sum(bF[b,k,,]);
      }
  }
  value <- c(value, sum(bFs) / ((K - 1) * 100));
  type <- c(type, "bootstrap error");
  
  # correlation
  value <- c(value, max(cor(resultForSave[[1]]@sampleSignatureDistribution) - 2 * diag(rep(1, K))));
  type <- c(type, "max. corr.");
}

type <- factor(type, levels = c("log-likelihood", "bootstrap error", "max. corr."));
UTUC_stat <- data.frame(value = value, type = type, sigNum = sigNum);


ggplot(UTUC_stat, aes(x = sigNum, y = value)) +
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
        strip.text = element_text(face="bold", size=rel(1.2)));

ggsave("../../result/simulation_stat.eps", width=15, height=5, units = "in");


