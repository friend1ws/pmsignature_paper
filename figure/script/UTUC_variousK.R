# goal: generating several figures from the results of PMSignature on Hoang et al. (2013) data.


.pardefault <- par(no.readonly = TRUE);
##########
# functions

library(pmsignature);
library(RColorBrewer);

PMSSigPlot <- function(K) {
  
  fileName <- paste("../../../UTUC/result/Param_ind/", as.character(K), ".Rdata", sep="");
  load(fileName);

  mat <- matrix(c(1:4, rep(K, 4)), 2, 4, byrow = TRUE);
  
  if (K <= 4) {
    for (k in K:4) {
      mat[1, k] <- 0;
    }
  }
  
  layout(mat);
  tempMar <- par("mar");
  par(mar= 0.4 * tempMar);
  for (k in 1:(K-1)) {
    visPMSignature(resultForSave[[1]], k);
  }
  
  par(mar= 1.3 * tempMar);
  Q <- resultForSave[[1]]@sampleSignatureDistribution;
  d <- dist(Q);
  h <- hclust(d);
  barplot(t(Q)[,h$order], col=brewer.pal(max(K, 3), "Set2"), xlab = "cancer genomes", ylab = "membership parameters", cex.axis = 1.5, cex.lab = 1.5);
  
  par(.pardefault);
}


##########

PMSSigPlot(2);
outputName <- "../../result/UTUC_signature_K2.eps";
dev.copy2eps(file=outputName, height = 7, width = 15, pointsize = 9);
par(.pardefault);

PMSSigPlot(3);
outputName <- "../../result/UTUC_signature_K3.eps";
dev.copy2eps(file=outputName, height = 7, width = 15, pointsize = 9);
par(.pardefault);

PMSSigPlot(4);
outputName <- "../../result/UTUC_signature_K4.eps";
dev.copy2eps(file=outputName, height = 7, width = 15, pointsize = 9);
par(.pardefault);

PMSSigPlot(5);
outputName <- "../../result/UTUC_signature_K5.eps";
dev.copy2eps(file=outputName, height = 7, width = 15, pointsize = 9);
par(.pardefault);




fileName <- paste("../../../UTUC/result/Param_ind/", as.character(4), ".Rdata", sep="");
load(fileName);

Q <- resultForSave[[1]]@sampleSignatureDistribution;
F <- resultForSave[[1]]@signatureFeatureDistribution;

for (k in 1:3) {
  # AA
  if (F[k,1,4] > 0.8 & F[k,3,2] > 0.5 & F[k,4,3] > 0.5 & F[k,6,2] > 0.6) {
    AA_ind <- k;
  } else if (F[k,1,4] > 0.6 & F[k,6,2] > 0.6) {
    AA_like_ind <- k;
  }
}

AAmem <- data.frame(AA_like = Q[,AA_ind], AA = Q[,AA_like_ind]);
ggplot(AAmem, aes(x=AA, y=AA_like)) + geom_point(size=3, col="blue") +
  theme(legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        axis.text.x = element_text(size=rel(1)),
        axis.title.x = element_text(size=rel(1)),
        axis.text.y = element_text(size=rel(1)),
        axis.title.y = element_text(size=rel(1)))


ggsave("../../result/UTUC_AA_AAlike_cor.eps", width=5, height=5);







