# goal: generating several figures from the results of PMSignature on Hoang et al. (2013) data.

# 2015 3/1 Y.S
# The size of legend is influenced by the plot device size in RStudio.
# I trid to fix this probelm, but I couldn't. 
# When make the figures, please enlarge the device so that the legend sizes in the figures are natural...

.pardefault <- par(no.readonly = TRUE);
##########
# functions

library(pmsignature);
library(RColorBrewer);

PMSSigPlot <- function(K) {
  
  fileName <- paste("../..//analysis/UCUT/result/Param_ind/", as.character(K), ".Rdata", sep="");
  load(fileName);

  # prepare matrix for drawing multiple figures
  mat <- matrix(c(1:3, rep(K, 3)), 2, 3, byrow = TRUE);
  if (K <= 3) {
    for (k in K:3) {
      mat[1, k] <- 0;
    }
  }
  layout(mat);
  
  # memorize the current parameter and change the margin slightly
  tempMar <- par("mar");
  par(mar = 0.4 * tempMar);
 
  
  F <- resultForSave[[1]]@signatureFeatureDistribution;
  
  # check the index corresponding each signature types
  AA_ind <- NULL;
  AA_like_ind <- NULL;
  APOBEC_ind <- NULL;
  for (k in 1:(K - 1)) {
    # AA
    if (F[k,1,4] > 0.8 & F[k,3,2] > 0.5 & F[k,4,3] > 0.5 & F[k,6,2] > 0.6) {
      AA_ind <- k;
    } else if (sum(F[k,1,1:3]) > 0.8 & F[k,3,4] > 0.7) {
      APOBEC_ind <- k;
    } else if (F[k,1,4] > 0.6 & F[k,6,2] > 0.6) {
      AA_like_ind <- k;
    }
  }
  
  
  ##########
  # here, draw the figure for each signature and set the index and color vector 
  # for the later membership barplot figrues
  
  mainSigInd <- c();
  remainSigInd <- 1:(K - 1);
  Q_col <- c();
  Q_legend <- c();

  if (!is.null(AA_ind)) {
    visPMSignature(resultForSave[[1]], AA_ind);
    mainSigInd <- c(mainSigInd, AA_ind);
    remainSigInd <- remainSigInd[-which(remainSigInd == AA_ind)];
    Q_col <- c(Q_col, brewer.pal(4, "Set2")[1]);
    Q_legend <- c(Q_legend, "AA");
  }
  if (!is.null(AA_like_ind)) {
    visPMSignature(resultForSave[[1]], AA_like_ind);
    mainSigInd <- c(mainSigInd, AA_like_ind);
    remainSigInd <- remainSigInd[-which(remainSigInd == AA_like_ind)];
    Q_col <- c(Q_col, brewer.pal(4, "Set2")[2]);
    Q_legend <- c(Q_legend, "AA_like");
    
  }  
  if (!is.null(APOBEC_ind)) {
    visPMSignature(resultForSave[[1]], APOBEC_ind);
    mainSigInd <- c(mainSigInd, APOBEC_ind);
    remainSigInd <- remainSigInd[-which(remainSigInd == APOBEC_ind)];
    Q_col <- c(Q_col, brewer.pal(4, "Set2")[3]);
    Q_legend <- c(Q_legend, "APOBEC");
  }  
  
  # lastly, draw the signature for the remaining (non AA, AA_like, APOBEC) signatures
  for (k in remainSigInd) {
    visPMSignature(resultForSave[[1]], k);
  }
  ##########
  
  # for background signatures 
  remainSigInd <- c(remainSigInd, K);
  Q_col <- c(Q_col, brewer.pal(4, "Set2")[4], K);
  Q_legend <- c(Q_legend, "BG");
  
  par(mar= 1 * tempMar);
  par(mar = c(5.1, 6.1, 4.1, 24.1), xpd = TRUE);

  # draw the membership parameter barplot
  Q <- resultForSave[[1]]@sampleSignatureDistribution[,c(mainSigInd, remainSigInd)];
  barplot(t(Q), col=Q_col, xlab = "cancer genomes", ylab = "membership parameters", 
          cex.axis = 1.5, cex.lab = 1.5);
  legend("topright", inset=c(-0.12, 0), fill = Q_col, legend = Q_legend, text.width = 3, y.intersp = 1.2);
  
  par(.pardefault);
}


##########

PMSSigPlot(2);
outputName <- "../result/SText2/UCUT_signature_K2.eps";
dev.copy2eps(file=outputName, height = 6, width = 15, pointsize = 9);
par(.pardefault);
dev.off();

PMSSigPlot(3);
outputName <- "../result/SText2/UCUT_signature_K3.eps";
dev.copy2eps(file=outputName, height = 6, width = 15, pointsize = 9);
par(.pardefault);
dev.off();

PMSSigPlot(4);
outputName <- "../result/SText2/UCUT_signature_K4.eps";
dev.copy2eps(file=outputName, height = 6, width = 15, pointsize = 9);
par(.pardefault);
dev.off();

# PMSSigPlot(5);
# outputName <- "../../supp/UTUC_signature_K5.eps";
# dev.copy2eps(file=outputName, height = 4, width = 15, pointsize = 9);
# par(.pardefault);



##########
# draw the figure showing the relationships between the AA and AA_like signature

fileName <- paste("../../analysis/UCUT/result/Param_ind/", as.character(4), ".Rdata", sep="");
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
  theme(legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        axis.text.x = element_text(size=rel(1.2)),
        axis.title.x = element_text(size=rel(1.2)),
        axis.text.y = element_text(size=rel(1.2)),
        axis.title.y = element_text(size=rel(1.2)))


ggsave("../result/SText2/UCUT_AA_AAlike_cor.eps", width=8, height=8);







