
PMSSummaryPlot <- function(type) {
  
  Ls <- rep(0, 5);
  Bs <- rep(0, 5);
  Cs <- rep(0, 5);
  
  for(K in 2:6) {
    fileName <- paste("example/", type, "/", K, "_1_1.Rdata", sep="");
    load(fileName);
    Ls[(K - 1)] <- lret[[1]][[3]];
    Cs[(K - 1)] <- max(cor(lret[[1]][[2]]) - 2 * diag(rep(1, K)));
    
    bF <- lret[[2]][[1]];
    bFs <- matrix(0, 100, K - 1);
    for(b in 1:100) {
      for(k in 1:(K - 1)) {
        bFs[b, k] <- sum(bF[b,k,,]);
      }
      Bs[(K - 1)] <- sum(bFs) / ((K - 1) * 100);
    } 
    
  }
  
  par(mfrow=c(1, 3));
  plot(2:6, Ls, xlab="#signature", ylab="log-likelihood", type="l");
  grid();
  plot(2:6, Bs, xlab="#signature", ylab="bootstrap-error", type="l");
  grid();
  plot(2:6, Cs, xlab="#signature", ylab="maximum cor. efficient", type="l");
  grid(); 
  
}

#' export eps width=1000, height=400
PMSSummaryPlot("hoang_ind5_dir"); # 8 * 4 inches


