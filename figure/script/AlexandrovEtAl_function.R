
##########
# functions for drawing figures for the Alexandrov et al analysis data


visPMS_ind <- function(vF, numBases, baseCol = NA, trDir, charSize = 1.2, scale = 0) {
  
  if (is.na(baseCol)) {
    gg_color_hue6 <- hcl(h = seq(15, 375, length = 7), l=65, c=100)[1:6]
    baseCol <- c(gg_color_hue6[3], gg_color_hue6[5], gg_color_hue6[2], gg_color_hue6[1], gg_color_hue6[4], gg_color_hue6[6]);
  }
  
  centerBase <- (1 + numBases) / 2;
  
  v1 <- vF[1,1:6];
  V2 <- vF[2:(numBases),1:4];
  A <- matrix(0, numBases, 4);
  B <- matrix(0, 4, 4);
  
  if (trDir == TRUE) {
    v3 <- vF[(numBases + 1),1:2];
  }
  
  for (l in 1:numBases) {
    if (l < centerBase) {
      A[l, ] <- V2[l, ];
    } else if (l > centerBase) {
      A[l, ] <- V2[l - 1, ];
    }
  }
  A[centerBase,2] <- sum(v1[1:3]);
  A[centerBase,4] <- sum(v1[4:6]);
  
  B[2, c(1, 3, 4)] <- v1[1:3] / sum(v1[1:3]);
  B[4, c(1, 2, 3)] <- v1[4:6] / sum(v1[4:6]);
  
  num2base <- c("A", "C", "G", "T");
  
  # frame();
  # plot.window(xlim=c(-0.05, 1.25 * numBases + 0.05), ylim=c(-0.25, 3.25));
  
  sizes <- 0.5 * (2 - apply(A, MARGIN = 1, FUN = function(x) {-sum(x * log2(x), na.rm = TRUE)}));
  sizes <- sizes ** scale
  
  startx <- 0;
  for(l in 1:numBases) {
    
    for(w in 1:4) {
      endx <- startx + A[l,w]
      polygon(c(startx, endx, endx, startx), c(0, 0, sizes[l], sizes[l]), col = baseCol[w], border=F);
      if (endx - startx > 1 / 4 & sizes[l] > 0.5) {
        text(0.5 * (endx + startx), 0.5 * sizes[l], num2base[w], col="white", cex=charSize)
      }
      startx <- endx;
    }
    startx <- startx + 0.25;
  }
  
  startx <- (centerBase - 1) * 1.25;
  for (w in 1:4) {
    starty <- 2;
    endx <- startx + A[centerBase,w];
    for(ww in 1:4) {
      endy <- starty + B[w,ww];
      polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[ww], border=F);
      if ((endy - starty > 1 / 4) & (endx - startx > 1 / 4)) {
        text(0.5 * (endx + startx), 0.5 * (endy + starty), num2base[ww], col="white", cex=charSize)
      }
      starty <- endy;
    }
    startx <- endx
    starty <- endy;
  }
  
  ##########
  # draw arrow
  xs <- c(1 /3, 2 / 3, 2 / 3, 5 / 6, 1 / 2, 1 / 6, 1 / 3, 1 / 3) + (centerBase - 1) * 1.25;
  ys <- c(1 / 4, 1 / 4, 1 / 2, 1 / 2, 3 / 4, 1 / 2, 1 / 2, 1 / 4) + 1;
  
  polygon(xs, ys, col=8, border=F);
  ##########
  
  ##########
  if (trDir == TRUE) {
    # draw direction bias
    # startx <- (numBases - 1) * 1.25 + 0.5;
    # endx <- (numBases - 1) * 1.25 + 0.75;
    # starty <- 1.9;
    # endy <- starty + v3[1];
    # polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[5], border=F);
    
    # if (endy - starty > 1 / 8) {
    #   text(0.5 * (startx + endx), 0.5 * (starty + endy), "+", col="white", cex=1.2)
    # }
    # starty <- endy;
    # endy <- 2.9;
    # polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[6], border=F);
    # if (endy - starty > 1 / 8) {
    #   text(0.5 * (startx + endx), 0.5 * (starty + endy), "-", col="white", cex=1.2)
    # } 
    
    # draw direction bias
    startx <- (numBases - 1) * 1.25 + 0.24;
    endx <- (numBases - 1) * 1.25 + 0.49;
    starty <- 2;
    endy <- starty + v3[1];
    polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[5], border=F);
    if (endy - starty > 1 / 8) {
      text(0.5 * (startx + endx), 0.5 * (starty + endy), "+", col="white", cex=charSize)
    }
    
    startx <- (numBases - 1) * 1.25 + 0.51;
    endx <- (numBases - 1) * 1.25 + 0.76;
    starty <- 2;
    endy <- starty + v3[2];
    polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[6], border=F);
    if (endy - starty > 1 / 8) {
      text(0.5 * (startx + endx), 0.5 * (starty + endy), "-", col="white", cex=charSize)
    }
    
  }
  ##########
  
}



# convert the mutation signature to one dimensional vector
convertFMatrixToVector <- function(Fmat, fdim) {
  
  M <- prod(fdim);
  Fvec <- rep(1, M);
  
  temp1 <- 1;
  temp2 <- 1;
  for (i in 1:length(fdim)) {
    temp1 <- temp1 * fdim[i];
    divInd <- (1:M - 1) %% temp1 + 1;
    for (j in 1:fdim[i]) {
      targetInd <- divInd > temp2 * (j - 1) & divInd <= temp2 * j
      Fvec[targetInd] <- Fvec[targetInd] * Fmat[i,j]
    }
    temp2 <- temp2 * fdim[i];
  }
  
  return(Fvec)
}

# get the mutation signature of the specified cancer type and index
getF <- function(type, K, index, trDir = FALSE) {
  
  if (trDir ==  TRUE) {
    inputName <- paste("../../AlexandrovEtAl/result/Param_ind5_dir/", type, ".", as.character(K), ".Rdata", sep="");
  } else {
    inputName <- paste("../../AlexandrovEtAl/result/Param_ind5/", type, ".", as.character(K), ".Rdata", sep="");    
  }
  
  load(inputName);
  F <- resultForSave[[1]]@signatureFeatureDistribution;
  
  return(F[index,,]);
  
}


# get the bootstrap error 
getBoot <- function(type, K, index, trDir = FALSE) {
  
  if (trDir ==  TRUE) {
    inputName <- paste("../../AlexandrovEtAl/result/Param_ind5_dir/", type, ".", as.character(K), ".Rdata", sep="");
  } else {
    inputName <- paste("../../AlexandrovEtAl/result/Param_ind5/", type, ".", as.character(K), ".Rdata", sep="");    
  }
  
  load(inputName);
  B <- resultForSave[[2]][[1]][,index,,];
  
  return(B);
  
}



# get the mutation signatures of the specified cancer type
getSigMat <- function(type, K, trDir = FALSE) {
  
  if (trDir ==  TRUE) {
    inputName <- paste("../../AlexandrovEtAl/result/Param_ind5_dir/", type, ".", as.character(K), ".Rdata", sep="");
  } else {
    inputName <- paste("../../AlexandrovEtAl/result/Param_ind5/", type, ".", as.character(K), ".Rdata", sep="");    
  }
  
  load(inputName);
  F <- resultForSave[[1]]@signatureFeatureDistribution;
  fdim <- resultForSave[[1]]@possibleFeatures
    
  resFmat <- matrix(0, K - 1, prod(fdim));
  for (k in 1:(K-1)) {
    tempFvec <- convertFMatrixToVector(F[k,,], fdim);
    Fvec <- tempFvec / sqrt(sum( tempFvec^2));
    resFmat[k,] <- Fvec;
  }
  
  rownames(resFmat) <- paste(type, 1:(K-1), sep="_");
  
  return(resFmat);
  
}


innerProd <- function(a, b) {
  na <- a / sqrt(sum(a^2));
  nb <- b / sqrt(sum(b^2));
  
  return(sum(na * nb));
}




##########
# get the order of the merged mutation signature

sortSignature <- function(Fs) {
  
  # change the order of the signature
  remainInd <- 1:length(Fs);
  newOrder <- c();

  for (k in remainInd) {
    if (Fs[[k]][1,1] > 0.6) {
      newOrder <- c(newOrder, k);
      remainInd <- remainInd[-which(remainInd == k)];
    }
  }

  for (k in remainInd) {
    if (Fs[[k]][1,2] > 0.6) {
      newOrder <- c(newOrder, k);
      remainInd <- remainInd[-which(remainInd == k)];
    }
  }

  for (k in remainInd) {
    if (Fs[[k]][1,3] > 0.6) {
      newOrder <- c(newOrder, k);
      remainInd <- remainInd[-which(remainInd == k)];
    }
  }

  for (k in remainInd) {
    if (sum(Fs[[k]][1,1:3]) > 0.75) {
      newOrder <- c(newOrder, k);
      remainInd <- remainInd[-which(remainInd == k)];
    }
  }


  for (k in remainInd) {
    if (sum(Fs[[k]][1,4]) > 0.6) {
      newOrder <- c(newOrder, k);
      remainInd <- remainInd[-which(remainInd == k)];
    }
  }

  for (k in remainInd) {
    if (sum(Fs[[k]][1,5]) > 0.6) {
      newOrder <- c(newOrder, k);
      remainInd <- remainInd[-which(remainInd == k)];
    }
  }

  for (k in remainInd) {
    if (sum(Fs[[k]][1,6]) > 0.6) {
      newOrder <- c(newOrder, k);
      remainInd <- remainInd[-which(remainInd == k)];
    }
  }


  for (k in remainInd) {
    if (sum(Fs[[k]][1,4:6]) > 0.75) {
      newOrder <- c(newOrder, k);
      remainInd <- remainInd[-which(remainInd == k)];
    }
  }

  return(c(newOrder, remainInd))
  
}


