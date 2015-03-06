
library(Matrix);
library(corrplot);
library(pmsignature);

.pardefault <- par(no.readonly = TRUE);

##########
# functions


visPMS_ind <- function(vF, numBases, baseCol = NA, trDir, charSize = 1.2) {
  
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
  
  frame();
  plot.window(xlim=c(-0.25, 1.25 * numBases + 0.25), ylim=c(-0.25, 3.25));
  
  startx <- 0;
  for(l in 1:numBases) {
    
    for(w in 1:4) {
      endx <- startx + A[l,w]
      polygon(c(startx, endx, endx, startx), c(0, 0, 1, 1), col = baseCol[w], border=F);
      if (endx - startx > 1 / 4) {
        text(0.5 * (endx + startx), 0.5, num2base[w], col="white", cex=charSize)
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
getF <- function(type, K, fdim, ind) {
  
  inputName <- paste("../../../AlexandrovEtAl/result/Param_ind5/", type, ".", as.character(K), ".Rdata", sep="");
  
  load(inputName);
  F <- Param@signatureFeatureDistribution;
  
  return(F[ind,,]);
  
}


# get the mutation signatures of the specified cancer type
getSigMat <- function(type, K, fdim) {
  
  inputName <- paste("../../../AlexandrovEtAl/result/Param_ind5/", type, ".", as.character(K), ".Rdata", sep="");
  
  load(inputName);
  F <- Param@signatureFeatureDistribution;
  
  resFmat <- matrix(0, K - 1, prod(fdim));
  for (k in 1:(K-1)) {
    tempFvec <- convertFMatrixToVector(F[k,,], fdim);
    Fvec <- tempFvec / sqrt(sum( tempFvec^2));
    resFmat[k,] <- Fvec;
  }
  
  rownames(resFmat) <- paste(type, 1:(K-1), sep="_");
  
  return(resFmat);
  
}

##########
#' Here, we focus on the case where substitution patterns 
#' and two 5' and 3' bases are considered
fdim <- c(6, 4, 4, 4, 4);
SigMat <- matrix(0, 0, prod(fdim));
# eps <- 2 - 2 * cos(pi / 5); # this value can be tuned.. where is the best value??
eps <- 0.6

#' create the matrix showing the relationships between cancer type 
#' and the specified numer of mutation signatures
type2sigNum <- read.table("../../data/AlexandrovEtAl_sigNum.txt", sep=",", header=FALSE);

#' reading and converting the mutation signatures for each cancer type
for (i in 1:nrow(type2sigNum)) {
  SigMat <- rbind(SigMat, getSigMat(type2sigNum[i,1], type2sigNum[i,2], fdim));
}


innerProd <- function(a, b) {
  na <- a / sqrt(sum(a^2));
  nb <- b / sqrt(sum(b^2));
  
  return(sum(na * nb));
}


SigV <- rownames(SigMat);

#' calculate cosine distances among mutation signatures
#' and list up the pairs having smaller distances than the specified threshould
# closeSig <- Matrix(0, nrow(SigMat), nrow(SigMat));
closeSig <- diag(nrow(SigMat));
dist <- c();
for (i in 1:(nrow(SigMat) - 1)) {
  for (j in (i+1):nrow(SigMat)) {
    d <- sqrt(sum( (SigMat[i,] - SigMat[j,])^2 ));
    # d <- 2 - 2 * innerProd(SigMat[i,], SigMat[j,]);
    # dist <- c(dist, d); # for checking the histgram of cosine distances
    if ( d < eps) {
      closeSig[i,j] <- 1; 
    }
  }
}


#' clustering the similar mutation signatures
EXPCLOSESIG <- expm(closeSig);
removeList <- rep(1, nrow(SigMat));
curInd <- 1;
typeVec <- c();

Fs <- list();
for (i in 1:length(removeList)) {
  
  samples <- c();
  # indList <- which(closeSig[i,] > 0);
  indList <- which(EXPCLOSESIG[i,] > 0 & removeList == 1);
  
  if (removeList[i] != 0 & length(indList) > 0) {
    
    F <- matrix(0, 5, 6);
    
    for (j in 1:length(indList)) {
      
      jj <- indList[j];
      typeName <- as.character(strsplit(SigV[jj], split="_")[[1]][1]);
      typeSigNum <- type2sigNum[type2sigNum == typeName, 2];
      sigInd <- as.integer(strsplit(SigV[jj], split="_")[[1]][2]);
      
      F <- F + getF(typeName, typeSigNum, fdim, sigInd) / length(indList);
      samples <- c(samples, SigV[jj]);
      removeList[jj] <- 0;
      
    }
    
    
    Fs[[curInd]] <- F;
    
    # outputName <- paste("sig_center/signature_", curInd, ".png", sep="");
    # png(outputName, pointsize = 16);
    # plot.new();
    # visPMS_ind5(F[1,1:6], F[2:5,1:4]);
    # dev.off();
    
    typeVec <- c(typeVec, paste0(samples, collapse=",") ); 
    
    curInd <- curInd + 1;
    
    
  }
  
}


#' generate the figure of mutation signatures

par(mar=c(1, 0, 2, 0));
par(mfrow=c(ceiling(length(Fs) / 4), 4));
# myOrder <- c(7, 11, 12, 15, 21,  # C > A
#              1, 6, 10, 13, 18, 20, 23, # C > T
#              2, 8, 14, 22, #C > any
#              9, 16, 19, 24, #T > any
#              3, 4, 5, 17, 25, 26)

# for (i in myOrder) {
for (i in 1:length(Fs)) {
  visPMS_ind(Fs[[i]], numBases = 5, trDir = FALSE, charSize = 1);
  # visPMS_ind5_mod3(Fs[[i]], 1 / 100000);
  mtext(paste("signature", i),
        outer = FALSE,      # 作図領域の外の余白に書く
        side = 3,          # 上の余白に書く
        cex = 1,         # 字の大きさ
        line = 0.2,          # 外に向かって 0.5行離れたところに書く．
        col = "black")        
}


outputName <- "../../result/AlexandrovEtAl_mergedSignature.eps";
dev.copy2eps(file=outputName, height = ceiling(length(Fs) / 4) * 2.7, width = 15, pointsize = 18);
par(.pardefault);


##########
#' comparing with the signatures observed in the Alexandrov et al. Nature 2013
nature2013_sig_raw <- read.table("../../data/AlexandrovEtAl_signatures.txt", header=T, sep="\t");
nature2013_sig <- t(nature2013_sig_raw[,4:30]);
colnames(nature2013_sig) <- nature2013_sig_raw[,3];

eps_diff <- 0.02;
corSigVec <- c();
for (i in 1:length(Fs)) {
  tempDiff <- matrix(convertFMatrixToVector(Fs[[i]][c(4, 3, 1),], c(4, 4, 6)), nrow(nature2013_sig), 96,byrow=T) - nature2013_sig
  tempSigVec <- rownames(nature2013_sig)[rowSums(tempDiff^2) < eps_diff];
  corSigVec <- c(corSigVec, paste0(tempSigVec, collapse=",") );
}
##########

#########
#' generate the table showing the relationships between signatures and cancer types
sig2type <- matrix(0, length(typeVec), length(type2sigNum[,1]));
for (i in 1:length(typeVec)) {
  types <- strsplit(typeVec[i], split=",")[[1]];
  for (j in 1:length(types)) {
    tempType <- substring(types[j], 1, nchar(types[j]) - 2);
    sig2type[i, which(tempType == type2sigNum[,1])] <- 1;
  }
}

rownames(sig2type) <- paste("signature", 1:length(typeVec));
colnames(sig2type) <- type2sigNum[,1];

myCol <- colorRampPalette(c("#F8F8FF", "#F8F8FF", "#F8F8FF", "#6B8E23"));

#' Export Width 800, Height 800
corrplot(sig2type, is.corr=FALSE, bg="#F8F8FF", col=myCol(200), tl.col="black", tl.cex=0.7, tl.srt=90, cl.pos="n");

outputName <- "../../result/corrplot.eps"
dev.copy2eps(file=outputName, height = 12, width = 12, pointsize = 12);
par(.pardefault);

##########

# generate figures for each signature group

#' Export eps, width=800, height=600 (sig2) 

dir.create("../../result/sig_group/");
for (i in 1:length(typeVec)) {
  
  # outputName <- paste("sig_group/signature_", i, ".png", sep="");
  # png(outputName, height = ceiling(length(types) / 3) * 300, width = 1200, pointsize = 32);
  
  types <- strsplit(typeVec[i], split=",")[[1]];
  par(mar=c(1, 0, 2, 0));
  par(mfrow=c(ceiling(length(types) / 4), 4));
  # plot.new();
  for (j in 1:length(types)) {
    
    typeName <- as.character(strsplit(types[j], split="_")[[1]][1]);
    typeSigNum <- type2sigNum[type2sigNum == typeName, 2];
    sigInd <- as.integer(strsplit(types[j], split="_")[[1]][2]);
    F <- getF(typeName, typeSigNum, fdim, sigInd);
    
    visPMS_ind(F, numBases = 5, trDir = FALSE, charSize = 1);
    mtext(types[j],
          outer = FALSE,      # 作図領域の外の余白に書く
          side = 3,          # 上の余白に書く
          cex = 1,         # 字の大きさ
          line = 0.3,          # 外に向かって 0.5行離れたところに書く．
          col = "black")        
    
    # if (j < length(types)) {
    #   frame();
    # }
    
  }
  
  # dev.off();
  outputName <- paste("../../result/sig_group/signature_", i, ".eps", sep="");
  dev.copy2eps(file=outputName, height = ceiling(length(types) / 4) * 2.7, width = 15, pointsize = 18);
  
}
par(.pardefault);


##########
# two 5' bases from the mutated site for APOBEC mutation signature
types <- strsplit(typeVec[2], split=",")[[1]];

vtype <- c();
vbase <- c();
vint <- c();
twoFivePrimeSig <- data.frame();
for (j in 1:length(types)) {
  
  typeName <- as.character(strsplit(types[j], split="_")[[1]][1]);
  typeSigNum <- type2sigNum[type2sigNum == typeName, 2];
  sigInd <- as.integer(strsplit(types[j], split="_")[[1]][2]);
  F <- getF(typeName, typeSigNum, fdim, sigInd);
  
  vtype <- c(vtype, rep(types[j], 4));
  vbase <- c(vbase, c("A", "C", "G", "T"));
  vint <- c(vint, F[2,1:4]);
  
}

twoFivePrimeSig <- data.frame(type =vtype, base = vbase, intensity = vint);


ggplot(twoFivePrimeSig, aes(x=type, y=intensity, fill=base)) +
  geom_bar(stat="identity") +
  # scale_fill_brewer(palette="Set2") + 
  theme_bw() +
  theme(text = element_text(size=20)) +
  theme(axis.text.x= element_text(angle=60,hjust=1));

ggsave("../../result/APOBEC_two5prime.eps", width=8, height=8, units="in");



