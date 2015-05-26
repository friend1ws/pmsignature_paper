# Script for generating figures in SText 1 (simulation study)
# First, we need the results ("../../analysis/UCUT/result"),

if (!file.exists("../result/Fig2")) {
  dir.create("../result/Fig2")
}

if (!file.exists("../result/Fig4")) {
  dir.create("../result/Fig4")
}

if (!file.exists("../result/Fig5")) {
  dir.create("../result/Fig5")
}

if (!file.exists("../result/Fig6")) {
  dir.create("../result/Fig6")
}

if (!file.exists("../result/SFig2")) {
  dir.create("../result/SFig2")
}


##########
library(Matrix);
library(corrplot);
library(pmsignature);

.pardefault <- par(no.readonly = TRUE);

# functions
source("utils.R")


#' Here, we focus on the case where substitution patterns 
#' and two 5' and 3' bases are considered
fdim <- c(6, 4, 4, 4, 4, 2);
SigMat <- matrix(0, 0, prod(fdim));
# eps <- 2 - 2 * cos(pi / 5); # this value can be tuned.. where is the best value??
eps <- 0.6

#' create the matrix showing the relationships between cancer type 
#' and the specified numer of mutation signatures
type2sigNum <- read.table("../data/AlexandrovEtAl_sigNum.txt", sep=",", header=FALSE);

#' reading and converting the mutation signatures for each cancer type
for (i in 1:nrow(type2sigNum)) {
  SigMat <- rbind(SigMat, getSigMat(type2sigNum[i,1], type2sigNum[i,2], trDir = TRUE));
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
    
    F <- matrix(0, length(fdim), max(fdim));
    
    for (j in 1:length(indList)) {
      
      jj <- indList[j];
      typeName <- as.character(strsplit(SigV[jj], split="_")[[1]][1]);
      typeSigNum <- type2sigNum[type2sigNum == typeName, 2];
      sigInd <- as.integer(strsplit(SigV[jj], split="_")[[1]][2]);
      
      F <- F + getF(typeName, typeSigNum, sigInd, trDir = TRUE) / length(indList);
      samples <- c(samples, SigV[jj]);
      removeList[jj] <- 0;
      
    }
    Fs[[curInd]] <- F;
  
    typeVec <- c(typeVec, paste0(samples, collapse=",") ); 
    curInd <- curInd + 1;
    
  }
  
}


# obtain new order
newOrder <- sortSignature(Fs)

##########
#' generate the figure of mutation signatures

par(mar=c(0, 0, 0, 0));
par(bg = rgb(0.9, 0.9, 0.9));
par(xaxs = "i", yaxs = "i");
par(mfrow=c(ceiling(length(Fs) / 4), 4));

curSig <- 1;
for (i in newOrder) {
  # for (i in 1:length(Fs)) {
  plot.new();
  plot.window(xlim=c(-0.3, 6.3), ylim=c(-0.3, 3.8));
  polygon(c(-0.25, 6.25, 6.25, -0.25), c(-0.25, -0.25, 3.75, 3.75), col = "white", border = FALSE);
  visPMS_ind(Fs[[i]], numBases = 5, trDir = TRUE, charSize = 1);
  # visPMS_ind5_mod3(Fs[[i]], 1 / 100000);
  mtext(paste("signature", curSig),
        outer = FALSE,      # draw outside the figure regions 
        side = 3,          # draw on upper margin 
        cex = 1,         # character size
        line = -1.5,          # draw at 0.5 line toward outside 
        col = "black")    
  curSig <- curSig + 1;
  
}

i <- max(newOrder);
edgeFlag <- 1;
if (i < ceiling(length(Fs) / 4) * 4) {
  i <- i + 1;
  for (j in i:(ceiling(length(Fs) / 4) * 4)) {
    plot.new();
    plot.window(xlim=c(-0.3, 6.3), ylim=c(-0.3, 3.8));
    if (edgeFlag == 1) {
      polygon(c(-0.25, 6.3, 6.3, -0.25), c(-0.3, -0.3, 3.75, 3.75), col = "white", border = FALSE);
      edgeFlag <- 0;
    } else {
      polygon(c(-0.3, 6.3, 6.3, -0.3), c(-0.3, -0.3, 3.75, 3.75), col = "white", border = FALSE);      
    }
  }
}


outputName <- "../result/Fig4/AlexandrovEtAl_mergedSignature.eps";
dev.copy2eps(file=outputName, height = ceiling(length(Fs) / 4) * 2.5, width = 16, pointsize = 18);
par(.pardefault);



##########
#' comparing with the signatures observed in the Alexandrov et al. Nature 2013
nature2013_sig_raw <- read.table("../data/AlexandrovEtAl_signatures.txt", header=T, sep="\t");
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

sig2type <- matrix(0, 0, length(type2sigNum[,1]));
sigNames <- c();
curSig <- 1;

for (i in newOrder) {
  # for (i in 1:length(typeVec)) {
  types <- strsplit(typeVec[i], split=",")[[1]];
  tempSig2type <- rep(0, length(type2sigNum[,1]));
  for (j in 1:length(types)) {
    tempType <- substring(types[j], 1, nchar(types[j]) - 2);
    # sig2type[i, which(tempType == type2sigNum[,1])] <- 1;
    tempSig2type[which(tempType == type2sigNum[,1])] <- 1;
  }
  sig2type <- rbind(sig2type, tempSig2type);
  sigNames <- c(sigNames, paste("signature", curSig));
  curSig <- curSig + 1;
}

# rownames(sig2type) <- paste("signature", 1:length(typeVec));
rownames(sig2type) <- sigNames;
colnames(sig2type) <- type2sigNum[,1];

myCol <- colorRampPalette(c("#F8F8FF", "#F8F8FF", "#F8F8FF", "#6B8E23"));

#' Export Width 800, Height 800
corrplot(sig2type, is.corr=FALSE, bg="#F8F8FF", col=myCol(200), tl.col="black", tl.cex=1.0, tl.srt=90, cl.pos="n");

outputName <- "../result/Fig5/corrplot.eps"
dev.copy2eps(file=outputName, height = 15, width = 15, pointsize = 12);
par(.pardefault);

##########

##########


# get the index for important signature
for (k in 1:length(Fs)) {
  # APOBEC
  if (Fs[[k]][1,2] > 0.35 & Fs[[k]][1,3] > 0.4 & Fs[[k]][3,4] > 0.85) {
    APOBEC_ind <- k;
  }
  
  if (Fs[[k]][1,1] > 0.6 & Fs[[k]][3,4] > 0.6 & Fs[[k]][4,4] > 0.6) {
    POLE1_ind <- k;
  }
  
  if (Fs[[k]][1,3] > 0.7 & Fs[[k]][3,4] > 0.7 & Fs[[k]][4,3] > 0.6) {
    POLE2_ind <- k;
  }
  
  if (Fs[[k]][1,1] > 0.4 & min(Fs[[k]][2,1:4]) > 0.15  & min(Fs[[k]][3,1:4]) > 0.15 & min(Fs[[k]][4,1:4]) > 0.15 & Fs[[k]][6,2] > 0.6) {
    LUNG1_ind <- k;
  }
 
  if (Fs[[k]][1,1] > 0.4 & max(Fs[[k]][2,1:4]) < 0.4  & max(Fs[[k]][3,1:4]) < 0.4 & max(Fs[[k]][4,1:4]) > 0.45 & Fs[[k]][6,2] > 0.6) {
    LUNG2_ind <- k;
  }
  
  if (Fs[[k]][1,3] > 0.8 & Fs[[k]][3,2] > 0.3 & Fs[[k]][3,4] > 0.35 & Fs[[k]][4,2] > 0.3) {
    UV_ind <- k;
  }
  
  if (Fs[[k]][1,3] > 0.8 & Fs[[k]][3,3] > 0.7 & Fs[[k]][4,3] > 0.5) {
    LMST_ind <- k;
  }
}


LUNG_ind <- paste(LUNG1_ind, LUNG2_ind, sep = ",")

importantSigs <- c(APOBEC_ind, POLE1_ind, POLE2_ind, LUNG_ind, UV_ind, LMST_ind);
names(importantSigs) <- c("APOBEC", "POLE1", "POLE2", "LUNG", "UV", "LMST");

# generate figures for each signature group

#' Export eps, width=800, height=600 (sig2) 

for (i in 1:length(importantSigs)) {
  
  # outputName <- paste("sig_group/signature_", i, ".png", sep="");
  # png(outputName, height = ceiling(length(types) / 3) * 300, width = 1200, pointsize = 32);
  
  par(mar=c(0, 0, 0, 0));
  par(bg = rgb(0.9, 0.9, 0.9));
  par(xaxs = "i", yaxs = "i");
  
  types <- c();
  inds <- as.numeric(strsplit(importantSigs[i], split = ",")[[1]]);
  for (iii in inds) {
    types <- c(types, strsplit(typeVec[iii], split=",")[[1]]);
  }
  
  if (length(types) > 4) {
    par(mfrow=c(ceiling(length(types) / 4), 4));
  } else {
    par(mfrow = c(1, length(types)));
  }
  
  
  # plot.new();
  for (j in 1:length(types)) {
    
    typeName <- as.character(strsplit(types[j], split="_")[[1]][1]);
    typeSigNum <- type2sigNum[type2sigNum == typeName, 2];
    sigInd <- as.integer(strsplit(types[j], split="_")[[1]][2]);
    F <- getF(typeName, typeSigNum, sigInd, trDir = TRUE);
    
    plot.new();
    plot.window(xlim=c(-0.3, 6.3), ylim=c(-0.3, 3.8));
    polygon(c(-0.25, 6.25, 6.25, -0.25), c(-0.25, -0.25, 3.75, 3.75), col = "white", border = FALSE);
    
    # the setting of charSize is temporary... why the size changes by the number of rows.....
    visPMS_ind(F, numBases = 5, trDir = TRUE, charSize = 0.35 + 0.05 * ceiling(length(types) / 4) + 0.1666 * min(length(types), 4));
    mtext(types[j],
          outer = FALSE,
          side = 3,
          cex = 1,
          line = -1.5,
          col = "black")        
    
  }
  
  j <- length(types);
  edgeFlag <- 1;
  if (length(types) > 4 & j < ceiling(length(types) / 4) * 4) {
    j <- j + 1;
    for (jj in j:(ceiling(length(types) / 4) * 4)) {
      plot.new();
      plot.window(xlim=c(-0.3, 6.3), ylim=c(-0.3, 3.8));
      if (edgeFlag == 1) {
        polygon(c(-0.25, 6.3, 6.3, -0.25), c(-0.3, -0.3, 3.75, 3.75), col = "white", border = FALSE);
        edgeFlag <- 0;
      } else {
        polygon(c(-0.3, 6.3, 6.3, -0.3), c(-0.3, -0.3, 3.75, 3.75), col = "white", border = FALSE);      
      }
    }
  }
  
  
  
  outputName <- paste("../result/SFig2/signatureList_", names(importantSigs)[i], ".eps", sep="");
  
  if (length(types) > 4) {
    dev.copy2eps(file=outputName, height = ceiling(length(types) / 4) * 2.4, width = 15);
  } else {
    dev.copy2eps(file=outputName, height = 2.4, width = 3.75 * length(types));    
  }
  
  
}
par(.pardefault);


###
# for method overview figure

for (k in 1:length(Fs)) {
  if (Fs[[k]][1,5] > 0.6 & Fs[[k]][3,1] > 0.25 & Fs[[k]][6,2] > 0.7) {
    TCR_ind <- k;
  }
}
  
  
importantSigs2 <- c(LUNG1_ind, UV_ind, TCR_ind);
names(importantSigs2) <- c("LUNG1", "UV", "TCR");

par(mar=c(0, 0, 0, 0));
par(xaxs = "i", yaxs = "i");

for (i in 1:length(importantSigs2)) {
  

  F <- Fs[[as.numeric(importantSigs2[i])]];
    
  plot.new();
  plot.window(xlim=c(-0.3, 6.3), ylim=c(-0.3, 3.8));
  
  # the setting of charSize is temporary... why the size changes by the number of rows.....
  visPMS_ind(F, numBases = 5, trDir = TRUE, charSize = 0);
  outputName <- paste("../result/Fig2/signatureList_", names(importantSigs2)[i], "_nonBase.eps", sep="");
  dev.copy2eps(file=outputName, height = 2.4, width = 3.75, pointsize = 1e-10);
  
  plot.new();
  plot.window(xlim=c(-0.3, 6.3), ylim=c(-0.3, 3.8));

  # the setting of charSize is temporary... why the size changes by the number of rows.....
  visPMS_ind(F, numBases = 5, trDir = TRUE, charSize = 0.6 + 0.1333);
  outputName <- paste("../result/Fig2/signatureList_", names(importantSigs2)[i], ".eps", sep="");
  dev.copy2eps(file=outputName, height = 2.4, width = 3.75, pointsize = 1e-10);
}

par(.pardefault);



###

##########

##########

# two 5' bases from the mutated site for several signature

for (i in 1:length(importantSigs)) {
  
  types <- c();
  inds <- as.numeric(strsplit(importantSigs[i], split = ",")[[1]]);
  for (iii in inds) {
    types <- c(types, strsplit(typeVec[iii], split=",")[[1]]);
  }

  vtype <- c();
  vbase <- c();
  vint <- c();
  vsd <- c();
  twoFivePrimeSig <- data.frame();
  for (j in 1:length(types)) {
  
    typeName <- as.character(strsplit(types[j], split="_")[[1]][1]);
    typeSigNum <- type2sigNum[type2sigNum == typeName, 2];
    sigInd <- as.integer(strsplit(types[j], split="_")[[1]][2]);
    F <- getF(typeName, typeSigNum, sigInd, trDir = TRUE);
    Boot <- getBoot(typeName, typeSigNum, sigInd, trDir = TRUE);
  
    vtype <- c(vtype, rep(types[j], 4));
    vbase <- c(vbase, c("A", "C", "G", "T"));
    vint <- c(vint, F[2,1:4]);
    vsd <- c(vsd, sqrt(apply(Boot, MARGIN = c(2, 3), mean)[2,1:4]))
  
  }

  twoFivePrimeSig <- data.frame(type =vtype, base = vbase, intensity = vint, sd = vsd);


  ggplot(twoFivePrimeSig, aes(x=type, y=intensity, fill=base)) +
    geom_bar(stat="identity", position = position_dodge()) +
    geom_errorbar(
      aes(ymin = intensity - sd, ymax = intensity + sd), 
      position = position_dodge(0.9), width = 0.15) + 
    # scale_fill_brewer(palette="Set2") + 
    theme_bw() +
    theme(text = element_text(size=20)) +
    theme(axis.text.x= element_text(angle=45,hjust=1));

  outputName <- paste("../result/Fig6/", names(importantSigs)[i], "_two5prime.eps", sep="");
  
  ggsave(outputName, width=length(types) * 1.0 + 2.5, height = 6, units="in");

  
}


# ##########
# # don't change the scale

# #' generate the figure of mutation signatures

# par(mar=c(0, 0, 0, 0));
# par(bg = rgb(0.9, 0.9, 0.9));
# par(xaxs = "i", yaxs = "i");
# par(mfrow=c(ceiling(length(Fs) / 4), 4));

# curSig <- 1;
# for (i in newOrder) {
#   # for (i in 1:length(Fs)) {
#   plot.new();
#   plot.window(xlim=c(-0.3, 6.3), ylim=c(-0.3, 3.8));
#   polygon(c(-0.25, 6.25, 6.25, -0.25), c(-0.25, -0.25, 3.75, 3.75), col = "white", border = FALSE);
#   visPMS_ind(Fs[[i]], numBases = 5, trDir = TRUE, charSize = 1, scale = FALSE);
#   # visPMS_ind5_mod3(Fs[[i]], 1 / 100000);
#   mtext(paste("signature", curSig),
#         outer = FALSE,
#         side = 3,
#         cex = 1,
#         line = -1.5,
#         col = "black")    
#   curSig <- curSig + 1;
#   
# }

# i <- max(newOrder);
# edgeFlag <- 1;
# if (i < ceiling(length(Fs) / 4) * 4) {
#   i <- i + 1;
#   for (j in i:(ceiling(length(Fs) / 4) * 4)) {
#     plot.new();
#     plot.window(xlim=c(-0.3, 6.3), ylim=c(-0.3, 3.8));
#     if (edgeFlag == 1) {
#       polygon(c(-0.25, 6.3, 6.3, -0.25), c(-0.3, -0.3, 3.75, 3.75), col = "white", border = FALSE);
#       edgeFlag <- 0;
#     } else {
#       polygon(c(-0.3, 6.3, 6.3, -0.3), c(-0.3, -0.3, 3.75, 3.75), col = "white", border = FALSE);      
#     }
#   }
# }


# outputName <- "../../supp/AlexandrovEtAl_mergedSignature_unchanged.eps";
# dev.copy2eps(file=outputName, height = ceiling(length(Fs) / 4) * 2.5, width = 16, pointsize = 18);
# par(.pardefault);



