
library(Matrix);
library(corrplot);
library(pmsignature);

.pardefault <- par(no.readonly = TRUE);

# functions
source("AlexandrovEtAl_function.R")


##########
#' Here, we focus on the case where substitution patterns 
#' and two 5' and 3' bases are considered
fdim <- c(6, 4, 4, 4, 4);
SigMat <- matrix(0, 0, prod(fdim));
# eps <- 2 - 2 * cos(pi / 5); # this value can be tuned.. where is the best value??
eps <- 0.6

#' create the matrix showing the relationships between cancer type 
#' and the specified numer of mutation signatures
type2sigNum <- read.table("../data/AlexandrovEtAl_sigNum.txt", sep=",", header=FALSE);

#' reading and converting the mutation signatures for each cancer type
for (i in 1:nrow(type2sigNum)) {
  SigMat <- rbind(SigMat, getSigMat(type2sigNum[i,1], type2sigNum[i,2], fdim));
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


####################
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

newOrder <- c(newOrder, remainInd)
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
  visPMS_ind(Fs[[i]], numBases = 5, trDir = FALSE, charSize = 1);
  # visPMS_ind5_mod3(Fs[[i]], 1 / 100000);
  mtext(paste("signature", curSig),
        outer = FALSE,      # 作図領域の外の余白に書く
        side = 3,          # 上の余白に書く
        cex = 1,         # 字の大きさ
        line = -1.5,          # 外に向かって 0.5行離れたところに書く．
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


outputName <- "../../manuscript/AlexandrovEtAl_mergedSignature.eps";
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
corrplot(sig2type, is.corr=FALSE, bg="#F8F8FF", col=myCol(200), tl.col="black", tl.cex=1.2, tl.srt=90, cl.pos="n");

outputName <- "../../manuscript/corrplot.eps"
dev.copy2eps(file=outputName, height = 15, width = 15, pointsize = 18);
par(.pardefault);

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

  if (Fs[[k]][1,3] > 0.8 & Fs[[k]][3,2] > 0.3 & Fs[[k]][3,4] > 0.35 & Fs[[k]][4,2] > 0.3) {
    UV_ind <- k;
  }
  
  if (Fs[[k]][1,3] > 0.8 & Fs[[k]][3,3] > 0.7 & Fs[[k]][4,3] > 0.5) {
    LMST_ind <- k;
  }
}


importantSigs <- c(APOBEC_ind, POLE1_ind, POLE2_ind, UV_ind, LMST_ind);
names(importantSigs) <- c("APOBEC", "POLE1", "POLE2", "UV", "LMST");

# generate figures for each signature group

#' Export eps, width=800, height=600 (sig2) 

for (i in 1:length(importantSigs)) {
  
  # outputName <- paste("sig_group/signature_", i, ".png", sep="");
  # png(outputName, height = ceiling(length(types) / 3) * 300, width = 1200, pointsize = 32);
  
  par(mar=c(0, 0, 0, 0));
  par(bg = rgb(0.9, 0.9, 0.9));
  par(xaxs = "i", yaxs = "i");
  types <- strsplit(typeVec[importantSigs[i]], split=",")[[1]];
  
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
    F <- getF(typeName, typeSigNum, fdim, sigInd);

    plot.new();
    plot.window(xlim=c(-0.3, 6.3), ylim=c(-0.3, 3.8));
    polygon(c(-0.25, 6.25, 6.25, -0.25), c(-0.25, -0.25, 3.75, 3.75), col = "white", border = FALSE);
    
    # the setting of charSize is temporary... why the size changes by the number of rows.....
    visPMS_ind(F, numBases = 5, trDir = FALSE, charSize = 0.6 + 0.1333 * ceiling(length(types) / 4));
    mtext(types[j],
          outer = FALSE,      # 作図領域の外の余白に書く
          side = 3,          # 上の余白に書く
          cex = 1,         # 字の大きさ
          line = -1.5,          # 外に向かって 0.5行離れたところに書く．
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
  
  
  
  outputName <- paste("../../manuscript/signatureList_", names(importantSigs)[i], ".eps", sep="");
  
  if (length(types) > 4) {
    dev.copy2eps(file=outputName, height = ceiling(length(types) / 4) * 2.4, width = 15);
  } else {
    dev.copy2eps(file=outputName, height = 2.4, width = 3.75 * length(types));    
  }
  
  
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

ggsave("../../manuscript/APOBEC_two5prime.eps", width=8, height=5, units="in");





##########
# additional script....

# generate figures for each signature group

#' Export eps, width=800, height=600 (sig2) 

dir.create("../../supp/sig_group/");
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





####################
# change the scale
# This might be not necessary for the paper
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
  visPMS_ind(Fs[[i]], numBases = 5, trDir = FALSE, charSize = 1, scale = 0.333);
  # visPMS_ind5_mod3(Fs[[i]], 1 / 100000);
  mtext(paste("signature", curSig),
        outer = FALSE,      # 作図領域の外の余白に書く
        side = 3,          # 上の余白に書く
        cex = 1,         # 字の大きさ
        line = -1.5,          # 外に向かって 0.5行離れたところに書く．
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


outputName <- "AlexandrovEtAl_mergedSignature_scale0.333.eps";
dev.copy2eps(file=outputName, height = ceiling(length(Fs) / 4) * 2.5, width = 16, pointsize = 18);
par(.pardefault);



