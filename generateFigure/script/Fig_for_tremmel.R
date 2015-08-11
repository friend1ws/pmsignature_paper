# Script for generating data for tremmel.


if (!file.exists("../result/Fig_tremmel")) {
  dir.create("../result/Fig_tremmel")
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



for (i in newOrder) {
  write.table(Fs[[i]], file = paste("../result/Fig_tremmel/sig_", i, ".txt", sep = ""), sep="\t", col.names = FALSE, row.names = FALSE)
}
