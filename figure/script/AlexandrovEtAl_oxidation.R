library(pmsignature);

.pardefault <- par(no.readonly = TRUE);

##########
# function for checking the index of the oxidative signature
checkOxidationInd <- function(Param) {
  
  maxK <- slot(Param, "signatureNum") - as.integer(slot(Param, "isBackGround"));
  
  sigInd <- 0;
  for (k in 1:maxK) {
    
    F <- Param@signatureFeatureDistribution;
    if (F[k,1,1] > 0.6 & F[k,2,2] > 0.5 & F[k,3,2] > 0.5) {
      sigInd <- k;
    }
    
  }
  
  return(sigInd);
  
}
##########



#' create the matrix showing the relationships between cancer type 
#' and the specified numer of mutation signatures
type2sigNum <- read.table("../data/AlexandrovEtAl_sigNum.txt", sep=",", header=FALSE);

oxidation_types <- c("Kidney-Clear-Cell", "Lung-Adeno", "Melanoma");

for (i in 1:length(oxidation_types)) {
  
  K <- type2sigNum[type2sigNum[,1] == oxidation_types[i], 2];
  inputMPdata <- paste("../../AlexandrovEtAl/result/MPFormat/", oxidation_types[i], ".mp.txt.gz", sep = "");
  inputRdata <- paste("../../AlexandrovEtAl/result//Param_ind5_dir/", oxidation_types[i], ".", K, ".Rdata", sep="");
        
  G <- readMPFile(inputMPdata, numBases = 5, trDir = TRUE);
  load(inputRdata);
  visMembership(G, resultForSave[[1]], toSample = 100, colourBrewer = "Set2");
  ggsave(paste("../../supp/", oxidation_types[i], "_oxidation_membership.eps", sep=""), height = 4, width = 15, pointsize = 9);

  oxInd <- checkOxidationInd(resultForSave[[1]])
  
  tempMar <- par("mar");
  par(mar = 0.4 * tempMar);
  
  visPMSignature(resultForSave[[1]], oxInd, charSize = 0.7);
  dev.copy2eps(file=paste("../../supp/", oxidation_types[i], "_oxidation_signature.eps", sep=""), height = 2.4, width = 3.75);    

  par(.pardefault);
  
}

