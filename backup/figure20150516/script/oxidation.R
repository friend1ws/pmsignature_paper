library(pmsignature);


#' create the matrix showing the relationships between cancer type 
#' and the specified numer of mutation signatures
type2sigNum <- read.table("../data/AlexandrovEtAl_sigNum.txt", sep=",", header=FALSE);

oxidation_types <- c("Kidney-Clear-Cell", "Lung-Adeno", "Melanoma");

for (i in 1:length(oxidation_types)) {
  
  K <- type2sigNum[type2sigNum[,1] == oxidation_types[i], 2];
  inputMPdata <- paste("../../AlexandrovEtAl/result/MPFormat/", oxidation_types[i], ".mp.txt.gz", sep = "");
  inputRdata <- paste("../../AlexandrovEtAl/result//Param_ind5_dir/", oxidation_types[i], ".", K, ".Rdata", sep="");
        
  G <- readMPFile(inputMPdata, numBases = 5, trDir = TRUE);
  visMembership(G, resultForSave[[1]], toSample = 100, colourBrewer = "Set2");
  ggsave("../../supp/", oxidation_types[i], "_oxidation_membership.eps", height = 4, width = 15, pointsize = 9);

  oxInd <- checkOxidationInd(resultForSave[[1]])
  visPMSignature(resultForSave[[1]], oxInd);
  dev.copy2eps(file="../../supp/", oxidation_types[i], "_oxidation_signature.eps", height = 2.4, width = 3.75);    

}


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

