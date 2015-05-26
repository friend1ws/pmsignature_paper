library(pmsignature);

##########
# function for checking the index of the POLE signatures
checkPoleInd <- function(Param) {
  
  maxK <- slot(Param, "signatureNum") - as.integer(slot(Param, "isBackGround"));
  
  sigInd1 <- 0;
  sigInd2 <- 0;
  for (k in 1:maxK) {
    
    F <- Param@signatureFeatureDistribution;

    if (F[k,1,1] > 0.6 & F[k,3,4] > 0.6 & F[k,4,4] > 0.6) {
      sigInd1 <- k;
    }
    
    if (F[k,1,3] > 0.6 & F[k,3,4] > 0.6 & F[k,4,3] > 0.6) {
      sigInd2 <- k;
    }
    
  }
  
  return(c(sigInd1, sigInd2));
  
}
###################



type2sigNum <- read.table("../data/AlexandrovEtAl_sigNum.txt", sep=",", header=FALSE);
POLE_sample <- read.table("../data/POLE_sample.txt", header=FALSE, stringsAsFactors = FALSE)[[1]];

POLE_types <- c("Colorectum", "Uterus");

mem_pole <- data.frame(membership_sig1 = c(), membership_sig8 = c(), type = c())

##########
# cololrectum
K <- type2sigNum[type2sigNum[,1] == "Colorectum", 2];
inputRdata <- paste("../../AlexandrovEtAl/result//Param_ind5_dir/", "Colorectum", ".", K, ".Rdata", sep="");
load(inputRdata);
poleInd_colo <- checkPoleInd(resultForSave[[1]]);
samples_colo <- unlist(lapply(strsplit(resultForSave[[1]]@sampleList, split="-"), function(x) {paste(x[1], x[2], x[3], sep="-")}));
mem_colo <- resultForSave[[1]]@sampleSignatureDistribution[samples_colo %in% POLE_sample ,poleInd_colo];
mem_pole <- cbind(mem_pole, mem_colo[,1], mem_colo[,2], rep("colorectum", nrow(mem_colo)));


##########
# uterus
K <- type2sigNum[type2sigNum[,1] == "Uterus", 2];

inputMPdata <- paste("../../AlexandrovEtAl/result/MPFormat/", "Uterus", ".mp.txt.gz", sep = "");
G <- readMPFile(inputMPdata, numBases = 5, trDir = TRUE);

inputRdata <- paste("../../AlexandrovEtAl/result//Param_ind5_dir/", "Uterus", ".", K, ".Rdata", sep="");
load(inputRdata);
poleInd_ute <- checkPoleInd(resultForSave[[1]]);
samples_ute <- unlist(lapply(strsplit(resultForSave[[1]]@sampleList, split="-"), function(x) {paste(x[1], x[2], x[3], sep="-")}));
mem_ute <- resultForSave[[1]]@sampleSignatureDistribution[samples_ute %in% POLE_sample ,poleInd_ute];



mem_pole <- data.frame(membership_sig1 = c(mem_colo[, 1], mem_ute[, 1]), 
                       membership_sig8 = c(mem_colo[, 2], mem_ute[, 2]),
                       type = c(rep("colorectum", nrow(mem_colo)), rep("uterus", nrow(mem_ute))))


gplot(mem_pole, aes(x = membership_sig1, y = membership_sig8, colour = type)) + geom_point() + xlim(0, 1) + ylim(0, 1);





