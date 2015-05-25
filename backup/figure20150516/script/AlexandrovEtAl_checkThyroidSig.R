
type2sigNum <- read.table("../data/AlexandrovEtAl_sigNum.txt", sep=",", header=FALSE);

matchMutNum <- matrix(0, length(type2sigNum[,1]), 2);
rownames(matchMutNum) <- type2sigNum[,1];

for (i in 1:length(type2sigNum[,1])) {
  
  type <-  type2sigNum[i, 1];
  print(type);
  inputFile <- paste("../../AlexandrovEtAl/result/MPFormat/", type, ".mp.txt.gz", sep="");
  G <- readMPFile(inputFile, numBases = 5, trDir = TRUE);

  FVList <- G@featureVectorList;

  mutInds_plus <- which(FVList[1,] == 5 & FVList[2, ] == 2 & FVList[3,] %in% c(1, 2) & FVList[6,] == 1);
  mutInds_minus <- which(FVList[1,] == 5 & FVList[2, ] == 2 & FVList[3,] %in% c(1, 2) & FVList[6,] == 2);

  countData <- G@countData;
  countData[3, countData[1, ] %in% mutInds_plus];
  countData[3, countData[1, ] %in% mutInds_minus];

  matchMutNum[i, 1] <- sum( countData[3, countData[1, ] %in% mutInds_plus] );
  matchMutNum[i, 2] <- sum( countData[3, countData[1, ] %in% mutInds_minus] );
  
}



