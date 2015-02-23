# initialization

library(pmsignature);


inputFile <- commandArgs()[5];
outputFile <- commandArgs()[6];
sigNum <- as.integer(commandArgs()[7]);
type <- commandArgs()[8];

# reading the inputFile (mutation feature format)

cat(type);
if (type == "ind") {
    G <- readMFVFile(inputFile, numBases = 5, type="independent", trDir=TRUE);
} else {
    G <- readMFVFile(inputFile, numBases = 5, type="full", trDir=TRUE);
}

##########
# estimation and bootstraping
BG_prob <- readBGFile(G);
Param <- getPMSignature(G, K = sigNum, BG = BG_prob);
Param_boot <- bootPMSignature(G, Param0 = Param, bootNum = 100, BG = BG_prob);
##########

# save the result
resultForSave <- list(Param, Param_boot);
save(resultForSave, file= outputFile);
 

