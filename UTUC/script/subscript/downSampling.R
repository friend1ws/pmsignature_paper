# initialization

library(pmsignature);


inputFile <- commandArgs()[5];
ratio <- as.numeric(commandArgs()[6]);
sigNum <- as.integer(commandArgs()[7]);
outputDir <- commandArgs()[8];
repNum <- commandArgs()[9];
type <- commandArgs()[10];

# reading the inputFile (mutation feature format)

cat(type);
if (type == "ind") {
    G <- readMFVFile(inputFile, numBases = 5, type="independent", trDir=TRUE);
} else {
    G <- readMFVFile(inputFile, numBases = 5, type="full", trDir=TRUE);
}

##########
# downsampling
for (kkk in 1:repNum) {

    outputFile <- paste(outputDir, "/", sigNum, "_", ratio, "_", kkk, ".Rdata", sep="");

    G_RLE <- list(lengths = G@countData[3,], values = 1:length(G@countData[3,]));
    class(G_RLE) <- "rle";

    G_RLE$values <- 1:length(G@countData[3,]);
    G_RLE$lengths <- G@countData[3,];

    countSample <- table(sample(inverse.rle(G_RLE), sum(G@countData[3,]) * ratio, replace=FALSE));

    G_down <- G;
    G_down@countData[3,as.integer(names(countSample))] <- countSample;

    ##########
    # estimation and bootstraping
    BG_prob <- readBGFile(G_down);
    Param <- getPMSignature(G_down, K = sigNum, BG = BG_prob);
    Param_boot <- bootPMSignature(G_down, Param0 = Param, bootNum = 100, BG = BG_prob);
    ##########

    # save the result
    resultForSave <- list(Param, Param_boot);
    save(resultForSave, file= outputFile);
 
}

