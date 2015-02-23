# initialization

library(pmsignature);


inputFile <- commandArgs()[5];
ratio <- as.numeric(commandArgs()[6]);
sigNum <- as.integer(commandArgs()[7]);
outputDir <- commandArgs()[8];
repNum <- as.integer(commandArgs()[9]);
type <- commandArgs()[10];

# reading the inputFile (mutation feature format)

if (type == "ind") {
    G <- readMFVFile(inputFile, numBases = 5, type="independent", trDir=TRUE);
} else {
    G <- readMFVFile(inputFile, numBases = 5, type="full", trDir=TRUE);
}

BG_prob <- readBGFile(G);
Param_GOLD <- getPMSignature(G, K = sigNum, BG = BG_prob);

# check each obtained signature (AA or APOBEC?)
F_true <- Param_GOLD@signatureFeatureDistribution;
if (type == "ind") {
    if (sum(F_true[1,1,1:3]) > 0.75 & F_true[2,1,4] > 0.75) {
        header <- c("APOBEC", "AA");
    } else if (sum(F_true[2,1,1:3]) > 0.75 & F_true[1,1,4] > 0.75) {
        header <- c("AA", "APOBEC");
    } else {
        stop("either AA or APOBEC signal is not detected");
    }
} else {
    if (sum(F_true[1,1,(256 * 0 + 1):(256 * 3)]) + sum(F_true[1,1,(256 * 6 + 1):(256 * 9)]) > 0.6 & sum(F_true[2,1,(256 * 3 + 1):(256 * 4)]) + sum(F_true[2,1,(256 * 9 + 1):(256 * 10)]) > 0.6) {
        header <- c("APOBEC", "AA");
    } else if (sum(F_true[2,1,(256 * 0 + 1):(256 * 3)]) + sum(F_true[2,1,(256 * 6 + 1):(256 * 9)]) > 0.6 & sum(F_true[1,1,(256 * 3 + 1):(256 * 4)]) + sum(F_true[1,1,(256 * 9 + 1):(256 * 10)]) > 0.6) {
        header <- c("AA", "APOBEC");
    } else {
        stop("either AA or APOBEC signal is not detected");
    }
}

 
Dists <- matrix(0, repNum, sigNum - 1);
##########
# downsampling
for (kkk in 1:repNum) {

    G_RLE <- list(lengths = G@countData[3,], values = 1:length(G@countData[3,]));
    class(G_RLE) <- "rle";

    G_RLE$values <- 1:length(G@countData[3,]);
    G_RLE$lengths <- G@countData[3,];

    countSample <- table(sample(inverse.rle(G_RLE), sum(G@countData[3,]) * ratio, replace=FALSE));

    G_down <- G;
    G_down@countData[3,] <- 0;
    G_down@countData[3,as.integer(names(countSample))] <- countSample;
    print(sum(G_down@countData[3,]));

    ##########
    # estimation and bootstraping
    BG_prob <- readBGFile(G_down);
    Param_down <- getPMSignature(G_down, K = sigNum, BG = BG_prob);
    F_down <- Param_down@signatureFeatureDistribution;
    ##########
    Dists[kkk, ] <- getCosDistance(F_true, F_down, Param_down@possibleFeatures);
}


colnames(Dists) <- header;
outputFile <- paste(outputDir, "/", sigNum, "_", ratio, ".txt", sep="");
write.table(Dists[,c("AA", "APOBEC")], file=outputFile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t");



