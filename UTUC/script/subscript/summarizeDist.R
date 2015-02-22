source("../lib/myRScripts/cosineDist.R");
inputPrefix <- commandArgs()[5];
trueFile <- commandArgs()[6];
outputFile <- commandArgs()[7];
type <- commandArgs()[8];

if (type == "ind3_dir") {
    fdim <- c(6, 4, 4, 2);
} else if (type == "ind5_dir") {
    fdim <- c(6, 4, 4, 4, 4, 2);
} else if (type == "full3_dir") {
    fdim <- c(192);
} else if (type == "full5_dir") {
    fdim <- c(3072);
}

load(trueFile);
F_true <- lret[[1]][[1]];


Dists <- matrix(0, 100, 2);
for (i in 1:100) {
    load(paste(inputPrefix, "_", i, ".Rdata", sep=""));
    F_est <- lret[[1]][[1]];
    Dists[i,] <- getCosDistance(F_true, F_est, fdim);
}

write.table(Dists, file=outputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t");

