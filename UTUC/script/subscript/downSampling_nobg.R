# initialization

source("PMSfunction.R");

adata <- read.table(commandArgs()[5], sep="\t", header=FALSE, fill=TRUE);
ratio <- as.numeric(commandArgs()[6]);


N <- adata[1,1];
fdim <- as.integer(adata[2,1:adata[1,2]]);
M <- prod(fdim);
L <- length(fdim);


sampleIDs <- as.integer(adata[3:nrow(adata),1]);
mutFeatures <- adata[3:nrow(adata),2:(length(fdim) + 1)];


G <- matrix(0, N, M);
if (length(fdim) > 1) {
    mutIDs <- rowSums(matrix(c(1, cumprod(fdim)[1:(L-1)]), length(sampleIDs), L, byrow=TRUE) * (mutFeatures - 1)) + 1;
} else {
    mutIDs <- rowSums(matrix(1, length(sampleIDs), L, byrow=TRUE) * (mutFeatures - 1)) + 1;
}

for (n in 1:N) {
    sMutIDs <- mutIDs[sampleIDs == n];
    spMutIDs <- sMutIDs[order(runif(length(sMutIDs), 0, 1))[1:ceiling(length(sMutIDs) * ratio)]];
    for (s in 1:length(spMutIDs)) {
        G[n, spMutIDs[s]] = G[n, spMutIDs[s]] + 1;
    }
}

# for (i in 1:length(mutIDs)) {
#     G[sampleIDs[i], mutIDs[i]] = G[sampleIDs[i], mutIDs[i]] + 1;
# }

# for (i in 1:nrow(adata)) {
#   mutID <- adata[i,2] + (adata[i,3] - 1) * 6 + (adata[i,4] - 1) * 24 + (adata[i,5] - 1) * 96 + (adata[i,6] - 1) * 384;
#   G[adata[i,1], mutID] <- G[adata[i,1], mutID] + 1;
# }



bdata <- read.table(commandArgs()[7], sep="\t");
BG0 <- bdata[,2];



K <- as.integer(commandArgs()[8]);

ret <- getPMSignature(G, K, fdim, FALSE, BG0, 10);
F0 <- ret[[1]];
Q0 <- ret[[2]];


# visPMS_ind5(F0[1,1,], F0[1,2:5,1:4]);


bret <- bootPMSignature(G, K, fdim, FALSE, BG0, F0, Q0, 100);

lret <- list(ret, bret);
save(lret, file=commandArgs()[9]); 
