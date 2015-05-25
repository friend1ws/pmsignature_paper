#! /home/yshira/local/bin/R

library(pmsignature);

inputFile <- commandArgs()[5];
outputFile <- commandArgs()[6];
sigNum <- as.numeric(commandArgs()[7]);
trDirFlag <- as.logical(commandArgs()[8]);
trialNum <- as.numeric(commandArgs()[9]);

G <- readMPFile(inputFile, numBases = 5, trDir = trDirFlag);

BG_prob <- readBGFile(G);
Param <- getPMSignature(G, K = sigNum , BG = BG_prob, numInit = trialNum);

Boot <- bootPMSignature(G, Param0 = Param, bootNum = 100, BG = BG_prob);

save(list(Param, Boot), file=outputFile); 

