
####################
# figure 1

# dummy signature estimation for generating class
inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt", package="pmsignature");
G <- readMPFile(inputFile, numBases = 3, type = "full");
Param <- getPMSignature(G, K = 2);

# reading the signature from Alexandrov et al. Nature 2013
nature2013_sig_raw <- read.table("../../data/AlexandrovEtAl_signatures.txt", header=T, sep="\t");
nature2013_sig <- t(nature2013_sig_raw[,4:30]);
colnames(nature2013_sig) <- nature2013_sig_raw[,3];

Param@signatureFeatureDistribution[1,1,] <- nature2013_sig[3,];
Param@signatureFeatureDistribution[2,1,] <- nature2013_sig[11,];

visPMSignature(Param, 1)
ggsave("../../result/example96_APOBEC.eps", width=12, height=3, units="in");

visPMSignature(Param, 2)
ggsave("../../result/example96_POLE.eps", width=12, height=3, units="in");
####################


####################
# figure 2

inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt", package="pmsignature");
G <- readMPFile(inputFile, numBases = 5, trDir = TRUE);
Param <- getPMSignature(G, K = 1);

Param@signatureFeatureDistribution[1,1,1:6] <- c(0.201, 0.080, 0.314, 0.000, 0.405, 0.000);
Param@signatureFeatureDistribution[1,2,1:4] <- c(0.146, 0.257, 0.128, 0.469);
Param@signatureFeatureDistribution[1,3,1:4] <- c(0.183, 0.389, 0.331, 0.097);
Param@signatureFeatureDistribution[1,4,1:4] <- c(0.161, 0.107, 0.013, 0.719);
Param@signatureFeatureDistribution[1,5,1:4] <- c(0.019, 0.265, 0.366, 0.350);
Param@signatureFeatureDistribution[1,6,1:2] <- c(0.375, 0.625);
  
visPMSignature(Param, 1);

outputName <- "../../result/example_signature.eps";
dev.copy2eps(file=outputName, height = 4, width = 8, pointsize = 18);
par(.pardefault);


