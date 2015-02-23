library(pmsignature);

load("../../../UTUC/result/Param_full/3.Rdata");
Param <- resultForSave[[1]];

if (sum(Param@signatureFeatureDistribution[1,1,(256 * 0 + 1):(256 * 3)]) + sum(Param@signatureFeatureDistribution[1,1,(256 * 6 + 1):(256 * 9)]) > 0.6 & sum(Param@signatureFeatureDistribution[2,1,(256 * 3 + 1):(256 * 4)]) + sum(Param@signatureFeatureDistribution[2,1,(256 * 9 + 1):(256 * 10)]) > 0.6) {
  visPMSignature(Param, 1);
  ggsave("../../result/UTUC_APOBEC_full.eps", height = 4, width = 6);
  visPMSignature(Param, 2);
  ggsave("../../result/UTUC_AA_full.eps", height = 4, width = 6); 
} else {
  visPMSignature(Param, 2);
  ggsave("../../result/UTUC_APOBEC_full.eps", height = 4, width = 6);
  visPMSignature(Param, 1);
  ggsave("../../result/UTUC_AA_full.eps", height = 4, width = 6);   
}
  
load("../../../UTUC/result/Param_ind/3.Rdata");
Param <- resultForSave[[1]];

if (sum(Param@signatureFeatureDistribution[1,1,1:3]) > 0.75 & Param@signatureFeatureDistribution[2,1,4] > 0.75) {
  visPMSignature(Param, 1, charSize = 0.8);
  dev.copy2eps(file="../../result/UTUC_APOBEC_ind.eps", height = 4, width = 6);
  visPMSignature(Param, 2, charSize = 0.8);
  dev.copy2eps(file="../../result/UTUC_AA_ind.eps", height = 4, width = 6);
} else {
  visPMSignature(Param, 2, charSize = 0.8);
  dev.copy2eps(file="../../result/UTUC_APOBEC_ind.eps", height = 4, width = 6);
  visPMSignature(Param, 1, charSize = 0.8);
  dev.copy2eps(file="../../result/UTUC_AA_ind.eps", height = 4, width = 6);
} 
