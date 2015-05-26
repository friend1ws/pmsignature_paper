library(pmsignature);

source("utils.R");

load("../../analysis/UCUT/result/Param_full/3.Rdata");
Param_full <- resultForSave[[1]];

if (sum(Param_full@signatureFeatureDistribution[1,1,(256 * 0 + 1):(256 * 3)]) + sum(Param_full@signatureFeatureDistribution[1,1,(256 * 6 + 1):(256 * 9)]) > 0.6 & sum(Param_full@signatureFeatureDistribution[2,1,(256 * 3 + 1):(256 * 4)]) + sum(Param_full@signatureFeatureDistribution[2,1,(256 * 9 + 1):(256 * 10)]) > 0.6) {
  visPMSignature(Param_full, 1);
  ggsave("../result/Fig3/UCUT_APOBEC_full.eps", height = 5, width = 7.5);
  visPMSignature(Param_full, 2);
  ggsave("../result/Fig3/UCUT_AA_full.eps", height = 5, width = 7.5); 
} else {
  visPMSignature(Param_full, 2);
  ggsave("../result/Fig3/UCUT_APOBEC_full.eps", height = 5, width = 7.5);
  visPMSignature(Param_full, 1);
  ggsave("../result/Fig3/UCUT_AA_full.eps", height = 5, width = 7.5);   
}
  
load("../../analysis/UCUT/result/Param_ind/3.Rdata");
Param_ind <- resultForSave[[1]];

if (sum(Param_ind@signatureFeatureDistribution[1,1,1:3]) > 0.75 & Param_ind@signatureFeatureDistribution[2,1,4] > 0.75) {
  visPMSignature(Param_ind, 1, charSize = 1, isScale = TRUE, alpha = 2);
  dev.copy2eps(file="../result/Fig3/UCUT_APOBEC_ind.eps", height = 5, width = 7.5);
  visPMSignature(Param_ind, 2, charSize = 1, isScale = TRUE, alpha = 2);
  dev.copy2eps(file="../result/Fig3/UCUT_AA_ind.eps", height = 5, width = 7.5);
} else {
  visPMSignature(Param_ind, 2, charSize = 1, isScale = TRUE, alpha = 2);
  dev.copy2eps(file="../result/Fig3/UCUT_APOBEC_ind.eps", height = 5, width = 7.5);
  visPMSignature(Param_ind, 1, charSize = 1, isScale = TRUE, alpha = 2);
  dev.copy2eps(file="../result/Fig3/UCUT_AA_ind.eps", height = 5, width = 7.5);
} 



