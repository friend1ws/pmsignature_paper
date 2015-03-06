library(pmsignature);

load("../../UTUC/result/Param_full/3.Rdata");
Param_full <- resultForSave[[1]];

if (sum(Param_full@signatureFeatureDistribution[1,1,(256 * 0 + 1):(256 * 3)]) + sum(Param_full@signatureFeatureDistribution[1,1,(256 * 6 + 1):(256 * 9)]) > 0.6 & sum(Param_full@signatureFeatureDistribution[2,1,(256 * 3 + 1):(256 * 4)]) + sum(Param_full@signatureFeatureDistribution[2,1,(256 * 9 + 1):(256 * 10)]) > 0.6) {
  visPMSignature(Param_full, 1);
  ggsave("../../manuscript/UTUC_APOBEC_full.eps", height = 5, width = 7.5);
  visPMSignature(Param_full, 2);
  ggsave("../../manuscript/UTUC_AA_full.eps", height = 5, width = 7.5); 
} else {
  visPMSignature(Param_full, 2);
  ggsave("../../manuscript/UTUC_APOBEC_full.eps", height = 5, width = 7.5);
  visPMSignature(Param_full, 1);
  ggsave("../../manuscript/UTUC_AA_full.eps", height = 5, width = 7.5);   
}
  
load("../../UTUC/result/Param_ind/3.Rdata");
Param_ind <- resultForSave[[1]];

if (sum(Param_ind@signatureFeatureDistribution[1,1,1:3]) > 0.75 & Param_ind@signatureFeatureDistribution[2,1,4] > 0.75) {
  visPMSignature(Param_ind, 1, charSize = 1);
  dev.copy2eps(file="../../manuscript/UTUC_APOBEC_ind.eps", height = 5, width = 7.5);
  visPMSignature(Param_ind, 2, charSize = 1);
  dev.copy2eps(file="../../manuscript/UTUC_AA_ind.eps", height = 5, width = 7.5);
} else {
  visPMSignature(Param_ind, 2, charSize = 1);
  dev.copy2eps(file="../../manuscript/UTUC_APOBEC_ind.eps", height = 5, width = 7.5);
  visPMSignature(Param_ind, 1, charSize = 1);
  dev.copy2eps(file="../../manuscript/UTUC_AA_ind.eps", height = 5, width = 7.5);
} 




# 2 base 5' to the mutated site
## ind model 

if (sum(Param_ind@signatureFeatureDistribution[1,1,1:3]) > 0.75 & Param_ind@signatureFeatureDistribution[2,1,4] > 0.75) {
  baseTwoFivePrime <- data.frame(model=rep("ind", 4), base=c("A", "C", "G", "T"), prob = Param_ind@signatureFeatureDistribution[1,2,1:4]);
} else {
  baseTwoFivePrime <- data.frame(model=rep("ind", 4), base=c("A", "C", "G", "T"), prob = Param_ind@signatureFeatureDistribution[2,2,1:4]);
}



## full model

if (sum(Param_full@signatureFeatureDistribution[1,1,(256 * 0 + 1):(256 * 3)]) + sum(Param_full@signatureFeatureDistribution[1,1,(256 * 6 + 1):(256 * 9)]) > 0.6 & sum(Param_full@signatureFeatureDistribution[2,1,(256 * 3 + 1):(256 * 4)]) + sum(Param_full@signatureFeatureDistribution[2,1,(256 * 9 + 1):(256 * 10)]) > 0.6) {
  tempA <- sum(Param_full@signatureFeatureDistribution[1,1,0:3071 %% 16 >= 0 & 0:3071 %% 16 < 4]);
  tempC <- sum(Param_full@signatureFeatureDistribution[1,1,0:3071 %% 16 >= 4 & 0:3071 %% 16 < 8]);
  tempG <- sum(Param_full@signatureFeatureDistribution[1,1,0:3071 %% 16 >= 8 & 0:3071 %% 16 < 12]);
  tempT <- sum(Param_full@signatureFeatureDistribution[1,1,0:3071 %% 16 >= 12 & 0:3071 %% 16 < 16]);
  temp <- data.frame(model=rep("full", 4), base=c("A", "C", "G", "T"), prob=c(tempA, tempC, tempG, tempT));
  baseTwoFivePrime <- rbind(baseTwoFivePrime, temp);
} else {
  tempA <- sum(Param_full@signatureFeatureDistribution[2,1,0:3071 %% 16 >= 0 & 0:3071 %% 16 < 4]);
  tempC <- sum(Param_full@signatureFeatureDistribution[2,1,0:3071 %% 16 >= 4 & 0:3071 %% 16 < 8]);
  tempG <- sum(Param_full@signatureFeatureDistribution[2,1,0:3071 %% 16 >= 8 & 0:3071 %% 16 < 12]);
  tempT <- sum(Param_full@signatureFeatureDistribution[2,1,0:3071 %% 16 >= 12 & 0:3071 %% 16 < 16]);
  temp <- data.frame(model=rep("full", 4), base=c("A", "C", "G", "T"), prob=c(tempA, tempC, tempG, tempT));
  baseTwoFivePrime <- rbind(baseTwoFivePrime, temp);
}



ggplot(baseTwoFivePrime, aes(x=model, y=prob, fill=base)) +
  geom_bar(stat="identity") +
  theme(legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.text.y = element_text(size=rel(1.5)),
        axis.title.y = element_text(size=rel(1.5)))

ggsave("../../supp/UTUC_APOBEC_twoFivePrime.eps", width=4, height=8);
