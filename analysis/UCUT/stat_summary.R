
# strand direction information
## ind model

load("summary/hoang_ind5_dir_1_1.Rdata");

strandInfo <- data.frame(sig=c("AA", "APOBEC"), type=c("ind", "ind"), strand=c("plus", "plus"), prob=c(lret[[1]][[1]][1,6,1], lret[[1]][[1]][2,6,1]));

temp <- data.frame(sig=c("AA", "APOBEC"), type=c("ind", "ind"), strand=c("minus", "minus"), prob=c(lret[[1]][[1]][1,6,2], lret[[1]][[1]][2,6,2]));
strandInfo <- rbind(strandInfo, temp);



load("summary/hoang_full5_dir_1_1.Rdata");

plus_AA_full <- sum(lret[[1]][[1]][1,,][seq(1, 3072, 2)]);
plus_APOBEC_full <- sum(lret[[1]][[1]][2,,][seq(1, 3072, 2)]);

temp <- data.frame(sig=c("AA", "APOBEC"), type=c("full", "full"), strand=c("plus", "plus"), prob=c(plus_AA_full, plus_APOBEC_full));
strandInfo <- rbind(strandInfo, temp);

minus_AA_full <- sum(lret[[1]][[1]][1,,][seq(2, 3072, 2)]);
minus_APOBEC_full <- sum(lret[[1]][[1]][2,,][seq(2, 3072, 2)]);

temp <- data.frame(sig=c("AA", "APOBEC"), type=c("full", "full"), strand=c("minus", "minus"), prob=c(minus_AA_full, minus_APOBEC_full));
strandInfo <- rbind(strandInfo, temp);


ggplot(strandInfo %.% filter(sig=="AA"), aes(x=type, y=prob, fill=strand)) + geom_bar(stat="identity");
> ggplot(strandInfo %.% filter(sig=="APOBEC"), aes(x=type, y=prob, fill=strand)) + geom_bar(stat="identity");

# 2 base 5' to the mutated site
## ind model 

load("summary/hoang_ind5_dir_1_1.Rdata");

baseTwoFivePrime <- data.frame(sig=rep("AA", 4), type=rep("ind", 4), base=c("A", "C", "G", "T"), prob=lret[[1]][[1]][1,2,1:4]);

temp <- data.frame(sig=rep("APOBEC", 4), type=rep("ind", 4), base=c("A", "C", "G", "T"), prob=lret[[1]][[1]][2,2,1:4]);
baseTwoFivePrime <- rbind(baseTwoFivePrime, temp);

## full model

load("summary/hoang_full5_dir_1_1.Rdata");

tempA <- sum(lret[[1]][[1]][1,1,0:3071 %% 32 >= 0 & 0:3071 %% 32 < 8]);
tempC <- sum(lret[[1]][[1]][1,1,0:3071 %% 32 >= 8 & 0:3071 %% 32 < 16]);
tempG <- sum(lret[[1]][[1]][1,1,0:3071 %% 32 >= 16 & 0:3071 %% 32 < 24]);
tempT <- sum(lret[[1]][[1]][1,1,0:3071 %% 32 >= 24 & 0:3071 %% 32 < 32]);

temp <- data.frame(sig=rep("AA", 4), type=rep("full", 4), base=c("A", "C", "G", "T"), prob=c(tempA, tempC, tempG, tempT));
baseTwoFivePrime <- rbind(baseTwoFivePrime, temp);


tempA <- sum(lret[[1]][[1]][2,1,0:3071 %% 32 >= 0 & 0:3071 %% 32 < 8]);
tempC <- sum(lret[[1]][[1]][2,1,0:3071 %% 32 >= 8 & 0:3071 %% 32 < 16]);
tempG <- sum(lret[[1]][[1]][2,1,0:3071 %% 32 >= 16 & 0:3071 %% 32 < 24]);
tempT <- sum(lret[[1]][[1]][2,1,0:3071 %% 32 >= 24 & 0:3071 %% 32 < 32]);

temp <- data.frame(sig=rep("APOBEC", 4), type=rep("full", 4), base=c("A", "C", "G", "T"), prob=c(tempA, tempC, tempG, tempT));
baseTwoFivePrime <- rbind(baseTwoFivePrime, temp);


ggplot(baseTwoFivePrime %.% filter(sig=="APOBEC"), aes(x=type, y=prob, fill=base)) +
  + geom_bar(stat="identity");


ggplot(baseTwoFivePrime %.% filter(sig=="APOBEC"), aes(x=type, y=prob, fill=base)) +
  geom_bar(stat="identity") +
  theme(legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.text.y = element_text(size=rel(1.5)),
        axis.title.y = element_text(size=rel(1.5)))

ggsave("APOBEC_baseTwoFivePrime.eps", width=8, height=8);

###############################

vis_full5 <- function(Fvec = rep(1 / 1536, 1536)) {
  
  X <- data.frame(probability = Fvec);
  X$subtype <- factor(c(rep("C>A", 256), rep("C>G", 256), rep("C>T", 256), rep("T>A", 256), rep("T>C", 256), rep("T>G", 256)));
  X$flank <- rep(as.vector(1:256), 6);
  
  ggplot(X, aes(x=flank, y=probability, fill=subtype)) +
    geom_bar(stat="identity", position="identity") + 
    facet_grid(. ~ subtype) + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=rel(1.2)),
          axis.title.y = element_text(size=rel(1.2)),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.text= element_text(face="bold", size=rel(1.2))) +
    guides(fill=FALSE);
  
}



load("example//breast_sanger_full5/4_1_1.Rdata")

tempA <- sum(lret[[1]][[1]][1,1,0:1535 %% 16 >= 0 & 0:1535 %% 16 < 4]);
tempC <- sum(lret[[1]][[1]][1,1,0:1535 %% 16 >= 4 & 0:1535 %% 16 < 8]);
tempG <- sum(lret[[1]][[1]][1,1,0:1535 %% 16 >= 8 & 0:1535 %% 16 < 12]);
tempT <- sum(lret[[1]][[1]][1,1,0:1535 %% 16 >= 12 & 0:1535 %% 16 < 16]);

temp <- data.frame(sig=rep("AA", 4), type=rep("full", 4), base=c("A", "C", "G", "T"), prob=c(tempA, tempC, tempG, tempT));
baseTwoFivePrime <- rbind(baseTwoFivePrime, temp);


tempA <- sum(lret[[1]][[1]][2,1,0:3071 %% 32 >= 0 & 0:3071 %% 32 < 8]);
tempC <- sum(lret[[1]][[1]][2,1,0:3071 %% 32 >= 8 & 0:3071 %% 32 < 16]);
tempG <- sum(lret[[1]][[1]][2,1,0:3071 %% 32 >= 16 & 0:3071 %% 32 < 24]);
tempT <- sum(lret[[1]][[1]][2,1,0:3071 %% 32 >= 24 & 0:3071 %% 32 < 32]);



