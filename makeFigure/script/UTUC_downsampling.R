library(ggplot2);
library(dplyr);

ratio <- c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5);

Dist <- rep(0, 200 * 6 * 2);
Type <- rep(0, 200 * 6 * 2);
Ratio <- rep(0, 200 * 6 * 2);
Sig <- rep(0, 200 * 6 * 2);


# full5: 1: AA 2: APOBEC
Type[1:1200] <- "full";
for (i in 1:6) {
  Inds <- 1:100 + 200 * (i - 1);
  a <- read.table(paste("../../UTUC/result/cosineDist_full/3_", ratio[i], ".txt", sep=""), sep="\t", header=TRUE);
  Ratio[Inds] <- ratio[i];
  Ratio[Inds + 100] <-ratio[i];                       
  Dist[Inds] <- a[,"AA"];
  Dist[Inds + 100] <- a[,"APOBEC"];
  Sig[Inds] <- "AA";
  Sig[Inds + 100] <- "APOBEC";
}

# ind5: 1: AA 2: APOBEC
Type[1201:2400] <- "ind";
for (i in 1:6) {
  Inds <- 1:100 + 200 * (i - 1) + 1200;
  a <- read.table(paste("../../UTUC/result/cosineDist_ind/3_", ratio[i], ".txt", sep=""), sep="\t", header=TRUE);
  Ratio[Inds] <- ratio[i];
  Ratio[Inds + 100] <-ratio[i];                       
  Dist[Inds] <- a[,"AA"];
  Dist[Inds + 100] <- a[,"APOBEC"];
  Sig[Inds] <- "AA";
  Sig[Inds + 100] <- "APOBEC";
}


UTUC.downsampling <- data.frame(type=Type, ratio=Ratio, distance=Dist, signature=Sig);
summary_UTUC.downsampling <- UTUC.downsampling %>% group_by(type, ratio, signature) %>% summarise(cosine_similality = mean(distance), seDist = sd(distance));

summary.AA <- summary_UTUC.downsampling %.% filter(signature=="AA");
summary.AA$seDist_max <- sapply(summary.AA$cosine_similality + summary.AA$seDist, min, 1) - summary.AA$cosine_similality;
summary.AA$seDist_min <- - sapply(summary.AA$cosine_similality - summary.AA$seDist, max, 0) + summary.AA$cosine_similality;

ggplot(summary.AA, aes(x=ratio, y=cosine_similality, colour=type, group=type, fill=type)) +
  geom_line() +
  geom_point(size=4, shape=21) +
  geom_errorbar(aes(ymin=cosine_similality - seDist_min, ymax=cosine_similality + seDist_max), width=0.01) + 
  ggtitle("AA") + 
  ylim(0, 1) +
  theme_bw() +
  theme(plot.title = element_text(size=rel(2)),
        axis.title = element_text(size = rel(2)), 
        axis.text = element_text(size = rel(1.5)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        legend.position = c(1, 0), legend.justification=c(1, 0),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.5)));

ggsave("../../manuscript/UTUC_downsampling_AA.eps", width=7.5, height=5.5, units = "in");


summary.APOBEC <- summary_UTUC.downsampling %>% filter(signature=="APOBEC");
summary.APOBEC$seDist_max <- sapply(summary.APOBEC$cosine_similality + summary.APOBEC$seDist, min, 1) - summary.APOBEC$cosine_similality;
summary.APOBEC$seDist_min <- - sapply(summary.APOBEC$cosine_similality - summary.APOBEC$seDist, max, 0) + summary.APOBEC$cosine_similality;


ggplot(summary.APOBEC, aes(x=ratio, y=cosine_similality, colour=type, group=type, fill=type)) +
  geom_line() +
  geom_point(size=4, shape=21) +
  geom_errorbar(aes(ymin=cosine_similality - seDist_min, ymax=cosine_similality + seDist_max), width=0.01) + 
  ggtitle("APOBEC") + 
  ylim(0, 1) +
  theme_bw() +
  theme(plot.title = element_text(size=rel(2)),
        axis.title = element_text(size = rel(2)), 
        axis.text = element_text(size = rel(1.5)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        legend.position = c(1, 0), legend.justification=c(1, 0),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.5)));

ggsave("../../manuscript/UTUC_downsampling_APOBEC.eps", width=7.5, height=5.5, units = "in");



