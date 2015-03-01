library(ggplot2);
library(dplyr);


mutNum <- c(10, 25, 50, 100, 250, 500, 1000);
sampleNum <- c(10, 25, 50, 100);
alpha <- c(0.5, 1, 2);
gamma <- c(0.5, 1, 2);

MutNum <- rep(0, 25200);
SampleNum <- rep(0, 25200);
Alpha <- rep(0, 25200);
Gamma <- rep(0, 25200);
Dist <- rep(0, 25200);

Ind <- 1:100;

for (i1 in 1:7) {
  for (i2 in 1:4) {
    for (i3 in 1:3) {
      for (i4 in 1:3) {
        a <- read.table(paste("../../simulation/result/cosineDist/", mutNum[i1], "_", sampleNum[i2], "_", alpha[i3], "_", gamma[i4], ".txt", sep=""), sep="\t", header=FALSE);

        Dist[Ind] <- rowSums(a) / 4;
        MutNum[Ind] <- mutNum[i1];
        SampleNum[Ind] <- sampleNum[i2];
        Alpha[Ind] <- alpha[i3];
        Gamma[Ind] <- gamma[i4];
        
        Ind <- Ind + 100;
      }
    }
  }
}



simulation.result <- data.frame(mutNum=MutNum, sampleNum=SampleNum, alpha=Alpha, gamma=Gamma, distance=Dist);
summary_simulation.result <- simulation.result  %.% group_by(mutNum, sampleNum, alpha, gamma) %.% summarise(cosine_similality = mean(distance), seDist = sd(distance));
summary_simulation.result$seDist2 <- sapply(summary_simulation.result$cosine_similality + summary_simulation.result$seDist, min, 1) - summary_simulation.result$cosine_similality;


summary_simulation.result.1_1 <- summary_simulation.result %.% filter(alpha == 1 & gamma == 1);
summary_simulation.result.1_1$seDist2 <- sapply(summary_simulation.result.1_1$cosine_similality + summary_simulation.result.1_1$seDist, min, 1) - summary_simulation.result.1_1$cosine_similality;


pd <- position_dodge(0.2);
ggplot(summary_simulation.result, aes(x=factor(mutNum), y=cosine_similality, ymax=1, ymin = 0.3, fill=factor(sampleNum), colour=factor(sampleNum), group=factor(sampleNum))) +
  geom_line(position=pd) +
  geom_point(position=pd, size=2, shape=21) +
  geom_errorbar(aes(ymin=cosine_similality - seDist, ymax=cosine_similality + seDist2), width=0.01, position=pd) + 
  facet_grid(gamma ~ alpha, labeller = label_both) +
  xlab("#mutation") +
  ylab("mean cosine similarity") +
  labs(fill="#sample", colour="#sample") +
  theme_bw() +
  theme(plot.title = element_text(size=rel(2)),
        axis.title = element_text(size = rel(1.2)), 
        axis.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        strip.text= element_text(face="bold", size=rel(1.2)));

ggsave("../../manuscript/simulation_result.eps", width=15, height=6);

