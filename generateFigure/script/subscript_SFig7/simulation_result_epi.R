library(ggplot2);
library(dplyr);


featNum <- c(0, 1, 2, 3, 4, 5, 6, 7, 8);
sampleNum <- c(10, 25, 50, 100);
alpha <- c(0.5, 1, 2);
gamma <- c(0.5, 1, 2);

FeatNum <- rep(0, 21600);
SampleNum <- rep(0, 21600);
Alpha <- rep(0, 21600);
Gamma <- rep(0, 21600);
Dist <- rep(0, 21600);

Ind <- 1:100;

for (i1 in 1:6) {
  for (i2 in 1:4) {
    for (i3 in 1:3) {
      for (i4 in 1:3) {
        a <- read.table(paste("../../analysis/simulation/result/cosineDist_epi/", featNum[i1], "_", sampleNum[i2], "_", alpha[i3], "_", gamma[i4], ".txt", sep=""), sep="\t", header=FALSE);
        
        Dist[Ind] <- rowSums(a) / 4;
        FeatNum[Ind] <- featNum[i1];
        SampleNum[Ind] <- sampleNum[i2];
        Alpha[Ind] <- alpha[i3];
        Gamma[Ind] <- gamma[i4];
        
        Ind <- Ind + 100;
      }
    }
  }
}


simulation.result <- data.frame(featNum=FeatNum, distance=Dist);
summary_simulation.result <- simulation.result  %>% group_by(featNum) %>% summarise(cosine_similality = mean(distance), seDist = sd(distance));
summary_simulation.result$seDist2 <- sapply(summary_simulation.result$cosine_similality + summary_simulation.result$seDist, min, 1) - summary_simulation.result$cosine_similality;

simulation.result <- data.frame(featNum=FeatNum, sampleNum=SampleNum, alpha=Alpha, gamma=Gamma, distance=Dist);
summary_simulation.result <- simulation.result  %>% group_by(featNum, sampleNum, alpha, gamma) %>% summarise(cosine_similality = mean(distance), seDist = sd(distance));
summary_simulation.result$seDist2 <- sapply(summary_simulation.result$cosine_similality + summary_simulation.result$seDist, min, 1) - summary_simulation.result$cosine_similality;



pd <- position_dodge(0.2);
ggplot(summary_simulation.result, aes(x=factor(featNum), y=cosine_similality, ymax=1, ymin = 0.7, fill=factor(sampleNum), colour=factor(sampleNum), group=factor(sampleNum))) +
  geom_line(position=pd) +
  geom_point(position=pd, size=2, shape=21) +
  geom_errorbar(aes(ymin=cosine_similality - seDist, ymax=cosine_similality + seDist2), width=0.01, position=pd) + 
  facet_grid(gamma ~ alpha, labeller = label_both) +
  xlab("#additional features") +
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

ggsave("../result/SFig7/simulation_result_epi.eps", width=12, height=6);


