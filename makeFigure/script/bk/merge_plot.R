
library(pmsignature);

.pardefault <- par(no.readonly = TRUE);
####################

load("/Users/friend1ws/Desktop/all.20.Rdata");
Fs <- Param@signatureFeatureDistribution;

# change the order of the signature
remainInd <- 1:dim(Fs)[1];
newOrder <- c();

for (k in remainInd) {
  if (Fs[k,1,1] > 0.6) {
    newOrder <- c(newOrder, k);
    remainInd <- remainInd[-which(remainInd == k)];
  }
}

for (k in remainInd) {
  if (Fs[k,1,2] > 0.6) {
    newOrder <- c(newOrder, k);
    remainInd <- remainInd[-which(remainInd == k)];
  }
}

for (k in remainInd) {
  if (Fs[k,1,3] > 0.6) {
    newOrder <- c(newOrder, k);
    remainInd <- remainInd[-which(remainInd == k)];
  }
}

for (k in remainInd) {
  if (sum(Fs[k,1,1:3]) > 0.75) {
    newOrder <- c(newOrder, k);
    remainInd <- remainInd[-which(remainInd == k)];
  }
}


for (k in remainInd) {
  if (sum(Fs[k,1,4]) > 0.6) {
    newOrder <- c(newOrder, k);
    remainInd <- remainInd[-which(remainInd == k)];
  }
}

for (k in remainInd) {
  if (sum(Fs[k,1,5]) > 0.6) {
    newOrder <- c(newOrder, k);
    remainInd <- remainInd[-which(remainInd == k)];
  }
}

for (k in remainInd) {
  if (sum(Fs[k,1,6]) > 0.6) {
    newOrder <- c(newOrder, k);
    remainInd <- remainInd[-which(remainInd == k)];
  }
}


for (k in remainInd) {
  if (sum(Fs[k,1,4:6]) > 0.75) {
    newOrder <- c(newOrder, k);
    remainInd <- remainInd[-which(remainInd == k)];
  }
}

newOrder <- c(newOrder, remainInd)
##########


#' generate the figure of mutation signatures

par(mar=c(0, 0, 0, 0));
par(bg = rgb(0.9, 0.9, 0.9));
par(xaxs = "i", yaxs = "i");
par(mfrow=c(ceiling(dim(Fs)[1] / 4), 4));

curSig <- 1;
# for (i in newOrder) {
for (i in 1:dim(Fs)[1]) {
  plot.new();
  plot.window(xlim=c(-0.3, 6.3), ylim=c(-0.3, 3.8));
  polygon(c(-0.25, 6.25, 6.25, -0.25), c(-0.25, -0.25, 3.75, 3.75), col = "white", border = FALSE);
  visPMS_ind(Fs[i,,], numBases = 5, trDir = FALSE, charSize = 1);
  # visPMS_ind5_mod3(Fs[[i]], 1 / 100000);
  mtext(paste("signature", curSig),
        outer = FALSE,      # 作図領域の外の余白に書く
        side = 3,          # 上の余白に書く
        cex = 1,         # 字の大きさ
        line = -1.5,          # 外に向かって 0.5行離れたところに書く．
        col = "black")    
  curSig <- curSig + 1;
  
}

##########
sample2type_temp <- read.table("../data/AlexandrovEtAl_sample2type.txt", sep="\t", header = FALSE)
sample2type <- sample2type_temp[,2];
names(sample2type) <- sample2type_temp[,1];
types <- as.character(sort(unique(sample2type_temp[,2])));

Q <- Param@sampleSignatureDistribution;
sampleList <- Param@sampleList;

colvec = c("orange",
           "blue",
           "yellow",
           "pink",
           "green",
           "purple",
           "red",
           "lightgreen",
           "darkblue",
           "mediumpurple",
           "lightyellow",
           "brown",
           "lightblue",
           "olivedrab",
           "palevioletred",
           "seagreen",
           "yellowgreen",
           "slateblue",
           "turquoise",
           "gray",
           "darkgreen",
           "lightgray",
           "red2",
           "lightblue2",
           "tan")


allQ <- matrix(0, dim(Q)[2], 0);

cumNum <- c();
tempNum <- 0;
for (t in types) {
  tempQ <- Q[which(sample2type[sampleList] == t), ];
  
  tempNum <- tempNum + length(which(sample2type[sampleList] == t));
  cumNum <- c(cumNum, tempNum);
  
  d <- dist(tempQ);
  h <- hclust(d);
  allQ <- cbind(allQ, t(tempQ)[,h$order]);
}

curMar <- par("mar");
par(mar = c(7.1, 7.1, 7.1, 7.1));
barplot(allQ, col=colvec, border = NA, space = 0);
legend("topright", fill = colvec[1:dim(Q)[2]], legend = 1:dim(Q)[2], box.lty = 0, text.width = 3, xpd = TRUE);



for (i in 1:length(cumNum)) {
  abline(v = cumNum[i], col = "black", lty = 4);
  
  if (i > 2) {
    text((cumNum[i] + cumNum[i - 1]) / 2, -0.02, labels = types[i], adj = c(0,0), srt = -45, xpd = TRUE);
  } else {
    text(cumNum[i] / 2, -0.02, labels = types[i], adj = c(0,0), srt = -45, xpd = TRUE)
  }
}


