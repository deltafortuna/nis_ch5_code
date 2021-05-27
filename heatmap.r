library(ggplot2)
library(cowplot)
library(reshape2)

heatmap <- function(datafile="sumstats") {
  d <- read.table(file = datafile, header = T);
  dd <- split(d, d$stat)
  size <- dim(dd$K)
  numcol = size[2] - 2
  kmat <- dd$K[,3:size[2]]
  kmat2 <- matrix(0, nrow = max(kmat)+5, ncol = numcol)
  for (i in 1:size[1]) {
    for (j in 1:numcol) {
      kmat2[kmat[i,j],j] <- kmat2[kmat[i,j],j] + 1
    }
  }
  kmat2.melted <- melt(kmat2)
  ggplot(kmat2.melted, aes(x = Var2, y = Var1, fill = value)) + geom_tile() + coord_equal() + scale_fill_gradient(low="white", high="black")
}
