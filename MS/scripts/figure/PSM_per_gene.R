#!/usr/bin/env Rscript
library(stringr)

args <- commandArgs(TRUE)

f_in <- args[1]
f_pdf <- args[2]
f_tsv <- args[3]

data <- read.table(file = f_in, header = T, sep = "\t")


### plot
pdf(file = f_pdf, width = 8, height = 6)

plot(density(log2(data$PSMs_per_gene)), col = "#CCCCCC", lwd = 2, xlim = c(0, 8), ylim = c(0, 1.5), axes = T, xlab = "PSMs per gene (log2)", ylab = "Density", main = "")
#hist(nums, breaks = 10, )

axis(1, seq(-1, 8, 1), c("", seq(0, 8 ,1)), lwd = 2, cex.axis=1.1)
#axis(2, seq(0, 16, 2), seq(0, 16, 2), lwd = 2, cex.axis=1.1)

dev.off()

### record source data
#df <- data.frame("PSMs_per_gene" = nums)
#write.table(df, file = f_tsv, row.names = F, sep = "\t")
