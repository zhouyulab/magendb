#!/usr/bin/env Rscript
library(stringr)

args <- commandArgs(TRUE)

f_list <- args[1:10]
f_pdf <- args[11]
f_tsv <- args[12]

nums <- vector()
for (f in f_list) {
    data = read.table(file = f, header = TRUE)
    PSMs <- strsplit(as.character(data$PSMs), split = "[;]")
    PSMs <- lapply(PSMs, as.numeric)
    PSMs <- as.vector(unlist(lapply(PSMs, sum)))
    nums <- c(nums, PSMs)
}

### count 1 - 10 peps gene
counts <- vector()
for (i in 1:9) {
    pep_num <- length(nums[nums == i])
    counts <- c(counts, pep_num)
}
pep_num <-length(nums[nums >= 10])
counts <- c(counts, pep_num)
print(counts)

### plot
pdf(file = f_pdf, width = 8, height = 6)

#plot(density(log2(nums)), col = "#CCCCCC", lwd = 2, xlim = c(0, 8), ylim = c(0, 1.5), axes = T, xlab = "PSMs per gene (log2)", ylab = "Density", main = "")
hist(nums, breaks = 10, )


axis(1, seq(-1, 8, 1), c("", seq(0, 8 ,1)), lwd = 2, cex.axis=1.1)
#axis(2, seq(0, 16, 2), seq(0, 16, 2), lwd = 2, cex.axis=1.1)

dev.off()

### record source data
df <- data.frame("PSMs_per_gene" = nums)
write.table(df, file = f_tsv, row.names = F, sep = "\t")
