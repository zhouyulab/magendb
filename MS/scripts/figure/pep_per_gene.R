#!/usr/bin/env Rscript
library(stringr)

args <- commandArgs(TRUE)

f_in <- args[1]
f_pdf <- args[2]
f_tsv <- args[3]

data = read.table(file = f_in, header = TRUE)

### count 1 - 10 peps gene
counts <- vector()
for (i in 1:9) {
    pep_num <- nrow(data[data$peptide_per_gene == i,])
    print(pep_num)
    counts <- c(counts, pep_num)
}

pep_num <- nrow(data[data$peptide_per_gene >= 10, ])
print(pep_num)
counts <- c(counts, pep_num)
counts_log2 <- log2(counts)

### barplot
pdf(file = f_pdf, width = 8, height = 6)

barplot(counts)
#barplot(counts_log2, col = "#CCCCCC", beside=F, horiz=F, axes=F, xlab = "Peptides per gene", ylab = "Gene number (log2)", xlim = c(0, 10*1.2-0.5), ylim = c(0, 16))
#axis(1, seq(0, 10, 1)*1.2-0.5, c("", as.character(seq(1, 9, 1)), ">=10"), lwd = 2, cex.axis=1.1)
#axis(2, seq(0, 16, 2), seq(0, 16, 2), lwd = 2, cex.axis=1.1)

dev.off()

### record source data
df <- data.frame("Peptides_per_gene" = c(as.character(seq(1, 9, 1)), ">=10"), "Gene_number_log2" = counts_log2)
write.table(df, file = f_tsv, row.names = F, sep = "\t")
