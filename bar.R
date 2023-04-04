library(rhierbaps)

setwd("C:/Users/lauri/Desktop/FINAL PROJ")

#fasta.file.name <- "small_single.fa"
fasta.file.name <- "big_single.fa"
snp.matrix <- load_fasta(fasta.file.name)


SNP <- data.matrix(snp.matrix)
SNP [SNP == 'a'] <- 'A'
SNP [SNP == 'g'] <- 'G'
SNP [SNP == 'c'] <- 'C'
SNP [SNP == 't'] <- 'T'

#Nucleotide Frequency Among 315 High Confidence Variant Sites in Coding Regions
nucFreq <- barplot(table(SNP),
        yaxp=c(0, 160000, 8),
        ylim = c(0, 160000),
        xlab="Nucleotide",
        ylab="Count",
        border="#BBACF2",
        col="#ACC0F2",
        cex.axis = 1.5,
        cex.names = 1.5,
        cex.lab = 1.5
)

text(nucFreq, table(SNP), labels = table(SNP), pos = 3, cex = 1.5)

# dpi = 300
# ggsave(filename="bar_small.png", width=8500, height=7500, dpi = dpi, scale = 0.45, limitsize = FALSE, units = "px")