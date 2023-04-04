library(rhierbaps)
library(ggtree)



setwd("C:/Users/lauri/Desktop/FINAL PROJ")
#par(mar=c(5, 6, 6, 5), xpd = TRUE)

#fasta.file.name <- "small_single.fa"
fasta.file.name <- "big_single.fa"
snp.matrix <- load_fasta(fasta.file.name)
snp.matrix [snp.matrix == 'a'] <- 'A'
snp.matrix [snp.matrix == 'g'] <- 'G'
snp.matrix [snp.matrix == 'c'] <- 'C'
snp.matrix [snp.matrix == 't'] <- 'T'
WT_minus <- 0
WT_plus <- 0
MutS_minus <- 0
MutS_plus <- 0

metadata <- read.csv("metadata1.csv")
selection <- metadata[39]
cipro <- c("WT cipro -", "WT cipro +", "MutS cipro -", "MutS cipro +")


for(x in 1:4){
  counts <- c(0,0,0,0,0,0,0,0,0,0,0,0) # count transversion pairs
  snp_x <- snp.matrix
  delete_counter <- 0
  noSNP <- 0
  for(y in 1:nrow(snp.matrix)){
    if(selection[rownames(snp.matrix)[y],] != cipro[x]){
      snp_x <- snp_x[-(y-delete_counter),]
      delete_counter <- delete_counter + 1
    }
  }
  A <- colCounts(snp_x, value = 'A')
  G <- colCounts(snp_x, value = 'G')
  T <- colCounts(snp_x, value = 'T')
  C <- colCounts(snp_x, value = 'C')
  for(i in 1:1536){
    # A majority
    if(A[i] > G[i] && A[i] > T[i] && A[i] > C[i]){
      # G minority
      if(T[i] == 0 && C[i] == 0){
        counts[1] <- counts[1] + G[i]
        if(G[i] == 0){
          noSNP <- noSNP + 1
        }
      }
      # T minority
      else if(G[i] == 0 && C[i] == 0){
        counts[6] <- counts[6] + T[i]
        if(T[i] == 0){
          noSNP <- noSNP + 1
        }
      }
      # C minority
      else if(G[i] == 0 && T[i] == 0){
        counts[5] <- counts[5] + C[i]
        if(C[i] == 0){
          noSNP <- noSNP + 1
        }
      }
    }
    # G majority
    else if(G[i] > A[i] && G[i] > T[i] && G[i] > C[i]){
      # A minority
      if(T[i] == 0 && C[i] == 0){
        counts[2] <- counts[2] + A[i]
        if(A[i] == 0){
          noSNP <- noSNP + 1
        }
      }
      # T minority
      else if(A[i] == 0 && C[i] == 0){
        counts[10] <- counts[10] + T[i]
        if(T[i] == 0){
          noSNP <- noSNP + 1
        }
      }
      # C minority
      else if(A[i] == 0 && T[i] == 0){
        counts[9] <- counts[9] + C[i]
        if(C[i] == 0){
          noSNP <- noSNP + 1
        }
      }
    }
    # T majority
    else if(T[i] > A[i] && T[i] > G[i] && T[i] > C[i]){
      # A minority
      if(G[i] == 0 && C[i] == 0){
        counts[11] <- counts[11] + A[i]
        if(A[i] == 0){
          noSNP <- noSNP + 1
        }
      }
      # G minority
      else if(A[i] == 0 && C[i] == 0){
        counts[12] <- counts[12] + G[i]
        if(G[i] == 0){
          noSNP <- noSNP + 1
        }
      }
      # C minority
      else if(A[i] == 0 && G[i] == 0){
        counts[4] <- counts[4] + C[i]
        if(C[i] == 0){
          noSNP <- noSNP + 1
        }
      }
    }
    # C majority
    else if(C[i] > A[i] && C[i] > G[i] && C[i] > T[i]){
      # A minority
      if(G[i] == 0 && T[i] == 0){
        counts[7] <- counts[7] + A[i]
        if(A[i] == 0){
          noSNP <- noSNP + 1
        }
      }
      # G minority
      else if(A[i] == 0 && T[i] == 0){
        counts[8] <- counts[8] + G[i]
        if(G[i] == 0){
          noSNP <- noSNP + 1
        }
      }
      # T minority
      else if(A[i] == 0 && G[i] == 0){
        counts[3] <- counts[3] + T[i]
        if(T[i] == 0){
          noSNP <- noSNP + 1
        }
      }
    }
    else{
      print("NA")
    }
  }
  if(x == 1){
    WT_minus <- counts
  }
  else if(x == 2){
    WT_plus <- counts
  }
  else if(x == 3){
    MutS_minus <- counts
  }
  else if(x == 4){
    MutS_plus <- counts
  }
  print(counts)
  print(noSNP)
}

A_G <- c(WT_minus[1], WT_plus[1], MutS_minus[1], MutS_plus[1])
G_A <- c(WT_minus[2], WT_plus[2], MutS_minus[2], MutS_plus[2])
C_T <- c(WT_minus[3], WT_plus[3], MutS_minus[3], MutS_plus[3])
T_C <- c(WT_minus[4], WT_plus[4], MutS_minus[4], MutS_plus[4])
A_C <- c(WT_minus[5], WT_plus[5], MutS_minus[5], MutS_plus[5])
A_T <- c(WT_minus[6], WT_plus[6], MutS_minus[6], MutS_plus[6])
C_A <- c(WT_minus[7], WT_plus[7], MutS_minus[7], MutS_plus[7])
C_G <- c(WT_minus[8], WT_plus[8], MutS_minus[8], MutS_plus[8])
G_C <- c(WT_minus[9], WT_plus[9], MutS_minus[9], MutS_plus[9])
G_T <- c(WT_minus[10], WT_plus[10], MutS_minus[10], MutS_plus[10])
T_A <- c(WT_minus[11], WT_plus[11], MutS_minus[11], MutS_plus[11])
T_G <- c(WT_minus[12], WT_plus[12], MutS_minus[12], MutS_plus[12])
multibar <- data.frame(A_G, G_A, C_T, T_C, A_C, A_T, C_A, C_G, G_C, G_T, T_A, T_G)


TRN <- barplot(as.matrix(multibar),
        # setting y label only
        # because x-label will be our
        # barplots name
        ylab="Count",
        yaxp=c(0, 350, 14),
        ylim = c(0, 360),
        xlim = c(0, 70),
        names=c("A - G", "G - A", "C - T", "T - C", "A - C", "A - T", "C - A", "C - G", "G - C", "G - T", "T - A", "T - G"),
        col = c("#196F3D", "#A9DFBF", "#5B2C6F", "#D2B4DE"),
        # to plot the bars separately http://127.0.0.1:28349/graphics/plot_zoom_png?width=1904&height=981
        beside=TRUE,
        # legend.text = TRUE,
        # args.legend = list(x = "right", inset = c(-0.15, 0), title = "test", adj[0, 1])
        cex.axis = 1.5,
        cex.names = 1.5,
        cex.lab = 1.5
)



legend("right",
       fill=c("#196F3D", "#A9DFBF", "#5B2C6F", "#D2B4DE"),
       c("WT cipro -", "WT cipro +", "MutS cipro -", "MutS cipro +"),
       bty="n",
       title = "Strain/Selection",
       title.adj = 0.15,
       x.intersp = 0.1,
       y.intersp = 0.6,
       inset = c(-0.2, -0.25),
       cex = 1.5
)

# for(x in 1:12){
#   text(TRN, labels = multibar[x], pos = 3)
# }
#text(40.5, 210 , label = "Strain/Selection",cex = 1.5)
text(TRN + 0.3, as.matrix(multibar), labels = as.matrix(multibar), pos = 3, cex = 1.5, srt = 90, offset = 1.5)