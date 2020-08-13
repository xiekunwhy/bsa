###########input################################
#replications <- 10000
#popStruc <- "ril"
#bulkSize <- c(50,60)
#interv <- c(95, 99)
#depth <- 10:500
#fvalue <- 0
#outfile <- "path"
###### prepare
interv <- sort(interv)
ninter <- length(interv)
interv_low <- c()
interv_up <- c()
thresh <- vector()
#inthead <- c("Delta_Low", "Delta_Up", "A_Index_Low", "A_Index_Up", "B_Index_Low", "B_Index_Up")
inthead <- c("Delta_Low", "Delta_Up", "A_Index_Up", "B_Index_Up")

for (int in interv) {
  if(int < 50){
    cat("ERROR: confidence interval can not < 50\n")
    stop()
  }
  if(int >= 100){
    cat("ERROR: confidence interval can not >= 100\n")
    stop()
  }
  tmp_low <- (100 - int)/100/2
  tmp_up <- 1 - tmp_low
  interv_low <- c(interv_low, tmp_low)
  interv_up <- c(interv_up, tmp_up)
  thresh <- c(thresh, paste0(inthead, int))
}
print(interv_low)
print(interv_up)
thresh <- c("depth", "b1size", "b2size", thresh)
## simulation allele freq
## modified from qtlseqr
simulateAlleleFreq <- function(n, pop = "ril") {
  if (pop == "f2") {
    mean(sample(x = c(0, 0.5, 1), size = n, prob = c(1, 2, 1), replace = TRUE))
  } else if(pop == "ril") {
    mean(sample(x = c(0, 1), size = n, prob = c(1, 1), replace = TRUE))
  } else if(pop == "bc"){
    mean(sample(x = c(0, 0.5), size = n, prob = c(1, 1), replace = TRUE))
  }
}
## snp index calculation
## modified from qtlseqr
simulateSNPindex <- function(depth,
                             altFreq1,
                             altFreq2,
                             filter = fvalue) {
  SNPindex_H <- rbinom(n = length(altFreq1), size = depth, prob = altFreq1) / depth
  SNPindex_L <- rbinom(n = length(altFreq2), size = depth, prob = altFreq2) / depth
  keep_index <- which(SNPindex_H >= filter | SNPindex_L >= filter)
  SNPindex_H <- SNPindex_H[keep_index]
  SNPindex_L <- SNPindex_L[keep_index]
  deltaSNP <- SNPindex_L - SNPindex_H
  return(list("H" = SNPindex_H, "L" = SNPindex_L, "D" = deltaSNP))
}
## do simulations
sink(file = outfile)
cat(thresh, sep = "\t")
cat("\n")
for (d in depth) {
  tmp_freq1 <- replicate(n = replications, simulateAlleleFreq(n = bulkSize[1], pop = popStruc))
  tmp_freq2 <- replicate(n = replications, simulateAlleleFreq(n = bulkSize[2], pop = popStruc))
  sim_results <- simulateSNPindex(depth = d, altFreq1 = tmp_freq1, altFreq2 = tmp_freq2)
  SNPindex_H <- sort(sim_results$H)
  SNPindex_L <- sort(sim_results$L)
  Delta <- sort(sim_results$D)
  ndelta <- length(Delta)
  tmp_results <- c(d, bulkSize[1], bulkSize[2])
  for (i in 1:ninter) {
    int_low <- interv_low[i]
    int_up <- interv_up[i]
    
    delta_low <- Delta[floor(int_low * ndelta)]
    delta_up <- Delta[ceiling(int_up * ndelta)]
    
    #index_A_low <- SNPindex_H[floor(int_low * ndelta)]
    index_A_up <- SNPindex_H[ceiling(int_up * ndelta)]
    
    #index_B_low <- SNPindex_L[floor(int_low * ndelta)]
    index_B_up <- SNPindex_L[ceiling(int_up * ndelta)]
    
    #tmp_results <- c(tmp_results, delta_low, delta_up, index_A_low, index_A_up, index_B_low, index_B_up)
    tmp_results <- c(tmp_results, delta_low, delta_up, index_A_up, index_B_up)
  }
  cat(tmp_results, sep = "\t")
  cat("\n")
}
sink()
