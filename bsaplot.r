bsaplot <- function(
  key,
  chrom,
  chrmax,
  chrmin = 1,
  linef = NULL,
  linep = c(1, 2),
  linev = c(3, 4),
  linec = c("black", "blue", "red", "blue", "red"),
  pointf = NULL,
  pointp = c(1, 2),
  pointv = 3,
  pointc = c("#33A02C", "#FF7F00"),
  hline = NULL,
  hcolor = "grey",
  hlty = 1,
  ymin = -1,
  ymax = 1,
  xlab = "Chromosome",
  ylab = expression(Delta(SNPindex)),
  cexa = 1.5,
  cexl = 2,
  cexp = 0.4,
  pch = 19,
  lwd = 3,
  las = 1,
  bandc = "grey90",
  alt = F,
  width = 15,
  height = 6,
  res = 600
){
  library(gtools)
  multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
    do.call(order, c(
      lapply(list(...), function(l){
        if(is.character(l)){
          factor(l, levels=mixedsort(unique(l)))
        } else {
          l
        }
      }),
      list(na.last = na.last, decreasing = decreasing)
    ))
  }
  ## global deal
  chrom <- as.character(chrom)
  mindf <- data.frame(chr = chrom, pos = chrmin, stringsAsFactors = F)
  maxdf <- data.frame(chr = chrom, pos = chrmax, stringsAsFactors = F)
  minmax <- rbind(mindf, maxdf)
  minmax <- minmax[multi.mixedorder(minmax[,1], minmax[,2]), ]
  chrom <- mixedsort(chrom)
  mins <- NULL
  maxs <- NULL
  print("pass1: global")
  ## line deal
  if(!is.null(linef)){
    linet <- read.table(linef, header = T, stringsAsFactors = F, check.names = F)
    linet <- linet[,c(linep, linev)]
    linet[, 1] <- as.character(linet[, 1])
    linet <- linet[linet[, 1] %in% chrom, ]
    lminmax <- minmax
    for(i in 3:ncol(linet)){
      lminmax[,i] <- NA
    }
    names(lminmax) <- names(linet)
    linet <- rbind(linet, lminmax)
    linet <- linet[multi.mixedorder(linet[,1], linet[,2]),]
    mins <- min(linet[, 3:ncol(linet)], na.rm = T)
    maxs <- max(linet[, 3:ncol(linet)], na.rm = T)
    class(linet) <- c("scanone", "data.frame")
  }
  print("pass2: line data")
  ## point deal
  if(!is.null(pointf)){
    pointc <- rep(pointc, ceiling(length(chrom)/length(pointc)))[1:length(chrom)]
    chrcolor <- data.frame(chrom, pointc, stringsAsFactors = F)
    pointt <- read.table(pointf, header = T, stringsAsFactors = F, check.names = F)
    pointt <- pointt[,c(pointp, pointv)]
    pointt[, 1] <- as.character(pointt[, 1])
    pointt <- pointt[pointt[, 1] %in% chrom,]
    pminmax <- minmax
    for(i in 3:ncol(pointt)){
      pminmax[,i] <- NA
    }
    names(pminmax) <- names(pointt)
    names(pminmax) <- names(pointt)
    pointt <- rbind(pointt, pminmax)
    pointt <- pointt[multi.mixedorder(pointt[,1], pointt[,2]),]
    if(is.null(mins)){
    	mins <- min(pointt[,3], na.rm = T)
    }else{
    	mins <- min(min(pointt[,3], na.rm = T), mins)
    }
    if(is.null(maxs)){
    	maxs <- max(pointt[,3], na.rm = T)
    }else{
    	maxs <- max(max(pointt[,3], na.rm = T), maxs)
    }
    class(pointt) <- c("scanone", "data.frame")
  }
  print("pass3: point data")
  ## plot
  if(is.null(ymin)){
    ymin <- mins
  }
  if(is.null(ymax)){
    ymax <- maxs
  }
  minmax$lod <- NA
  class(minmax) <- c("scanone", "data.frame")

  # png
  png(filename = paste0(key, ".plplot.png"), width = width, height = height, units = "in", res = res)
  dev.control("enable")
  # global parameters
  chara_len <- max(nchar(chrom))
  marx <- 5
  if(chara_len > 6){
    marx <- marx * chara_len/5
  }
  par(mgp = c(2.5, 0.8, 0), mai = c(0.8, 1, 0.2, 0.2),
      mar = c(marx,5,0.8,0.8), cex.axis = cexa, cex.lab = cexl)
  # empty plot
  plot.scanone(x = minmax, type="n", ylim = c(ymin, ymax), bty = "n",
               bandcol = bandc, incl.markers = F, yaxs="i", alternate.chrid = alt,
               xlab = xlab, ylab = ylab, las = las)
  box()
  # point plot
  if(!is.null(pointf)){
    nchr <- length(chrcolor[,2])
    if(nchr > 1){
      for (i in 1:nchr) {
        ptmp <- pointt
        ptmp[ptmp[, 1] != chrcolor[i, 1], ][, 3] <- NA
        plot.scanone(ptmp, add = T, lodcolumn = 1, col = chrcolor[i, 2], type = "p", pch = pch, cex = cexp)
      }
    }else{
      ptmp <- pointt
      plot.scanone(ptmp, add = T, lodcolumn = 1, col = chrcolor[,2], type = "p", pch = pch, cex = cexp)
    }
    print("pass4: point plot")
  }
  # hline
  if(!is.null(hline)){
    abline(h = hline, lty = hlty, col = hcolor, lwd = lwd)
  }
  # line plot
  if(!is.null(linef)){
    for (l in 1:(ncol(linet) - 2)) {
      plot.scanone(linet, add=T, lodcolumn = l, col = linec[l], type = "l",
                   pch = pch, lwd = lwd, cex = cexp)
    }
    print("pass5: line plot")
  }
  # box
  #box()
  # copy to pdf
  dev.copy(pdf, file = paste0(key, ".plplot.pdf"), width = width, height = height)
  dev.off(which = dev.cur())
  # close all image device
  for (i in 1:length(dev.list())) {
    dev.off()
  }
}
