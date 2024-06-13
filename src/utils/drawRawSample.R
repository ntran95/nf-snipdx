drawRawSample <- function(pl, sample_id, plotFolderMain, smoothWind=501) {
  # columns
  # include: TRUE/FALSE
  # File2A, File2R

  chroms <- paste0('chr', c(1:22, 'X'))
  chrMaxPos <- tapply(pl$Position, pl$Chromosome, max)[chroms]

  chrShift <- cumsum(c(0, chrMaxPos[1:(length(chrMaxPos) - 1)]))
  names(chrShift) <- names(chrMaxPos)

  pl$af <- pl$File2A / (pl$File2A + pl$File2R)

  # count approximate number of heterozygotes
  nHets <- sum(pl$af > 0.05 & pl$af > 0.95, na.rm = TRUE)
  noSnpsTotal <- nrow(pl)
  noSnps100 <- sum((pl$File2A + pl$File2R) > 100, na.rm = TRUE)

  pl$Position.cum <- pl$Position + chrShift[pl$Chromosome]

  pl$col <- 'gray'
  pl$col[pl$include] <- 'black'

  plot_file <- file.path(plotFolderMain, paste0(sample_id, '.rawsummary.png'))
  message("... creating: ", plot_file)
  png(plot_file, width = 1600, height = 700)
  layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))

  # genome-wide plots
  plot(pl$Position.cum,
       pl$af,
       pch = 19,
       main = paste(sample_id,
                    'noSnps:', noSnpsTotal,
                    'noSnps>100x:',  noSnps100,
                    'noHets:', nHets) ,
       ylim = c(0, 1.1),
       cex = 1,
       cex.lab = 2,
       cex.axis = 2,
       cex.main = 2,
       ylab = 'allele fraction',
       xlab = 'position',
       col = pl$col)

  # draw chromosomes
  for (chr in chroms) {
    abline(v = chrShift[chr], col = 'violet')
    text(chrShift[chr] + 2e7,
         1.05,
         labels = chr,
         col = 'violet',
         pos = 3,
         srt = 90)
  }

  pl.inc <- subset(pl, include == TRUE)

  plot(pl.inc$Position.cum,
       RcppRoll::roll_mean(pl.inc$logR, smoothWind, fill = NA),
       pch = 19,
       main = 'LogR (501 rolling median)' ,
       ylim = c(-3, 3),
       cex = 1,
       cex.lab = 2,
       cex.axis = 2,
       cex.main = 2,
       ylab = 'LogR',
       xlab = 'position',
       col = pl.inc$col)

  abline(h = mean(pl.inc$logR), col = 'gray')

  # draw chromosomes
  for (chr in chroms) {
    abline(v = chrShift[chr], col = 'violet')
    text(chrShift[chr] + 2e7,
         1.05,
         labels = chr,
         col = 'violet',
         pos = 3,
         srt = 90)
  }
  dev.off()
}
