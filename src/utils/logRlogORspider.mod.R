logRlogORspider.mod <- function (cncf, dipLogR = 0, nfrac = 0.005, xlim=NULL, ylim=NULL)
{
    rho <- seq(0, 0.95, by = 0.01)
    nrho <- length(rho)
    logACR <- logCNR <- matrix(0, nrho, 19)
    l <- 1
    logCNR[, l] <- log2(2 * (1 - rho) + 1 * rho) - 1
    logACR[, l] <- log(1/(1 - rho))
    for (i in 2:7) {
        for (j in 0:floor(i/2)) {
            l <- l + 1
            logCNR[, l] <- log2(2 * (1 - rho) + i * rho) - 1
            logACR[, l] <- log(1 - rho + (i - j) * rho) - log(1 -
                rho + j * rho)
        }
    }
    if (is.null(xlim) && is.null(ylim)) {
    plot(c(-0.95, 1.8), c(0, 5.2), type = "n", xlab = "Expected(logR - dipLogR)",
        ylab = " Expected(|logOR|)")
    } else {
    plot(c(-0.95, 1.8), c(0, 5.2), type = "n", xlab = "Expected(logR - dipLogR)",
        ylab = " Expected(|logOR|)", xlim=xlim, ylim=ylim)
    }

    l <- 1
    i <- 1
    j <- 0
    linecols <- c("black", "cyan3", "green3", "blue")
    lines(logCNR[, l], logACR[, l], lty = 1, col = j + 1, lwd = 1.25)
    text(logCNR[nrho, l] + 0.03, logACR[nrho, l], paste(i, j,
        sep = "-"), cex = 0.65)
    for (i in 2:7) {
        for (j in 0:floor(i/2)) {
            l <- l + 1
            lines(logCNR[, l], logACR[, l], lty = i - 1, col = linecols[j +
                1], lwd = 1.25)
            text(logCNR[nrho, l] + 0.03, logACR[nrho, l], paste(i,
                j, sep = "-"), cex = 0.65)
        }
    }
    nsnps <- sum(cncf$num.mark)
    nhets <- sum(cncf$nhet)
    ii <- cncf$num.mark > nfrac * nsnps & cncf$nhet > nfrac *
        nhets
    cex <- 0.3 + 2.7 * (cncf$num.mark[ii]/sum(0.1 * cncf$num.mark[ii]))
    chrcol <- rainbow(24)
    points(cncf$cnlr.median[ii] - dipLogR, sqrt(abs(cncf$mafR[ii])),
        cex = cex, pch = 10, col = chrcol[cncf$chrom[ii]], lwd = 1.5)

    text(cncf$cnlr.median[ii] - dipLogR, sqrt(abs(cncf$mafR[ii])),
         labels=cncf$chrom[ii],
        pos=3, col = chrcol[cncf$chrom[ii]], lwd = 1.5, cex=0.5)

    legend(-1, 5.25, paste("chr", c(1:22, "X"), sep = ""), ncol = 4,
        pch = 10, col = chrcol[1:23], cex = 0.65)
}
