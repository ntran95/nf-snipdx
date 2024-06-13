normalizeWrtPrimer <- function(pl2, plot_fn) {
    pl2.clean <- subset(pl2, !is.na(closestPrimerIsSnpTiled) & !is.na(primerMinDistance))
    pl2.clean <- pl2.clean[order(pl2.clean$primerMinDistance),]

    pl2.clean.gene <- subset(pl2.clean, closestPrimerIsSnpTiled==FALSE)
    out.lowess.gene <- lowess(x=pl2.clean.gene$primerMinDistance,
                         y=pl2.clean.gene$r,
                         f=0.1,  delta=0)
    pl2.clean.gene$r.lowess <- out.lowess.gene$y

    pl2.clean.tiled <- subset(pl2.clean, closestPrimerIsSnpTiled==TRUE)
    out.lowess.tiled <- lowess(x=pl2.clean.tiled$primerMinDistance,
                         y=pl2.clean.tiled$r,
                         f=0.1,  delta=0)
    pl2.clean.tiled$r.lowess <- out.lowess.tiled$y

    pl2.clean.with.pred <- rbind(pl2.clean.gene, pl2.clean.tiled)

    pl2.clean.with.pred$r.norm <- pl2.clean.with.pred$r/pl2.clean.with.pred$r.lowess

    #loessMod03 <- loess(xR , data=pl2.clean)
    #smoothed03 <- predict(loessMod03)
    #


    #smoothingSpline <-  smooth.spline(pl2.clean$primerMinDistance, pl2.clean$logR, spar=0.35)

    png(plot_fn, width=2000, height=2000, res=200)
    par(mfrow=c(2,2))
    plot(x=pl2.clean.with.pred$primerMinDistance, y=pl2.clean.with.pred$r,  type="p", pch=19, cex=0.2, main="Raw data", xlab="Distance to primer", ylab="R")
    points(out.lowess.gene, col='green', pch=19, cex=0.2)
    points(out.lowess.tiled, col='salmon', pch=19, cex=0.2)

    plot(out.lowess.gene, col='green', pch=19, cex=0.2)
    points(out.lowess.tiled, col='salmon', pch=19, cex=0.2)

    plot(x=pl2.clean.with.pred$primerMinDistance, y=pl2.clean.with.pred$r.norm,  type="p", pch=19, cex=0.2, main="Raw data", xlab="Distance to primer", ylab="R")
    l2 <- lowess(x=pl2.clean.with.pred$primerMinDistance,
                         y=pl2.clean.with.pred$r.norm,
                         f=0.1,  delta=0)
    points(l2, col='blue', pch=19, cex=0.2)

    plot(pl2.clean.with.pred$File2R, pl2.clean.with.pred$File2R/pl2.clean.with.pred$r.lowess, col='blue', pch=19, cex=0.2, xlab='before normalization', ylab='after normalization')

    dev.off()
    pl2.clean.with.pred$File2R <- round(pl2.clean.with.pred$File2R/pl2.clean.with.pred$r.lowess)
    pl2.clean.with.pred$File2A <- round(pl2.clean.with.pred$File2A/pl2.clean.with.pred$r.lowess)
    pl2.clean.with.pred <- pl2.clean.with.pred[order(as.integer(as.character(pl2.clean.with.pred$Chromosome)), pl2.clean.with.pred$Position),]
    return(pl2.clean.with.pred)
}