getSamplePrimerDistance <- function(pl, minCov=100, plot_fn, specifyCol=NULL, flagCol=NULL) {
# expecting the following columns:
# File2R
# File2A
# primerMinDistance

    png(plot_fn, width=1600, height=1200, res=200)
    if (is.null(specifyCol)) {
        smoothingSpline <-  smooth.spline(pl$primerMinDistance, pl$File2R + pl$File2A, spar=0.35)
    } else {
        smoothingSpline <-  smooth.spline(pl$primerMinDistance, pl[,specifyCol], spar=0.35)
    }

    if (is.null(specifyCol)) {
    plot(smoothingSpline,xlim=c(0,500), ylim=c(0,2500),
     xlab='distance from primer (3), bp', ylab='sequencing coverage (x)',
     pch=19,cex.lab=2, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col='black')
    } else {
        plot(smoothingSpline,xlim=c(0,500),
         xlab='distance from primer (3), bp', ylab='sequencing coverage (x)',
         pch=19,cex.lab=2, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col='black')
    }
    abline(h=minCov)

    if (!is.null(flagCol)) {
        pl.true <- pl[!is.na(pl[,flagCol]) & pl[,flagCol]==TRUE ,]
        smoothingSpline.t <-  smooth.spline(pl.true$primerMinDistance, pl.true[,specifyCol], spar=0.5)
        lines(smoothingSpline.t, col='green')

        pl.false <- pl[!is.na(pl[,flagCol]) & pl[,flagCol]==FALSE ,]
        smoothingSpline.f <-  smooth.spline(pl.false$primerMinDistance, pl.false[,specifyCol], spar=0.5)
        lines(smoothingSpline.f,
          col='salmon')
    }


    dev.off()

    if (all(smoothingSpline$y[1:400]>minCov)) {
        max.dist <- 400
    } else {
        first.below.thresh <- min(which(smoothingSpline$y[1:400]<minCov))
        max.dist <- smoothingSpline$x[first.below.thresh]
    }

    return(max.dist)

}