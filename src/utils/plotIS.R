plotIS <- function(fn, cluster.centers, cluser.df, draw=TRUE, plot_fn = NULL) {

    sample.id <- gsub('.insert.stats','',basename(fn))

    insert.sizes <- read.table(fn)
    insert.sizes.nz <- subset(insert.sizes, V1>0)
    insert.sizes.nz$cumsum <- cumsum(insert.sizes.nz$V2)

    mean.is <- sum(insert.sizes.nz$V1 * insert.sizes.nz$V2)/sum(insert.sizes.nz$V2)
    meadian.is <- min(which((insert.sizes.nz$cumsum>(sum(insert.sizes.nz$V2)/2))>0))
    first.qt.is <- min(which((insert.sizes.nz$cumsum>(sum(insert.sizes.nz$V2)/4))>0))
    third.qt.is <- min(which((insert.sizes.nz$cumsum>(sum(insert.sizes.nz$V2)*3/4))>0))
    qt1.is <- min(which((insert.sizes.nz$cumsum>(sum(insert.sizes.nz$V2)/10))>0))
    qt9.is <- min(which((insert.sizes.nz$cumsum>(sum(insert.sizes.nz$V2)*9/10))>0))

    v <- insert.sizes.nz$V2
    names(v) <- insert.sizes.nz$V1
    if (draw & !is.null(plot_fn)) {
        png (plot_fn,width=1600, height=2000, res = 200)
        par(mfrow=c(2,1))
    }
    if (draw) {

        plot(v/sum(v,na.rm=TRUE), type='l', ylim=c(0,0.01), lwd=3, col='black', main=sample.id, xlab='insert length [bp]', ylab='density')
        abline(v=first.qt.is, col='red')
        abline(v=meadian.is, col='red')
        abline(v=third.qt.is, col='red')
        abline(v=qt1.is, col='red')
        abline(v=qt9.is, col='red')
        abline(v=mean.is, col='blue')

        for (i in 1:nrow(cluster.centers)) {
            if (i==1) {
                plot(cluster.centers[i,], type='l', ylim=c(0,0.01), lwd=3, col=cluser.df[rownames(cluster.centers)[i],'col'] , xlab='insert length [bp]', ylab='density')
            } else {
                lines(cluster.centers[i,], lwd=2,col=cluser.df[rownames(cluster.centers)[i],'col'])
            }
        }

        lines(v/sum(v,na.rm=TRUE), lwd=2,col='black')

        legend('topright', legend=cluser.df$ascCode,
        col=cluser.df$col, lty=1, lwd=3,
        text.font=4)


    }
    if (draw & !is.null(plot_fn)) {
        dev.off()
    }

    r <- data.frame(sample.id=sample.id,
                     mean.is=mean.is,
                      meadian.is=meadian.is,
                      first.qt.is=first.qt.is,
                      third.qt.is=third.qt.is,
                      qt1.is=qt1.is,
                      qt9.is=qt9.is
                     )


    v <- insert.sizes.nz$V2[1:500]

    rl <- list()
    rl[['r']] <- r
    rl[['v']] <- v

    return(rl)
}