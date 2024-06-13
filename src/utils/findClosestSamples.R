findClosestSamples <- function(sample.v,
                                  sample.vecs.df,
                                  full.sample.id,
                                  normal.sample.df,
                                  plot_fn=NULL,
                                  top.n=10) {



    euclidean <- function(a, b) sqrt(sum((a - b)^2))


    sample.v[is.na(sample.v)] <- 0

    sample.dists <- vector()
    for (si in 1:nrow(normal.sample.df)) {

        nid <- as.character(normal.sample.df$Sample.Sequencing.Name[si])
        if (nid!=full.sample.id) {
            norm.v <- sample.vecs.df[as.character(nid),]
            d <- euclidean(sample.v/sum(sample.v, na.rm=TRUE),
                       norm.v/sum(norm.v, na.rm=TRUE))
            sample.dists[si] <- d
        } else {
            sample.dists[si] <- NA
        }
    }
    normal.sample.df$distance.to.x <- sample.dists
    normal.sample.df <- normal.sample.df[order(normal.sample.df$distance.to.x),]
    names(sample.dists) <- normal.sample.df$Sample.Sequencing.Name
    sample.dists <- sample.dists[names(sample.dists) != full.sample.id]
    sample.dists <- sort(sample.dists)

    normal.sample.df$isTop <- FALSE
    normal.sample.df[1:min(top.n, nrow(normal.sample.df)),'isTop'] <- TRUE
    top.normal.sample.df <- normal.sample.df[1:min(top.n, nrow(normal.sample.df)),]

    if (!is.null(plot_fn)) {
        png (plot_fn,width=1600, height=2000, res = 200)
        par(mfrow=c(2,1))

        barplot(normal.sample.df$distance.to.x , col='gray', border=NA, main='samples closest by insert size distribution' )
        abline(v=top.n+0.5)

        for(ni in 1:top.n) {

            nid <- as.character(top.normal.sample.df$Sample.Sequencing.Name[ni])
            norm.v <- sample.vecs.df[nid,]
            norm.v[is.na(norm.v)] <- 0
            if (ni==1) {
                plot(norm.v/sum(norm.v,na.rm=TRUE), type='l', ylim=c(0,0.01), lwd=1, col='gray',
                     main=full.sample.id, xlab='insert length [bp]', ylab='density')
            } else {
                lines(norm.v/sum(norm.v,na.rm=TRUE), lwd=1, col='gray')
            }
            lines(sample.v/sum(sample.v,na.rm=TRUE), col='black', lwd=2)

        }



        dev.off()


    }

    return(normal.sample.df)
}