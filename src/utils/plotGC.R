plotGC <- function(gc_fn, gc_csv_fn, xx) {

    png(gc_fn, width=4000, height=2000, res = 200)
    par(mfrow=c(1,2))
    # Raw LogR
    xx$jointseg$raw.logr <- log2(xx$jointseg$rCountT /  xx$jointseg$rCountN)
    plot(xx$jointseg$gcpct, xx$jointseg$raw.logr , xlab='GC', ylab='raw LogR', main=opt$sample_id, ylim=c(-2,2), pch=19)
    out.lowess.tiled <- lowess(x=xx$jointseg$gcpct,
                             y=xx$jointseg$raw.logr,
                             f=0.05,  delta=0)
    points(out.lowess.tiled, col='orange', pch=19, cex=0.2)

    # Normalized LogR
    plot(xx$jointseg$gcpct, xx$jointseg$cnlr, xlab='GC', ylab='normalized LogR', main=opt$sample_id, ylim=c(-2,2), pch=19)

    out.lowess.tiled <- lowess(x=xx$jointseg$gcpct,
                             y=xx$jointseg$cnlr,
                             f=0.05,  delta=0)
    points(out.lowess.tiled, col='green', pch=19, cex=0.2)
    dev.off()


    gc.measure <- cor(xx$jointseg$raw.logr, xx$jointseg$cnlr, method = "pearson")
    cat(gc.measure, file=gc_csv_fn)

}
