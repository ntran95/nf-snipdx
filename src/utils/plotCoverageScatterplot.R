plotCoverageScatterplot <- function(pl2, plot_fn=NULL) {

    if (!is.null(plot_fn)) {
        png(plot_fn, width=1000, height=1000, res=200)
        plot(pl2$File1R + pl2$File1A, pl2$File2R + pl2$File2A, xlab='PON coverage (X)', ylab='Sample coverage', pch=10, cex=0.5,
             col=factor(pl2$Chromosome))
        dev.off()
    }

}
