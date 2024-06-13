plotSampleModified <- function(x, emfit=NULL, clustered=FALSE, plot.type=c("em","naive","both","none"), sname=NULL) {
    def.par <- par(no.readonly = TRUE) # save default, for resetting...
    # plot.type
    plot.type <- match.arg(plot.type)
    # layout of multi panel figure
    if (plot.type=="none") layout(matrix(1:2, ncol=1))
    if (plot.type=="em") layout(matrix(rep(1:4, c(9,9,6,1)), ncol=1))
    if (plot.type=="naive") layout(matrix(rep(1:4, c(9,9,6,1)), ncol=1))
    if (plot.type=="both") layout(matrix(rep(1:6, c(9,9,6,1,6,1)), ncol=1))
    par(mar=c(0.25,3,0.25,1), mgp=c(1.75, 0.6, 0), oma=c(3,0,1.25,0))
    # raw data used for joint segmentation
    jseg <- x$jointseg
    # chromosome boundaries
    chrbdry <- which(diff(jseg$chrom) != 0)
    if (missing(emfit)) {
        out <- x$out
        if (plot.type=="em" | plot.type=="both") {
            warning("emfit is missing; plot.type set to naive")
            plot.type <- "naive"
        }
    } else {
        out <- emfit$cncf
        # add the naive tcn, lcn and cf to out
        out$tcn <- x$out$tcn
        out$lcn <- x$out$lcn
        out$cf <- x$out$cf
    }
    # determine which of the cnlr.median & mafR to show
    if (clustered) {
        cnlr.median <- out$cnlr.median.clust
        mafR <- out$mafR.clust
        mafR[is.na(mafR)] <- out$mafR[is.na(mafR)]
    } else {
        cnlr.median <- out$cnlr.median
        mafR <- out$mafR
    }
    mafR <- abs(mafR)
    # chromosome colors
    chrcol <- 1+rep(out$chrom-2*floor(out$chrom/2), out$num.mark)
    nn <- cumsum(table(jseg$chrom[is.finite(jseg$cnlr)]))
    segbdry <- cumsum(c(0,out$num.mark))
    segstart <- segbdry[-length(segbdry)]
    segend <- segbdry[-1]
    # plot the logR data and segment medians

    cnlr.bottom <- quantile(jseg$cnlr[is.finite(jseg$cnlr)],0.01)
    cnlr.ceiling <-  quantile(jseg$cnlr[is.finite(jseg$cnlr)],0.99)

    plot(jseg$cnlr[is.finite(jseg$cnlr)], pch=".", cex=2, col = c("grey","lightblue","azure4","slateblue")[chrcol], ylab="log-ratio", xaxt="n",
         ylim=c(min(cnlr.bottom-0.1, min(cnlr.median)-0.5, -1), max(cnlr.ceiling+0.1,max(cnlr.median)+0.5,1)))

    abline(v=chrbdry, lwd=0.25)
    abline(h=median(jseg$cnlr, na.rm=TRUE), col="green2")
    abline(h = x$dipLogR, col = "magenta4")
    segments(segstart, cnlr.median, segend, cnlr.median, lwd=1.75, col=2)
    # plot the logOR data and mafR
    plot(jseg$valor[is.finite(jseg$cnlr)], pch=".", cex=2.5, col = c("grey","lightblue","azure4","slateblue")[chrcol], ylab="log-odds-ratio", ylim=c(-4,4), xaxt="n")
    abline(v=chrbdry, lwd=0.25)
    segments(segstart, sqrt(mafR), segend, sqrt(mafR), lwd=1.75, col=2)
    segments(segstart, -sqrt(mafR), segend, -sqrt(mafR), lwd=1.75, col=2)
    # naive copy number and cellular faction pieces
    cfpalette <- c(colorRampPalette(c("white", "steelblue"))(10),"bisque2")
    if (plot.type=="naive" | plot.type=="both") {
        # plot the estimated copy numbers and cf
        out$tcn[out$tcn > 10] <- 9 + log10(out$tcn[out$tcn > 10])
        ii <- which(out$lcn > 5)
        if (length(ii)>0) out$lcn[ii] <- 5 + log10(out$lcn[ii])
        plot(c(0,length(jseg$cnlr)), c(0,max(out$tcn)), type="n", ylab="copy number (nv)", xaxt="n")
        abline(v=chrbdry, lwd=0.25)
        segments(segstart, out$lcn, segend, out$lcn, lwd=1.75, col=2)
        segments(segstart, out$tcn, segend, out$tcn, lwd=1.75, col=1)
        # add the cf
        plot(c(0,length(jseg$cnlr)), 0:1, type="n", ylab="", xaxt="n", yaxt="n")
        mtext("cf-nv", side=2, at=0.5, line=0.3, las=2, cex=0.75)
        cfcol <- cfpalette[round(10*out$cf+0.501)]
        rect(segstart, 0, segend, 1, col=cfcol, border=NA)
    }
    # EM copy number and cellular faction pieces
    if (plot.type=="em" | plot.type=="both") {
        # plot the estimated copy numbers and cf
        out$tcn.em[out$tcn.em > 10] <- 9 + log10(out$tcn.em[out$tcn.em > 10])
        ii <- which(out$lcn.em > 5)
        if (length(ii)>0) out$lcn.em[ii] <- 5 + log10(out$lcn.em[ii])
        plot(c(0,length(jseg$cnlr)), c(0,max(out$tcn.em)), type="n", ylab="copy number (em)", xaxt="n")
        abline(v=chrbdry, lwd=0.25)
        segments(segstart, out$lcn.em, segend, out$lcn.em, lwd=1.75, col=2)
        segments(segstart, out$tcn.em, segend, out$tcn.em, lwd=1.75, col=1)
        # add the cf
        plot(c(0,length(jseg$cnlr)), 0:1, type="n", ylab="", xaxt="n", yaxt="n")
        mtext("cf-em", side=2, at=0.5, line=0.2, las=2, cex=0.75)
        cfcol <- cfpalette[round(10*out$cf.em+0.501)]
        rect(segstart, 0, segend, 1, col=cfcol, border=NA)
    }

    # now add the chromosome ticks on x-axis
    chromlevels <- x$chromlevels
    # just make sure chromlevels actually exists
    if (is.null(chromlevels)) chromlevels <- 1:length(nn)
    axis(labels=chromlevels, side=1, at=(nn+c(0,nn[-length(nn)]))/2, cex=0.65)
    mtext(side=1, line=1.75, "Chromosome", cex=0.8)
    if (!missing(sname)) mtext(sname, side=3, line=0, outer=TRUE, cex=0.8)
    par(def.par)  #- reset to default
}