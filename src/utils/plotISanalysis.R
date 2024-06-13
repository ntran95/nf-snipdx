  
plotISanalysis <- function(xx, pon.df, sample.is.fn, is.plot.fn) {

    jointseg.annot <- xx$jointseg
    jointseg.annot$sampleloci <- paste0('chr',jointseg.annot$chrom,' ', jointseg.annot$maploc)
    
    jointseg.annot$pon.is <- pon.df[jointseg.annot$sampleloci, paste0('is.',class.used)]

    jointseg.annot$snipdx.gene <- pon.df[jointseg.annot$sampleloci, 'snipdx.gene']
    jointseg.annot$target.type <- pon.df[jointseg.annot$sampleloci, 'target.type']
    jointseg.annot$closestPrimer <- pon.df[jointseg.annot$sampleloci, 'closestPrimer']


    sample.is <- read.csv(sample.is.fn)
    rownames(sample.is) <- paste(sample.is$Chromosome, sample.is$Position)
    jointseg.annot$sample.is <- sample.is[jointseg.annot$sampleloci ,4]


    jointseg.annot$pch=1
    jointseg.annot$pch[jointseg.annot$target.type=='exon']=19
    jointseg.annot$col='darkgray'
    jointseg.annot$col[jointseg.annot$snipdx.gene=='RNASEH2B']='green'
    jointseg.annot$col[jointseg.annot$snipdx.gene=='SETD2']='red'

    png(is.plot.fn, width=3000, height=2000, res = 100)
    par(mfrow=c(2,3))
    par(mar = c(5, 5, 5, 5))
    
    plot(jointseg.annot$rCountN, jointseg.annot$rCountT, col=jointseg.annot$col, pch=jointseg.annot$pch, main='Coverage', 
           xlab=paste0('PON ',class.used, ' coverage (X)'), ylab='Sample coverage (X)', cex.lab=2, cex.axis=1.5)
    points(subset(jointseg.annot, col!='darkgray')$rCountN, subset(jointseg.annot, col!='darkgray')$rCountT, col=subset(jointseg.annot, col!='darkgray')$col, pch=subset(jointseg.annot, col!='darkgray')$pch)
    lines(c(0,100000), c(0,100000), col='gray')

    plot(jointseg.annot$pon.is, jointseg.annot$sample.is, col=jointseg.annot$col, pch=jointseg.annot$pch, xlim=c(100,300), ylim=c(100,300), 
           main='Insert size length', xlab=paste0('PON ',class.used, ' insert size (bp)'), ylab='Sample insert size (bp)', cex.lab=2, cex.axis=1.5)
    points(subset(jointseg.annot, col!='darkgray')$pon.is, subset(jointseg.annot, col!='darkgray')$sample.is, 
           col=subset(jointseg.annot, col!='darkgray')$col, pch=subset(jointseg.annot, col!='darkgray')$pch)
    lines(c(0,100000), c(0,100000), col='gray')

    plot(jointseg.annot$gcpct,  jointseg.annot$rCountT/jointseg.annot$rCountN, col=jointseg.annot$col, pch=jointseg.annot$pch,
        main='GC bias', xlab='GC content', ylab='Coverage Tumor / Normal', cex.lab=2, cex.axis=1.5)
    points(subset(jointseg.annot, col!='darkgray')$gcpct,  subset(jointseg.annot, col!='darkgray')$rCountT/subset(jointseg.annot, col!='darkgray')$rCountN, , 
           col=subset(jointseg.annot, col!='darkgray')$col, pch=subset(jointseg.annot, col!='darkgray')$pch)
    lines(c(0,100000), c(0,100000), col='gray')



    plot(jointseg.annot$sample.is,  jointseg.annot$cnlr, col=jointseg.annot$col, pch=jointseg.annot$pch, xlim=c(0,300),
          xlab='Sample insert size (bp)', ylab='Sample LogR', cex.lab=2, cex.axis=1.5)
    points(subset(jointseg.annot, col!='darkgray')$sample.is,  
           subset(jointseg.annot, col!='darkgray')$cnlr, , col=subset(jointseg.annot, col!='darkgray')$col, pch=subset(jointseg.annot, col!='darkgray')$pch)
    lines(c(0,100000), c(0,100000), col='gray')
        out.lowess.tiled <- lowess(x=subset(jointseg.annot, target.type=='SNP')$sample.is, 
                                 y=subset(jointseg.annot, target.type=='SNP')$cnlr, 
                                 f=0.2,  delta=0)
        points(out.lowess.tiled, col='orange', pch=19, cex=0.2)

        out.lowess.tiled <- lowess(x=subset(jointseg.annot, target.type!='SNP')$sample.is, 
                                 y=subset(jointseg.annot, target.type!='SNP')$cnlr, 
                                 f=0.2,  delta=0)
        points(out.lowess.tiled, col='purple', pch=19, cex=0.2)





    plot(jointseg.annot$pon.is,  jointseg.annot$cnlr, , col=jointseg.annot$col, pch=jointseg.annot$pch, xlim=c(0,300),
        xlab=paste0('PON ',class.used, ' insert size (bp)'), ylab='Sample LogR', cex.lab=2, cex.axis=1.5)
    points(subset(jointseg.annot, col!='darkgray')$pon.is,  
           subset(jointseg.annot, col!='darkgray')$cnlr, , col=subset(jointseg.annot, col!='darkgray')$col, pch=subset(jointseg.annot, col!='darkgray')$pch)
    lines(c(0,100000), c(0,100000), col='gray')
        out.lowess.tiled <- lowess(x=subset(jointseg.annot, target.type=='SNP')$pon.is, 
                                 y=subset(jointseg.annot, target.type=='SNP')$cnlr, 
                                 f=0.2,  delta=0)
        points(out.lowess.tiled, col='orange', pch=19, cex=0.2)
        out.lowess.tiled <- lowess(x=subset(jointseg.annot, target.type!='SNP')$pon.is, 
                                 y=subset(jointseg.annot, target.type!='SNP')$cnlr, 
                                 f=0.2,  delta=0)
        points(out.lowess.tiled, col='purple', pch=19, cex=0.2)




    plot.xlims <- quantile(jointseg.annot$sample.is/jointseg.annot$pon.is, c(0.001, 0.999), na.rm=TRUE) + c(-0.1, 0.1)
    plot.ylims <- quantile(jointseg.annot$cnlr, c(0.001, 0.999), na.rm=TRUE) + c(-0.1, 0.1)
     
    plot(jointseg.annot$sample.is/jointseg.annot$pon.is,  jointseg.annot$cnlr, , col=jointseg.annot$col, pch=jointseg.annot$pch, xlim= plot.xlims, ylim=plot.ylims,
        xlab=paste0('Insert size bias: Sample / PON'), ylab='Sample LogR', cex.lab=2, cex.axis=1.5)
    points(subset(jointseg.annot, col!='darkgray')$sample.is/subset(jointseg.annot, col!='darkgray')$pon.is,  
           subset(jointseg.annot, col!='darkgray')$cnlr, , col=subset(jointseg.annot, col!='darkgray')$col, pch=subset(jointseg.annot, col!='darkgray')$pch)
    lines(c(0,100000), c(0,100000), col='gray')

        out.lowess.tiled <- lowess(x=subset(jointseg.annot, target.type=='SNP')$sample.is/subset(jointseg.annot, target.type=='SNP')$pon.is, 
                                 y=subset(jointseg.annot, target.type=='SNP')$cnlr, 
                                 f=0.2,  delta=0)
        points(out.lowess.tiled, col='orange', pch=19, cex=0.2)
        out.lowess.tiled <- lowess(x=subset(jointseg.annot, target.type!='SNP')$sample.is/subset(jointseg.annot, target.type!='SNP')$pon.is, 
                                 y=subset(jointseg.annot, target.type!='SNP')$cnlr, 
                                 f=0.2,  delta=0)
        points(out.lowess.tiled, col='purple', pch=19, cex=0.2)

    dev.off()

#    png(paste0(plotFolderMain, opt$sample_id, '.is.barplots.png'), width=2000, height=2000, res = 100)
#    par(mfrow=c(1,2))
    
#    boxplot(sample.is ~ target.type, data=jointseg.annot)
#        boxplot(sample.is ~ target.type + snipdx.gene, data=jointseg.annot, horizontal=TRUE, las=2)
#    dev.off()
}