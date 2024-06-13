extract.ei.bias <- function (ei.bias.plot.fn, ei.bias.table.fn, xx, opt) {
    xx$jointseg$snp.id <- paste0('chr', xx$jointseg$chrom, ' ', xx$jointseg$maploc)
    xx$jointseg$target.type <- pon.df[xx$jointseg$snp.id, 'target.type']
    xx$jointseg$snipdx.gene <- pon.df[xx$jointseg$snp.id, 'snipdx.gene']


    png(ei.bias.plot.fn, width=1400, height=4000, res = 200)
    par(mar=c(2,10,2,2))
    boxplot(cnlr ~ target.type + snipdx.gene, data=xx$jointseg, col=c('darkblue', 'orange'), ylim=c(-2,2), main=opt$sample_id, horizontal=TRUE, las=2, ylab='')
    abline(v=0)
    dev.off()

    am <- aggregate(xx$jointseg[, c('cnlr')], list(xx$jointseg$snipdx.gene, xx$jointseg$target.type), median ,na.rm=TRUE)
    sample.matrix <- dcast(am, Group.1 ~ Group.2, value.var='x')
    rownames(sample.matrix) <- sample.matrix$Group.1
    sample.matrix$Group.1 <- NULL
    write.csv(sample.matrix, file=ei.bias.table.fn, quote=FALSE)
    }

