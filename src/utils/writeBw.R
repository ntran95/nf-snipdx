writeBw <- function(bw.fn, opt, xx) {
for (chri in 1:22) {
    if (chri==1) {
        cat(paste0('track type=wiggle_0 name="norm.',opt$sample_id,'LogR" description="normalized LogR" visibility=full autoScale=on viewLimits=-1:1 color=50,150,255 yLineMark=0 yLineOnOff=on priority=10'),file=bw.fn,sep="\n")
        cat(paste0('variableStep chrom=chr',chri , ' span=1'),file=bw.fn,sep="\n", append=TRUE)
    } else {
        cat(paste0('track type=wiggle_0 name="',opt$sample_id,'LogR" description="normalized LogR" visibility=full autoScale=on viewLimits=-1:1 color=50,150,255 yLineMark=0 yLineOnOff=on priority=10'),file=bw.fn,sep="\n", append=TRUE)
        cat(paste0('variableStep chrom=chr',chri , ' span=1'),file=bw.fn,sep="\n", append=TRUE)
    }

    table.chr1 <- subset(xx$jointseg ,chrom==as.character(chri))

    table.chr1 <- subset(table.chr1, !is.infinite(table.chr1$cnlr) & !is.na(table.chr1$cnlr) )
    #table.chr1 <- subset(table.chr1, !is.infinite(table.chr1$logR) & !is.na(table.chr1$logR) )
    bw.table <- table.chr1[,c('maploc', 'cnlr')]

    write.table(bw.table, file=bw.fn, sep=' ', row.names=FALSE, col.names=FALSE, append=TRUE)
}


}