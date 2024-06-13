normalize_by_is <- function(pl2, sample.is, pon.df, is_plot_fn) {
    
        pl2$r <- ((pl2$File2R + pl2$File2A)/ median(pl2$File2R + pl2$File2A, na.rm=TRUE)) / ((pl2$File1R + pl2$File1A) / median(pl2$File1R + pl2$File1A, na.rm=TRUE))
        pl2$is.sample.class <- pon.df[paste0('chr', pl2$Chromosome, ' ',pl2$Position ), paste0('is.',class.used)]
        pl2$target.type <- pon.df[paste0('chr', pl2$Chromosome, ' ',pl2$Position ), 'target.type']
        pl2$sample.is <- sample.is[paste0('chr', pl2$Chromosome, ' ',pl2$Position) ,4]
        pl2$is.ratio <- pl2$sample.is/pl2$is.sample.class

        
        pl2.snps <- subset(pl2,target.type=='SNP')
        pl2.snps.lowess <- lowess(x=pl2.snps$is.ratio, 
                         y=pl2.snps$r , 
                         f=0.1,  delta=0)
                 
        jj <- match(pl2.snps$is.ratio, pl2.snps.lowess$x)
        pl2.snps$r.lowess <-pl2.snps.lowess$y[jj]
        
        pl2.exons <- subset(pl2,target.type!='SNP')
        pl2.exons.lowess <- lowess(x=pl2.exons$is.ratio, 
                         y=pl2.exons$r , 
                         f=0.1,  delta=0)
        
        jj <- match(pl2.exons$is.ratio, pl2.exons.lowess$x)
        pl2.exons$r.lowess <-pl2.exons.lowess$y[jj]
        

        pl2 <- rbind(pl2.snps,pl2.exons )
        pl2 <- pl2[order(as.integer(as.character(pl2$Chromosome)), pl2$Position),]
        
        pl2$File2R.old <- pl2$File2R 
        pl2$File2A.old <-  pl2$File2A 
        pl2$r.old <-  pl2$r 
        
        png(is_plot_fn, width=1000, height=1000, res = 100)
        par(mfrow=c(2,2))
        plot(pl2$sample.is/pl2$is.sample.class, pl2$r,  pch=19)
            points(pl2.snps.lowess, col='orange', pch=19, cex=0.2)
            points(pl2.exons.lowess, col='purple', pch=19, cex=0.2)

        pl2$File2R <- round(pl2$File2R/pl2$r.lowess)
        pl2$File2A <- round(pl2$File2A/pl2$r.lowess)
        pl2$r <- ((pl2$File2R + pl2$File2A)/ median(pl2$File2R + pl2$File2A, na.rm=TRUE)) / ((pl2$File1R + pl2$File1A) / median(pl2$File1R + pl2$File1A, na.rm=TRUE))
        pl2.snps <- subset(pl2,target.type=='SNP')
        pl2.exons <- subset(pl2,target.type!='SNP')
        pl2.snps.lowess2 <- lowess(x=pl2.snps$is.ratio, 
                         y=pl2.snps$r , 
                         f=0.1,  delta=0)
        pl2.exons.lowess2 <- lowess(x=pl2.exons$is.ratio, 
                         y=pl2.exons$r , 
                         f=0.1,  delta=0)       
        
        
        plot(pl2$is.ratio, pl2$r,  pch=19)
            points(pl2.snps.lowess2, col='orange', pch=19, cex=0.2)
            points(pl2.exons.lowess2, col='purple', pch=19, cex=0.2)

        plot(pl2$File2R.old + pl2$File2A.old,  pl2$File2R + pl2$File2A ,  pch=19)
        plot(pl2$r.old,  pl2$r ,  pch=19)
        
        dev.off()  
    
    return(pl2)
}