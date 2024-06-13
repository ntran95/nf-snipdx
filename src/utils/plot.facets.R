plot.facets <- function(facets.fn = NULL,
                        sample.id = 'sample',
                        plotFolder = '.',
                        facets.profile = NULL,
                        ylim2 = NULL,
                        ylim3 = NULL,
                        ylab2 = 'log-odds-ratio',
                        format = 'pdf',
                        onlyAutosomal = FALSE,
                        sample.label = NULL,
                        majorCol = 'tcn.em') {

    # needs the following columns
    # loc.start
    # loc.end
    # chrom
    # cnlr.median
    # mafR
    # tcn.em
    # lcn.em

    if (is.null(sample.label)) {
      sample.label <- sample.id
    }

    # prepare chromosome lengths
    chromSizes <- vector()
    chromSizes['chr1'] <- 249250621
    chromSizes['chr2'] <- 243199373
    chromSizes['chr3'] <- 198022430
    chromSizes['chr4'] <- 191154276
    chromSizes['chr5'] <- 180915260
    chromSizes['chr6'] <- 171115067
    chromSizes['chr7'] <- 159138663
    chromSizes['chr8'] <- 146364022
    chromSizes['chr9'] <- 141213431
    chromSizes['chr10'] <- 135534747
    chromSizes['chr11'] <- 135006516
    chromSizes['chr12'] <- 133851895
    chromSizes['chr13'] <- 115169878
    chromSizes['chr14'] <- 107349540
    chromSizes['chr15'] <- 102531392
    chromSizes['chr16'] <- 90354753
    chromSizes['chr17'] <- 81195210
    chromSizes['chr18'] <- 78077248
    chromSizes['chr20'] <- 63025520
    chromSizes['chr19'] <- 59128983
    chromSizes['chr22'] <- 51304566
    chromSizes['chr21'] <- 48129895

    if (!onlyAutosomal) {
      chromSizes['chr23'] <- 155270560 # X
      chromSizes['chr24'] <- 59373566 # Y
    }

    cumSizes = cumsum(c(0, chromSizes[1:length(chromSizes) - 1]))
    names(cumSizes) <- names(chromSizes)
    midPoints <- rowMeans(cbind(cumSizes, cumsum(chromSizes)))
    names(midPoints) <- names(chromSizes)

    # load Facets profile
    if (!is.null(facets.fn)) {
      facets.profile <- read.table(facets.fn, sep = '\t',
                                   header = TRUE, stringsAsFactors = FALSE)
    }

    facets.profile$loc.start.abs <-
      facets.profile$loc.start + cumSizes[paste0('chr', facets.profile$chrom)]
    facets.profile$loc.end.abs <-
      facets.profile$loc.end + cumSizes[paste0('chr', facets.profile$chrom)]

    if (onlyAutosomal) {
      facets.profile <- subset(facets.profile, chrom %in% as.character(1:22))
    }

    # absolute x coordinates
    facets.profile$loc.start.abs <-
      facets.profile$loc.start + cumSizes[paste0('chr', facets.profile$chrom)]
    facets.profile$loc.end.abs <-
      facets.profile$loc.end + cumSizes[paste0('chr', facets.profile$chrom)]

    if (format == 'pdf') {
      plot_file <- file.path(plotFolder, paste0(sample.id, '.pdf'))
      pdf(plot_file, width = 20, heigh = 8)
    } else if (format == 'png') {
      plot_file <- file.path(plotFolder, paste0(sample.id, '.png'))
      png(plot_file, width = 2000, height = 1300, res = 200)
    }

    message("... creating: ", plot_file)

    #options(repr.plot.width=20, repr.plot.height=8)

    #par(mfrow = c(3, 1),     # 2x1 layout
    #oma = c(0, 1, 0, 0), # two rows of text at the outer left and bottom margin
    #mar = c(2, 4, 3, 0), # space for one row of text at ticks and to separate plots
    #    mar = c(0.5,5,5,0.5) # space for one row of text at ticks and to separate plots
    #mgp = c(2, 1, 0)    # axis label at 2 rows distance, tick labels at 1 row
    #)

    #par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)

    par(mar = c(0.5, 5, 5, 0.5), mfrow = c(3, 1), cex = 0.4, cex.main = 3, cex.axis = 2)
    par(cex.lab = 2)
    par(cex.main = 2)

    if (is.null(ylim3)) {
      ylim3 <- c(-0.5, max(c(5, max(facets.profile$tcn.em))))
    }

    plot(NULL, ylab = 'copy number', xlab = '',
         xlim = c(0, sum(chromSizes)), ylim = ylim3,
         xaxt = 'n', main = sample.label)

    segments(facets.profile$loc.start.abs,
             facets.profile[, majorCol] + 0.05,
             facets.profile$loc.end.abs,
             facets.profile[, majorCol] + 0.05,
             col = 'darkred', lwd = 3)

    segments(facets.profile$loc.start.abs,
             facets.profile$lcn.em,
             facets.profile$loc.end.abs,
             facets.profile$lcn.em,
             col = 'darkgreen', lwd = 3)

    abline(v = 0, col = 'gray')
    for (c in names(chromSizes)) {
      abline(v = chromSizes[c] + cumSizes[c], col = 'gray')
    }

    text(midPoints, -0.8, substr(names(midPoints), 4, 100), pos = 3, cex = 2)

    #browser()

    plot(NULL, ylab = 'log-ratio', xlab = '',
      xlim = c(0, sum(chromSizes)), ylim = c(-2, 2),
      xaxt = 'n', main = '')

    segments(facets.profile$loc.start.abs,
             facets.profile$cnlr.median,
             facets.profile$loc.end.abs,
             facets.profile$cnlr.median,
             col = 'red', lwd = 3)

    abline(v = 0, col = 'gray')
    for (c in names(chromSizes)) {
      abline(v = c(0, chromSizes[c] + cumSizes[c]), col = 'gray')
    }

    #par(
    #    mar = c(0.5,5,5,0.5)
    #    #mar = c(2, 4, 0, 0)
    #)
    #
    if (is.null(ylim2)) {
      ylim2 <- c(-2, 2)
    }

    plot(NULL, ylab = ylab2, xlab = '', xlim = c(0, sum(chromSizes)), ylim = ylim2, xaxt = 'n')
    segments(facets.profile$loc.start.abs,
             facets.profile$mafR,
             facets.profile$loc.end.abs,
             facets.profile$mafR,
             col = 'red', lwd = 3)

    abline(v = 0, col = 'gray')
    for (c in names(chromSizes)) {
      abline(v = chromSizes[c] + cumSizes[c], col = 'gray')
    }

    if (format == 'pdf' | format == 'png') {
      dev.off()
    }
}
