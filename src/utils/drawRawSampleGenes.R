drawRawSampleGenes <- function(pl, sample_id,
                               plotFolderMain,
                               panel.genes,
                               facets.cncf,
                               oo = NULL,
                               margin = 5e6,
                               allChroms = TRUE) {

    # pl: the data frame of SNPs in Facets format
    # sample_id: ID of smaple, for naming files
    # plotFolderMain folder where plot should be saved
    # panel.genes: data frame of panel genes, with chromosome coordinates
    # facets.cncf: facets output - integer copy number (major and minor) for chromosomal regions (not just genes)
    # oo: Result of Facets segmentation
    # margin: how much margin to draw around genes. Recommend 1e5 to show gene in detail, 5e6 to show more of chromosomal context
    # allChroms: if TRUE, show one page with all genes. If FALSE, show I gene per page

    pl$af <- pl$File2A / (pl$File2A + pl$File2R)

    # count approximate number of heterozygotes
    nHets <- sum(pl$af > 0.05 & pl$af > 0.95, na.rm = TRUE)
    noSnpsTotal <- nrow(pl)
    noSnps100 <- sum((pl$File2A + pl$File2R) > 100, na.rm = TRUE)

    # calculate 'absolute' x axis, combining results from all chromosomes
    chroms <- paste0('chr', c(1:22, 'X'))
    chrMaxPos <- tapply(pl$Position, pl$Chromosome, max)[chroms]
    chrShift <- cumsum(c(0, chrMaxPos[1:(length(chrMaxPos) - 1)]))
    names(chrShift) <- names(chrMaxPos)
    pl$Position.cum <- pl$Position + chrShift[pl$Chromosome]

    # included SNPs are black.
    pl$col <- 'gray'
    pl$col[pl$include] <- 'black'

    # prepare plot files to write to
    # if all chromosome summary
    allChroms.str <- 'gene.per.page'
    if (allChroms) {
        plot_file <- file.path(plotFolderMain, paste0(sample_id, '.rawsummary.genes.', margin / 1e6, 'Mb.all.chroms.png'))
        png(plot_file, width = 3500, height = 1100)
        layout(matrix(1:(3 * 26), 3, 26, byrow = FALSE))
        par(mar = c(4, 5, 2, 0))
        # if gene per page summary
    } else {
        plot_file <- file.path(plotFolderMain, paste0(sample_id, '.rawsummary.genes.', margin / 1e6, 'Mb.gene.per.page.pdf'))
        pdf(plot_file, width = 12, height = 11)
        par(mfrow = c(3, 1))
    }

    message("... creating: ", plot_file)

    # mean LogR across all SNPs
    baselineLogR <- mean(subset(pl, include == TRUE)$logR, na.rm = TRUE)
    # get maximum copy number of a sample
    maxCn <- max(facets.cncf$tcn.em)

    # loop over all genes
    for (gi in 1:nrow(panel.genes)) {
        # get gene names and coordinates
        g_n <- panel.genes$gene.name[gi]
        g_chrom <- panel.genes$chromosome_name[gi]
        g_window_start <- panel.genes$transcript_start[gi] - margin
        g_window_end <- panel.genes$transcript_end[gi] + margin

        # get SNPs that overlap with the window to show
        pl_sel <- subset(pl, Chromosome == paste0('chr', g_chrom) &
                             Position > g_window_start &
                             Position < g_window_end)

        # plot the allele fraction
        if (gi == 1 || !allChroms) {
            plot(pl_sel$Position / 1e6,
                 pl_sel$af,
                 pch = 19,
                 main = g_n ,
                 ylim = c(0, 1.1),
                 cex = 1,
                 cex.lab = 2,
                 cex.axis = 2,
                 cex.main = 2,
                 xlab = '',
                 ylab = 'allele fraction',
                 col = pl_sel$col,
                 xlim = c(g_window_start / 1e6, g_window_end / 1e6),
                 xaxt = 'n')
        } else {
            plot(pl_sel$Position / 1e6,
                 pl_sel$af,
                 pch = 19,
                 xlab = '',
                 ylab = '',
                 main = g_n ,
                 ylim = c(0, 1.1),
                 cex = 1,
                 cex.lab = 2,
                 cex.axis = 2,
                 cex.main = 2,
                 col = pl_sel$col,
                 xlim = c(g_window_start / 1e6, g_window_end / 1e6),
                 xaxt = 'n',
                 yaxt = 'n')
        }

        # show boundaries of gene body
        abline(v = panel.genes$transcript_start[gi] / 1e6, col = 'gold1')
        abline(v = panel.genes$transcript_end[gi] / 1e6, col = 'gold1')

        pl.inc <- subset(pl_sel, include == TRUE)

        # Facets segmentation
        window.segments <- subset(facets.cncf,
                                  chrom == g_chrom &
                                      end > g_window_start &
                                      start < g_window_end)

        # plot read coverage in region
        # if Facets segmentation output is availble, use this as it is GC-corrected
        if (!is.null(oo)) {
            oo.jointseg.sel <- subset(oo$jointseg,
                                      chrom == g_chrom  &
                                          maploc > g_window_start &
                                          maploc < g_window_end)
            # for first chromosome, show all axis labels
            if (gi == 1 || !allChroms) {
                plot(oo.jointseg.sel$maploc / 1e6,
                     oo.jointseg.sel$cnlr,
                     pch = 19,
                     main = g_n ,
                     ylim = c(-3, 3),
                     cex = 1,
                     cex.lab = 2,
                     cex.axis = 2,
                     cex.main = 2,
                     ylab = 'LogR',
                     xlab = '',
                     col = 'darkgreen',
                     xlim = c(g_window_start / 1e6, g_window_end / 1e6))
            } else {
                # for further chromosomes, don't show y axis
                plot(oo.jointseg.sel$maploc / 1e6,
                     oo.jointseg.sel$cnlr,
                     pch = 19,
                     main = g_n ,
                     ylim = c(-3, 3),
                     cex = 1,
                     cex.lab = 2,
                     cex.axis = 2,
                     cex.main = 2,
                     xlab = '',
                     ylab = '',
                     col = 'darkgreen',
                     xlim = c(g_window_start / 1e6, g_window_end / 1e6),
                     yaxt = 'n',
                     xaxt = 'n')
            }

            # draw Faces LogR segments, on top, in red
            if (nrow(window.segments) > 0) {
                segments(window.segments$start / 1e6,
                         window.segments$cnlr.median,
                         window.segments$end / 1e6,
                         window.segments$cnlr.median,
                         col = 'red', lwd = 3)
            }
            abline(h = 0, col = 'gray')
        } else {
            # if Facets segmentation is not available, draw raw LogR data
            if (gi == 1 || !allChroms) {
                plot(pl.inc$Position / 1e6,
                     pl.inc$logR,
                     pch = 19,
                     main = g_n ,
                     ylim = c(-3, 3),
                     cex = 1,
                     cex.lab = 2,
                     cex.axis = 2,
                     cex.main = 2,
                     ylab = 'LogR',
                     xlab = '',
                     col = 'darkmagenta',
                     xlim = c(g_window_start / 1e6, g_window_end / 1e6),
                     xaxt = 'n')
            } else {
                plot(pl.inc$Position / 1e6,
                     pl.inc$logR,
                     pch = 19,
                     main = g_n ,
                     ylim = c(-3, 3),
                     cex = 1,
                     cex.lab = 2,
                     cex.axis = 2,
                     cex.main = 2,
                     xlab = '',
                     ylab = '',
                     col = 'darkmagenta',
                     xlim = c(g_window_start / 1e6, g_window_end / 1e6),
                     yaxt = 'n',
                     xaxt = 'n')
            }

            abline(h = mean(pl.inc$logR, na.rm = TRUE),
                   col = 'darkmagenta', lwd = 3)
            abline(h = baselineLogR, col = 'sandybrown', lwd = 3)
        }

        abline(v = panel.genes$transcript_start[gi] / 1e6, col = 'gold1')
        abline(v = panel.genes$transcript_end[gi] / 1e6, col = 'gold1')

        # draw region copy number, as estimated by Facets
        if (nrow(window.segments) > 0) {
            # on first gene, draw all axes
            if (gi == 1 || !allChroms) {
                plot(1,
                     type = "n",
                     ylab = "copy number",
                     main = g_n,
                     xlim = c(g_window_start / 1e6, g_window_end / 1e6),
                     ylim = c(0, maxCn + 1),
                     cex = 1,
                     cex.lab = 2,
                     cex.axis = 2,
                     cex.main = 2,
                     xlab = paste0('chr', g_chrom))

                segments(window.segments$start / 1e6,
                         window.segments$tcn.em,
                         window.segments$end / 1e6,
                         window.segments$tcn.em,
                         lwd = 5)

                segments(window.segments$start / 1e6,
                         window.segments$lcn.em,
                         window.segments$end / 1e6,
                         window.segments$lcn.em,
                         lwd = 5, col = 'darkred')
                abline(h = 0)

                # for further genes, don't show y axis labels
            } else {
                plot(1,
                     type = "n",
                     ylab = "",
                     main = g_n,
                     xlim = c(g_window_start / 1e6, g_window_end / 1e6),
                     ylim = c(0, maxCn + 1),
                     cex = 1,
                     cex.lab = 2,
                     cex.axis = 2,
                     cex.main = 2,
                     xlab = paste0('chr', g_chrom),
                     yaxt = 'n')

                segments(window.segments$start / 1e6,
                         window.segments$tcn.em,
                         window.segments$end / 1e6,
                         window.segments$tcn.em,
                         lwd = 5)

                segments(window.segments$start / 1e6,
                         window.segments$lcn.em,
                         window.segments$end / 1e6,
                         window.segments$lcn.em,
                         lwd = 5, col = 'darkred')
                abline(h = 0)
            }
            # make an empty plot if no segmentation is available
        } else {
            plot(1,
                 type = "n",
                 ylab = "",
                 main = g_n,
                 xlim = c(g_window_start / 1e6, g_window_end / 1e6),
                 ylim = c(0, maxCn + 1),
                 cex = 1,
                 cex.lab = 2,
                 cex.axis = 2,
                 cex.main = 2,
                 xlab = paste0('chr', g_chrom),
                 yaxt = 'n')
        } # end if

        # add gene boundaries
        abline(v = panel.genes$transcript_start[gi] / 1e6, col = 'gold1')
        abline(v = panel.genes$transcript_end[gi] / 1e6, col = 'gold1')

    } # end loop over genes

    dev.off()
} # end function
