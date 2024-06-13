#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def main(args):
    # read in files
    df_cnv = pd.read_table(args.cnv_table)
    df_snps = pd.read_csv(args.jointseg)
    df_gene_coord = pd.read_table(args.gene_coord_file)
    sample_id = args.cnv_table.split('.')[0]
    pdf_out = f'{sample_id}.chrom_plots.pdf'
    df_cnv['chr_string'] = 'chr' + df_cnv.chrom.astype(str)
    df_cnv['seg_string'] = 'seg' + df_cnv.seg.astype(str)
    df_snps['chr_string'] = 'chr' + df_snps.chrom.astype(str)
    df_cnv_seg = df_cnv[['chr_string', 'seg_string', 'chrom', 'cnlr.median', 'mafR', 'loc.start', 'loc.end', 'ploidy',
                         'purity', 'tcn.em', 'lcn.em']]
    # get baf for folded snp plot
    get_baf_func = np.vectorize(get_baf, otypes=[float])
    df_snps['baf'] = get_baf_func(df_snps['vafT'])
    # filter to keep heterozygous snps
    df_het_snps = df_snps[df_snps['het'] == 1]
    grouped_snps = df_snps[['chr_string', 'chrom', 'maploc', 'baf', 'cnlr']].groupby('chrom')
    segments = []
    with PdfPages(pdf_out) as pdf:
        # iterate through chromosomes
        for name, group in grouped_snps:
            name = 'chr' + str(name)
            outfile = f'{sample_id}.{name}.plot_w_segments.png'
            df_cnv_chrom = df_cnv_seg[df_cnv_seg['chr_string'] == name].reset_index(drop=True)
            df_het_snps_chrom = df_het_snps[df_het_snps['chr_string'] == name].reset_index(drop=True)
            # start and end of segmented region of the chromosome
            min_pos = df_cnv_chrom['loc.start'].iloc[0]
            max_pos = df_cnv_chrom['loc.end'].iloc[-1]
            lr_y_min = min(group['cnlr']) - 0.5 if min(group['cnlr']) < -2 else -2
            lr_y_max = max(group['cnlr']) + 0.5 if max(group['cnlr']) > 2 else 2
            df_chr_gene = df_gene_coord[df_gene_coord['chrom'] == name].reset_index(drop=True)
            # 2 subplots: lr top, baf bottom
            fig, axs = plt.subplots(2, sharex=True)
            sns.scatterplot(data=group, x='maploc', y='cnlr', s=8, ax=axs[0])
            sns.scatterplot(data=df_het_snps_chrom, x='maploc', y='baf', s=8, ax=axs[1])
            axs[0].set_xlim(min_pos, max_pos)
            axs[0].set_ylim(lr_y_min, lr_y_max)
            axs[1].set_xlim(min_pos, max_pos)
            axs[1].set_ylim(0, 0.5)
            # iterate through each segment
            for i, row in df_cnv_chrom.iterrows():
                start = row['loc.start']
                end = row['loc.end']
                # get heterozygous snps located within segment
                df_het_snps_seg = df_het_snps_chrom.query(f'maploc >= {start} & maploc <= {end}')
                # in some cases there will be no heterzygous snps on a segment
                if df_het_snps_seg.empty:
                    seg_median_baf = None
                # use median baf to set location of segment line
                else:
                    seg_median_baf = df_het_snps_seg['baf'].median()
                seg_genes_snipdx = []
                seg_genes_other = []
                # add segments to the plot
                axs[0].hlines(y=row['cnlr.median'], xmin=row['loc.start'], xmax=row['loc.end'], color='red')
                if seg_median_baf:
                    axs[1].hlines(y=seg_median_baf, xmin=row['loc.start'], xmax=row['loc.end'], color='red')
                seg_length = row['loc.end'] - row['loc.start'] + 1
                # find any snipdx genes or other genes of interest on the segment
                if df_chr_gene.empty:
                    seg_genes_snipdx_str = ''
                    seg_genes_other_str = ''
                else:
                    for j, row2 in df_chr_gene.iterrows():
                        if ((row['loc.start'] < row2['start'] < row['loc.end'] or
                             row['loc.start'] < row2['end'] < row['loc.end']) or
                                (row2['start'] >= row['loc.start'] and row2['end'] <= row['loc.end'])):
                            if row2['on_snipdx'] == 'yes':
                                seg_genes_snipdx.append(row2['gene'])
                            else:
                                seg_genes_other.append(row2['gene'])
                    seg_genes_snipdx_str = ','.join(seg_genes_snipdx)
                    seg_genes_other_str = ','.join(seg_genes_other)
                # save segment data for output dataframe
                segments.append([row['chr_string'], row['loc.start'], row['loc.end'], row['cnlr.median'], row['mafR'],
                                 seg_median_baf, seg_length, round(seg_length / 1000000, 1), round(row['ploidy'], 2),
                                 round(row['purity'], 2), row['tcn.em'], row['lcn.em'], seg_genes_snipdx_str,
                                 seg_genes_other_str, int(row['chrom'])])
            # add vertical lines to plot for each gene located on chromosome
            if len(df_chr_gene) > 0:
                for k, row in df_chr_gene.iterrows():
                    # chr17 has the most genes, this is the number of colors to account for that, if additional genes
                    # added will need to add colors
                    colors = ['tab:green', 'tab:purple', 'tab:orange', 'tab:olive', 'tab:blue', 'tab:red', 'tab:brown']
                    axs[0].vlines(x=row['start'], ymin=lr_y_min, ymax=lr_y_max, color=colors[k], label=row['gene'])
                    axs[1].vlines(x=row['start'], ymin=0, ymax=0.5, color=colors[k])
                # format and place legend with gene names/colors
                axs[0].legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=len(df_chr_gene), columnspacing=0.75,
                              handletextpad=0.25)
            chrom_length_for_plot = (max_pos - min_pos)
            # set data ticks on x-axis for genomic positions at even intervals of segmented region
            tick_data_points = [min_pos, (chrom_length_for_plot * 0.25) + min_pos,
                                (chrom_length_for_plot * 0.5) + min_pos, (chrom_length_for_plot * 0.75) + min_pos,
                                max_pos]
            # use Mb
            tick_data_points_mb = [round(float(x) / 1000000, 2) for x in tick_data_points]
            axs[0].xaxis.set_visible(False)
            axs[1].set_xticks(tick_data_points)
            axs[1].ticklabel_format(axis='x', style='plain')
            axs[1].set_xticklabels([str(int(x)) for x in tick_data_points_mb], rotation=0, fontsize=6, ha='center')
            axs[1].set_xlabel('genomic position (Mb)')
            fig.suptitle(name)
            fig.get_figure()
            fig.savefig(outfile)
            pdf.savefig(fig)
            plt.close()
        # create output dataframe, sort by chromosome + position, write to .csv file
        df_seg_output = pd.DataFrame(segments, columns=['chrom', 'start', 'end', 'median_lr', 'mafR', 'median_baf',
                                                        'seg_length', 'seg_length_Mb', 'ploidy', 'purity', 'total_cn',
                                                        'minor_cn', 'snipdx_genes_on_segment', 'other_genes_on_segment',
                                                        'chrom_num'])
        df_seg_output = df_seg_output.sort_values(by=['chrom_num', 'start']).drop('chrom_num', axis=1)
        outfile_seg = f'{sample_id}.segments.csv'
        df_seg_output.to_csv(outfile_seg, index=False)


def get_baf(vaf):
    """ get baf for folded snp plot, range of 0-0.5 """
    if vaf < 0.5:
        return vaf
    else:
        return 1 - vaf


def parse_args():
    parser = argparse.ArgumentParser(description='create chromosome CN plot from facets data')
    parser.add_argument('cnv_table', help='cnv.table.txt file output by facets')
    parser.add_argument('jointseg', help='jointseg.csv file output by facets')
    parser.add_argument('gene_coord_file', help='tab-delimited file with snipdx gene coordinates')
    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
