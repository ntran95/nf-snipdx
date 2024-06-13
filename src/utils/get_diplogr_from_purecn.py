#!/usr/bin/env python3

import argparse
import pandas as pd


def main(args):
    outfile = f'{args.tumor_sample_name}_purecn_diplogr.txt'
    df_seg = pd.read_csv(args.seg_file)
    df_seg_diploid = df_seg[df_seg['C'] == 2]
    mean_diplogr = df_seg_diploid['seg.mean'].mean()
    with open(outfile, 'w') as fh_out:
        fh_out.write(str(mean_diplogr) + '\n')


def parse_args():
    parser = argparse.ArgumentParser(description='get average diplogr from purecn segment file')
    parser.add_argument('seg_file', help='{sample}_loh.csv file output from purecn')
    parser.add_argument('tumor_sample_name', help='tumor sample name')
    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
