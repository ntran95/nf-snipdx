import argparse
import statistics
import pandas as pd
import pysam
import os
import sys

if __name__ == "__main__":

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                 usage=__doc__)
    ap.add_argument('--targets_csv', required=True,
                    help='target SNP coordinates')
    
    ap.add_argument('--bam', required=True,
                    help='bam file path')
    
    ap.add_argument('--output_fp', required=True,
                    help='output file path')
    
    args = ap.parse_args()
    
    
    targets=pd.read_csv(args.targets_csv)
    
    samfile = pysam.AlignmentFile(args.bam, "rb")
    for i, row in targets.iterrows():
        tlen_list=[]
        for read in samfile.fetch(row[0], row[1]-1, row[1]):
            if (abs(read.tlen)<1e4):
                tlen_list.append(abs(read.tlen))
        if len(tlen_list)>0:
            targets.at[i,'is_mean'] = statistics.mean(tlen_list)
    samfile.close()

    
    targets.to_csv(args.output_fp, index=False)
    
