from argparse import ArgumentParser
import gzip
import numpy as np
import os
import pandas as pd
from scipy import stats


def get_chromsize_dict(chromsize_path):
    chromsize_dict = {}
    with open(chromsize_path, "r") as chrom_link:
        for chrom_line in chrom_link:
            chrom, chromlen = chrom_line.rstrip().split("\t")
            chromsize_dict[chrom] = int(chromlen)
    return chromsize_dict


def get_uint_ar(uint_path, kmer):
    quant_ar = np.zeros(1000)
    if os.path.exists(uint_path):
        with gzip.open(uint_path, "rb") as uint_link:
            uint_ar = np.frombuffer(uint_link.read(), dtype=np.uint8)
        less_than_kmer = uint_ar <= kmer
        not_zero = uint_ar != 0
        uniquely_mappable = np.array(
            np.logical_and(less_than_kmer, not_zero), dtype=int)
        LEN_AR = len(uniquely_mappable)
        poses_start, = np.where(uniquely_mappable == 1)
        ar_quant = np.zeros(len(uniquely_mappable), dtype=float)
        for ind_st in poses_start:
            ind_end = ind_st + kmer
            if ind_end > LEN_AR:
                ind_end = LEN_AR
            ar_quant[ind_st:ind_end] = ar_quant[ind_st:ind_end] + 1.0
        print("Created multi-read mappability array for k{}".format(kmer))
    return quant_ar


def get_bed_ar(bed_path, chrom, chromsize):
    skip_rows = 0
    bed_ar = np.zeros(chromsize, dtype=float)
    if os.path.exists(bed_path):
        with gzip.open(bed_path, 'rb') as bed_link:
            header_bed = bed_link.readline()
            if len(header_bed.rstrip().split('\t')) < 2:
                skip_rows = 1
                print("Assuming header exists in bed")
        bed_df = pd.read_csv(bed_path, sep="\t",
                             compression="gzip", header=None,
                             skip_rows=skip_rows)
        bed_df = bed_df[bed_df.iloc[:, 0] == chrom, ]
        bed_ar = np.zeros(chromsize, dtype=float)
        for i in range(bed_df.shape[0]):
            bed_ar[bed_df.iloc[i, 1] - 1:bed_df.iloc[i, 2]] = bed_df.iloc[i, 3]
        print("Successfully made uniquely mappable bed file array")
    return bed_ar


def get_bedg_ar(bedgraph_path, chrom, chromsize):
    skip_rows = 0
    bedg_ar = np.zeros(chromsize, dtype=float)
    if os.path.exists(bedgraph_path):
        with gzip.open(bedgraph_path, 'rb') as bedg_link:
            header_bedg = bedg_link.readline()
            if len(header_bedg.rstrip().split('\t')) < 3:
                skip_rows = 1
                print("Assuming header exists in bedGraph")
        bedg_df = pd.read_csv(bedgraph_path, sep="\t",
                              compression="gzip", header=None,
                              skip_rows=skip_rows)
        bedg_df = bedg_df[bedg_df.iloc[:, 0] == chrom, ]
        bedg_ar = np.zeros(chromsize, dtype=float)
        idx_ones, = np.where(bedg_df.iloc[:, 2] - bedg_df.iloc[:, 1] == 1)
        idx_others, = np.where(bedg_df.iloc[:, 2] - bedg_df.iloc[:, 1] > 1)
        bedg_ar[bedg_df.iloc[idx_ones, 1] - 1] = bedg_df.iloc[idx_ones, 3]
        for each_idx in idx_others:
            st = bedg_df.iloc[each_idx, 1]
            end = bedg_df.iloc[each_idx, 2]
            val = bedg_df.iloc[each_idx, 3]
            bedg_ar[st - 1:end] = val
        print("Successfully converted bedGraph to numpy array")
    return bedg_ar


def get_wig_ar(wig_path, chrom, chromsize):
    wig_ar = np.zeros(chromsize, dtype=float)
    start = 1
    if os.path.exists(wig_path):
        with gzip.open(wig_path, "rb") as wig_link:
            for wig_line in wig_link:
                if 'fixedStep' in wig_line:
                    start_str = wig_line.rstrip().split(" ")[-2]
                    new_start = int(start_str.replace("start=", ""))
                    start = new_start
                else:
                    ad_val = float(wig_line.rstrip())
                    wig_ar[(start - 1)] = ad_val
                    start = start + 1
    return wig_ar


def compare_bedg_uint(bedg_ar, uint_ar, out_link, chrom):
    if sum(bedg_ar) == 0 or sum(uint_ar) == 0:
        print("Not comparing bedGraph and uint")
        unique_both = 0
        total_uint = len(uint_ar)
        prsn_r = "NA"
    else:
        unique_both = len(
            np.where(np.logical_and(
                bedg_ar > 0, uint_ar > 0))[0])
        total_uint = len(np.where(uint_ar > 0)[0])
        prsn_r, pval = stats.pearsonr(bedg_ar, uint_ar)
    ad_list = ["bedGraph", "uint", str(prsn_r),
               str(total_uint/unique_both),
               str(prsn_r), chrom]
    out_link.write("\t".join(ad_list) + "\n")


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Checks if output files of umap "
        "are consistent with each other.")
    parser.add_argument(
        "out_path",
        help="Path to output file for saving results")
    parser.add_argument(
        "kmer",
        type=int,
        help="Kmer size (integer)")
    parser.add_argument(
        "chromsize_path",
        help="Path to chromosome size file")
    parser.add_argument(
        "chrom",
        help="Name of the chromosome")
    parser.add_argument(
        "--bedGraph_path",
        default="NA",
        help="Optional -- path to gzipped bedGraph")
    parser.add_argument(
        "--uint_path",
        default="NA",
        help="OPtional -- path to gzipped unique.uint8 file")
    parser.add_argument(
        "--bed_path",
        default="NA",
        help="Optional -- path to gzipped bed file of uniquely "
        "mappable regions")
    parser.add_argument(
        "--wig_path",
        default="NA",
        help="Optional -- path to gzipped wiggle file")
    args = parser.parse_args()
    chromsize_dict = get_chromsize_dict(args.chromsize_path)
    chromsize = chromsize_dict[args.chrom]
    out_link = gzip.open(args.out_path, "wb")
    out_link.write("File 1\tFile 2\tPearson\tOverlap\tChromosome\n")
    bedg_ar = get_bedg_ar(args.bedGraph_path, args.chrom, chromsize)
    uint_ar = get_uint_ar(args.uint_path, args.kmer)
    bed_ar = get_bed_ar(args.bed_path, args.chrom, chromsize)
    wig_ar = get_wig_ar(args.wig_path, args.chrom, chromsize)
    compare_bedg_uint(bedg_ar, uint_ar, out_link, args.chrom)
    compare_bedg_bed(bedg_ar, bed_ar, out_link, args.chrom)
    compare_bedg_wig(bedg_ar, args.wig_path, out_link, args.chrom)
    compare_wig_uint(args.wig_path, uinr_ar, out_link, args.chrom)
    compare_wig_bed(args.wig_path, bed_ar, out_link, args.chrom)
