from argparse import ArgumentParser
import gzip
import numpy as np
import os
import pandas as pd
import re


def get_args():
    parser = ArgumentParser(
        description="This script expects an input "
        "gzipped BED file, a directory with umap/bismap "
        "produced files names as <chr>.uint8."
        "unique.gz. It will add one column to the bed "
        "file which would determine the fraction of "
        "uniquely mappable reads that passover this "
        "region. The first column added to input is "
        "single-read mappability, which is the percent "
        "of the region that is uniquely mappable. "
        "The second added column is the multi-read"
        " mappability; the probability that a randonly "
        "selected read overlapping with that region "
        "is uniquely mappable. If the input BED file "
        "has 1 bp intervals, you can specify -SingleNucleotide "
        "and -wigdir to use a faster algorithm.")
    parser.add_argument(
        "bed_path",
        help="Path to gzipped bed file")
    parser.add_argument(
        "out_path",
        help="Path to gzipped output file")
    parser.add_argument(
        "umap_dir",
        help="Path to directory with uint8 binary "
        "files produced by Umap/Bismap.")
    parser.add_argument(
        "kmer",
        type=int,
        help="The read length for defining single-read "
        "and multi-read mappability.")
    parser.add_argument(
        "-SingleNucleotide",
        action="store_true",
        help="If specified, assumes each region is only one "
        "nucleotide. You must specify -wig as well.")
    parser.add_argument(
        "-wigdir",
        default="",
        help="Path to directory with <kmer>_<chrom>."
        "MultiTrackMappability.wg.gz files.")
    args = parser.parse_args()
    if args.SingleNucleotide:
        if not os.path.exists(args.wigdir):
            raise ValueError("-wigdir didn't exist!")
    out_list = [args.bed_path, args.out_path,
                args.umap_dir, args.kmer,
                args.wigdir]
    return out_list


def load_wig(wig_path, chrom_size):
    '''
    Translated a fixed step wig file of a chromosome
    with step size 1 to a numpy array.
    Will be used internally by BedMappability.
    Unexpected behaviour if wig file has multiple
    chromosomes.
    '''
    wig_ar = np.zeros(chrom_size, dtype=float)
    pos = 0
    with gzip.open(wig_path, "rb") as wig_link:
        for wig_line in wig_link:
            if "start" in wig_line:
                pos_str = wig_line.rstrip().split(" ")[2]
                pos = int(pos_str.replace("start=", "")) - 1
            else:
                value = float(wig_line.rstrip())
                wig_ar[pos] = value
                pos = pos + 1
    print("Loaded {}".format(wig_path))
    return wig_ar


def bed_handler(bed_path):
    '''
    Saves a gzipped BED file to a pandas
    dataframe.
    '''
    no_header = True
    with gzip.open(bed_path, "rb") as bed_link:
        header_line = bed_link.readline()
        header_list = header_line.rstrip().split("\t")
        if len(header_list) < 3:
            no_header = False
        elif "chr" not in header_list[0]:
            no_header = False
    if no_header:
        bed_df = pd.read_csv(
            bed_path,
            sep="\t",
            header=None,
            compression="gzip")
    else:
        bed_df = pd.read_csv(
            bed_path,
            skiprows=1,
            sep="\t",
            header=None,
            compression="gzip")
    column_names = list(bed_df.columns)
    column_names[:3] = ["Chromosome", "Start", "End"]
    bed_df.columns = column_names
    return bed_df


def subsetter(my_list, regex):
    '''
    Finds elements of a list matching a regex
    '''
    out_list = []
    for each in my_list:
        RE = re.search(regex, each)
        if RE:
            out_list.append(each)
    return out_list


class BedMappability:
    '''
    This class provides methods for identifying
    single-read and multi-read mappability of
    an input BED file. it uses local functions
    such as bed_handler and get_wig if necessary.
    '''
    def __init__(self, bed_path, umap_dir,
                 out_path, kmer, wigdir="NA"):
        '''

        :param bed_path: Input gzipped BED path
        :param umap_dir: Umap/Bismap unique uint directory
        :param out_path: Path to output gzipped BED
        :param kmer: Integer kmer size
        :param wigdir: Directory of multi-read mappability
            files of different chromosomes. Optional.
        '''
        self.bed_path = bed_path
        self.umap_dir = umap_dir
        self.out_path = out_path
        self.kmer = kmer
        self.kmer_str = "k{}".format(kmer)
        self.wigdir = wigdir

    def annotate_bed(self):
        '''
        Main motor method for getting single-read
        and multi-read mappability of BedMappability object.

        Returns:
            Pandas dataframe BED file with mappability info
        '''
        bed_df = bed_handler(self.bed_path)
        all_chrs = list(pd.unique(bed_df["Chromosome"]))
        all_chrs.sort()
        bed_out_list = []
        for each_chr in all_chrs:
            print("Started working on {}".format(each_chr))
            temp_df = bed_df[bed_df["Chromosome"] == each_chr]
            mapped_df = self.get_df_mappability(temp_df, each_chr)
            if len(mapped_df) > 0:
                bed_out_list.append(mapped_df)
        out_df = pd.concat(bed_out_list)
        return out_df

    def get_df_mappability(self, temp_df, cur_chr):
        '''
        Internally used method for finding the single-read
        and multi-read mappability of a bed file for a given
        chromosome.
        '''
        dists = pd.unique(temp_df.iloc[:, 2] - temp_df.iloc[:, 1])
        FIXED_DIST_1 = False
        if len(dists) < 2:
            if 1 in dists or 0 in dists:
                FIXED_DIST_1 = True
        temp_df["SingleRead.Mappability"] = 0
        temp_df["MultiRead.Mappability"] = 0
        umap_paths = subsetter(os.listdir(self.umap_dir),
                               "uint8.unique.gz")
        umap_path = subsetter(umap_paths, "{}.uint".format(cur_chr))
        if len(umap_path) == 0:
            temp_df = []
        else:
            umap_path = umap_path[0]
            uint_path = "{}/{}".format(
                self.umap_dir, umap_path)
            uint_link = gzip.open(uint_path, "rb")
            uint_ar = np.frombuffer(uint_link.read(), dtype=np.uint8)
            uint_link.close()
            print("Using {}".format(umap_path))
            kmer = self.kmer
            wig_path = "{}/{}_{}.MultiTrackMappability.wg.gz".format(
                self.wigdir, self.kmer_str, cur_chr)
            if os.path.exists(wig_path) and FIXED_DIST_1:
                wig_ar = load_wig(wig_path, len(uint_ar))
                temp_df["MultiRead.Mappability"] = \
                    wig_ar[temp_df.iloc[:, 1] - 1]
                temp_df["SingleRead.Mappability"] = \
                    temp_df["MultiRead.Mappability"]
                ind_unique, = np.where(temp_df["SingleRead.Mappability"] > 0)
                temp_df["SingleRead.Mappability"].iloc[ind_unique] = 1
            else:
                less_than_kmer = uint_ar <= kmer
                not_zero = uint_ar != 0
                uniquely_mappable = np.array(
                    np.logical_and(less_than_kmer, not_zero),
                    dtype=int)
                for i in range(temp_df.shape[0]):
                    chr, start, end = temp_df.iloc[i, :3]
                    st_vec = start - kmer
                    if st_vec < 0:
                        st_vec = 0
                    region_vec = np.zeros(
                        len(range(st_vec, end + 1)),
                        dtype=float)
                    single_vec = np.zeros(
                        len(range(st_vec, end + 1)),
                        dtype=float)
                    idxs_unique, = np.where(
                        [each_idx == 1 for each_idx in
                         uniquely_mappable[st_vec:(end + 1)]])
                    for idx_unique in idxs_unique:
                        end_idx = idx_unique + kmer
                        if end_idx > len(region_vec):
                            end_idx = len(region_vec)
                        region_vec[idx_unique:end_idx] = 1
                        single_vec[idx_unique:end_idx] = \
                            single_vec[idx_unique:end_idx] + 1
                    unique_map_percent = float(
                        sum(single_vec[kmer:(kmer + end - start + 1)] > 0)) /\
                        float(len(single_vec[kmer:(kmer + end - start + 1)]))
                    map_prob = float(len(idxs_unique)) / float(len(region_vec))
                    temp_df.iloc[i, -1] = map_prob
                    temp_df.iloc[i, -2] = unique_map_percent
        return temp_df

    def write_bed(self, mapped_df):
        '''
        Use this method to write the output of
        BedMappability.annotate_bed() to a gzipped file,
        '''
        out_link = gzip.open(self.out_path, "wb")
        mapped_df.to_csv(out_link, sep="\t", index=False, header=False)
        out_link.close()
        print("Successfully created {}".format(self.out_path))


if __name__ == "__main__":
    bed_path, out_path, umap_dir, kmer, wigdir = get_args()
    BedMapObj = BedMappability(bed_path, umap_dir,
                               out_path, kmer, wigdir)
    mapped_df = BedMapObj.annotate_bed()
    BedMapObj.write_bed(mapped_df)
