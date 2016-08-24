import argparse
from datetime import datetime
import gzip
import numpy as np
import os
import re
import warnings


def subset_list(my_list, regex, reg_neg=False):
    out_list = []
    for item in my_list:
        RE = re.search(regex, item)
        if RE and not reg_neg:
            out_list.append(item)
        elif not RE and reg_neg:
            out_list.append(item)
    return out_list


class Int8Handler:
    def __init__(self, in_dir, out_dir, C2T, G2A, chrsize_path):
        """Infers important directory structure from in_dir

        The Int8Handler class requires in_dir argument which
        is a directory with subfolders named as k<integer>
        that have <chrom>.<kmer>.uint8.unique.gz files.
        The Int8Handler.write_beds method creates one
        browser extensible data (BED) file for each of the
        k<integer> subfolders.

        :param in_dir: Directory with <chrom>.uint8.unique.gz files
        :param out_dir: Directory for saving output BED files
        """
        self.in_dir = in_dir
        self.C2T = C2T
        self.G2A = G2A
        self.chrsize_path = chrsize_path
        self.in_prev_dir = "/".join(in_dir.split("/")[:-1])
        self.fix = ".uint8.unique.gz"
        self.uint8s = subset_list(os.listdir(in_dir), self.fix)
        self.prev_dir = "/".join(self.in_dir.split("/")[:-1])
        self.kmers = subset_list(next(os.walk(self.in_prev_dir))[1],
                                 "^k")
        self.out_dir = out_dir
        self.type = self.find_type()
        self.chrsize_dict = self.make_chrsize_dict()
        self.chroms = [each_chr for each_chr in
                       self.chrsize_dict.keys() if
                       "RC" not in each_chr]
        self.chroms.sort()

    def find_type(self):
        if self.C2T:
            self.type = "Bismap.Forward"
        elif self.G2A:
            self.type = "Bismap.Reverse"
        else:
            self.type = "Umap"
        return self.type

    def write_beds(self, out_label, kmer_cur,
                   WriteUnique=False):
        """Convert uint8 files to BED

        A method of class Int8Handler, it reads uint8 files of
        different chromosomes and saves them into a BED file.
        This is the only method in this class that needs
        to be called manually by the user.

        Args:
            out_label: Will be used in output names: <kmer>.<out_label>.bed.gz
            kmer_cur: k-mer size in format of k<integer>
            WriteUnique: Defaults to False. Will merge and save uint files.
                Use if working with -Bismap and have reverse
                complemented chromosomes.

        Returns:
            Saved ths output to a gzipped BED file

        Raises:
            Warning: If no uniquely mappable read is found in
                any of the uint8 arrays.
        """
        kmer = int(kmer_cur.replace("k", ""))
        STRAND = "+"
        if self.G2A:
            STRAND = "-"
        in_dir = self.in_dir
        all_chrs = self.chroms
        all_chrs.sort()
        for cur_chr in self.chroms:
            other_chr = cur_chr + "_RC"
            uint_path = "{}/{}{}".format(
                in_dir, cur_chr, self.fix)
            uniquely_mappable = self.load_uint_ar(
                uint_path, kmer, cur_chr)
            if self.type != "Umap":
                uint_path_other = in_dir +\
                    "/" + uint_path.replace(cur_chr, other_chr)
                uniquely_mappable_other = self.load_uint_ar(
                    uint_path_other, kmer, other_chr)

                # Remove the last k nucleotides of the reverse chromosome
                uniquely_mappable_other = uniquely_mappable_other[:-kmer]

                # Reverse the chromosome (Hence it's a reverse complement)
                uniquely_mappable_other = uniquely_mappable_other[::-1]

                # Remove the last k nucleotides of the forward chromosome
                uniquely_mappable = uniquely_mappable[:-kmer]

                # Merge to find uniquely mappable reads in both strands
                uniquely_mappable = np.array(
                    np.logical_and(
                        uniquely_mappable_other, uniquely_mappable),
                    dtype=int)
                if WriteUnique:
                    unique_ar = np.array(uniquely_mappable, dtype=np.uint8)
                    out_dir_uint = "{}/k{}".format(self.out_dir, kmer)
                    if not os.path.exists(out_dir_uint):
                        os.makedirs(out_dir_uint)
                    out_path_uint = "{}/{}.k{}.uint8.unique.gz".format(
                        (out_dir_uint, cur_chr, kmer))
                    out_link_uint = gzip.open(out_path_uint, "wb")
                    out_link_uint.write(unique_ar.tobytes())
                    out_link_uint.close()

            bed_kmer_pos = self.get_bed6(
                uniquely_mappable, kmer, STRAND, uint_path, cur_chr)

            # Write the BED6 to a gzipped tsv file
            out_name = "{}/{}.{}.bed.gz".format(
                self.out_dir, kmer_cur, out_label)
            if cur_chr == all_chrs[0]:
                header = self.make_header(kmer_cur, type)
                out_link = gzip.open(out_name, "wb")
                out_link.write(header)
            else:
                out_link = gzip.open(out_name, "ab")
            if len(bed_kmer_pos) > 0:
                print "Found %d regions in %s" %\
                    (bed_kmer_pos.shape[0], cur_chr)
                for bed_line in bed_kmer_pos:
                    line_out = [str(val) for val in bed_line]
                    line_out = "\t".join(line_out) + "\n"
                    out_link.write(line_out)
                print "Created data of %s at %s" %\
                    (cur_chr, str(datetime.now()))
            out_link.close()

    def get_bed6(self, uniquely_mappable, kmer,
                 STRAND, uint_path, cur_chr):
        """Make BED6 from a binary vector

        Converts a binary vector with the same length
        as the chromosome to a BED6 file. Each 1
        entry in the binary file indicates that the
        k-mer starting at that position and ending at
        k nucleotides downstream is uniquely mappable.
        Thus the BED6 would show any region in the genome
        that is uniquely mappable by at least one k-mer.

        :param uniquely_mappable: numpy binary array (0 or 1 values)
        :param kmer: Integer scalar showing read length (k-mer)
        :param STRAND: Strand that will be saved to the BED6
        :param uint_path: Path of the uint8 array that was used
        :param cur_chr: Chromosome name the the data is from

        :returns: A numpy string array with BED6 information
        """
        unimap_diff = np.diff(uniquely_mappable)
        poses_start, = np.where(unimap_diff == 1)
        poses_end, = np.where(unimap_diff == -1)
        if len(poses_start) != len(poses_end):
            if len(poses_start) > len(poses_end):
                poses_end = np.append(poses_end, [len(uniquely_mappable)])
            else:
                poses_start = np.append([0], poses_start)
        elif uniquely_mappable[0] == 1:
            poses_start = np.append([0], poses_start)
            poses_end = np.append(poses_end, [len(uniquely_mappable)])
        if len(poses_start) == 0:
            warnings.warn(
                "Found no uniquely mappable reads for {}!".format(
                    uint_path))
            bed_kmer_pos = []
        else:
            chr_length = self.chrsize_dict[cur_chr]
            bed_kmer_pos = np.empty(shape=(len(poses_start), 6),
                                    dtype="S16")
            ind_high = []
            for ind_st in range(len(poses_start)):
                if poses_end[ind_st] + kmer - 1 > chr_length:
                    ind_high.append(ind_st)
            bed_kmer_pos = np.array(
                [[cur_chr, str(poses_start[i] + 1),
                  str(poses_end[i] + kmer - 1),
                  "k" + str(kmer),
                  1, STRAND] for i in range(len(poses_start))],
                dtype="S64")
            for each_ind in ind_high:
                if int(bed_kmer_pos[each_ind][2]) > chr_length:
                    bed_kmer_pos[each_ind][2] = str(chr_length)
        return bed_kmer_pos

    def load_uint_ar(self, uint_path, kmer, cur_chr):
        """Loads a gzipped unsigned 8-bit integer as a numpy array
        """
        uint_link = gzip.open(uint_path, "rb")
        uint_ar = np.frombuffer(uint_link.read(), dtype=np.uint8)
        uint_link.close()
        print "Processing {} for {} {}".format(
            uint_path, kmer, cur_chr)
        less_than_kmer = uint_ar <= kmer
        not_zero = uint_ar != 0
        uniquely_mappable = np.array(
            np.logical_and(less_than_kmer, not_zero), dtype=int)
        return uniquely_mappable

    def make_header(self, kmer, type):
        """Created BED6 header

        :param str kmer: k<integer>
        :param type: One of Bismap.Forward, Bismap.Reverse or Umap
        """
        dict_col = {"Bisulfite.Converted": "240,40,80",
                    "Umap": "80,40,240",
                    "Bismap.Forward": "220,20,80",
                    "Bismap.Reverse": "80,20,220"}
        header = 'track name="{} {}"'.format(type, kmer) +\
            'description="Single-read mappability with {}-mers'.format(
                kmer) +\
            'color=%s \n'.format(dict_col.get(type, "40,40,240"))
        return header

    def make_chrsize_dict(self):
        """Creates dictionary from 2-column chromosome size file
        """
        chrsize_dict = {}
        with open(self.chrsize_path, "r") as chrsize_link:
            for each_line in chrsize_link:
                chr = each_line.rstrip("\n").split("\t")[0]
                length = int(each_line.rstrip("\n").split("\t")[1])
                chrsize_dict[chr] = length
        return chrsize_dict

    def bin_arr_to_wig(self, uniquely_mappable, kmer):
        """Converts a uniquely a binary array to mult-read mappability

        Unique mappability array is a binary array where each
        value of 1 indicates that the k-mer starting at that
        position is uniquely mappable. This function generates
        an array of the same length and saves the multi-read mappability
        of each position in that array. Multi-read mappability is the
        probability of finding a uniquely mappable k-mer among all of
        the k-mers that overlap with a given position.

        :param uniquely_mappable: A binary numpy array
        :param kmer: An integer defining read length

        :returns: Multi-read mappability array.
        """
        LEN_AR = len(uniquely_mappable)
        poses_start, = np.where(uniquely_mappable == 1)
        ar_quant = np.zeros(len(uniquely_mappable), dtype=float)
        for ind_st in poses_start:
            ind_end = ind_st + kmer
            if ind_end > LEN_AR:
                ind_end = LEN_AR
            ar_quant[ind_st:ind_end] = ar_quant[ind_st:ind_end] + 1.0
        print("Finished generating multi-read mappability date at {}".format(
            str(datetime.now())))
        del poses_start
        ar_quant = ar_quant / float(kmer)
        return ar_quant

    def write_as_wig(self, uint_path, out_path, kmer, chrom):
        """unsigned 8-bit integer array file to wiggle

        For a given numeric unsigned 8-bit integer vector that
        is generated by Umap, this method save the wiggle file
        which is specific for one chromosome over one read length.

        :param uint_path: Path to a gzipped numeric unsigned 8-bit array
        :param out_path: Gzipped path for saving wiggle
        :param kmer: Integer defining read length
        :param chrom: Chromosome that the uint_path is specific to
        """
        uniquely_mappable = self.load_uint_ar(uint_path, kmer, chrom)
        ar_quant = self.bin_arr_to_wig(uniquely_mappable, kmer)
        out_link = gzip.open(out_path, "wb")
        unimap_diff = np.diff(np.array(ar_quant > 0, dtype=int))
        poses_start, = np.where(unimap_diff == 1)
        poses_end, = np.where(unimap_diff == -1)
        if len(poses_start) != len(poses_end):
            if len(poses_start) > len(poses_end):
                poses_end = np.append(poses_end, [len(uniquely_mappable)])
            else:
                poses_start = np.append([0], poses_start)
        elif uniquely_mappable[0] == 1:
            poses_start = np.append([0], poses_start)
            poses_end = np.append(poses_end, [len(uniquely_mappable)])
        for ind_st in range(len(poses_start)):
            pos_st = poses_start[ind_st] + 1
            pos_end = poses_end[ind_st] + 1
            start_line = "fixedStep chrom={} start={} step=1 span=1\n".format(
                chrom, pos_st)
            out_link.write(start_line)
            for each_pos in range(pos_st, pos_end):
                out_link.write(str(ar_quant[each_pos]) + "\n")
        print(
            "Finished saving the data to {} at {}".format(
                out_path, str(datetime.now())))
        out_link.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts\
    unionized uint8 outputs of umap/ubismap to bed files")
    parser.add_argument(
        "in_dir",
        help="folder with <chrom>.uint8.unique.gz files")
    parser.add_argument(
        "out_dir",
        help="Folder for writing the output files")
    parser.add_argument(
        "out_label",
        help="File names would be kmer.<out_label>.bed.gz")
    parser.add_argument(
        "-C2T",
        action="store_true",
        help="If using converted genomes specify -C2T or -G2A")
    parser.add_argument(
        "-G2A",
        action="store_true",
        help="If using converted genomes specify -C2T or -G2A")
    parser.add_argument(
        "-chrsize_path",
        default="../../chrsize.tsv",
        help="Path to a 2 column file of chromosome and length. "
        "By default it goes to ../../chrsize.tsv from out_dir")
    parser.add_argument(
        "-WriteUnique",
        action="store_true",
        help="If -Bismap is true and want to store the merged "
        "uint file, specify this option")
    parser.add_argument(
        "-wiggle",
        action="store_true",
        help="If specified, will generate wiggle files "
        "for each chromosome. Make sure to specify -job_id "
        "or run in job array for parallel computation.")
    parser.add_argument(
        "-bed",
        action="store_true",
        help="If specified, will generate bed files that specify "
        "all of the regions in the genome that are uniquely mappable "
        "by each of the k-mers")
    parser.add_argument(
        "-job_id",
        type=int,
        default=0,
        help="If not using job array, specify this index "
        "which will be used for selecting the chromosomes")
    parser.add_argument(
        "-var_id",
        default="SGE_TASK_ID",
        help="Environmental variable for finding chromosome indices")
    args = parser.parse_args()
    if not args.bed and not args.wiggle:
        raise ValueError("Please specify only one of -bed or -wiggle")
    if args.bed and args.wiggle:
        raise ValueError("Please specify at least one of -bed or -wiggle")
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    FileHandler = Int8Handler(args.in_dir,
                              args.out_dir, args.C2T, args.G2A,
                              args.chrsize_path)
    job_id = args.job_id
    PARALLEL = True
    if job_id == 0:
        job_id = os.environ[args.var_id]
        if job_id == "":
            PARALLEL = False
    kmers = FileHandler.kmers
    # for kmer in kmers:
    for kmer in kmers:
        if args.bed:
            print("Creating BED file")
            FileHandler.write_beds(args.out_label,
                                   kmer, args.WriteUnique)
        elif PARALLEL and args.wiggle:
            chrom = FileHandler.chroms[int(job_id) - 1]
            print(
                "Saving Wiggle using job id {}, selected {}".format(
                    job_id, chrom))
            out_path = "{}/{}.{}.{}.MultiReadMappability.wg.gz".format(
                args.out_dir, args.out_label, chrom, kmer)
            uint_path = "{}/{}.uint8.unique.gz".format(
                args.in_dir, chrom)
            kmer_num = int(kmer.replace("k", ""))
            FileHandler.write_as_wig(uint_path, out_path, kmer_num, chrom)
        elif args.wiggle:
            print("Creating wiggles for all chromosomes consequently")
            print("This may take long...")
            for chrom in FileHandler.chroms:
                out_path = "{}/{}.{}.{}.MultiReadMappability.wg.gz".format(
                    args.out_dir, args.out_label, chrom, kmer)
                uint_path = "{}/{}.uint8.unique.gz".format(
                    args.in_dir, chrom)
                kmer_num = int(kmer.replace("k", ""))
                FileHandler.write_as_wig(uint_path, out_path, kmer_num, chrom)
                print(
                    "Tired of waiting? Try parallel processing "
                    "by specifying -job_id")
