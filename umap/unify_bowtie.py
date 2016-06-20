from argparse import ArgumentParser
import gzip
import numpy as np
import os
import pandas as pd
import re


def subset_list(list_items, regex):
    out_list = []
    for each_item in list_items:
        RE = re.search(regex, each_item)
        if RE:
            out_list.append(each_item)
    return out_list


class UnifyBowtie:
    def __init__(self, bowtie_outdir, chrsize_path, job_id):
        """Merges bowtie.gz outputs of run_bowtie

        Based on outputs of run_bowtie that are saved in format
        of <chr>.<kmer>.<job_id>.bowtie.gz, this class uses
        a variable ID to select a particular chromosome and
        merge the data in all of the different bowtie.gz files
        of that chromosome.

        :param bowtie_outdir: Directory with <chr>.<kmer>.<job_id>.bowtie.gz
        :param chrsize_path: Path to 2-column file: <chr>\t<size>\n...
        :param int job_id: A 0-based index to select chromosome
            based on chrsize_path

        :returns: Saves the output to bowtie_outdir/
            <chr>.<kmer>.uint8.unique.gz file
        """
        self.bowtie_outdir = bowtie_outdir
        self.chrsize_path = chrsize_path
        self.ind_chr = job_id
        self.chr_dict = self.make_chr_dict()
        self.bowtie_to_unique()

    def make_chr_dict(self):
        """Makes a dictionary using self.chrsize_path"""
        chr_dict = {}
        with open(self.chrsize_path, "r") as chrsize_link:
            for chrsize_line in chrsize_link:
                chrom, size = chrsize_line.rstrip().split("\t")
                chr_dict[chrom] = int(size)
        return chr_dict

    def get_mapped_positions(self, bowtie_path):
        """Finds mapped regions in bowtie output

        In a gzipped bowtie output with perfect matches,
        filters the results to those in the forward strand
        and saves them to an array.

        Why results filtered for forward strand?
        The results uniquely mapped
        to the reverse strand are rare and a mistake of the
        aligning algorithm because all of the k-mers are
        generated from the forward strand. If match happens
        on reverse strand, it means that the match is not the
        only unique match (reads were generated from forward
        strand).
        """
        bowtie_df = pd.read_csv(
            bowtie_path, sep="\t",
            compression="gzip",
            names=["Ind", "Strand", "Chr", "Start"])
        bowtie_df = bowtie_df[bowtie_df["Strand"] == "+"]
        ind_ar = bowtie_df.iloc[:, 3]
        return ind_ar

    def get_other_chr_name(self, chr_name):
        """Finds name of the reverse complement chromosome

        In Bismap, we generate a set of reverse complemented
        chromosomes to account for bisulfite conversion before
        reverse complementation. This function finds those
        chromosomes"""
        new_name = chr_name
        if "RC" in chr_name:
            len_paths = len(
                subset_list(
                    os.listdir(self.bowtie_outdir), chr_name))
            if len_paths == 0:
                new_name = chr_name + "_RC"
        return new_name

    def bowtie_to_unique(self):
        """Wrapper method of UnifyBowtie

        Uses information of bowtie.gz files to find mappability
        of a given chromosome. Is automatically called by
        UnifyBowtie."""
        KMER = int(self.bowtie_outdir.split("/")[-1].split("k")[-1])
        all_chrs = self.chr_dict.keys()
        all_chrs.sort()
        chrom = all_chrs[self.ind_chr]
        size = self.chr_dict[chrom]
        new_chr_name = self.get_other_chr_name(chrom)
        chr_paths = ["{}/{}".format(self.bowtie_outdir, bowtie_path)
                     for bowtie_path in subset_list(
                     os.listdir(self.bowtie_outdir),
                     "{}\.".format(new_chr_name))]
        bowtie_paths = subset_list(chr_paths, ".bowtie.gz")
        unique_ar = np.zeros(size, dtype=np.uint8)
        for bowtie_path in bowtie_paths:
            mapped_indices = self.get_mapped_positions(bowtie_path)
            for st_index in mapped_indices:
                end_index = st_index + KMER
                unique_ar[st_index:end_index] = 1
            print("Done with {}".format(bowtie_path))
        out_path = "{}/{}.k{}.uint8.unique.gz".format(
            self.bowtie_outdir, chrom, KMER)
        out_link = gzip.open(out_path, "wb")
        # np.save(out_link, unique_ar)
        out_link.write(unique_ar.tobytes())
        # unique_ar.tofile(out_link)
        out_link.close()
        print("Saved {}".format(out_path))
        print("Exiting successfully")


if __name__ == "__main__":
    parser = ArgumentParser(
        description="")
    parser.add_argument(
        "bowtie_outdir",
        help="Directory containing bowtie output files")
    parser.add_argument(
        "chrsize_path",
        help="A file containing the order of chromosome names\
        to consider (one chromosome name per line)")
    parser.add_argument(
        "-job_id",
        default=0,
        type=int,
        help="If not using a cluster for submitting jobs, "
        "specify the job_id by integer ranging from 1 to "
        "total number of chromosomes in chrsize_path")
    parser.add_argument(
        "-var_id",
        default="SGE_TASK_ID",
        help="HPC variable name for job ID (1-based index)")
    args = parser.parse_args()
    job_id = args.job_id
    if job_id == 0:
        job_id = int(os.environ[args.var_id]) - 1
    UnifyBowtie(args.bowtie_outdir, args.chrsize_path, job_id)
