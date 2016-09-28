from argparse import ArgumentParser
from datetime import datetime
import os
import re
import subprocess


def subset_list(list_items, regex):
    out_list = []
    for each_item in list_items:
        RE = re.search(regex, each_item)
        if RE:
            out_list.append(each_item)
    return out_list


class BowtieWrapper:
    def __init__(self, kmer_dir, bowtie_dir,
                 index_dir, index_name,
                 job_id, Bismap):
        """Runs Bowtie one <chrom>.<kmer>.<jobid>.kmer.gz

        Using the job_id, this function identifies
        one kmer.gz file and runs bowtie on that and saves
        the output to <chrom>.<kmer>.<jobid>.bowtie.gz.

        :param kmer_dir: Directory with <chrom>.<kmer>.<jobid>.kmer.gz files
        :param bowtie_dir: Directory with Bowtie 1.1.0 executable files.
        :param index_dir: Directory with Bowtie index
        :param index_name: Name used for generating Bowtie index files
        :param int job_id: will be used for finding kmer.gz file
        :param bool Bismap: Run bowtie with --norc

        :returns: Saves the output to a file in the same directory as kmer_dir
        """
        self.kmer_dir = kmer_dir
        self.bowtie_dir = bowtie_dir
        self.index_dir = index_dir
        self.index_name = index_name
        self.job_id = job_id
        self.Bismap = Bismap
        self.execute_bowtie_command()

    def execute_bowtie_command(self):
        """The only method of BowtieWrapper

        Will be executed automatically by BowtieWrapper

        :raises ValueError: If job_id is out of expected range
        """
        if self.Bismap:
            rev_comp = " --norc "
        else:
            rev_comp = ""
        kmer_names = ["{}/{}".format(self.kmer_dir, each_kmer) for each_kmer
                      in subset_list(os.listdir(self.kmer_dir), ".kmer.gz$")]
        kmer_names.sort()
        LongIndex = False
        short_ind_path = "{}/{}.1.ebwtl".format(
            self.index_dir, self.index_name)
        if os.path.exists(short_ind_path):
            LongIndex = True
            print("Switching to use of long index")
        if job_id <= len(kmer_names):
            try:
                kmer_file = kmer_names[job_id]
            except:
                raise ValueError(
                    "{} does not exist. Time: {}".format(
                        job_id, str(datetime.now())))
            print("processing Kmer File {}".format(kmer_file))
            # kmer_name = kmer_dir.split("/")[-1]
            kmer_path = "{}/{}".format(self.kmer_dir, kmer_file.split("/")[-1])
            bowtie_out_path = kmer_path.replace(".kmer.gz", ".bowtie.gz")
            first_part_of_command = "gunzip -c {} | {}/bowtie ".format(
                kmer_path, self.bowtie_dir)
            if LongIndex:
                first_part_of_command = first_part_of_command +\
                    "--large-index "
            bowtiecmd = first_part_of_command +\
                "{}/{} ".format(self.index_dir, self.index_name) +\
                "-v 0 -k 1 -m 1 {}--mm ".format(rev_comp) +\
                "-r --refidx --suppress 5,6,7,8 - " +\
                "| gzip -c > {}".format(bowtie_out_path)
            subprocess.call(bowtiecmd, shell=True)
            print("Executing {}".format(bowtiecmd))
        else:
            print("The length of files was {} but the index was {}".format(
                len(kmer_names), job_id))


def check_genome(index_dir, index_name):
    all_idx_paths = [
        each_path for each_path in
        os.listdir(index_dir) if index_name in each_path]
    if len(all_idx_paths) < 6:
        raise ValueError("Index does not exist")


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Umap wrapper for running bowtie "
        "on individual k-mer files.")
    parser.add_argument(
        "kmer_dir",
        help="Directory containing the .kmer files")
    parser.add_argument(
        "bowtie_dir",
        help="Directory containing bowtie executable")
    parser.add_argument(
        "index_dir",
        help="Directory containing bowtie index")
    parser.add_argument(
        "index_name",
        help="prefix name of bowtie index")
    parser.add_argument(
        "-Bismap",
        action="store_true",
        help="Run bowtie with --norc")
    parser.add_argument(
        "-var_id",
        default="SGE_TASK_ID",
        help="HPC environmental variable for JOB ID")
    parser.add_argument(
        "-job_id",
        type=int,
        default=0,
        help="1-based index for selecting a k-mer file")
    args = parser.parse_args()
    check_genome(args.index_dir, args.index_name)
    job_id = args.job_id
    if job_id == 0:
        job_id = int(os.environ[args.var_id]) - 1
    BowtieWrapper(args.kmer_dir, args.bowtie_dir,
                  args.index_dir, args.index_name, job_id,
                  args.Bismap)
