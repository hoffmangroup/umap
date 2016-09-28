from argparse import ArgumentParser
import gzip
import numpy as np
import os
import pandas as pd


def get_args():
    parser = ArgumentParser(
        description="Creates list of unique Kmers "
        "for each 1Mb of the genome. Relies on an environmental "
        "variable such as SGE_TASK_ID to identify the correct "
        "megabase of the genome for this purpose. "
        "If not using job arrays, specify -job_id manually.")
    parser.add_argument(
        "chromsize_path",
        help="Path to 2 column tsv file where first column is "
        "chromosome name and second column is chromosome size")
    parser.add_argument(
        "out_dir",
        help="Path to the directory for creating "
        "<chromosome>.<Megabase>.<kmer>.kmer.gz files")
    parser.add_argument(
        "chr_dir",
        help="Path to directory with <chromosome>.fasta files.")
    parser.add_argument(
        "idx_path",
        help="Path to the 4 column file with the following columns: "
        "Index | Chromosome | Start | End. This file will be used "
        "for identifying the chunk of the chromosome.")
    parser.add_argument(
        "--kmer",
        default="infer",
        help="The software would infer it based on the "
        "name of the 'out_dir'. If it is set and "
        "contradicts the 'out_dir', a subfolder "
        "under out_dir will be created named 'kmer' and "
        "out_dir will be changed to that.")
    parser.add_argument(
        "--job_id",
        default=0,
        type=int,
        help="If not submitted in job array, would require this "
        "parameter to be set. (1-based index)")
    parser.add_argument(
        "--var_id",
        default="SGE_TASK_ID",
        help="The variable name that the script would use "
        "for identifying the job id. By default: SGE_TASK_ID.")
    args = parser.parse_args()
    out_dir = args.out_dir
    inferred_kmer = args.out_dir.split("/")[-1]
    kmer = args.kmer
    job_id = args.job_id
    if args.kmer == "infer":
        kmer = inferred_kmer
    elif args.kmer != inferred_kmer:
        kmer = args.kmer
        out_dir = "{}/{}".format(out_dir, kmer)
    if job_id == 0:
        job_id = int(os.environ[args.var_id]) - 1
    out_list = [args.chromsize_path, out_dir,
                kmer, job_id, args.chr_dir,
                args.idx_path]
    return out_list


class GetKmers:
    def __init__(self, out_dir, kmer, job_id,
                 chr_dir, chromsize_path, idx_path):
        """Creates all the possible k-mers for part of the genome.

        Used a referece file to find the appropriate chromosome,
        start and end position. Passes through the fasta file
        of the chromosome and generates all of the possible k-mers.

        Args:
        :param out_dir: Directory for saving <chrom>.<jobid>.kmer.gz files
        :param str kmer: k-mer string such as 'k24'
        :param int job_id: Reference ID used for finding chrom, start and end
        :param chr_dir: Path to directory with chromosome fasta files
        :param chromsize_path: Path to 2 column file of chrom\tsize\n
        :param idx_path: Path to 4 column file of index\tchrom\t\st\tend\n

        :returns: An object with methods such as get_step_fasta(),
            get_seq_ar(), write_kmers() and write_regions().

        :raises ValueError: if expected chromosome path does not exist
        """
        self.out_dir = out_dir
        self.kmer = kmer
        self.job_id = job_id
        self.chromsize_path = chromsize_path
        self.chr_dir = chr_dir
        self.idx_path = idx_path
        self.chrom, self.start, self.end = self.get_region()
        self.chrom_path = "{}/{}.fasta".format(
            self.chr_dir, self.chrom)
        if not os.path.exists(self.chrom_path):
            raise ValueError(
                "{} does not exist".format(self.chrom_path))
        elif not os.path.exists(self.idx_path):
            raise ValueError(
                "{} does not exist".format(self.idx_path))

    def get_region(self):
        """Find chromosomal chunk

        Using job_id and chromsize_path, finds chromosome,
        start and end. Looks if chrsize_index.tsv exists.
        If it doesn't, it will create chrsize_index.tsv.
        Otherwise will use it for finding that information.

        :returns: chrom, start, end

        :raises ValueError: If job_id is out of expected range
        """
        job_id = self.job_id
        idx_path = self.idx_path
        idx_df = pd.read_csv(idx_path, sep="\t", index_col=0)
        dict_idx = idx_df.to_dict()
        try:
            chrom, start, end = [
                dict_idx[each_line][job_id] for
                each_line in ["Chromosome", "Start", "End"]]
        except:
            raise ValueError(
                "{} Job id is larger than available indices".format(
                    job_id))
        return chrom, int(start), int(end)

    def get_step_fasta(self):
        """Finds number of nucleotides at each line of FASTA


        :raises ValueError: If top 10 FASTA lines have
            sequences of varying length
        """
        chrom_path = self.chrom_path
        line_nums = 10
        len_lines = []
        with open(chrom_path, "r") as chr_link:
            for line_num in range(line_nums):
                chr_line = chr_link.readline().rstrip()
                if ">" not in chr_line:
                    len_line = len(chr_line)
                    if len_line not in len_lines:
                        if chr_line != "":
                            len_lines.append(len_line)
        if len(len_lines) == 1:
            return(len_lines[0])
        else:
            raise ValueError(
                "Top 10 FASTA lines have sequences of varying length")

    def get_seq_ar(self, fasta_fix):
        """Gets sequence as numpy array

        get_seq_ar creates a numpy array matching sequence
        of the region specified by chromosome, start and
        end. However it starts the sequence with kmer - 1
        less than the specified start to account for the
        missing kmers from the previous chunk of chromosome.
        Each element is one line of a fasta file that may be
        cut to match the exact start/end positions.

        Args:
        :param int fasta_fix: Number of nucleotides per fasta line.

        :returns: Numpy array with each member being one line of FASTA
        """
        start = self.start
        end, chrom_path = self.end, self.chrom_path
        kmer = int(self.kmer.replace("k", ""))
        start = start - kmer
        if start < 0:
            start = 0
        st_line = int(np.floor(start / fasta_fix))
        st_line_pos = start % fasta_fix
        end_line = int(np.floor(end / fasta_fix))
        end_line_pos = (end % fasta_fix) + 1
        seq_ar = np.empty(int(end_line - st_line + 1),
                          dtype="|S{}".format(fasta_fix + 1))
        ind_seq = 0
        STORE_LINE = False
        with open(chrom_path, "r") as chr_link:
            line_ind = 0
            for each_line in chr_link:
                if line_ind == st_line:
                    STORE_LINE = True
                    each_line = each_line.rstrip()[st_line_pos:]
                if STORE_LINE:
                    if line_ind == end_line:
                        STORE_LINE = False
                        each_line = each_line.rstrip()[:end_line_pos]
                    seq_ar[ind_seq] = each_line.rstrip()
                    ind_seq = ind_seq + 1
                line_ind = line_ind + 1
        if "chr" in seq_ar[0]:
            seq_ar = seq_ar[1:]
        return seq_ar

    def write_kmers(self, cur_seq, kmer_int, out_link):
        """Creates k-mers and saves them to the file

        Using a numpy array of sequences, generates and
        saves all of the possible k-mers.

        Args:
        :param cur_seq: A string of sequences
        :param kmer_int: An integer denoting the k-mer size
        :param out_link: File link for writing k-mers

        :returns: May modify cur_seq by deducing the k-mers that
            are written to the out_link file object.
        """
        if len(cur_seq) < kmer_int:
            return cur_seq
        elif len(cur_seq) == kmer_int:
            out_link.write(cur_seq + "\n")
            cur_seq = cur_seq[1:]
            return cur_seq
        else:
            for i in range(len(cur_seq) - kmer_int):
                j = i + kmer_int
                if j > len(cur_seq):
                    j = len(cur_seq)
                seq_out = cur_seq[i:j]
                if "N" not in seq_out:
                    out_link.write(seq_out + "\n")
            cur_seq = cur_seq[(i + 1):]
            return cur_seq

    def write_region(self):
        """Wrapper of GetKmers to get and write k-mers

        write_region method executes all the other methods
        of the GetKmers class for writing out all of the
        possible k-mers from chunk of a FASTA file.
        """
        out_dir, kmer, job_id = self.out_dir, self.kmer, self.job_id
        chrom = self.chrom
        start, end = self.start, self.end
        out_path = "{}/{}.{}.{}.kmer.gz".format(
            out_dir, chrom, job_id, kmer.replace("k", ""))
        out_link = gzip.open(out_path, "wb")
        fasta_fix = self.get_step_fasta()
        kmer_int = int(kmer.replace("k", ""))
        seq_ar = self.get_seq_ar(fasta_fix)
        cur_seq = ""
        for each_seq in seq_ar:
            cur_seq = cur_seq + each_seq
            cur_seq = self.write_kmers(cur_seq, kmer_int, out_link)
        print("Created all sequences for {}:{}-{}".format(
            chrom, start, end))
        out_link.close()


if __name__ == "__main__":
    chromsize_path, out_dir, kmer, job_id, chr_dir, idx_path = get_args()
    GetKmerObj = GetKmers(
        out_dir, kmer, job_id,
        chr_dir, chromsize_path,
        idx_path)
    GetKmerObj.write_region()
