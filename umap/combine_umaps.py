from argparse import ArgumentParser
import gzip
import numpy as np
import os


def get_args():
    parser = ArgumentParser(
        description="Combines mappability uint8 vectors of "
        "several kmers into 1 uint8 vector per chromosome. "
        "It requires a directory with subfolders names as "
        "k<read length>. This script requires a number to infer "
        "chromosome. If not specifying -job_id, it will "
        "identify the chromosome using -var_id environmental "
        "varibale which by default is SGE_TASK_ID.")
    parser.add_argument(
        "kmer_dir",
        help="Directory with subfolders "
        "named as k<read length>)")
    parser.add_argument(
        "chrsize_path",
        help="Path to 2 column tsv file with first column "
        "as chromosome and second column as its size. Will "
        "be used to identify order of the chromosomes.")
    parser.add_argument(
        "-out_dir",
        default="infer",
        help="If not specified, a subfolder "
        "will be created in kmer_dir names as "
        "globalmap_k<smallestkmer>tok<largestkmer>")
    parser.add_argument(
        "-job_id",
        default=0,
        type=int,
        help="1-based index for finding chromosome from -chrsize_path. "
        "If not specified, it will user -var_id to "
        "infer the chromosome for combining mappabilitiy of "
        "different kmers.")
    parser.add_argument(
        "-var_id",
        default="SGE_TASK_ID",
        help="If -job_id is not specified, job_id will be inferred "
        "from environmental variable -var_id.")
    parser.add_argument(
        "-kmer_dir_2",
        default="NA",
        help="Specify to merge kmers of two different directories "
        "by logical operation AND.")
    args = parser.parse_args()
    job_id = args.job_id
    out_dir = args.out_dir
    kmers = [each_kmer for each_kmer in next(os.walk(args.kmer_dir))[1]
             if "k" == each_kmer[0]]
    if len(kmers) < 1:
        raise ValueError("{} lacks kmer folders".format(args.kmer_dir))
    if args.job_id == 0:
        job_id = int(os.environ[args.var_id]) - 1
    else:
        job_id = job_id - 1
    if args.out_dir == "infer":
        kmer_ints = [int(each_kmer.replace("k", "")) for each_kmer in kmers]
        out_dir = "{}/globalmap_k{}tok{}".format(
            args.kmer_dir, min(kmer_ints), max(kmer_ints))
    try:
        out_list = [args.kmer_dir, args.chrsize_path,
                    out_dir, job_id, args.kmer_dir_2]
    except:
        raise ValueError(
            "chrsize_path or job_id were invalid.")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    return out_list


class CombineUmaps:
    def __init__(self, kmer_dir, chrsize_path,
                 out_dir, job_id, kmer_dir_2):
        """Combine <chr>.<kmer>.uint8.unique.gz files

        The kmer_dir has subfolders named as k<integer> with
        <chr>.<kmer>.uint8.unique.gz files inside. This script
        will use the job_id to use a particular chromosome and
        merge the uint8.gz files across all of the different
        kmers into one uint8.gz. Note: All methods would
        run consequently when creating the class instance.

        Args:
            kmer_dir: Directory with k<integer> subfolders
            chrsize_path: 2-column tsv file of <chrom>\t<size>\n...
            out_dir: Output files will be saved to this folder
            job_id: 0-based index to find chromosome from chrsize_path
            kmer_dir_2: If using Bismap and want to merge C2T and G2A data

        Returns
            Saves the output to outdir/<chrom>.uint8.unique.gz files
        """
        self.kmer_dir = kmer_dir
        self.kmer_dir_2 = kmer_dir_2
        self.chrsize_path = chrsize_path
        self.out_dir = out_dir
        self.job_id = job_id
        self.chrom, self.size = self.get_chrom_size()
        self.kmers = [each_kmer for each_kmer in
                      next(os.walk(self.kmer_dir))[1]
                      if "k" == each_kmer[0]]
        combined_ar = self.combine_uints()
        self.write_ar(combined_ar)

    def get_chrom_size(self):
        """Finds chrom and size from self.chrsize_path"""
        with open(self.chrsize_path, "r") as chrsize_link:
            ind_chr = 0
            for chrsize_line in chrsize_link:
                if ind_chr == self.job_id:
                    chromosome = chrsize_line.rstrip().split("\t")[0]
                    size = int(chrsize_line.rstrip().split("\t")[1])
                ind_chr = ind_chr + 1
        return chromosome, size

    def combine_uints(self):
        """Merged different kmer arrays of one chromosome
        """
        MergeKmers = False
        if os.path.exists(self.kmer_dir_2):
            print(
                "Limit mappability to regions {} {} and {}".format(
                    "that are unique in both",
                    kmer_dir, kmer_dir_2))
            MergeKmers = True
        combined_ar = np.zeros(self.size, dtype=np.uint8)
        kmer_nums = [int(kmer.replace("k", "")) for kmer in self.kmers]
        kmer_nums.sort()
        for kmer_num in kmer_nums:
            kmer = "k{}".format(kmer_num)
            full_kmer_path = "{}/{}/{}.{}.uint8.unique.gz".format(
                self.kmer_dir, kmer, self.chrom, kmer)
            if not os.path.exists(full_kmer_path):
                raise ValueError("{} does not exist".format(full_kmer_path))
            kmer_link = gzip.open(full_kmer_path, "rb")
            kmer_ar = np.frombuffer(kmer_link.read(), dtype=np.uint8)
            if len(kmer_ar) != self.size:
                kmer_ar = np.append(kmer_ar, np.zeros(kmer_num))
            kmer_link.close()
            index_comb_0 = combined_ar == 0
            index_kmer = kmer_ar != 0
            index_adkmer = np.logical_and(index_comb_0, index_kmer)
            if MergeKmers:
                full_kmer_path_2 = "{}/{}/{}.{}.uint8.unique.gz".format(
                    self.kmer_dir_2, kmer, self.chrom, kmer)
                kmer_link_2 = gzip.open(full_kmer_path_2, "rb")
                kmer_ar_2 = np.frombuffer(kmer_link_2.read(), dtype=np.uint8)
                if len(kmer_ar_2) != self.size:
                    kmer_ar_2 = np.append(kmer_ar_2, np.zeros(kmer_num))
                kmer_link_2.close()
                index_kmer_2 = kmer_ar_2 != 0
                index_adkmer_2 = np.logical_and(index_comb_0, index_kmer_2)
                index_adkmer = np.logical_and(index_adkmer, index_adkmer_2)
            combined_ar[index_adkmer] = kmer_num
            print(
                "Added information of {} for {}".format(
                    kmer, self.chrom))
        return combined_ar

    def write_ar(self, combined_ar):
        """Writes merged array to unsigned 8bit integer file.

        Used self.out_dir and self.chrom.

        Args:
            combined_ar: Can be any numpy array
        """
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        out_path = "{}/{}.uint8.unique.gz".format(
            self.out_dir, self.chrom)
        out_link = gzip.open(out_path, "wb")
        out_link.write(combined_ar.tobytes())
        out_link.close()


if __name__ == "__main__":
    kmer_dir, chrsize_path, out_dir, job_id, kmer_dir_2 = get_args()
    CombineUmaps(kmer_dir, chrsize_path, out_dir, job_id, kmer_dir_2)
