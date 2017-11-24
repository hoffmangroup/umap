from argparse import ArgumentParser
from datetime import datetime
import numpy as np
import os
from subprocess import call
import warnings


class FastaHandler:
    def __init__(self, in_path, out_path, chrsize_path,
                 chr_dir, complement=True, conversion="None"):
        """Reverse complement the genome for Bismap

        This class generates the necessary FASTA files
        for running Umap or Bismap. If using Bismap,
        a reverse complemented set of chromosomes for
        the specified conversion will be produced.

        :param in_path: Path to genome .fasta file
        :param out_path: Path to output .fasta file
        :param chrsize_path: Path to 2-column chrom\tsize\n file
        :param chr_dir: Path to directory for
            writing chromosome fasta files
        :param complement: If reverse complemented chromosomes
            must be created
        :param conversion: One of 'None', 'C2T' or 'G2A'
        """
        self.in_path = in_path
        self.out_path = out_path
        self.chrsize_path = chrsize_path
        self.chr_dir = chr_dir
        self.complement = complement
        self.conversion = conversion

    def handle_fasta(self):
        """Wrapper function of FastaHandler

        Writes individual fasta files for each chromosome
        and converts nucleotides or reverse complements
        sequences if necessary (for Bismap feature)
        """
        chrsize_dict = self.get_chrsize_dict()
        print("Reverse complementation: {}".format(self.complement))
        print("Nucleotide conversion: {}".format(self.conversion))
        if not os.path.exists(self.chr_dir):
            os.makedirs(self.chr_dir)
        chr_path = self.chr_dir + "/chr1.fasta"
        chr_link = open(chr_path, "w")
        if self.complement:
            chr_path_rc = chr_path.replace(".fasta", "_RC.fasta")
            chr_link_rc = open(chr_path_rc, "w")
            out_path_rc = self.out_path.replace(
                ".fasta", "_RC.fasta")
            if out_path_rc == self.out_path:
                out_path_rc = self.out_path.replace(
                    ".fa", "_RC.fasta")
            out_link_rc = open(out_path_rc, "w")
        out_link = open(self.out_path, "w")
        ad_fin = ""
        chr_len = 0
        if self.complement:
            complement_ar = np.zeros(10)
        with open(self.in_path, "r") as fas_link:
            for each_line in fas_link:
                ad_line = each_line.rstrip()
                if ">" in ad_line:
                    if " " in ad_line:
                        ad_line = ad_line.split(" ")[0]
                    if "\t" in ad_line:
                        ad_line = ad_line.split("\t")[0]
                    chr_link.close()
                    print ad_fin
                    chrom = ad_line[1:]
                    CHROM_EXISTS = chrom in chrsize_dict.keys()
                    chr_path = "{}/{}.fasta".format(
                        self.chr_dir, chrom)
                    chr_link = open(chr_path, "w")
                    if self.complement and CHROM_EXISTS:
                        if complement_ar[0] != 0:
                            self.write_reverse(complement_ar, chr_link_rc,
                                               out_link_rc, chr_len, chrom)
                        chr_len = chrsize_dict.get(chrom, 0)
                        if chr_len == 0:
                            print(
                                "Failed to find length of {}".format(
                                    chrom))
                            warnings.warn(
                                "Failed to find length of {}".format(
                                    chrom))
                        complement_ar = np.empty((chr_len), dtype="|S8")
                        start = 0
                        ad_line_rc = ad_line + "_RC"
                        chr_link_rc.close()
                        chr_link_rc = open(
                            chr_path.replace(".fasta", "_RC.fasta"), "w")
                        out_link_rc.write(ad_line_rc + "\n")
                        chr_link_rc.write(ad_line_rc + "\n")
                    print(
                        "{} started at {}".format(
                            ad_line, str(datetime.now())))
                elif CHROM_EXISTS:
                    ad_line = ad_line.upper()
                    if self.complement:
                        ad_line_rc = ad_line.replace("A", "1")
                        ad_line_rc = ad_line_rc.replace("T", "A")
                        ad_line_rc = ad_line_rc.replace("C", "2")
                        ad_line_rc = ad_line_rc.replace("G", "C")
                        ad_line_rc = ad_line_rc.replace("1", "T")
                        ad_line_rc = ad_line_rc.replace("2", "G")
                        # ad_line_rc = ad_line_rc[::-1]
                    if self.conversion == "C2T":
                        ad_line = ad_line.replace("C", "T")
                        ad_line_rc = ad_line_rc.replace("C", "T")
                    elif self.conversion == "G2A":
                        ad_line = ad_line.replace("G", "A")
                        ad_line_rc = ad_line_rc.replace("G", "A")
                ad_fin = ad_line + "\n"
                out_link.write(ad_fin)
                chr_link.write(ad_fin)
                if self.complement and CHROM_EXISTS:
                    ad_fin_rc = ad_line_rc + "\n"
                    if "chr" not in ad_fin_rc:
                        end_pos = start + len(ad_line_rc)
                        try:
                            complement_ar[start:end_pos] = list(ad_line_rc)
                        except:
                            print(
                                "{} to {} with length of {}".format(
                                    start, end_pos,
                                    len(range(start, end_pos))))
                            print(ad_line_rc)
                            warnings.warn(
                                "Failed to add sequences for {}".format(chrom))
                        start = start + len(ad_line_rc)
        out_link.close()
        chr_link.close()
        if self.complement:
            if CHROM_EXISTS:
                self.write_reverse(complement_ar, chr_link_rc,
                                   out_link_rc, chr_len, chrom)
                chr_link_rc.close()
            out_link_rc.close()
            call(" ".join(["cat", out_path_rc,
                           ">>", self.out_path]), shell=True)
            call(["mv", out_path_rc,
                  out_path_rc.replace(".fa", "_Deprecated_RC.fa")])

    def write_reverse(self, complement_ar, chr_link_rc,
                      out_link_rc, chr_len, chrom):
        complement_ar = complement_ar[::-1]
        print("Length of {} was {}".format(chrom, chr_len))
        for start in range(0, chr_len, 50):
            end = start + 50
            if end > chr_len:
                end = chr_len
            ad_seq = "".join(complement_ar[start:end])
            out_link_rc.write(ad_seq + "\n")
            chr_link_rc.write(ad_seq + "\n")

    def get_chrsize_dict(self):
        chrsize_dict = {}
        with open(self.chrsize_path, "r") as chrsize_link:
            for chrsize_line in chrsize_link:
                chrom, size = chrsize_line.rstrip().split("\t")
                chrsize_dict[chrom] = int(size)
        return chrsize_dict


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Reverse complements a fasta genome")
    parser.add_argument(
        "in_fasta",
        help="FASTA format genome")
    parser.add_argument(
        "out_fasta",
        help="Output FASTA full path")
    parser.add_argument(
        "chrsize_path",
        help="Path to 2 column file with first column "
        "being the chromosome name and the second "
        "column being the chromosome size")
    parser.add_argument(
        "chr_dir",
        help="Path to directory for writing individual "
        "fasta files for each chromosome")
    parser.add_argument(
        "--Complement",
        action="store_true",
        help="Create a double genome with both "
        "+ and - strand information")
    parser.add_argument(
        "-Conversion",
        help="Specify C2T or G2A")
    args = parser.parse_args()
    FastaObj = FastaHandler(args.in_fasta, args.out_fasta, args.chrsize_path,
                            args.chr_dir, args.Complement, args.Conversion)
    FastaObj.handle_fasta()
