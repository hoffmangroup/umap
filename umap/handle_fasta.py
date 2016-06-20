from argparse import ArgumentParser
from datetime import datetime
import numpy as np
import os
from subprocess import call


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
        chr_path = self.chr_dir + "/chr1.fa"
        chr_link = open(chr_path, "w")
        if self.complement:
            chr_path_rc = chr_path.replace(".fa", "_RC.fa")
            chr_link_rc = open(chr_path_rc, "w")
            out_link_rc = open(self.out_path.replace(".fa", "_RC.fa"), "w")
        out_link = open(self.out_path, "w")
        ad_fin = ""
        chr_len = 0
        if self.complement:
            complement_ar = np.zeros(10)
        with open(self.in_path, "r") as fas_link:
            for each_line in fas_link:
                ad_line = each_line.rstrip()
                if ">chr" in ad_line:
                    if " " in ad_line:
                        ad_line = ad_line.split(" ")[0]
                    if "\t" in ad_line:
                        ad_line = ad_line.split("\t")[0]
                    chr_link.close()
                    print ad_fin
                    chrom = ad_line[1:]
                    chr_path = "{}/{}.fa".format(
                        self.chr_dir, chrom)
                    chr_link = open(chr_path, "w")
                    if self.complement:
                        if complement_ar[0] != 0:
                            self.write_reverse(complement_ar, chr_link_rc,
                                               out_link_rc, chr_len, chrom)
                        chr_len = chrsize_dict[chrom]
                        complement_ar = np.empty((chr_len), dtype="|S8")
                        start = 0
                        ad_line_rc = ad_line + "_RC"
                        chr_link_rc.close()
                        chr_link_rc = open(
                            chr_path.replace(".fa", "_RC.fa"), "w")
                        out_link_rc.write(ad_line_rc + "\n")
                        chr_link_rc.write(ad_line_rc + "\n")
                    print "%s started at %s" % (ad_line, str(datetime.now()))
                else:
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
                if self.complement:
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
                            raise ValueError("Failed to add sequences")
                        start = start + len(ad_line_rc)
        out_link.close()
        chr_link.close()
        if self.complement:
            self.write_reverse(complement_ar, chr_link_rc,
                               out_link_rc, chr_len, chrom)
            out_link_rc.close()
            chr_link_rc.close()
            call(" ".join(["cat", self.out_path.replace(".fa", "_RC.fa"),
                           ">>", self.out_path]), shell=True)
            call(["mv", self.out_path.replace(".fa", "_RC.fa"),
                 self.out_path.replace(".fa", "_Deprecated_RC.fa")])

    def write_reverse(complement_ar, chr_link_rc,
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
