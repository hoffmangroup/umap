from argparse import ArgumentParser
from pathlib import Path
from sys import stdout
from typing import BinaryIO

from umap.util import optional_gzip_open


DEFAULT_KMER_SIZE = 100
FASTA_FILE_IGNORE_DELIMITERS = (b'>', b';')


def get_args():
    parser = ArgumentParser(
        description="Prints all unique kmers per line from all sequences in a"
        " given FASTA file to standard output")

    parser.add_argument(
        "--kmer-length", "-k",
        default=DEFAULT_KMER_SIZE,
        type=int,
        help="Length of kmers to generate. Default is {}".format(
            DEFAULT_KMER_SIZE))

    parser.add_argument(
        "--fasta-file", "-f",
        help="Filename of reference fasta file for kmer generation")

    args = parser.parse_args()

    kmer_length = args.kmer_length
    fasta_filename = args.fasta_file

    out_list = [kmer_length, fasta_filename]

    return out_list


def create_unique_kmer_list(fasta_file: BinaryIO,
                            kmer_length: int) -> set:
    # Create an empty set
    kmer_set = set()

    # Create a line-based buffer for the current set of nucleotides
    nucleotide_buffer = b''
    kmer_current_index = 0  # NB: index into the nucleotide_buffer

    # For each line in the fasta
    for line in fasta_file:
        # If it does not start with '>' or ';'
        if not line.startswith(FASTA_FILE_IGNORE_DELIMITERS):
            # Slide a kmer window across the current nucleotide list
            nucleotide_buffer += line.rstrip()

            # Record the current kmer window indices
            kmer_start_index = kmer_current_index
            kmer_end_index = kmer_current_index + kmer_length

            # While there are enough nucleotides in the buffer
            while kmer_end_index <= len(nucleotide_buffer):
                # Add the current sliding window kmer to the set
                kmer_set.add(
                    nucleotide_buffer[kmer_current_index:kmer_end_index])

                # Advance the current window slice indexes by one
                kmer_current_index += 1
                kmer_end_index += 1

            # After removing all possible kmers from the buffer
            # Calculate amount the sliding window moved by
            kmer_window_shift = kmer_current_index - kmer_start_index

            # Truncate the beginning of nucleotide buffer by the shift amount
            nucleotide_buffer = nucleotide_buffer[kmer_window_shift:]
            # And move back the current window index by the shift amount
            kmer_current_index -= kmer_window_shift
        # Otherwise assume that either of the delimiters are indicators of a
        # new sequence, notably including comments
        else:
            # Empty the nucleotide buffer
            nucleotide_buffer = b''
            # Reset the current window index
            kmer_current_index = 0

    return kmer_set


def print_kmers(kmer_length: int,
                fasta_filename: Path):
    # NB: We open the file in binary mode to get read-only bytes
    with optional_gzip_open(fasta_filename, "rb") as fasta_file:
        # Create a unique list of kmers for this file
        # NB: The fasta_file type is BinaryIO but mypy can't infer
        # from the "rb" parameter to optional_gzip_open
        unique_kmer_list = create_unique_kmer_list(fasta_file,  # type: ignore
                                                   kmer_length)

        for elem in unique_kmer_list:
            stdout.buffer.write(elem + b'\n')


def main():
    kmer_length, fasta_filename = get_args()

    print_kmers(kmer_length,
                Path(fasta_filename))


if __name__ == "__main__":
    main()
