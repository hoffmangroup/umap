from argparse import ArgumentParser
from random import choice, seed

NUCLEOTIDES_LIST = ["A", "T", "C", "G", "N"]
FASTA_LINE_SEQUENCE_LENGTH = 80


def get_args():
    parser = ArgumentParser(description="Generates a mock fasta file")

    parser.add_argument(
        "-n", "--name",
        help='"Chromosome" name. Also used for filename output.',
        required=True)

    parser.add_argument("-l", "--length",
                        help="Chromosome length", required=True)

    parser.add_argument("-s", "--seed",
                        help="Random seed")

    args = parser.parse_args()

    chromosome_name = args.name
    chromosome_length = int(args.length)

    # Set random seed
    seed(args.seed)

    return [chromosome_name, chromosome_length]


def generate_fasta(chromosome_name, length):
    current_length = 0
    current_sequence_line = ""

    output_filename = chromosome_name + ".fa"
    with open(output_filename, "w") as output_file:
        output_file.write(">" + chromosome_name + "\n")

        while current_length < length:
            current_sequence_line += choice(NUCLEOTIDES_LIST)

            if len(current_sequence_line) == FASTA_LINE_SEQUENCE_LENGTH:
                output_file.write(current_sequence_line + "\n")
                current_sequence_line = ""

            current_length += 1

        output_file.write(current_sequence_line + "\n")


def main():
    chromosome_name, length = get_args()
    generate_fasta(chromosome_name, length)


if __name__ == "__main__":
    main()
