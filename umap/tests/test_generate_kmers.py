from pathlib import Path
import unittest

from util import TEST_DATA_DIRECTORY

from umap.get_kmers import create_unique_kmer_list


class TestGenerateKmers(unittest.TestCase):

    def setUp(self):
        # File looks like:
        # CGCANCAGA
        # GCANCGNCG
        self.fasta_file = open(
                Path(TEST_DATA_DIRECTORY / 'chr2.fa'), "rb")
        self.genome_fasta_file = open(
                Path(TEST_DATA_DIRECTORY / 'genome.fa'), "rb")

    def test_unique_kmer_list(self):
        kmer_3_set = create_unique_kmer_list(self.fasta_file, 3)
        self.assertEqual(
            set([b'CGC', b'GCA', b'CAN', b'ANC', b'NCA', b'CAG', b'AGA',
                 b'GAG', b'AGC',
                 b'GCA', b'CAN', b'ANC', b'NCG', b'CGN', b'GNC', b'NCG',
                 ]),
            kmer_3_set
        )

        self.fasta_file.seek(0)
        kmer_12_set = create_unique_kmer_list(self.fasta_file, 12)
        self.assertEqual(
            set([b'CGCANCAGAGCA', b'GCANCAGAGCAN', b'CANCAGAGCANC',
                 b'ANCAGAGCANCG', b'NCAGAGCANCGN', b'CAGAGCANCGNC',
                 b'AGAGCANCGNCG']),
            kmer_12_set
        )

    def test_unique_kmers_mulitple_sequences(self):
        kmer_3_set = create_unique_kmer_list(self.genome_fasta_file, 3)
        self.assertEqual(
            set([  # chr1
                 b'AAA', b'AAT', b'ATT', b'TTT', b'TTA', b'TAT',
                 b'ATC', b'TCG', b'CGA', b'GAA',
                   # chr2
                 b'AAA', b'AAG', b'AGG',
                 b'GGG', b'GGC', b'GCC',
                 b'CCC'
                 ]),
            kmer_3_set
        )

        self.genome_fasta_file.seek(0)
        kmer_10_set = create_unique_kmer_list(self.genome_fasta_file, 10)
        self.assertEqual(
            set([  # chr1
                 b'AAAAATTTTT', b'AAAATTTTTA', b'AAATTTTTAT', b'AATTTTTATC',
                 b'ATTTTTATCG', b'TTTTTATCGA', b'TTTTATCGAA', b'TTTATCGAAT',
                 b'TTATCGAATC', b'TATCGAATCG', b'ATCGAATCGA',
                   # chr2
                 b'AAAAAAAAAA', b'AAAAAAAAAG', b'AAAAAAAAGG', b'AAAAAAAGGG',
                 b'AAAAAAGGGG', b'AAAAAGGGGG', b'AAAAGGGGGG', b'AAAGGGGGGG',
                 b'AAGGGGGGGG', b'AGGGGGGGGG', b'GGGGGGGGGG', b'GGGGGGGGGC',
                 b'GGGGGGGGCC', b'GGGGGGGCCC', b'GGGGGGCCCC', b'GGGGGCCCCC',
                 b'GGGGCCCCCC', b'GGGCCCCCCC', b'GGCCCCCCCC', b'GCCCCCCCCC',
                 b'CCCCCCCCCC',
                ]),
            kmer_10_set
        )

    def tearDown(self):
        self.fasta_file.close()
        self.genome_fasta_file.close()


if __name__ == '__main__':
    unittest.main()
