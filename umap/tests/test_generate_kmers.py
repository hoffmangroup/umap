from pathlib import Path
import unittest

from util import TEST_DATA_DIRECTORY

from umap.get_kmers import create_unique_kmer_list

class TestGenerateKmers(unittest.TestCase):

    def setUp(self):
        # File looks like:
        # CGCANCAGA
        # GCANCGNCG
        self.fasta_file = open(Path(TEST_DATA_DIRECTORY / 'chr2.fa'), "rb")

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

    def tearDown(self):
        self.fasta_file.close()


if __name__ == '__main__':
    unittest.main()
