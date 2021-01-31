import unittest
import os
from app.utils import fasta_to_oneline


class MyTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_fasta_to_oneline(self):
        fasta_to_oneline("test_chr1.fasta", "test_chr1.fasta.output")

        with open("test_chr1.fasta.output", 'r') as file:
            file_data = file.read()

        with open("test_chr1_oneLine.txt", 'r') as file:
            test_file_data = file.read()

        self.assertEqual(file_data, test_file_data)

    def tearDown(self):
        if os.path.exists("test_chr1.fasta.output"):
            os.remove("test_chr1.fasta.output")


if __name__ == '__main__':
    unittest.main()
