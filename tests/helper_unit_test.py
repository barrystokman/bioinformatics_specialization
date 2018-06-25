import unittest

import helper.dataset_reader as dataset_reader


class DatasetReaderTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_read_dataset(self):
        """
        Test main class ReadDataset, result should be an instance of that class
        """
        test_dataset = dataset_reader.ReadDataset('../tests/testdata/test_file.txt')
        self.assertIsInstance(test_dataset, dataset_reader.ReadDataset)

    def test_read_pattern_count(self):
        """
        Test PatternCountDataset subclass
        """
        test_dataset = dataset_reader.PatternCountDataset('../tests/testdata/' +
                                                          'pattern_count.txt')
        self.assertEqual(test_dataset.get_text(), 'GCGCG')
        self.assertEqual(test_dataset.get_pattern(), 'GCG')
        self.assertEqual(test_dataset.get_expected_result(), 2)

    def test_read_pattern_count_challenge(self):
        """
        """
        test_dataset = dataset_reader.PatternCountDataset('../tests/testdata/' +
                                                          'pattern_count_challenge.txt')
        self.assertEqual(test_dataset.get_text_challenge(), 'GCGCG')
        self.assertEqual(test_dataset.get_pattern_challenge(), 'GCG')

    def test_read_frequent_words(self):
        """
        """
        test_dataset = dataset_reader.FrequentWordsDataset('../tests/testdata/' +
                                                           'frequent_words.txt')
        self.assertEqual(test_dataset.get_text(), 'ACGTTGCATGTCGCATGATGCATGAGAGCT')
        self.assertEqual(test_dataset.get_k(), 4)
        self.assertEqual(test_dataset.get_expected_result(), 'CATG GCAT')

    def test_read_frequent_words_challenge(self):
        """
        """
        test_dataset = dataset_reader.FrequentWordsDataset('../tests/testdata/' +
                                                           'frequent_words_challenge.txt')
        self.assertEqual(test_dataset.get_text_challenge(), 'ACGTTGCATGTCGCATGATGCATGAGAGCT')
        self.assertEqual(test_dataset.get_k_challenge(), 4)

    def test_read_reverse_complement(self):
        """
        """
        test_dataset = dataset_reader.ReverseComplementDataset('../tests/testdata/' +
                                                               'reverse_complement.txt')
        self.assertEqual(test_dataset.get_pattern(), 'AAAACCCGGT')
        self.assertEqual(test_dataset.get_expected_result(), 'ACCGGGTTTT')

    def test_read_reverse_complement_challenge(self):
        """
        """
        test_dataset = dataset_reader.ReverseComplementDataset('../tests/testdata/' +
                                                               'reverse_complement_challenge.txt')
        self.assertEqual(test_dataset.get_pattern_challenge(), 'AAAACCCGGT')

    def test_read_pattern_matching(self):
        """
        """
        test_dataset = dataset_reader.PatternMatchingDataset('../tests/testdata/' +
                                                             'pattern_matching.txt')
        self.assertEqual(test_dataset.get_pattern(), 'ATAT')
        self.assertEqual(test_dataset.get_genome(), 'GATATATGCATATACTT')
        self.assertEqual(test_dataset.get_expected_result(), '1 3 9')

    def test_read_pattern_matching_challenge(self):
        """
        """
        test_dataset = dataset_reader.PatternMatchingDataset('../tests/testdata/' +
                                                             'pattern_matching_challenge.txt')
        self.assertEqual(test_dataset.get_pattern_challenge(), 'ATAT')
        self.assertEqual(test_dataset.get_genome_challenge(), 'GATATATGCATATACTT')

    def test_read_clump_finding(self):
        """
        """
        test_dataset = dataset_reader.ClumpFindingDataset('../tests/testdata/' +
                                                          'clump_finding.txt')
        self.assertEqual(test_dataset.get_genome(),
                         'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACAT' +
                         'TGTAA')
        self.assertEqual(test_dataset.get_k(), 5)
        self.assertEqual(test_dataset.get_l(), 50)
        self.assertEqual(test_dataset.get_t(), 4)
        self.assertEqual(test_dataset.get_expected_result(), 'CGACA GAAGA')

    def test_read_clump_finding_challenge(self):
        """
        """
        test_dataset = dataset_reader.ClumpFindingDataset('../tests/testdata/' +
                                                          'clump_finding_challenge.txt')
        self.assertEqual(test_dataset.get_genome_challenge(),
                         'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACAT' +
                         'TGTAA')
        self.assertEqual(test_dataset.get_k_challenge(), 5)
        self.assertEqual(test_dataset.get_l_challenge(), 50)
        self.assertEqual(test_dataset.get_t_challenge(), 4)

    def test_read_computing_frequencies(self):
        """
        """
        test_dataset = dataset_reader.ComputingFrequenciesDataset('../tests/testdata/' +
                                                                  'frequency_array.txt')
        self.assertEqual(test_dataset.get_text(), 'ACGCGGCTCTGAAA')
        self.assertEqual(test_dataset.get_k(), 2)
        self.assertEqual(test_dataset.get_expected_result(), '2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0')

    def test_read_computing_frequencies_challenge(self):
        """
        """
        test_dataset = dataset_reader.ComputingFrequenciesDataset('../tests/testdata/' +
                                                                  'frequency_array_challenge.txt')
        self.assertEqual(test_dataset.get_text_challenge(), 'ACGCGGCTCTGAAA')
        self.assertEqual(test_dataset.get_k_challenge(), 2)

    def test_read_pattern_to_number(self):
        """
        """
        test_dataset = dataset_reader.PatternToNumberDataset('../tests/testdata/' +
                                                             'pattern_to_number.txt')
        self.assertEqual(test_dataset.get_pattern(), 'CTTCTCACGTACAACAAAATC')
        self.assertEqual(test_dataset.get_expected_result(), 2161555804173)

    def test_read_pattern_to_number_challenge(self):
        """
        """
        test_dataset = dataset_reader.PatternToNumberDataset('../tests/testdata/' +
                                                             'pattern_to_number_challenge.txt')
        self.assertEqual(test_dataset.get_pattern_challenge(), 'ACATGCTTGTTTTGCTTCG')

    def test_read_number_to_pattern(self):
        """
        """
        test_dataset = dataset_reader.NumberToPatternDataset('../tests/testdata/' +
                                                             'number_to_pattern.txt')
        self.assertEqual(test_dataset.get_number(), 5353)
        self.assertEqual(test_dataset.get_k(), 7)
        self.assertEqual(test_dataset.get_expected_result(), 'CCATGGC')

    def test_read_number_to_pattern_challenge(self):
        """
        """
        test_dataset = dataset_reader.NumberToPatternDataset('../tests/testdata/' +
                                                             'number_to_pattern_challenge.txt')
        self.assertEqual(test_dataset.get_number_challenge(), 8718)
        self.assertEqual(test_dataset.get_k_challenge(), 8)

    def test_read_minimum_skew(self):
        """
        """
        test_dataset = dataset_reader.MinimumSkewDataset('../tests/testdata/' +
                                                         'minimum_skew.txt')
        self.assertEqual(test_dataset.get_genome(),
                         'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')
        self.assertEqual(test_dataset.get_expected_result(), '11 24')

    def test_read_minimum_skew_challenge(self):
        """
        """
        test_dataset = dataset_reader.MinimumSkewDataset('../tests/testdata/' +
                                                         'minimum_skew_challenge.txt')
        self.assertEqual(test_dataset.get_genome_challenge(),
                         'AGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGA')

    def test_read_hamming_distance(self):
        """
        """
        test_dataset = dataset_reader.HammingDistanceDataset('../tests/testdata/' +
                                                             'hamming_distance.txt')
        self.assertEqual(test_dataset.get_string1(), 'GGGCCGTTGGT')
        self.assertEqual(test_dataset.get_string2(), 'GGACCGTTGAC')
        self.assertEqual(test_dataset.get_expected_result(), 3)

    def test_read_hamming_distance_challenge(self):
        """
        """
        test_dataset = dataset_reader.HammingDistanceDataset('../tests/testdata/' +
                                                             'hamming_distance_challenge.txt')
        self.assertEqual(test_dataset.get_string1_challenge(), 'GGGCCGTTGGT')
        self.assertEqual(test_dataset.get_string2_challenge(), 'GGACCGTTGAC')

    def test_read_approx_matching(self):
        """
        """
        test_dataset = dataset_reader.ApproxMatchDataset('../tests/testdata/' +
                                                         'approximate_match.txt')
        self.assertEqual(test_dataset.get_pattern(), 'ATTCTGGA')
        self.assertEqual(test_dataset.get_text(),
                         'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT')
        self.assertEqual(test_dataset.get_d(), 3)
        self.assertEqual(test_dataset.get_expected_result(), '6 7 26 27')

    def test_read_approx_matching_challenge(self):
        """
        """
        test_dataset = dataset_reader.ApproxMatchDataset('../tests/testdata/' +
                                                         'approximate_match_challenge.txt')
        self.assertEqual(test_dataset.get_pattern_challenge(), 'ATTCTGGA')
        self.assertEqual(test_dataset.get_text_challenge(),
                         'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT')
        self.assertEqual(test_dataset.get_d_challenge(), 3)


if __name__ == "__main__":
    unittest.main()
