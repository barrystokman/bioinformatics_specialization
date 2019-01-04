import unittest
import cli.biocli_module as biocli
from click.testing import CliRunner


class BiocliTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_biocli(self):

        runner = CliRunner()
        result = runner.invoke(biocli.biocli)
        assert result.exit_code == 0

    def test_pattern_count_no_challenge(self):

        test_dataset = '../tests/testdata/pattern_count.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['pattern-count', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\n2\n'

    def test_pattern_count_challenge(self):

        test_dataset = '../tests/testdata/pattern_count_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'pattern-count', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\n2\n'

    def test_frequent_words_no_challenge(self):

        test_dataset = '../tests/testdata/frequent_words.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['frequent-words', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\nCATG GCAT\n'

    def test_frequent_words_challenge(self):

        test_dataset = '../tests/testdata/frequent_words_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'frequent-words', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\nCATG GCAT\n'

    def test_reverse_complement_no_challenge(self):

        test_dataset = '../tests/testdata/reverse_complement.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['reverse-complement', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\nACCGGGTTTT\n'

    def test_reverse_complement_challenge(self):

        test_dataset = '../tests/testdata/reverse_complement_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'reverse-complement', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\nACCGGGTTTT\n'

    def test_pattern_matching_no_challenge(self):

        test_dataset = '../tests/testdata/pattern_matching.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['pattern-matching', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'THIS FUNCTION IS ONLY AVAILABLE IN CODE CHALLENGE MODE!\n'

    def test_pattern_matching_challenge(self):

        test_dataset = '../tests/testdata/pattern_matching_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'pattern-matching', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\n1 3 9\n'

    def test_clump_finding_no_challenge(self):

        test_dataset = '../tests/testdata/clump_finding.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['clump-finding', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\nCGACA GAAGA\n'

    def test_clump_finding_challenge(self):

        test_dataset = '../tests/testdata/clump_finding_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'clump-finding', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\nCGACA GAAGA\n'

    def test_computing_frequencies_no_challenge(self):

        test_dataset = '../tests/testdata/frequency_array.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['computing-frequencies', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\n2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0\n'

    def test_computing_frequencies_challenge(self):

        test_dataset = '../tests/testdata/frequency_array_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'computing-frequencies', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\n2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0\n'

    def test_pattern_to_number_no_challenge(self):

        test_dataset = '../tests/testdata/pattern_to_number.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['pattern-to-number', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\n2161555804173\n'

    def test_pattern_to_number_challenge(self):

        test_dataset = '../tests/testdata/pattern_to_number_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'pattern-to-number', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\n21071133174\n'

    def test_number_to_pattern_no_challenge(self):

        test_dataset = '../tests/testdata/number_to_pattern.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['number-to-pattern', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\nCCATGGC\n'

    def test_number_to_pattern_challenge(self):

        test_dataset = '../tests/testdata/number_to_pattern_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'number-to-pattern', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\nAGAGAATG\n'

    def test_minimum_skew_no_challenge(self):

        test_dataset = '../tests/testdata/minimum_skew.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['minimum-skew', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\n11 24\n'

    def test_minimum_skew_challenge(self):

        test_dataset = '../tests/testdata/pattern_count_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'minimum-skew', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\n0 2 4\n'

    def test_hamming_distance_no_challenge(self):

        test_dataset = '../tests/testdata/hamming_distance.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['hamming-distance', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\n3\n'

    def test_hamming_distance_challenge(self):

        test_dataset = '../tests/testdata/hamming_distance_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'hamming-distance', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\n3\n'

    def test_approx_matching_no_challenge(self):

        test_dataset = '../tests/testdata/approximate_match.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['approx-matching', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\n6 7 26 27\n'

    def test_approx_matching_challenge(self):

        test_dataset = '../tests/testdata/approximate_match_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'approx-matching', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\n6 7 26 27\n'

    def test_approx_count_no_challenge(self):

        test_dataset = '../tests/testdata/approximate_pattern_count.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['approx-count', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\n4\n'

    def test_approx_count_challenge(self):

        test_dataset = '../tests/testdata/approximate_pattern_count_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'approx-count', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\n4\n'

    def test_neighbors_no_challenge(self):

        test_dataset = '../tests/testdata/neighbors.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['neighbors', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\nAAG ACA ACC ACG ACT AGG ATG CCG GCG TCG\n'

    def test_neighbors_challenge(self):

        test_dataset = '../tests/testdata/neighbors_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'neighbors', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\n\nAAG\nACA\nACC\nACG\nACT\nAGG\nATG\nCCG\nGCG\nTCG\n'

    def test_frequent_words_mismatches_no_challenge(self):

        test_dataset = '../tests/testdata/frequent_words_mismatch.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['frequent-words-mismatches', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\nGCACACAGAC GCGCACACAC\n'

    def test_frequent_words_mismatches_challenge(self):

        test_dataset = '../tests/testdata/frequent_words_mismatch_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'frequent-words-mismatches', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\nAAAAA\n'

    def test_frequent_words_mismatches_and_reverse_complement_no_challenge(self):

        test_dataset = '../tests/testdata/frequent_words_mismatch_complements.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['frequent-words-mismatches-and-reverse-complement', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of this function is:\nAGCGCCGCT AGCGGCGCT\n'

    def test_frequent_words_mismatches_and_reverse_complement_challenge(self):

        test_dataset = '../tests/testdata/frequent_words_mismatch_complements_challenge.txt'
        runner = CliRunner()
        result = runner.invoke(biocli.biocli,  ['--code-challenge', 'frequent-words-mismatches-and-reverse-complement', test_dataset])
        assert result.exit_code == 0
        assert result.output == 'The result of the Coding Challenge is:\nCCCCCC GGGGGG\n'
