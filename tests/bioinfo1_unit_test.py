import sys
import unittest

sys.path.insert(0, '/home/barry/proj/bioinformatics/bioinformatics_1')

import bioinformatics_1 as bioinfo1


class PatternCountTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_pattern_count_sample(self):
        """
        test off-by-one errors at start and end, and overlap in the middle
        """
        test_text = "GCGCG"
        test_pattern = "GCG"
        test_count = 2
        self.assertEqual(bioinfo1.pattern_count(test_text, test_pattern), test_count)

    def test_pattern_count_1(self):
        """
        tests simple counting
        """
        test_text = "ACGTACGTACGT"
        test_pattern = "CG"
        test_count = 3
        self.assertEqual(bioinfo1.pattern_count(test_text, test_pattern), test_count)

    def test_pattern_count_2(self):
        """
        checks for off-by-one errors at the beginning of the text
        no overlap, no pattern at the end
        """
        test_text = "AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAATTTTATTGAC" +\
                    "TTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAATAATTACAGAGTA" +\
                    "CACAACATCCAT"
        test_pattern = "AAA"
        test_count = 4
        self.assertEqual(bioinfo1.pattern_count(test_text, test_pattern), test_count)

    def test_pattern_count_3(self):
        """
        checks for off-by-one errors at the end of the text
        no overlap, no pattern at the beginning

        """
        test_text = "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTG" +\
                    "CATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGC" +\
                    "CGACTTT"
        test_pattern = "TTT"
        test_count = 4
        self.assertEqual(bioinfo1.pattern_count(test_text, test_pattern), test_count)

    def test_pattern_count_4(self):
        """
        tests counting occurences of Reverse Complement of Pattern
        this should not happen in this function, only perfect matches should be counted
        """
        test_text = "GGACTTACTGACGTACG"
        test_pattern = "ACT"
        test_count = 2
        self.assertEqual(bioinfo1.pattern_count(test_text, test_pattern), test_count)

    def test_pattern_count_5(self):
        """
        tests pattern overlap
        """
        test_text = "ATCCGATCCCATGCCCATG"
        test_pattern = "CC"
        test_count = 5
        self.assertEqual(bioinfo1.pattern_count(test_text, test_pattern), test_count)

    def test_pattern_count_6(self):
        """
        tests function on full data set
        """
        test_text = "CTGTTTTTGATCCATGATATGTTATCTCTCCGTCATCAGAAGAACAGTGACGGATCGCCCTCTC" +\
                    "TCTTGGTCAGGCGACCGTTTGCCATAATGCCCATGCTTTCCAGCCAGCTCTCAAACTCCGGTGA" +\
                    "CTCGCGCAGGTTGAGTA"
        test_pattern = "CTC"
        test_count = 9
        self.assertEqual(bioinfo1.pattern_count(test_text, test_pattern), test_count)


class FrequentWordsTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_frequent_words_sample(self):
        """
        """
        test_text = "GCGCGCGTTGCATGTCGCATGATGCATGAGAGCT"
        test_k = 4
        test_output = ['CATG', 'GCAT']
        self.assertEqual(sorted(bioinfo1.frequent_words(test_text, test_k)), test_output)

    def test_frequent_words_1(self):
        """
        tests if the first k-mer in string is counted
        if the first k-mer is not counted, the output is ACT CAC CCA CTT GGT
        """
        test_text = "TGGTAGCGACGTTGGTCCCGCCGCTTGAGAATCTGGATGAACATAAGCTCCCACTTGGCTTATT" +\
                    "CAGAGAACTGGTCAACACTTGTCTCTCCCAGCCAGGTCTGACCACCGGGCAACTTTTAGAGCAC" +\
                    "TATCGTGGTACAAATAATGCTGCCAC"
        test_k = 3
        test_output = ['TGG']
        self.assertEqual(sorted(bioinfo1.frequent_words(test_text, test_k)), test_output)

    def test_frequent_words_2(self):
        """
        tests if the last k-mer in string is counted
        if the last k-mer is not counted, the output is AACG AATA ACAA CAAC CTGG CTTT
        TTGC TTTG
        """
        test_text = "CAGTGGCAGATGACATTTTGCTGGTCGACTGGTTACAACAACGCCTGGGGCTTTTGAGCAACGA" +\
                    "GACTTTTCAATGTTGCACCGTTTGCTGCATGATATTGAAAACAATATCACCAAATAAATAACGC" +\
                    "CTTAGTAAGTAGCTTTT"
        test_k = 4
        test_output = ['TTTT']
        self.assertEqual(sorted(bioinfo1.frequent_words(test_text, test_k)), test_output)

    def test_frequent_words_3(self):
        """
        tests function for handling overlapping k-mers
        """
        test_text = "ATACAATTACAGTCTGGAACCGGATGAACTGGCCGCAGGTTAACAACAGAGTTGCCAGGCACTG" +\
                    "CCGCTGACCAGCAACAACAACAATGACTTTGACGCGAAGGGGATGGCATGAGCGAACTGATCGT" +\
                    "CAGCCGTCAGCAACGAGTATTGTTGCTGACCCTTAACAATCCCGCCGCACGTAATGCGCTAACT" +\
                    "AATGCCCTGCTG"
        test_k = 5
        test_output = ['AACAA']
        self.assertEqual(sorted(bioinfo1.frequent_words(test_text, test_k)), test_output)

    def test_frequent_words_4(self):
        """
        tests if function handles more than one most frequent k-mer correctly
        """
        test_text = "CCAGCGGGGGTTGATGCTCTGGGGGTCACAAGATTGCATTTTTATGGGGTTGCAAAAATGTTTT" +\
                    "TTACGGCAGATTCATTTAAAATGCCCACTGGCTGGAGACATAGCCCGGATGCGCGTCTTTTACA" +\
                    "ACGTATTGCGGGGTAAAATCGTAGATGTTTTAAAATAGGCGTAAC"
        test_k = 5
        test_output = ['AAAAT', 'GGGGT', 'TTTTA']
        self.assertEqual(sorted(bioinfo1.frequent_words(test_text, test_k)), test_output)


class FasterFrequentWordsTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_faster_frequent_words_sample(self):
        """
        """
        test_text = "GCGCGCGTTGCATGTCGCATGATGCATGAGAGCT"
        test_k = 4
        test_output = ['CATG', 'GCAT']
        self.assertEqual(sorted(bioinfo1.faster_frequent_words(test_text, test_k)), test_output)

    def test_faster_frequent_words_1(self):
        """
        tests if the first k-mer in string is counted
        if the first k-mer is not counted, the output is ACT CAC CCA CTT GGT
        """
        test_text = "TGGTAGCGACGTTGGTCCCGCCGCTTGAGAATCTGGATGAACATAAGCTCCCACTTGGCTTATT" +\
                    "CAGAGAACTGGTCAACACTTGTCTCTCCCAGCCAGGTCTGACCACCGGGCAACTTTTAGAGCAC" +\
                    "TATCGTGGTACAAATAATGCTGCCAC"
        test_k = 3
        test_output = ['TGG']
        self.assertEqual(sorted(bioinfo1.faster_frequent_words(test_text, test_k)), test_output)

    def test_faster_frequent_words_2(self):
        """
        tests if the last k-mer in string is counted
        if the last k-mer is not counted, the output is AACG AATA ACAA CAAC CTGG CTTT
        TTGC TTTG
        """
        test_text = "CAGTGGCAGATGACATTTTGCTGGTCGACTGGTTACAACAACGCCTGGGGCTTTTGAGCAACGA" +\
                    "GACTTTTCAATGTTGCACCGTTTGCTGCATGATATTGAAAACAATATCACCAAATAAATAACGC" +\
                    "CTTAGTAAGTAGCTTTT"
        test_k = 4
        test_output = ['TTTT']
        self.assertEqual(sorted(bioinfo1.faster_frequent_words(test_text, test_k)), test_output)

    def test_faster_frequent_words_3(self):
        """
        tests function for handling overlapping k-mers
        """
        test_text = "ATACAATTACAGTCTGGAACCGGATGAACTGGCCGCAGGTTAACAACAGAGTTGCCAGGCACTG" +\
                    "CCGCTGACCAGCAACAACAACAATGACTTTGACGCGAAGGGGATGGCATGAGCGAACTGATCGT" +\
                    "CAGCCGTCAGCAACGAGTATTGTTGCTGACCCTTAACAATCCCGCCGCACGTAATGCGCTAACT" +\
                    "AATGCCCTGCTG"
        test_k = 5
        test_output = ['AACAA']
        self.assertEqual(sorted(bioinfo1.faster_frequent_words(test_text, test_k)), test_output)

    def test_faster_frequent_words_4(self):
        """
        tests if function handles more than one most frequent k-mer correctly
        """
        test_text = "CCAGCGGGGGTTGATGCTCTGGGGGTCACAAGATTGCATTTTTATGGGGTTGCAAAAATGTTTT" +\
                    "TTACGGCAGATTCATTTAAAATGCCCACTGGCTGGAGACATAGCCCGGATGCGCGTCTTTTACA" +\
                    "ACGTATTGCGGGGTAAAATCGTAGATGTTTTAAAATAGGCGTAAC"
        test_k = 5
        test_output = ['AAAAT', 'GGGGT', 'TTTTA']
        self.assertEqual(sorted(bioinfo1.faster_frequent_words(test_text, test_k)), test_output)


class FrequentWordsBySortingTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_frequent_words_by_sorting_sample(self):
        """
        """
        test_text = "GCGCGCGTTGCATGTCGCATGATGCATGAGAGCT"
        test_k = 4
        test_output = ['CATG', 'GCAT']
        self.assertEqual(sorted(bioinfo1.frequent_words_by_sorting(test_text, test_k)), test_output)

    def test_frequent_words_by_sorting_1(self):
        """
        tests if the first k-mer in string is counted
        if the first k-mer is not counted, the output is ACT CAC CCA CTT GGT
        """
        test_text = "TGGTAGCGACGTTGGTCCCGCCGCTTGAGAATCTGGATGAACATAAGCTCCCACTTGGCTTATT" +\
                    "CAGAGAACTGGTCAACACTTGTCTCTCCCAGCCAGGTCTGACCACCGGGCAACTTTTAGAGCAC" +\
                    "TATCGTGGTACAAATAATGCTGCCAC"
        test_k = 3
        test_output = ['TGG']
        self.assertEqual(sorted(bioinfo1.frequent_words_by_sorting(test_text, test_k)), test_output)

    def test_frequent_words_by_sorting_2(self):
        """
        tests if the last k-mer in string is counted
        if the last k-mer is not counted, the output is AACG AATA ACAA CAAC CTGG CTTT
        TTGC TTTG
        """
        test_text = "CAGTGGCAGATGACATTTTGCTGGTCGACTGGTTACAACAACGCCTGGGGCTTTTGAGCAACGA" +\
                    "GACTTTTCAATGTTGCACCGTTTGCTGCATGATATTGAAAACAATATCACCAAATAAATAACGC" +\
                    "CTTAGTAAGTAGCTTTT"
        test_k = 4
        test_output = ['TTTT']
        self.assertEqual(sorted(bioinfo1.frequent_words_by_sorting(test_text, test_k)), test_output)

    def test_frequent_words_by_sorting_3(self):
        """
        tests function for handling overlapping k-mers
        """
        test_text = "ATACAATTACAGTCTGGAACCGGATGAACTGGCCGCAGGTTAACAACAGAGTTGCCAGGCACTG" +\
                    "CCGCTGACCAGCAACAACAACAATGACTTTGACGCGAAGGGGATGGCATGAGCGAACTGATCGT" +\
                    "CAGCCGTCAGCAACGAGTATTGTTGCTGACCCTTAACAATCCCGCCGCACGTAATGCGCTAACT" +\
                    "AATGCCCTGCTG"
        test_k = 5
        test_output = ['AACAA']
        self.assertEqual(sorted(bioinfo1.frequent_words_by_sorting(test_text, test_k)), test_output)

    def test_frequent_words_by_sorting_4(self):
        """
        tests if function handles more than one most frequent k-mer correctly
        """
        test_text = "CCAGCGGGGGTTGATGCTCTGGGGGTCACAAGATTGCATTTTTATGGGGTTGCAAAAATGTTTT" +\
                    "TTACGGCAGATTCATTTAAAATGCCCACTGGCTGGAGACATAGCCCGGATGCGCGTCTTTTACA" +\
                    "ACGTATTGCGGGGTAAAATCGTAGATGTTTTAAAATAGGCGTAAC"
        test_k = 5
        test_output = ['AAAAT', 'GGGGT', 'TTTTA']
        self.assertEqual(sorted(bioinfo1.frequent_words_by_sorting(test_text, test_k)), test_output)


class ReverseComplementTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_reverse_complement_sample(self):
        """
        """
        test_pattern = "AAAACCCGGT"
        test_output = "ACCGGGTTTT"
        self.assertEqual(bioinfo1.reverse_complement(test_pattern), test_output)

    def test_reverse_complement_1(self):
        """
        checks for both complement and reverse
        """
        test_pattern = "ACACAC"
        test_output = "GTGTGT"
        self.assertEqual(bioinfo1.reverse_complement(test_pattern), test_output)


class PatternMatchingProblemTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_pattern_matching_problem_sample(self):
        """
        """
        test_pattern = "ATAT"
        test_genome = "GATATATGCATATACTT"
        test_output = [1, 3, 9]
        self.assertEqual(bioinfo1.pattern_matching_problem(test_pattern, test_genome), test_output)

    def test_pattern_matching_problem_1(self):
        """
        This dataset checks if your code is written correctly but is also taking into account
        reverse complements, which we are not yet doing. Even though the reverse complement of
        “ACAC” (which is “GTGT”) occurs in Genome, we only want to count occurrences of “ACAC”
        specifically, which only occurs at index
        """
        test_pattern = "ACAC"
        test_genome = "TTTTACACTTTTTTGTGTAAAAA"
        test_output = [4]
        self.assertEqual(bioinfo1.pattern_matching_problem(test_pattern, test_genome), test_output)

    def test_pattern_matching_problem_2(self):
        """
        This dataset checks for off-by-one errors at the beginning of Genome. Notice that “AAA”
        occurs at the very beginning of Genome, so if you were to miss the first kmer of Genome,
        your code would output the following:
        """
        test_pattern = "AAA"
        test_genome = "AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAATTTTATTGAC" +\
                      "TTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAATAATTACAGAGTA" +\
                      "CACAACATCCAT"
        test_output = [0, 46, 51, 74]
        self.assertEqual(bioinfo1.pattern_matching_problem(test_pattern, test_genome), test_output)

    def test_pattern_matching_problem_3(self):
        """
        This dataset checks for off-by-one errors at the endof Genome. Notice that “TTT” occurs
        at the very end of Genome, so if you were to miss the last kmer of Genome, your code would
        output the following:
        """
        test_pattern = "TTT"
        test_genome = "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTG" +\
                      "CATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGC" +\
                      "CGACTTT"
        test_output = [88, 92, 98, 132]
        self.assertEqual(bioinfo1.pattern_matching_problem(test_pattern, test_genome), test_output)

    def test_pattern_matching_problem_4(self):
        """
        This test dataset checks if your code correctly handles cases where instances of Pattern
        overlap in Genome. In this case, if you did not count overlaps, you would only find the
        first and last instances of ATA (ATATATA and ATATATA). However, there is indeed a third
        occurrence, where the other two overlap (ATATATA).
        """
        test_pattern = "ATA"
        test_genome = "ATATATA"
        test_output = [0, 2, 4]
        self.assertEqual(bioinfo1.pattern_matching_problem(test_pattern, test_genome), test_output)


class ClumpFindingProblemTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_clump_finding_problem_sample(self):
        """
        """
        test_genome = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGG" +\
                      "AACATTGTAA"
        test_k = 5
        test_l = 50
        test_t = 4
        test_output = ['CGACA', 'GAAGA']
        self.assertEqual(sorted(bioinfo1.clump_finding_problem(test_genome, test_k, test_l,
                                                               test_t)), test_output)

    def test_clump_finding_problem_1(self):
        """
        This dataset makes sure that your code only counts kmers that fall COMPLETELY
        within a given L-window. For example, take the 4-window starting at index 4
        (AAAA*CGTC*GAAAAA). One might think that the 2-mer “CG” occurs twice in this window
        since the first letter of the second occurrence happens at the very end of the window.
        However, since the second occurrence of “CG” does not fall entirely in this 4-window,
        it does not count. Thus, the only result is “AA”.
        """
        test_genome = "AAAACGTCGAAAAA"
        test_k = 2
        test_l = 4
        test_t = 2
        test_output = ['AA']
        self.assertEqual(sorted(bioinfo1.clump_finding_problem(test_genome, test_k, test_l,
                                                               test_t)), test_output)

    def test_clump_finding_problem_2(self):
        """
        This dataset checks if your code has an off-by-one error when checking kmers within an
        L-window. Notice that, for each 1-mer (A, C, G, and T), there are 3 nucleotides between the
        first and second occurrence. In other words, each nucleotide occurs twice in a specific
        5-window: once at the beginning of the 5-window, and once at the end: *`ACGT`A*CGT,
        A*`CGTA`C*GT, AC*`GTAC`G*T, and ACG*`TACG'T*.
        """
        test_genome = "ACGTACGT"
        test_k = 1
        test_l = 5
        test_t = 2
        test_output = ['A', 'C', 'G', 'T']
        self.assertEqual(sorted(bioinfo1.clump_finding_problem(test_genome, test_k, test_l,
                                                               test_t)), test_output)

    def test_clump_finding_problem_3(self):
        """
        This dataset checks if your code is correctly handling overlapping kmers. For example,
        “ATA” forms a (5, 2)-clump in CCCATATACCC (CCC*`ATA`TA*CCC and CCC*AT`ATA`*CCC).
        """
        test_genome = "CCACGCGGTGTACGCTGCAAAAAGCCTTGCTGAATCAAATAAGGTTCCAGCACATCCTCAATGG" +\
                      "TTTCACGTTCTTCGCCAATGGCTGCCGCCAGGTTATCCAGACCTACAGGTCCACCAAAGAACTT" +\
                      "ATCGATTACCGCCAGCAACAATTTGCGGTCCATATAATCGAAACCTTCAGCATCGACATTCAAC" +\
                      "ATATCCAGCG"
        test_k = 3
        test_l = 25
        test_t = 3
        test_output = ['AAA', 'CAG', 'CAT', 'CCA', 'GCC', 'TTC']
        self.assertEqual(sorted(bioinfo1.clump_finding_problem(test_genome, test_k, test_l,
                                                               test_t)), test_output)


class ClumpFindingTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_clump_finding_sample(self):
        """
        """
        test_genome = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGG" +\
                      "AACATTGTAA"
        test_k = 5
        test_l = 50
        test_t = 4
        test_output = ['CGACA', 'GAAGA']
        self.assertEqual(sorted(bioinfo1.clump_finding(test_genome, test_k, test_l,
                                                       test_t)), test_output)

    def test_clump_finding_1(self):
        """
        This dataset makes sure that your code only counts kmers that fall COMPLETELY
        within a given L-window. For example, take the 4-window starting at index 4
        (AAAA*CGTC*GAAAAA). One might think that the 2-mer “CG” occurs twice in this window
        since the first letter of the second occurrence happens at the very end of the window.
        However, since the second occurrence of “CG” does not fall entirely in this 4-window,
        it does not count. Thus, the only result is “AA”.
        """
        test_genome = "AAAACGTCGAAAAA"
        test_k = 2
        test_l = 4
        test_t = 2
        test_output = ['AA']
        self.assertEqual(sorted(bioinfo1.clump_finding(test_genome, test_k, test_l,
                                                       test_t)), test_output)

    def test_clump_finding_2(self):
        """
        This dataset checks if your code has an off-by-one error when checking kmers within an
        L-window. Notice that, for each 1-mer (A, C, G, and T), there are 3 nucleotides between the
        first and second occurrence. In other words, each nucleotide occurs twice in a specific
        5-window: once at the beginning of the 5-window, and once at the end: *`ACGT`A*CGT,
        A*`CGTA`C*GT, AC*`GTAC`G*T, and ACG*`TACG'T*.
        """
        test_genome = "ACGTACGT"
        test_k = 1
        test_l = 5
        test_t = 2
        test_output = ['A', 'C', 'G', 'T']
        self.assertEqual(sorted(bioinfo1.clump_finding(test_genome, test_k, test_l,
                                                       test_t)), test_output)

    def test_clump_finding_3(self):
        """
        This dataset checks if your code is correctly handling overlapping kmers. For example,
        “ATA” forms a (5, 2)-clump in CCCATATACCC (CCC*`ATA`TA*CCC and CCC*AT`ATA`*CCC).
        """
        test_genome = "CCACGCGGTGTACGCTGCAAAAAGCCTTGCTGAATCAAATAAGGTTCCAGCACATCCTCAATGG" +\
                      "TTTCACGTTCTTCGCCAATGGCTGCCGCCAGGTTATCCAGACCTACAGGTCCACCAAAGAACTT" +\
                      "ATCGATTACCGCCAGCAACAATTTGCGGTCCATATAATCGAAACCTTCAGCATCGACATTCAAC" +\
                      "ATATCCAGCG"
        test_k = 3
        test_l = 25
        test_t = 3
        test_output = ['AAA', 'CAG', 'CAT', 'CCA', 'GCC', 'TTC']
        self.assertEqual(sorted(bioinfo1.clump_finding(test_genome, test_k, test_l,
                                                       test_t)), test_output)


class BetterClumpFindingTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_better_clump_finding_sample(self):
        """
        """
        test_genome = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGG" +\
                      "AACATTGTAA"
        test_k = 5
        test_l = 50
        test_t = 4
        test_output = ['CGACA', 'GAAGA']
        self.assertEqual(sorted(bioinfo1.better_clump_finding(test_genome, test_k, test_l,
                                                              test_t)), test_output)

    def test_better_clump_finding_1(self):
        """
        This dataset makes sure that your code only counts kmers that fall COMPLETELY
        within a given L-window. For example, take the 4-window starting at index 4
        (AAAA*CGTC*GAAAAA). One might think that the 2-mer “CG” occurs twice in this window
        since the first letter of the second occurrence happens at the very end of the window.
        However, since the second occurrence of “CG” does not fall entirely in this 4-window,
        it does not count. Thus, the only result is “AA”.
        """
        test_genome = "AAAACGTCGAAAAA"
        test_k = 2
        test_l = 4
        test_t = 2
        test_output = ['AA']
        self.assertEqual(sorted(bioinfo1.better_clump_finding(test_genome, test_k, test_l,
                                                              test_t)), test_output)

    def test_better_clump_finding_2(self):
        """
        This dataset checks if your code has an off-by-one error when checking kmers within an
        L-window. Notice that, for each 1-mer (A, C, G, and T), there are 3 nucleotides between the
        first and second occurrence. In other words, each nucleotide occurs twice in a specific
        5-window: once at the beginning of the 5-window, and once at the end: *`ACGT`A*CGT,
        A*`CGTA`C*GT, AC*`GTAC`G*T, and ACG*`TACG'T*.
        """
        test_genome = "ACGTACGT"
        test_k = 1
        test_l = 5
        test_t = 2
        test_output = ['A', 'C', 'G', 'T']
        self.assertEqual(sorted(bioinfo1.better_clump_finding(test_genome, test_k, test_l,
                                                              test_t)), test_output)

    def test_better_clump_finding_3(self):
        """
        This dataset checks if your code is correctly handling overlapping kmers. For example,
        “ATA” forms a (5, 2)-clump in CCCATATACCC (CCC*`ATA`TA*CCC and CCC*AT`ATA`*CCC).
        """
        test_genome = "CCACGCGGTGTACGCTGCAAAAAGCCTTGCTGAATCAAATAAGGTTCCAGCACATCCTCAATGG" +\
                      "TTTCACGTTCTTCGCCAATGGCTGCCGCCAGGTTATCCAGACCTACAGGTCCACCAAAGAACTT" +\
                      "ATCGATTACCGCCAGCAACAATTTGCGGTCCATATAATCGAAACCTTCAGCATCGACATTCAAC" +\
                      "ATATCCAGCG"
        test_k = 3
        test_l = 25
        test_t = 3
        test_output = ['AAA', 'CAG', 'CAT', 'CCA', 'GCC', 'TTC']
        self.assertEqual(sorted(bioinfo1.better_clump_finding(test_genome, test_k, test_l,
                                                              test_t)), test_output)


class ComputingFrequenciesTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_computing_frequencies_sample(self):
        """
        """
        test_text = "ACGCGGCTCTGAAA"
        test_k = 2
        test_output = [2, 1, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 1, 1, 0]
        self.assertEqual(bioinfo1.computing_frequencies(test_text, test_k), test_output)

    def test_computing_frequencies_1(self):
        """
        This dataset checks if you have an off-by-one error at the end of Text (i.e. you are not
        counting the last kmer in Text). There are three instances of AA (*AA*AAC, A*AA*AC, and
        AA*AA*C), but there is one instance of AC at the end (AAA*AC*).
        """
        test_text = "AAAAC"
        test_k = 2
        test_output = [3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(bioinfo1.computing_frequencies(test_text, test_k), test_output)

    def test_computing_frequencies_2(self):
        """
        This dataset checks if you have an off-by-one error at the beginning of Text (i.e. you are
        not counting the first kmer in Text). There are two instances of AA (TT*AA*A and TTA*AA*), but
        there is one instance of TTA (*TTA*AA) and one instance of TAA (T*TAA*A)
        """
        test_text = "TTAAA"
        test_k = 2
        test_output = [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1]
        self.assertEqual(bioinfo1.computing_frequencies(test_text, test_k), test_output)

    def test_computing_frequencies_3(self):
        """
        This dataset checks if your code actually increments each count, or if your code instead
        just sets the count equal to one each time. In other words, this dataset checks if your code is
        doing something like array[kmer] = 1 instead of array[kmer] += 1
        """
        test_text = "AAA"
        test_k = 2
        test_output = [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(bioinfo1.computing_frequencies(test_text, test_k), test_output)


if __name__ == "__main__":
        unittest.main()
