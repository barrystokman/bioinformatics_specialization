import unittest

import bioinformatics_1.functions as bioinfo1


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
        not counting the first kmer in Text). There are two instances of AA (TT*AA*A and TTA*AA*),
        but there is one instance of TTA (*TTA*AA) and one instance of TAA (T*TAA*A)
        """
        test_text = "TTAAA"
        test_k = 2
        test_output = [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1]
        self.assertEqual(bioinfo1.computing_frequencies(test_text, test_k), test_output)

    def test_computing_frequencies_3(self):
        """
        This dataset checks if your code actually increments each count, or if your code instead
        just sets the count equal to one each time. In other words, this dataset checks if your code
        is doing something like array[kmer] = 1 instead of array[kmer] += 1
        """
        test_text = "AAA"
        test_k = 2
        test_output = [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(bioinfo1.computing_frequencies(test_text, test_k), test_output)


class MinimumSkewTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_minimum_skew_sample(self):
        """
        """
        test_genome = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
        test_output = [11, 24]
        self.assertEqual(bioinfo1.minimum_skew(test_genome), test_output)

    def test_minimum_skew_1(self):
        """
        This dataset checks if your code’s indexing is off. Specifically, it verifies that your
        code is not returning an index 1 too high (i.e. 4) or 1 too low (i.e. 2)
        """
        test_genome = 'ACCG'
        test_output = [3]
        self.assertEqual(bioinfo1.minimum_skew(test_genome), test_output)

    def test_minimum_skew_2(self):
        """
        This dataset checks to see if your code is missing the last symbol of Genome.
        """
        test_genome = 'ACCC'
        test_output = [4]
        self.assertEqual(bioinfo1.minimum_skew(test_genome), test_output)

    def test_minimum_skew_3(self):
        """
        This dataset makes sure you’re not accidentally finding the maximum skew instead of the
        minimum skew.
        """
        test_genome = 'CCGGGT'
        test_output = [2]
        self.assertEqual(bioinfo1.minimum_skew(test_genome), test_output)

    def test_minimum_skew_4(self):
        """
        First, this dataset checks if you are only finding 1 index (and not multiple indices).
        Then, it checks if you are using a delimiter to separate your indices (ideally a space
        character)
        """
        test_genome = 'CCGGCCGG'
        test_output = [2, 6]
        self.assertEqual(bioinfo1.minimum_skew(test_genome), test_output)


class HammingDistanceTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_hamming_distance_sample(self):
        """
        """
        test_string1 = 'GGGCCGTTGGT'
        test_string2 = 'GGACCGTTGAC'
        test_output = 3
        self.assertEqual(bioinfo1.hamming_distance(test_string1, test_string2), test_output)

    def test_hamming_distance_1(self):
        """
        This dataset checks if your code isn’t keeping count (i.e. returns ‘0’ when the answer is
        clearly nonzero) or if your code returns a negative value, which is impossible.
        """
        test_string1 = 'AAAA'
        test_string2 = 'TTTT'
        test_output = 4
        self.assertEqual(bioinfo1.hamming_distance(test_string1, test_string2), test_output)

    def test_hamming_distance_2(self):
        """
        This dataset checks if your code is finding Edit Distance (which would be 2) instead of
        Hamming Distance.
        """
        test_string1 = 'ACGTACGT'
        test_string2 = 'TACGTACG'
        test_output = 8
        self.assertEqual(bioinfo1.hamming_distance(test_string1, test_string2), test_output)

    def test_hamming_distance_3(self):
        """
        This dataset checks if your code is returning the number of matches (2) instead of the
        number of mismatches (6).
        """
        test_string1 = 'ACGTACGT'
        test_string2 = 'CCCCCCCC'
        test_output = 6
        self.assertEqual(bioinfo1.hamming_distance(test_string1, test_string2), test_output)

    def test_hamming_distance_4(self):
        """
        This dataset checks if your code works on a dataset where the two input strings have no
        matches.
        """
        test_string1 = 'ACGTACGT'
        test_string2 = 'TGCATGCA'
        test_output = 8
        self.assertEqual(bioinfo1.hamming_distance(test_string1, test_string2), test_output)

    def test_hamming_distance_5(self):
        """
        This dataset checks if you have an off­by­one error at the beginning (i.e. you are starting
        at the second character of the strings instead of the first character).
        """
        test_string1 = 'GATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACT'
        test_string2 = 'AATAGCAGCTTCTCAACTGGTTACCTCGTATGAGTAAATTAGGTCATTATTGACTCAGGTCACTAACGTCT'
        test_output = 15
        self.assertEqual(bioinfo1.hamming_distance(test_string1, test_string2), test_output)

    def test_hamming_distance_6(self):
        """
        This dataset checks if you have an off­by­one error at the end (i.e. you are ending at the
        second­to­last character of the strings instead of the last character).
        """
        test_string1 = 'AGAAACAGACCGCTATGTTCAACGATTTGTTTTATCTCGTCACCGGGATATTGCGGCCACTCATCGGTCAG' +\
                       'TTGATTACGCAGGGCGTAAATCGCCAGAATCAGGCTG'
        test_string2 = 'AGAAACCCACCGCTAAAAACAACGATTTGCGTAGTCAGGTCACCGGGATATTGCGGCCACTAAGGCCTTGG' +\
                       'ATGATTACGCAGAACGTATTGACCCAGAATCAGGCTC'
        test_output = 28
        self.assertEqual(bioinfo1.hamming_distance(test_string1, test_string2), test_output)


class ApproximatePatternMatchTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_approximate_pattern_match_sample(self):
        """

        """
        test_pattern = 'ATTCTGGA'
        test_text = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
        test_distance = 3
        test_output = [6, 7, 26, 27]
        self.assertEqual(bioinfo1.approx_pattern_match(test_pattern, test_text, test_distance),
                         test_output)

    def test_approximate_pattern_match_1(self):
        """
        This dataset checks if you are only counting instances where the number of mismatches is
        exactly equal to d (i.e. ignoring instances where mismatch < d).
        """
        test_pattern = 'AAA'
        test_text = 'TTTTTTAAATTTTAAATTTTTT'
        test_distance = 2
        test_output = [4, 5, 6, 7, 8, 11, 12, 13, 14, 15]
        self.assertEqual(bioinfo1.approx_pattern_match(test_pattern, test_text, test_distance),
                         test_output)

    def test_approximate_pattern_match_2(self):
        """
        This dataset checks if your code has an off-by-one error at the beginning of Text (i.e.
        your code is not checking the the leftmost substring of Text).
        """
        test_pattern = 'GAGCGCTGG'
        test_text = 'GAGCGCTGGGTTAACTCGCTACTTCCCGACGAGCGCTGTGGCGCAAATTGGCGATGAAACTGCAGAGAGAACTG' +\
                    'GTCATCCAACTGAATTCTCCCCGCTATCGCATTTTGATGCGCGCCGCGTCGATT'
        test_distance = 2
        test_output = [0, 30, 66]
        self.assertEqual(bioinfo1.approx_pattern_match(test_pattern, test_text, test_distance),
                         test_output)

    def test_approximate_pattern_match_3(self):
        """
        This dataset checks if your code has an off-by-one error at the end of Text (i.e. your
        code is not checking the the rightmost substring of Text).
        """
        test_pattern = 'AATCCTTTCA'
        test_text = 'CCAAATCCCCTCATGGCATGCATTCCCGCAGTATTTAATCCTTTCATTCTGCATATAAGTAGTGAAGGTATAGA' +\
                    'AACCCGTTCAAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGCGCCATAATCCAAACA'
        test_distance = 3
        test_output = [3, 36, 74, 137]
        self.assertEqual(bioinfo1.approx_pattern_match(test_pattern, test_text, test_distance),
                         test_output)

    def test_approximate_pattern_match_4(self):
        """
        This  dataset  checks  if  your  code  is  correctly  accounting  for  overlapping
        instances  of Pattern in Text.
        """
        test_pattern = 'CCGTCATCC'
        test_text = 'CCGTCATCCGTCATCCTCGCCACGTTGGCATGCATTCCGTCATCCCGTCAGGCATACTTCTGCATATAAGTACA' +\
                    'AACATCCGTCATGTCAAAGGGAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGC'
        test_distance = 3
        test_output = [0, 7, 36, 44, 48, 72, 79, 112]
        self.assertEqual(bioinfo1.approx_pattern_match(test_pattern, test_text, test_distance),
                         test_output)

    def test_approximate_pattern_match_5(self):
        """
        This  dataset  checks  if  you  are  only  counting  instances  of Pattern with  less
        than d mismatches (as opposed to instances of Pattern with less than or equal to d
        mismatches).
        """
        test_pattern = 'TTT'
        test_text = 'AAAAAA'
        test_distance = 3
        test_output = [0, 1, 2, 3]
        self.assertEqual(bioinfo1.approx_pattern_match(test_pattern, test_text, test_distance),
                         test_output)

    def test_approximate_pattern_match_6(self):
        """
        This dataset checks if your code works with input where d = 0 (i.e. only perfect matches
        are allowed).
        """
        test_pattern = 'CCA'
        test_text = 'CCACCT'
        test_distance = 0
        test_output = [0]
        self.assertEqual(bioinfo1.approx_pattern_match(test_pattern, test_text, test_distance),
                         test_output)


class ApproximatePatternCountTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_approximate_pattern_count_sample(self):
        """

        """
        test_pattern = 'GAGG'
        test_text = 'TTTAGAGCCTTCAGAGG'
        test_distance = 2
        test_output = 4
        self.assertEqual(bioinfo1.approx_pattern_count(test_pattern, test_text, test_distance),
                         test_output)

    def test_approximate_pattern_count_1(self):
        """

        """
        test_pattern = 'AA'
        test_text = 'AAA'
        test_distance = 0
        test_output = 2
        self.assertEqual(bioinfo1.approx_pattern_count(test_pattern, test_text, test_distance),
                         test_output)

    def test_approximate_pattern_count_2(self):
        """
        This dataset checks if your code is allowing occurrences with less d than mismatches (which
        it should). It is a common mistake to only allow occurrences with exactly d mismatches,
        where as we want all occurrences with less than or equal to d mismatches.
        """
        test_pattern = 'ATA'
        test_text = 'ATA'
        test_distance = 1
        test_output = 1
        self.assertEqual(bioinfo1.approx_pattern_count(test_pattern, test_text, test_distance),
                         test_output)


class NeighborsTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_neighbors_sample(self):
        """
        """
        test_pattern = 'ACG'
        test_distance = 1
        test_output = sorted(['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG', 'ACA', 'ACC', 'ACT', 'ACG'])
        self.assertEqual(sorted(bioinfo1.neighbors(test_pattern, test_distance)), test_output)

    def test_neighbors_extra(self):
        """
        """
        test_pattern = 'GGCCCAGAG'
        test_distance = 3
        test_output = sorted(['GGCACAGAT',
                              'GTCCTACAG',
                              'GTATCAGAG',
                              'GGATCTGAG',
                              'GTCCTTGAG',
                              'GGCACAGAG',
                              'GGCACAGAC',
                              'GGCACAGAA',
                              'GGCGCGTAG',
                              'TGCCGAGAA',
                              'GTCACAAAG',
                              'GGCTCCTAG',
                              'AGGCCAGGG',
                              'GCCGCACAG',
                              'AGAACAGAG',
                              'GGAGCGGAG',
                              'GGGCGAGAC',
                              'GGTACTGAG',
                              'GTCCCACCG',
                              'GTCCCCGAT',
                              'CGCACAGTG',
                              'GGCCCCCGG',
                              'GGGCGCGAG',
                              'TGCTCAGGG',
                              'GTCCCCGAC',
                              'GTCCCCGAA',
                              'GTCCCCGAG',
                              'GGCCAGAAG',
                              'TGTCCATAG',
                              'ATTCCAGAG',
                              'GGGCAAGAT',
                              'CGCCTATAG',
                              'GACCGTGAG',
                              'AGCCCTGTG',
                              'GGGCCGGCG',
                              'GGGCAAGAG',
                              'GGGCAAGAC',
                              'GGGCAAGAA',
                              'TGCCCTGAT',
                              'CTCCCATAG',
                              'CGTCTAGAG',
                              'GGACCTCAG',
                              'CACCCTGAG',
                              'CGCCCCGCG',
                              'GTCCCAGGC',
                              'TGCCCTGAG',
                              'GTCCCAGGA',
                              'TGCCCTGAA',
                              'GTCCCAGGG',
                              'TGCCCTGAC',
                              'GACTCAGGG',
                              'GGACCATTG',
                              'TCCCCAGCG',
                              'GGCCCCAAT',
                              'GGCCATGGG',
                              'GTAGCAGAG',
                              'GGCCGATCG',
                              'AGCCGACAG',
                              'GGATCAGGG',
                              'GGCGTACAG',
                              'AGCCGATAG',
                              'AGCCAAGGG',
                              'GGTCCTGTG',
                              'GCCCCATTG',
                              'CGGCCAAAG',
                              'GGCGCGCAG',
                              'GACGCAGCG',
                              'GGTCGAGTG',
                              'GGCGCAGAT',
                              'CGCCCAAGG',
                              'GGCCTGTAG',
                              'TGCGCAGAT',
                              'GTGCCTGAG',
                              'GGCGCAGAG',
                              'GACGTAGAG',
                              'GGCGCAGAA',
                              'GGCGCAGAC',
                              'TGCGCAGAG',
                              'TGCGCAGAC',
                              'GGCGCATGG',
                              'GGCAGGGAG',
                              'GGCCATCAG',
                              'ACTCCAGAG',
                              'GGAAGAGAG',
                              'GGCTACGAG',
                              'AGCTCAAAG',
                              'CGCCCATTG',
                              'GATACAGAG',
                              'GGCTCTAAG',
                              'TGCCCAACG',
                              'GGCTTAGAC',
                              'GGTCCGGAA',
                              'GGCTTAGAA',
                              'TGGCCGGAG',
                              'GGCTTAGAG',
                              'GGTCCGGAG',
                              'GCCCTACAG',
                              'GGTCCGGAT',
                              'CGCCCGGCG',
                              'GGCTTAGAT',
                              'GGGCCCTAG',
                              'GGCCTGGCG',
                              'GAGGCAGAG',
                              'GATCCAAAG',
                              'GAAGCAGAG',
                              'GGCCTTGGG',
                              'GTCTCATAG',
                              'GCTCTAGAG',
                              'GGTGCCGAG',
                              'GGCCCCTGG',
                              'GGCTCCGGG',
                              'GGCCTAATG',
                              'GGCCGGGAT',
                              'GGCTGAAAG',
                              'GACAGAGAG',
                              'GGTGAAGAG',
                              'GCCCCCCAG',
                              'GACCTAGGG',
                              'CTCCGAGAG',
                              'GACTGAGAG',
                              'TACCCTGAG',
                              'GGGCCAGGG',
                              'TGCCCAGGG',
                              'GGCGGAAAG',
                              'AGACCAGAT',
                              'TGCCCAGGC',
                              'GGCCCAGAA',
                              'GGCCGCAAG',
                              'TACCCACAG',
                              'GGCCCAGAG',
                              'AGACCAGAA',
                              'GGACCTTAG',
                              'AGACCAGAC',
                              'AGACCAGAG',
                              'GGCCCAGAT',
                              'GACCCTCAG',
                              'GATCCAGCG',
                              'TGCCCTCAG',
                              'GGACAAGAG',
                              'GGCCCGGTG',
                              'GGACAAGAA',
                              'GGCCCGGTA',
                              'GGACAAGAC',
                              'GGCCCGGTC',
                              'GGCATCGAG',
                              'GGCTCAACG',
                              'GGAGCCGAG',
                              'GGCCCGGTT',
                              'TGCCCGGGG',
                              'AGCACAGTG',
                              'GGACCGTAG',
                              'GAGCCGGAG',
                              'GAGCAAGAG',
                              'AGCCCCGGG',
                              'TGCCGACAG',
                              'GACCAAGTG',
                              'AGCGCCGAG',
                              'GGGGCAGCG',
                              'CGTCCAGAG',
                              'GTCCCTGAT',
                              'CGTCCAGAC',
                              'CGTCCAGAA',
                              'GGTCTAGGG',
                              'GACCCACAG',
                              'ATCCCAGAT',
                              'GTACCCGAG',
                              'GTCCCTGAG',
                              'CCCCCTGAG',
                              'CGTCCAGAT',
                              'GTCCCTGAC',
                              'GTCCCTGAA',
                              'ATCCCAGAG',
                              'AGCGCAGCG',
                              'GGTAAAGAG',
                              'ATCCCAGAA',
                              'GGCGCTGGG',
                              'GACTCATAG',
                              'GAGCCACAG',
                              'AGACGAGAG',
                              'CGCGCATAG',
                              'TGCCCGTAG',
                              'AGGTCAGAG',
                              'AGCCCGCAG',
                              'GGCCTCAAG',
                              'GGGCCCGCG',
                              'GCACCATAG',
                              'AGCCCAGGT',
                              'CGACCAGTG',
                              'GGCTCAGAA',
                              'GGCTCAGAC',
                              'GGTTCAGCG',
                              'GGCTCAGAT',
                              'AGCCCAGGG',
                              'GTCACCGAG',
                              'AGCCCAGGA',
                              'AGCCCAGGC',
                              'AGCCCTGGG',
                              'GGCCGAACG',
                              'CATCCAGAG',
                              'TCCCCAAAG',
                              'GGACGAGGG',
                              'GGCGCAGCT',
                              'GTCCCGCAG',
                              'GACACAGAA',
                              'GACACAGAC',
                              'GACACAGAG',
                              'GGCACAGCG',
                              'AGCCCCCAG',
                              'GACACAGAT',
                              'GGCCGGGGG',
                              'GGACCAGGT',
                              'GCCCAACAG',
                              'TGGTCAGAG',
                              'GGGCTACAG',
                              'GCACCAGCG',
                              'GGCTCATCG',
                              'GGACCAGTT',
                              'GGACCAGGG',
                              'GGACCAGGA',
                              'GGACCAGGC',
                              'TGTCCACAG',
                              'GTGCGAGAG',
                              'GGCCCTATG',
                              'GACGCATAG',
                              'TGCCTACAG',
                              'CGGACAGAG',
                              'GGGCCAAGG',
                              'AGCTCAGCG',
                              'CGCCAACAG',
                              'GGCACTCAG',
                              'GACGGAGAG',
                              'GTCCCCGTG',
                              'GTCGCTGAG',
                              'GGACCCGGG',
                              'TGCACTGAG',
                              'GCCCCACCG',
                              'GGCCCCCTG',
                              'GTACCAGGG',
                              'GAATCAGAG',
                              'GGTCGAGGG',
                              'GGCCGACAA',
                              'GGCCGACAC',
                              'GGCATAAAG',
                              'GGCCGACAG',
                              'GGCAGAAAG',
                              'ATCTCAGAG',
                              'AACCTAGAG',
                              'TGTTCAGAG',
                              'GGCCGACAT',
                              'GCCAAAGAG',
                              'GGACCATCG',
                              'GGACCAGTG',
                              'GTCCCAGGT',
                              'GGCAAAGGG',
                              'GGCCTAGAG',
                              'GGTACATAG',
                              'GGCCTAGAC',
                              'GGCTCCAAG',
                              'GCGCGAGAG',
                              'TGCCCAAAG',
                              'GGAGAAGAG',
                              'GCCCCAAAT',
                              'CGCCCAGTC',
                              'CGCCCAGTA',
                              'CGCCCAGTG',
                              'CGCCCTGCG',
                              'GCCCCAAAG',
                              'CGCCATGAG',
                              'GCCCCAAAC',
                              'GGACCAGTC',
                              'GCCTCCGAG',
                              'AGCCTAAAG',
                              'CGCCCAGTT',
                              'CTCCTAGAG',
                              'GGCCCGAAT',
                              'GGTCCTGAG',
                              'GGTCCTGAA',
                              'GGTCCTGAC',
                              'GCCCCATGG',
                              'GGCTCACTG',
                              'GGTCCTGAT',
                              'GGCCCGAAG',
                              'GGCCAACTG',
                              'GGCCCGAAA',
                              'GGCCCGAAC',
                              'GGTCCCAAG',
                              'CGACTAGAG',
                              'GGACTGGAG',
                              'TGTACAGAG',
                              'GGCCTAACG',
                              'GTCCGGGAG',
                              'TGACCAGTG',
                              'TGCATAGAG',
                              'TAGCCAGAG',
                              'CGCCGAGGG',
                              'GACCCCGGG',
                              'GTCTCACAG',
                              'GGAGCAGGG',
                              'AGCGCAGGG',
                              'GTCCAAGCG',
                              'GTGCAAGAG',
                              'AGCGAAGAG',
                              'GTTCCAGAT',
                              'AGTCTAGAG',
                              'GCTCCACAG',
                              'CCACCAGAG',
                              'GGTCCCGCG',
                              'CGCCGACAG',
                              'GGCTGACAG',
                              'TGCCCACTG',
                              'GTTCCAGAG',
                              'GTTCCAGAC',
                              'GTTCCAGAA',
                              'AGCCCCGTG',
                              'TTCCCAGAC',
                              'GGCCCGCCG',
                              'TTCCCAGAG',
                              'GTACCAAAG',
                              'AACCCAGCG',
                              'TTCCCAGAT',
                              'GACCAAGAG',
                              'GACCAAGAA',
                              'GACCAAGAC',
                              'GACCCAGGG',
                              'GACCCAGGA',
                              'GGGCCATTG',
                              'GGTCAAGGG',
                              'GTCCAAAAG',
                              'GACCCAGGT',
                              'AGCCTCGAG',
                              'GCAACAGAG',
                              'GACCAAGAT',
                              'TGGCCAGTG',
                              'GCCCGAGGG',
                              'AGTTCAGAG',
                              'AGCCCTCAG',
                              'GGAACAGCG',
                              'GTCCGATAG',
                              'GCGCCAGCG',
                              'GGACCGCAG',
                              'GGCCCTAGG',
                              'GGCCCCACG',
                              'CGACCAGCG',
                              'GGGTCAGGG',
                              'GGCCCTTGG',
                              'GGCCCCGGA',
                              'GGCCCCGGC',
                              'GGTACGGAG',
                              'GGCCCCGGG',
                              'CGCTCAGAC',
                              'TGCCCCGAT',
                              'GGCCCCAGG',
                              'GGCCCCGGT',
                              'GGTACAGTG',
                              'TGCCCCGAG',
                              'GGTCTAGTG',
                              'TGCCCCGAC',
                              'CGCCCAGCG',
                              'GGGCCGGAC',
                              'GGGCCGGAA',
                              'GGGCCGGAG',
                              'GCCGGAGAG',
                              'GGGCCGGAT',
                              'GGCAAGGAG',
                              'CGCCCGGAT',
                              'GGGGCAGTG',
                              'GGCACAGGG',
                              'GGCACAGGA',
                              'GGTCCATAT',
                              'GGCACAGGC',
                              'CGCCCGGAC',
                              'CGCCCGGAA',
                              'CGCCCGGAG',
                              'GGTCCAACG',
                              'GTCCCAAGG',
                              'GGTCCATAA',
                              'GGCACAGGT',
                              'GGTCCATAC',
                              'GGTCCATAG',
                              'GGCACCGGG',
                              'CGCCTAGAT',
                              'GGCCCTCGG',
                              'TACCCAGTG',
                              'AGTCCAGAA',
                              'GGCGTATAG',
                              'AGTCCAGAC',
                              'CGCCTAGAC',
                              'CGCCTAGAA',
                              'GTGCCATAG',
                              'CGCCTAGAG',
                              'GTCCCCGCG',
                              'AGTCCAGAT',
                              'GTCCCACAA',
                              'GTCCCACAC',
                              'GGACCCGTG',
                              'GTCCCACAG',
                              'GACCCTGTG',
                              'AGCCCAGTG',
                              'GACCCGAAG',
                              'AGCCCAGTA',
                              'AGCCCAGTC',
                              'GGCCGTCAG',
                              'GGCATTGAG',
                              'AGCCCAGTT',
                              'GGGCAAGGG',
                              'TGCCGCGAG',
                              'CGCGCTGAG',
                              'CGACCATAG',
                              'GGCACTGTG',
                              'GGCCAAGCT',
                              'GACTCAGAT',
                              'CGCCCCGAG',
                              'GGCCCAGCC',
                              'AGTGCAGAG',
                              'CGCCCCGAC',
                              'GTCACACAG',
                              'GCCACAGGG',
                              'GGCGCAGTG',
                              'GACTCAGAG',
                              'GGCCAAGCG',
                              'GGCCAAGCA',
                              'GACTCAGAC',
                              'GGCCAAGCC',
                              'GACTCAGAA',
                              'TGCCCTTAG',
                              'CGCCCCGAT',
                              'TGCGGAGAG',
                              'TCCCCAGAC',
                              'CGCCCACGG',
                              'GGCTCGGTG',
                              'GGCCGGGTG',
                              'GGATCAGAT',
                              'TGCCCTGCG',
                              'GTGCCGGAG',
                              'GGCCCGTGG',
                              'GTCCTCGAG',
                              'GGATCAGAG',
                              'GGATCAGAC',
                              'GGATCAGAA',
                              'AGCCAAGAT',
                              'GGACTAGCG',
                              'GTCCTAGAT',
                              'AGCCAAGAG',
                              'AGCCAAGAA',
                              'AGCCAAGAC',
                              'GCCCGCGAG',
                              'GGGCCAATG',
                              'GACGCAGAT',
                              'GGCACAGCC',
                              'AGACCATAG',
                              'GACGCAGAC',
                              'GACGCAGAA',
                              'GTCCCTTAG',
                              'GACGCAGAG',
                              'TGCACGGAG',
                              'GGTCCAGGT',
                              'TGACCAGGG',
                              'GGCGCATAT',
                              'GGGATAGAG',
                              'TGACAAGAG',
                              'GGTCCAGGG',
                              'GGTCCAGGA',
                              'GGTCCAGGC',
                              'GGCGCATAG',
                              'TCCACAGAG',
                              'GGCGCATAC',
                              'GGCGCATAA',
                              'GGCGCAGCC',
                              'GGCGCAGCA',
                              'GCCCAATAG',
                              'GGCGCAGCG',
                              'GGACAATAG',
                              'GACCCAATG',
                              'TGCCCACGG',
                              'GGCTTAAAG',
                              'GCCCCGTAG',
                              'CGGCCAGAT',
                              'AGCTTAGAG',
                              'CCCACAGAG',
                              'GGTCCGGGG',
                              'TGCCCAAAT',
                              'GGCTCATGG',
                              'TGCGCAGCG',
                              'CGGCCAGAA',
                              'CGGCCAGAC',
                              'CGGCCAGAG',
                              'TGCCCAAAA',
                              'GTCCCGTAG',
                              'TGCCCAAAC',
                              'GGACCAGTA',
                              'GGCCCGCTG',
                              'GGCCATTAG',
                              'GACACATAG',
                              'GAGCCAAAG',
                              'TCCTCAGAG',
                              'AACCCAGTG',
                              'GCTCCATAG',
                              'GGCCCATAT',
                              'AGCCCCGAG',
                              'TGCGCAAAG',
                              'AGCCCCGAA',
                              'GGCCCAGTT',
                              'AGCCCCGAC',
                              'GGCCCATAC',
                              'GGCACGCAG',
                              'GGCCCATAA',
                              'GGCCCATAG',
                              'GCCCCTTAG',
                              'AGCCCCGAT',
                              'ATCCGAGAG',
                              'AGCTCACAG',
                              'CGCGCAGTG',
                              'GGCCCAGCG',
                              'GGGCGACAG',
                              'CTCCCAGTG',
                              'TGCAAAGAG',
                              'GGCCCAGCA',
                              'AGACCAGCG',
                              'GCCCTATAG',
                              'CGCACAGCG',
                              'GATCCAGAC',
                              'GGATGAGAG',
                              'TACCAAGAG',
                              'GGACAAGCG',
                              'GGCAATGAG',
                              'GATCCAGAT',
                              'GTCGCATAG',
                              'CCCCCAGGG',
                              'GATCCAGAG',
                              'GATCCAGAA',
                              'GGCGCTAAG',
                              'ACCGCAGAG',
                              'GTGTCAGAG',
                              'AGCCCACGG',
                              'AGCCCTTAG',
                              'GGGACAGCG',
                              'GCCACATAG',
                              'GGCTATGAG',
                              'CGCGCAGCG',
                              'GGTCTAGAG',
                              'GGTCTAGAA',
                              'GGCGGAGAT',
                              'GGTCTAGAC',
                              'GGCCCCGTA',
                              'GGCCCCGTC',
                              'GGCCCCGTG',
                              'GGCGGAGAA',
                              'GACCCTAAG',
                              'GGCGGAGAC',
                              'ATCCCAGCG',
                              'GGTACAGCG',
                              'GGCGGAGAG',
                              'CGTCCAGCG',
                              'GGCCGCCAG',
                              'TACCCAAAG',
                              'TGCCCGCAG',
                              'GGCCTTTAG',
                              'GGCTCTGAC',
                              'GTCCCTGGG',
                              'GGCCCATTC',
                              'GCCCGGGAG',
                              'GGACTACAG',
                              'TGCGCAGAA',
                              'GGGTAAGAG',
                              'TGCTCATAG',
                              'GGGCCCGAG',
                              'CCCCAAGAG',
                              'GGGCCCGAA',
                              'CACCCAGAC',
                              'CACCCAGAA',
                              'CACCCAGAG',
                              'TGGCCAGAC',
                              'GGGCCCGAT',
                              'GGCACAATG',
                              'GTTACAGAG',
                              'CGCTCACAG',
                              'CACCCAGAT',
                              'GGCCGAAAG',
                              'GGAACTGAG',
                              'GGCCGAAAC',
                              'GGCGAGGAG',
                              'GAAACAGAG',
                              'TTCACAGAG',
                              'GGTTCAGAG',
                              'GGCTCTGCG',
                              'GGCCGAAAT',
                              'GGGCCACGG',
                              'GGTTCAGAT',
                              'GTCCCAGTG',
                              'GTCCCAGTC',
                              'GTCCCAGTA',
                              'GACACAGCG',
                              'GACCCTGCG',
                              'ATCCCAAAG',
                              'AGCCGCGAG',
                              'GTCCCAGTT',
                              'CGCCCACTG',
                              'GGGCAAGTG',
                              'GGTCCACTG',
                              'GGGTCGGAG',
                              'GCACCAGAC',
                              'GCACCAGAA',
                              'GCACCAGAG',
                              'GGCGAAGAA',
                              'GGCTCATAG',
                              'GGCTCATAC',
                              'GGCTCATAA',
                              'CGGCCATAG',
                              'GGCACTAAG',
                              'TCTCCAGAG',
                              'GGCTCTCAG',
                              'GGTTCAAAG',
                              'GCACCAGAT',
                              'GGTCCACAC',
                              'GGCTCATAT',
                              'GCGCCAGTG',
                              'GACCCACTG',
                              'GCCACCGAG',
                              'GGACCACGG',
                              'TGGCCAGGG',
                              'GGAGCAGAC',
                              'GCCCCCGCG',
                              'GGTACAAAG',
                              'CCCCCGGAG',
                              'GGACCCGAG',
                              'CTTCCAGAG',
                              'TCCCCACAG',
                              'GGACCCGAC',
                              'GGTCTAGAT',
                              'GGACCCGAT',
                              'GACCGAGGG',
                              'GGACCTAAG',
                              'GACCTACAG',
                              'GGTCGAGAG',
                              'GGTCGAGAC',
                              'ACCCCATAG',
                              'GGCCGACCG',
                              'GGCCCGACG',
                              'GGTCGAGAT',
                              'GGCAAAGAG',
                              'GGCCTAGCG',
                              'GGCAAAGAA',
                              'GGCCTAGCA',
                              'GGCAAAGAC',
                              'GGCCTAGCC',
                              'GGACCATAT',
                              'GGCACTGCG',
                              'TGTCCAAAG',
                              'GGCCTAGCT',
                              'GGCAAAGAT',
                              'GACACTGAG',
                              'GTCGGAGAG',
                              'GGACCATAA',
                              'GGACCATAC',
                              'GGTCCGGAC',
                              'TGCCTAAAG',
                              'CACCTAGAG',
                              'GGACCTGCG',
                              'GGCCACGGG',
                              'GCCCCATAG',
                              'CGCTGAGAG',
                              'GCCCCATAC',
                              'GCCCCATAA',
                              'GCCCCAAGG',
                              'GCCCCATAT',
                              'GCCTCTGAG',
                              'GGTCCTGCG',
                              'GGCTCCCAG',
                              'GCCCCGGCG',
                              'GGGCATGAG',
                              'GGCAGACAG',
                              'GGTGCGGAG',
                              'GGCCTAGTC',
                              'GCCCAAGCG',
                              'GGCAGTGAG',
                              'GTCACAGGG',
                              'GGCCTAAAA',
                              'GACCCACGG',
                              'GGCCTAAAC',
                              'CGCCAGGAG',
                              'GACCCCGAT',
                              'CCCCCCGAG',
                              'TGCCGAGTG',
                              'AGTACAGAG',
                              'AGCCCGAAG',
                              'GAGCCTGAG',
                              'GACCCCGAG',
                              'GCCATAGAG',
                              'GACCCCGAA',
                              'GACCCCGAC',
                              'TGAGCAGAG',
                              'GTCACATAG',
                              'CACACAGAG',
                              'TGGCGAGAG',
                              'GGCTCGGGG',
                              'GTTCCAGGG',
                              'TGTCCAGCG',
                              'GAACCAGAA',
                              'TACCCCGAG',
                              'GAACCAGAC',
                              'GGCCCGCAG',
                              'GAACCAGAG',
                              'GGTCCCGAG',
                              'AGCACAGGG',
                              'GGTCCCGAA',
                              'GGTCCCGAC',
                              'GAACCAGAT',
                              'GGCCCGCAT',
                              'GGTCCCGAT',
                              'GGCCTAGTT',
                              'GAGCCCGAG',
                              'GGGTTAGAG',
                              'CGCCAAGTG',
                              'TGCCAACAG',
                              'GGCCAAGGA',
                              'GTCTCAAAG',
                              'GTCCTAGGG',
                              'GGGCCCGAC',
                              'GACCAAGCG',
                              'GGGCACGAG',
                              'GGTCGAGCG',
                              'CGCTCAGGG',
                              'GGGCCTGAC',
                              'GGGCCTGAA',
                              'GGGCCTGAG',
                              'GGCGATGAG',
                              'GGGCCTGAT',
                              'GGCCCTTAT',
                              'CGCCAATAG',
                              'CGCCGAAAG',
                              'GGACGAGTG',
                              'GGCTTGGAG',
                              'GGCTAGGAG',
                              'GGCCCTGGT',
                              'TGGCTAGAG',
                              'CGACCAGAT',
                              'GGGCCACTG',
                              'GGCCCTTAG',
                              'GGCCCTTAA',
                              'GGCCCTTAC',
                              'CGACCAGAC',
                              'GGCCCTGGA',
                              'CGACCAGAA',
                              'GGCCCTGGC',
                              'CGACCAGAG',
                              'TTCCCAAAG',
                              'GGCCCTGGG',
                              'AGACCTGAG',
                              'GGTGCAGTG',
                              'ACCCCAGCG',
                              'AGCCACGAG',
                              'CGTCCAAAG',
                              'GACTTAGAG',
                              'GGCTGCGAG',
                              'GTTCTAGAG',
                              'AGCCCTAAG',
                              'TGCTCAGCG',
                              'GTCCAACAG',
                              'GCCCCAGTT',
                              'GCCCCTGAT',
                              'GCACGAGAG',
                              'GCCCCAGTC',
                              'GCCCCAGTA',
                              'GCCCCAGTG',
                              'GCCCCTGAC',
                              'GCCCCTGAA',
                              'CGTGCAGAG',
                              'AGGCCAGCG',
                              'TGCACAGCG',
                              'CGCCAAGCG',
                              'CGCCCAGGT',
                              'GGACCGAAG',
                              'GTGCTAGAG',
                              'GTCCCAAAG',
                              'GTCCCAAAC',
                              'GTCCCAAAA',
                              'GGTCCATCG',
                              'CGCCCTGTG',
                              'CGCCCAGGC',
                              'CGCCCAGGA',
                              'GGCGGTGAG',
                              'GGCACCGAT',
                              'CGGCGAGAG',
                              'GTCCCACGG',
                              'GGCACCGAG',
                              'GGCACCGAA',
                              'GGCACCGAC',
                              'GACCCATAT',
                              'GGCCGAAAA',
                              'GATCTAGAG',
                              'CGCCTAGCG',
                              'GGTTCAGAC',
                              'GACCCATAG',
                              'GACCCATAA',
                              'GACCCATAC',
                              'AGCTCATAG',
                              'GGTTCAGAA',
                              'GGCACATAT',
                              'TGTGCAGAG',
                              'GGCCGCTAG',
                              'GCGCCTGAG',
                              'GGCCTCTAG',
                              'TGACCATAG',
                              'GGCGCCAAG',
                              'GGCACATAG',
                              'GGCACATAA',
                              'GGCACATAC',
                              'GGCGGATAG',
                              'GCGGCAGAG',
                              'CGCCCGAAG',
                              'GGTCCACGG',
                              'AGTCCATAG',
                              'GTCCCAGCG',
                              'GGCCAAGAA',
                              'GACTCAGCG',
                              'GGCCAAGAG',
                              'GTCCCAGCC',
                              'GTCCCAGCA',
                              'GGATCAGCG',
                              'TCCCCATAG',
                              'GCATCAGAG',
                              'GGCCAAGAT',
                              'GACCCGGAT',
                              'AGCCGAGGG',
                              'GGTCCAAAA',
                              'GGTCCAAAC',
                              'GGACTAGAT',
                              'GCTCCCGAG',
                              'GACCCGGAA',
                              'GGCCCGTCG',
                              'GACCCGGAC',
                              'GACCCGGAG',
                              'AGCCAATAG',
                              'GGACTAGAC',
                              'GGTCCAAAT',
                              'GTGCCAGCG',
                              'GGACTAGAG',
                              'GAACTAGAG',
                              'AGCCAAGCG',
                              'GGACCACTG',
                              'GGCCTTCAG',
                              'GGTCAACAG',
                              'GGCCCAGAC',
                              'GGGCCGAAG',
                              'GACTCAAAG',
                              'GAGTCAGAG',
                              'GCCGCAGAA',
                              'GCCGCAGAC',
                              'GGCGCATCG',
                              'GCCGCAGAG',
                              'TGTCGAGAG',
                              'TGCCCAGGT',
                              'TGCCCACAT',
                              'GGCGTAGAT',
                              'AGCACTGAG',
                              'GGATCAAAG',
                              'GGTCCTGGG',
                              'TGCCCACAG',
                              'GGCGTAGAA',
                              'GGGTCAGTG',
                              'TGCCCACAC',
                              'TGCCCACAA',
                              'CGGCCAGCG',
                              'CGATCAGAG',
                              'TTCCGAGAG',
                              'GCGCTAGAG',
                              'CGCGCAGGG',
                              'CCCTCAGAG',
                              'GCCCGATAG',
                              'GACGCAAAG',
                              'GCTCAAGAG',
                              'AGCCCAAGG',
                              'GCCTCGGAG',
                              'GCCCCGGTG',
                              'AACTCAGAG',
                              'CTCCCAAAG',
                              'GGCCCTGCT',
                              'AGCCCCGCG',
                              'GGGACAGAC',
                              'GGGACAGAA',
                              'GGGACAGAG',
                              'GGACCTGAA',
                              'GGCTTACAG',
                              'GGTCCGCAG',
                              'GGCTTTGAG',
                              'GGACCTGAC',
                              'GGGACAGAT',
                              'GGTCATGAG',
                              'CGCACAGAT',
                              'GGCCTGGTG',
                              'GGCCCATCA',
                              'GGCCCATCC',
                              'GGGTGAGAG',
                              'GCCCTAGAT',
                              'GGCCCATCG',
                              'CGCACAGAC',
                              'GGCGCTCAG',
                              'CGCACAGAA',
                              'CGCACAGAG',
                              'GTCGCAGCG',
                              'GCCCTAGAC',
                              'GCCCTAGAA',
                              'GGCCCATCT',
                              'GCCCTAGAG',
                              'GGCTTATAG',
                              'GGGGCTGAG',
                              'AGCCCATGG',
                              'ACCCCAGTG',
                              'GGGCCAGGT',
                              'GGCCGCGCG',
                              'GAGCTAGAG',
                              'ATCCCACAG',
                              'TTGCCAGAG',
                              'AGCCCCTAG',
                              'GGGCCAGGC',
                              'GGGCCAGGA',
                              'AGCACGGAG',
                              'GTTCCAGTG',
                              'GGCACGAAG',
                              'GGCCCTGCC',
                              'GGTGCAGGG',
                              'GGCGGAGCG',
                              'GCCTCACAG',
                              'GGGGCAGGG',
                              'AGGCCAGTG',
                              'GGGCGAAAG',
                              'GGTACAGAA',
                              'GGTACAGAC',
                              'GCCTAAGAG',
                              'GGACAACAG',
                              'GGTACAGAG',
                              'GCCGTAGAG',
                              'GCACCCGAG',
                              'GGTACAGAT',
                              'CACCCAAAG',
                              'AACCCAGGG',
                              'GGGCTAGAG',
                              'GTACCAGTG',
                              'CGCCCAAAA',
                              'GGGCCCGGG',
                              'GAGCCAGGG',
                              'GGCGCACAT',
                              'GGACCAAGG',
                              'GGCGCACAA',
                              'GGCGCACAC',
                              'AATCCAGAG',
                              'GGCGCACAG',
                              'TGCTCAGTG',
                              'AGCCCTGCG',
                              'GGCCGAAGG',
                              'TGCGCACAG',
                              'GACCCTGAG',
                              'GACCCTGAA',
                              'GACCCTGAC',
                              'ACGCCAGAG',
                              'GCCCCACGG',
                              'GGCATGGAG',
                              'GACCCTGAT',
                              'TGTCCAGAA',
                              'GGCCGGGCG',
                              'TGCCGGGAG',
                              'TGCCCGAAG',
                              'GGTTCATAG',
                              'TGTCCAGAC',
                              'GCCCCTGTG',
                              'GGTCCAGCC',
                              'GGCTAAGTG',
                              'GTCAAAGAG',
                              'CCCCCATAG',
                              'GGCGCGGTG',
                              'ACCCCACAG',
                              'GGTCCGGTG',
                              'GGGCTAGAT',
                              'CGCTCATAG',
                              'GGTCCAGCG',
                              'GCCTCAAAG',
                              'CTCCCCGAG',
                              'GCCCCCGAT',
                              'GGCACACTG',
                              'GGCCTATGG',
                              'TGGCCAGAG',
                              'CGCTCAAAG',
                              'TGGCCAGAA',
                              'GCCCCCGAG',
                              'GCCCCCGAC',
                              'GCCCCCGAA',
                              'TGGCCAGAT',
                              'GGGCCGCAG',
                              'AGCCTGGAG',
                              'GTCACTGAG',
                              'TCCCCCGAG',
                              'GGCGCAACG',
                              'GGACCCGCG',
                              'TCCCAAGAG',
                              'TCCCGAGAG',
                              'GTCGCACAG',
                              'CGCCCAATG',
                              'GGCTAACAG',
                              'CGCGCCGAG',
                              'GCACCAAAG',
                              'CGTACAGAG',
                              'GACCCAAGG',
                              'GGCACTGAG',
                              'GTCCGAGAT',
                              'GGCACTGAA',
                              'GGCACTGAC',
                              'CTACCAGAG',
                              'GGCACTGAT',
                              'GTCCGAGAA',
                              'GTCCGAGAC',
                              'GGCATATAG',
                              'GTCCGAGAG',
                              'GGCAAAGCG',
                              'CGCCCATGG',
                              'GGTCCATTG',
                              'CGCCCTGGG',
                              'TGGGCAGAG',
                              'GGTCCCTAG',
                              'AACCCTGAG',
                              'GTGACAGAG',
                              'GGCCCAAAC',
                              'GTCCCACTG',
                              'GTCCAATAG',
                              'AGCCTTGAG',
                              'GGCCCAAAA',
                              'GGAGCATAG',
                              'AACGCAGAG',
                              'GGACTAGTG',
                              'AGCCCAATG',
                              'AGCGCATAG',
                              'TTCGCAGAG',
                              'GGACCTGAG',
                              'GACCTAAAG',
                              'AGCCGAGTG',
                              'GGCCTAAAG',
                              'GGACCTGAT',
                              'GCCCCATCG',
                              'GGCCATGCG',
                              'CGTTCAGAG',
                              'GCCCAAGAG',
                              'GGGCTAGAC',
                              'GGGCTAGAA',
                              'GCCCAAGAC',
                              'TGTCCAGAT',
                              'GCCCAAGAA',
                              'GCCCCGGAT',
                              'GCCCCGGAC',
                              'GCCCCGGAA',
                              'GCCCAAGAT',
                              'GCCCCGGAG',
                              'GGCCTGGGG',
                              'TGTCCAGAG',
                              'GGAGCAGCG',
                              'TGACGAGAG',
                              'GGCCGATTG',
                              'ATCCCAGAC',
                              'AGATCAGAG',
                              'GACCCCGCG',
                              'GGATTAGAG',
                              'CGCCCGTAG',
                              'GCCCGACAG',
                              'GGCGAACAG',
                              'GCCGCTGAG',
                              'AGCGTAGAG',
                              'GGCTCGGAG',
                              'GGCTCGGAC',
                              'GGCTCGGAA',
                              'GGCTCGGAT',
                              'GAACCAGCG',
                              'GGGCCAGTC',
                              'GGGCCAGTA',
                              'GGGCCAGTG',
                              'AGACCCGAG',
                              'GTTCAAGAG',
                              'TTCCAAGAG',
                              'TCCCCGGAG',
                              'GGGCCAGTT',
                              'GGGAGAGAG',
                              'GGCCTACAC',
                              'GGCCTACAA',
                              'GGCCTACAG',
                              'GGCCAATCG',
                              'GCCGCATAG',
                              'GGCCTACAT',
                              'CGCCTGGAG',
                              'AGCCCGGAA',
                              'AGCCCGGAC',
                              'CGCCGAGAT',
                              'AGCCCGGAG',
                              'GGGCCTGCG',
                              'GGTCGAGAA',
                              'GGCCACTAG',
                              'CGCCGAGAA',
                              'GGCCCTACG',
                              'CGCCGAGAC',
                              'AGCCCGGAT',
                              'CGCCGAGAG',
                              'GGGCGAGGG',
                              'GGCCTATTG',
                              'GGCACAAGG',
                              'GGTCCTCAG',
                              'GGCTCACGG',
                              'GGCCCTGAT',
                              'GGTTGAGAG',
                              'CGCTCAGAT',
                              'CGAACAGAG',
                              'GGACCAATG',
                              'GGCCCTGAA',
                              'ACCCCAGAC',
                              'GGCACCAAG',
                              'ACCCCAGAA',
                              'ACCCCAGAG',
                              'ATGCCAGAG',
                              'TGCCGAGAG',
                              'GGTATAGAG',
                              'GGCCCCGCT',
                              'GGCCCTTCG',
                              'ACCCCAGAT',
                              'GGCCCCGCG',
                              'GGCCCCGCA',
                              'TTCCCGGAG',
                              'GGCCCCGCC',
                              'CGCACTGAG',
                              'GAGCGAGAG',
                              'GTCTCAGAC',
                              'GGCCGAATG',
                              'GTCTCAGAG',
                              'GTCCACGAG',
                              'GGCCGATGG',
                              'TGCCTAGTG',
                              'TGCTCAGAT',
                              'GTCTCAGAT',
                              'GCCCCTGCG',
                              'TGCTCAGAC',
                              'TGCTCAGAA',
                              'TGCTCAGAG',
                              'CGCCAAGAG',
                              'CGCCCAGAT',
                              'CGCCAAGAC',
                              'AGCACATAG',
                              'CGCCAAGAA',
                              'CGCTCAGAA',
                              'TGCACAGAT',
                              'AGGCCAGAT',
                              'CGCCCAGAG',
                              'CGCCAAGAT',
                              'CGCCCAGAA',
                              'CGCCCAGAC',
                              'AGGCCAGAC',
                              'TGCACAGAA',
                              'AGGCCAGAA',
                              'TGCACAGAC',
                              'GCGCAAGAG',
                              'TGCACAGAG',
                              'GCCCTAGTG',
                              'GGCCCTAAG',
                              'GGAACAGGG',
                              'GGTCCACAG',
                              'GGCCCGGCA',
                              'GGCCCGGCC',
                              'CGCTCAGAG',
                              'GCGCCAGGG',
                              'GGCCCGGCG',
                              'GGCCCTAAA',
                              'GGCACCGCG',
                              'GGCTCAATG',
                              'GGCCCGGCT',
                              'GGGCCCGTG',
                              'CGTCCACAG',
                              'GAGCCAGTG',
                              'TGCCATGAG',
                              'GTCCCCGGG',
                              'GGGGCCGAG',
                              'GGCGAAAAG',
                              'TGCCAAAAG',
                              'TGCCCAATG',
                              'GTCTGAGAG',
                              'GGGTCAAAG',
                              'AGCCTAGCG',
                              'GGGCCGTAG',
                              'GGCCAGGGG',
                              'GGCTGGGAG',
                              'GGACCGGAA',
                              'GGACCGGAC',
                              'CTCACAGAG',
                              'GGACCGGAG',
                              'TGACCAAAG',
                              'CAACCAGAG',
                              'GTCCCAGAG',
                              'GTCCCAGAA',
                              'GTCCCAGAC',
                              'GTCCCAGAT',
                              'GGTCCACAT',
                              'TTACCAGAG',
                              'GTGCCAGAT',
                              'TGCCGTGAG',
                              'GGTCCAAGG',
                              'CACCCGGAG',
                              'GTCCCAACG',
                              'GACCCGGCG',
                              'GTGCCAGAA',
                              'CGCCTTGAG',
                              'GTGCCAGAC',
                              'GTGCCAGAG',
                              'GGACTAGGG',
                              'GGACCATAG',
                              'GGACGCGAG',
                              'CGCGCGGAG',
                              'GGCCCGTAG',
                              'GGCCCGTAA',
                              'GGCCCGTAC',
                              'GACGCGGAG',
                              'GGTCGAAAG',
                              'GGTGCTGAG',
                              'GGCCCGTAT',
                              'GGGGAAGAG',
                              'AGCCCAAAG',
                              'GGACGTGAG',
                              'AGCCCAAAA',
                              'GGCCAAAAT',
                              'AGCCCAAAC',
                              'GCCGCAGCG',
                              'GGAGCAGTG',
                              'GGCCAAAAA',
                              'AGCCCAAAT',
                              'GGCCAAAAC',
                              'GGCCAAAAG',
                              'GGACGAGCG',
                              'AGCCAACAG',
                              'GGACTATAG',
                              'AGCGCAGTG',
                              'AGCTAAGAG',
                              'GGCCGTGAC',
                              'GGCCGTGAA',
                              'GGCCGTGAG',
                              'GGCTAAGGG',
                              'CGCTTAGAG',
                              'GGCGCCGTG',
                              'GGCCCATTA',
                              'GCGCCGGAG',
                              'GTTCCATAG',
                              'GGCCGTGAT',
                              'CCCCCACAG',
                              'GGGCTTGAG',
                              'GGCCTTAAG',
                              'GGTCAAAAG',
                              'TTCCCATAG',
                              'GGCCTCGCG',
                              'GACCGAGTG',
                              'GTCTTAGAG',
                              'GTCCCATCG',
                              'CGACCTGAG',
                              'GAACCAGTG',
                              'GGGCCAGAG',
                              'GGGCCAGAA',
                              'GGGCCAGAC',
                              'CGTCAAGAG',
                              'GGGCCAGAT',
                              'TACCCGGAG',
                              'GGGGCGGAG',
                              'TTCCCACAG',
                              'GCCACTGAG',
                              'GACCTCGAG',
                              'GGCCCAAGA',
                              'GGCCCAAGC',
                              'GGCCCAAGG',
                              'GGCCCAGGT',
                              'GCCCTAGCG',
                              'GGCCCAGGC',
                              'GGAACAGTG',
                              'GGCCCAGGA',
                              'GGCCCAGGG',
                              'GGCCCAAGT',
                              'GACCCCGTG',
                              'TGCCACGAG',
                              'AACCCCGAG',
                              'CGACGAGAG',
                              'AGTCCGGAG',
                              'TGCCCACCG',
                              'GGCGTAGCG',
                              'GACGCACAG',
                              'GTTGCAGAG',
                              'GGCCGCGAG',
                              'AGCCCACCG',
                              'GGCCGCGAC',
                              'GGCCGCGAA',
                              'GGCCGCGAT',
                              'GACTCACAG',
                              'GGCCAACCG',
                              'GGGCCATGG',
                              'GGACTTGAG',
                              'GGCACAGTT',
                              'GCGCCCGAG',
                              'GGTCCGAAG',
                              'GCCCCAGGA',
                              'GTCCCCAAG',
                              'GGTGCAGAG',
                              'GGTGCAGAC',
                              'TTTCCAGAG',
                              'GCCACAGAC',
                              'GGCACAGTG',
                              'GGCACAGTA',
                              'GGCACAGTC',
                              'GGTGCAGAT',
                              'TACTCAGAG',
                              'CGCACAAAG',
                              'TGCCCATAA',
                              'GGTACAGGG',
                              'TGCCCATAC',
                              'TGCCCATAG',
                              'TCCGCAGAG',
                              'GGTCACGAG',
                              'GGCATAGCG',
                              'AGACCAGGG',
                              'AGCTGAGAG',
                              'TGCCCAAGG',
                              'CACGCAGAG',
                              'GTCGCAGAA',
                              'GTCGCAGAC',
                              'GTCGCAGAG',
                              'CGCCTAGTG',
                              'GTCGCAGAT',
                              'GACCCGTAG',
                              'GGCCCTGAG',
                              'AGCCCTGAT',
                              'GGACCAAAG',
                              'GGACCAAAA',
                              'GGCTTCGAG',
                              'TGCCCTAAG',
                              'GGGTCCGAG',
                              'GGTTCGGAG',
                              'GGTTAAGAG',
                              'GTCGAAGAG',
                              'GGCGCACCG',
                              'GGACCAAAT',
                              'AGCCCTGAG',
                              'TGCCTGGAG',
                              'AGCCCTGAA',
                              'AGCCCTGAC',
                              'GGGACAAAG',
                              'CGCCCTAAG',
                              'GGACAAAAG',
                              'GAACCGGAG',
                              'GGCCAGGTG',
                              'GGGCGAGAA',
                              'GGCCAGGAA',
                              'GACCGAGCG',
                              'GGGCAACAG',
                              'CACCCACAG',
                              'CGAGCAGAG',
                              'GGCTGTGAG',
                              'GGCCAGGAC',
                              'GGCTCAAAG',
                              'GACACAGGG',
                              'GGGCGAGAT',
                              'TGCCCGGAC',
                              'AGCACCGAG',
                              'GGCCGGGAC',
                              'GGCCGGGAA',
                              'GGCCGGGAG',
                              'GGCGCTTAG',
                              'GGCCCTGAC',
                              'GGTTCACAG',
                              'GGCCTAAAT',
                              'CGTCCATAG',
                              'TGGCCAAAG',
                              'GACGCTGAG',
                              'GGCCTATAG',
                              'TGGCCAGCG',
                              'GGCCTATAA',
                              'GGCCTATAC',
                              'TGCCTCGAG',
                              'GAACGAGAG',
                              'CACCCATAG',
                              'GGAACGGAG',
                              'GGCCTATAT',
                              'GTCCCATAA',
                              'GGCCAGGAT',
                              'GGACCACCG',
                              'GGGCAGGAG',
                              'GACCAACAG',
                              'GGCCACCAG',
                              'GGTCCACCG',
                              'GGCACGGAT',
                              'GGCGCAAAG',
                              'GCCCCCGGG',
                              'GGCGCAAAC',
                              'TCGCCAGAG',
                              'GGCGCAAAA',
                              'GGCACGGAA',
                              'GGCACGGAC',
                              'GGCACGGAG',
                              'GCCACAGTG',
                              'GGCGCAAAT',
                              'AACCAAGAG',
                              'CGGCCTGAG',
                              'GACAAAGAG',
                              'CTCCAAGAG',
                              'GGCAGAGGG',
                              'AACCCGGAG',
                              'GGCCGACGG',
                              'AACCCATAG',
                              'GTACCAGAC',
                              'GTACCAGAA',
                              'GTACCAGAG',
                              'GACCCAAAT',
                              'CGCCCATAT',
                              'GGCGGGGAG',
                              'GTCCGAGGG',
                              'GGAAAAGAG',
                              'GTACCAGAT',
                              'CGCCCATAG',
                              'GACCCAAAG',
                              'GACCCAAAA',
                              'CGCCCATAC',
                              'GACCCAAAC',
                              'CGCCCATAA',
                              'TGCACAGTG',
                              'GGCTCATTG',
                              'GGCTCAGTC',
                              'GTCGCGGAG',
                              'TCCCCTGAG',
                              'GGCCATGAT',
                              'GCTCCGGAG',
                              'GGCCATGAA',
                              'GGCCATGAC',
                              'TGCGTAGAG',
                              'GGCCATGAG',
                              'AGCACACAG',
                              'GTGCCAGTG',
                              'GTCTCAGAA',
                              'GCCCAAGGG',
                              'GGGCTAGCG',
                              'GCACCACAG',
                              'GGCCCATTT',
                              'GGCTCTGAT',
                              'GTTCCCGAG',
                              'GCACCTGAG',
                              'GGCCCAGTG',
                              'GACCCACCG',
                              'GGCCCAGTC',
                              'GGCCAGCAG',
                              'GGCCCAGTA',
                              'GGCTCTGAG',
                              'GGCACGGGG',
                              'GGCCCATTG',
                              'CTCGCAGAG',
                              'GGAACCGAG',
                              'GGACCATGG',
                              'GGCTCTGAA',
                              'GTACGAGAG',
                              'GGGCAAAAG',
                              'GGCGAATAG',
                              'GACACGGAG',
                              'GGAGCAGAT',
                              'CACCGAGAG',
                              'AGACCAGTG',
                              'GGAGCAGAA',
                              'GGCTCGAAG',
                              'GGAGCAGAG',
                              'AGCCCATTG',
                              'CCCCTAGAG',
                              'GGCTCGGCG',
                              'GGACCCCAG',
                              'GTCCTATAG',
                              'GCTTCAGAG',
                              'CGGTCAGAG',
                              'AGCCCACTG',
                              'GACCCTTAG',
                              'ACCCTAGAG',
                              'CACTCAGAG',
                              'GGGGCATAG',
                              'TGCGCCGAG',
                              'GGCACTTAG',
                              'TATCCAGAG',
                              'GGCAAACAG',
                              'GGCCTACCG',
                              'GCCCCGAAG',
                              'GTCCCTGTG',
                              'TGCGCAGTG',
                              'GGTAGAGAG',
                              'GTCCCGGTG',
                              'AGCCCGGCG',
                              'GTACTAGAG',
                              'CGGCTAGAG',
                              'GGCCAATAA',
                              'GTCCTAGCG',
                              'GGCCAATAC',
                              'GCCCGAAAG',
                              'GGCCAATAG',
                              'CGCCGAGCG',
                              'GGCACACCG',
                              'GGTCCTAAG',
                              'AGCCGTGAG',
                              'GGCCAATAT',
                              'CGCCGATAG',
                              'AGCGCAGAC',
                              'AGCGCAGAA',
                              'AGCGCAGAG',
                              'GGCACAAAT',
                              'GGAATAGAG',
                              'CGCCCCGAA',
                              'GGCCCTGCG',
                              'GGCCCTGCA',
                              'AGCGCAGAT',
                              'GACCAAAAG',
                              'GGCACAAAG',
                              'GGCACAAAA',
                              'GGCACAAAC',
                              'GGCTCAGAG',
                              'GGCCCCGAT',
                              'AGCTCCGAG',
                              'CTCCCGGAG',
                              'GGCCCCGAC',
                              'GGCCCCGAA',
                              'GGCCCCGAG',
                              'GGACCCGAA',
                              'GTCTCAGCG',
                              'GGCGAAGCG',
                              'GGCGCGGCG',
                              'GACCCCCAG',
                              'GTCACGGAG',
                              'TGCCGAGAT',
                              'TGCACAGGG',
                              'CGCCAAGGG',
                              'GGCCCGCAC',
                              'GGCCTCGAC',
                              'GGCCTGCAG',
                              'GGCCCGGAG',
                              'GTGGCAGAG',
                              'GGCCCGGAC',
                              'GGCCCGCAA',
                              'GGCCCGGAA',
                              'GGCCCTAAC',
                              'AGACCGGAG',
                              'GGCTCACAC',
                              'ACCCCAAAG',
                              'GGCCCTAAT',
                              'GATCCGGAG',
                              'GGCCCGGAT',
                              'CGCTAAGAG',
                              'CGCCTAGGG',
                              'CGTCCCGAG',
                              'AGGCCAGAG',
                              'GGCTCACAG',
                              'GACACAGTG',
                              'TGCCTTGAG',
                              'GACTCCGAG',
                              'ATCCCCGAG',
                              'AGCCTAGAA',
                              'AGCCTAGAC',
                              'CGTCGAGAG',
                              'AGCCTAGAG',
                              'TGCCGAGGG',
                              'TGCTCAAAG',
                              'GCCACAGAT',
                              'GACTAAGAG',
                              'AGCCTAGAT',
                              'ATCCCAGTG',
                              'GGACCGGCG',
                              'GCCACAGAA',
                              'GCCACAGAG',
                              'GCTCCTGAG',
                              'GGCCTCGAT',
                              'GGCCAGGAG',
                              'GCCTCAGAT',
                              'GGCTAAAAG',
                              'AGGCCAAAG',
                              'AGACAAGAG',
                              'ACCCCCGAG',
                              'CGCACATAG',
                              'GCCTCAGAA',
                              'GCCTCAGAC',
                              'GGCACCCAG',
                              'GCCTCAGAG',
                              'GGCTGAGCG',
                              'GCTCCAGTG',
                              'TGCCAAGGG',
                              'GGCTCAGTG',
                              'CGCGAAGAG',
                              'CGCCCAGCA',
                              'GGTCCATGG',
                              'GGCTCAGTA',
                              'GTGCCAGGG',
                              'GCCCCAATG',
                              'GACACCGAG',
                              'CGACCAAAG',
                              'CGGGCAGAG',
                              'GGCTCAGTT',
                              'CGCCCAGCT',
                              'CGCCGGGAG',
                              'GGGTCTGAG',
                              'GCTGCAGAG',
                              'CGCGCACAG',
                              'GGCCCCTCG',
                              'GGCCGACTG',
                              'ACCTCAGAG',
                              'TGACCAGAT',
                              'AGCCCAACG',
                              'GCGCCATAG',
                              'GGCCTCGAA',
                              'GGCGCCTAG',
                              'GGCCTCGAG',
                              'TGACCAGAG',
                              'TGACCAGAA',
                              'TGACCAGAC',
                              'GGCCAAACG',
                              'GGCCTAGAA',
                              'GGCCATGTG',
                              'AGCGGAGAG',
                              'GGCCCACGT',
                              'AGCCGAGCG',
                              'GGCCCACGG',
                              'GAACCTGAG',
                              'GGCCCACGC',
                              'GTACCATAG',
                              'GGCCCACGA',
                              'TTCCCTGAG',
                              'GGCTTAGGG',
                              'GGCCGTGCG',
                              'GGCTAAGAG',
                              'GGCTAAGAC',
                              'GGCTAAGAA',
                              'CAGCCAGAG',
                              'GGCTAAGAT',
                              'GGCCTAGAT',
                              'GGCCACGTG',
                              'TGCTTAGAG',
                              'CGGCCAGGG',
                              'GGTCGACAG',
                              'GGGCGATAG',
                              'GCCCCAGCG',
                              'GGCGTGGAG',
                              'GGGCCAGCC',
                              'GGGCCAGCA',
                              'GGGCCAGCG',
                              'TTCCCCGAG',
                              'GGGCCTGTG',
                              'GTCCCATAT',
                              'GACCTATAG',
                              'GGGCCAGCT',
                              'AGTCCTGAG',
                              'GTCCTAGTG',
                              'CGACCACAG',
                              'TGTCCAGTG',
                              'GGGGCACAG',
                              'GGCCCATGG',
                              'GGCCCATGA',
                              'GGCCCATGC',
                              'AGCCAAAAG',
                              'GGCCCATGT',
                              'CGCACGGAG',
                              'GGCAAAGTG',
                              'GCGTCAGAG',
                              'GGACGAGAC',
                              'GGACGAGAA',
                              'GACCTTGAG',
                              'GGCCTTGAT',
                              'AGACTAGAG',
                              'CCCCCAAAG',
                              'GGCCTTGAA',
                              'GTCCCTAAG',
                              'GGCCTTGAG',
                              'GGACGAGAT',
                              'GGCTCCGAC',
                              'GGGCCATAT',
                              'GGCCAACAT',
                              'GGCTCCGAG',
                              'GCCCTGGAG',
                              'TGTCCTGAG',
                              'GGCTCACCG',
                              'AGCCCACAT',
                              'GGCCAACAG',
                              'GGGCCATAG',
                              'CTGCCAGAG',
                              'GGCCAACAC',
                              'GGGCCATAC',
                              'GGCCAACAA',
                              'GGGCCATAA',
                              'AGCCCACAC',
                              'AGCCCACAA',
                              'AGCCCACAG',
                              'GGTGCAGCG',
                              'GCCCCAAAA',
                              'TGCCCAGAT',
                              'TGCTGAGAG',
                              'GCCCCAGAT',
                              'GCCACAGCG',
                              'GGGACATAG',
                              'TGCCCAGAC',
                              'TGCCCAGAA',
                              'TGCCCAGAG',
                              'GCCCCAGAG',
                              'GCCCCAGAA',
                              'CGACCGGAG',
                              'GCCCCAGAC',
                              'GGCGGAGGG',
                              'GTCCCGGGG',
                              'TGCCCATCG',
                              'TGCAGAGAG',
                              'GGCTCAGGT',
                              'GCCGCGGAG',
                              'GGCATAGAA',
                              'GGCATAGAC',
                              'GGTCCAATG',
                              'GGCATAGAG',
                              'GGCTCAGGC',
                              'CGCCTCGAG',
                              'GGCTCAGGG',
                              'GTCGCAGGG',
                              'CGCTCGGAG',
                              'GGACCGGAT',
                              'GAGCCAGCG',
                              'GGTCTAAAG',
                              'CGCCGTGAG',
                              'TGCCCAGCT',
                              'GGCGCCGCG',
                              'GGCACACGG',
                              'TGCCTAGGG',
                              'GGACCAACG',
                              'GACCGAGAC',
                              'GACCGAGAA',
                              'GGCCGAGGT',
                              'GACCGAGAG',
                              'GGCGCTGAT',
                              'GTACCACAG',
                              'GGGCGAGCG',
                              'CGCCCCGTG',
                              'GGCCGAGGG',
                              'GGCGCTGAC',
                              'GGCGCTGAA',
                              'GGCCGAGGC',
                              'GGCGCTGAG',
                              'GGCCGAGGA',
                              'GACCGAGAT',
                              'CGCACACAG',
                              'TGCGAAGAG',
                              'TGCACATAG',
                              'GTCCAGGAG',
                              'GTCTCAGTG',
                              'TGCCGATAG',
                              'GAACCAAAG',
                              'GTTCCGGAG',
                              'TGCTCTGAG',
                              'AGCCCGGTG',
                              'GGACCACAC',
                              'CGCCCTCAG',
                              'GGACCACAA',
                              'GGACCACAG',
                              'GGCCTATCG',
                              'GCCTCAGCG',
                              'GGCGCAAGG',
                              'GGTCCACAA',
                              'GACCTGGAG',
                              'GGCACGGCG',
                              'GGATCATAG',
                              'GGTGCAAAG',
                              'GGCACGTAG',
                              'TCCCCAGTG',
                              'GTCCCATAG',
                              'CGCCCCTAG',
                              'TTCTCAGAG',
                              'GGGACGGAG',
                              'GGTGGAGAG',
                              'ATCGCAGAG',
                              'GATCAAGAG',
                              'GTCCCATAC',
                              'TGCTCGGAG',
                              'AAGCCAGAG',
                              'GGCCGGAAG',
                              'GGACCGGTG',
                              'GGCCCACTA',
                              'GGCCCACTG',
                              'TGGCCACAG',
                              'GACGCCGAG',
                              'GCCCACGAG',
                              'CGCCCATCG',
                              'GACCCAACG',
                              'GGCAGCGAG',
                              'GGCCCACTT',
                              'GGGCTGGAG',
                              'GGCCTGAAG',
                              'GGCGTAGTG',
                              'AGGCCATAG',
                              'CGGCCCGAG',
                              'TACACAGAG',
                              'GTACCAGCG',
                              'AGCCCCAAG',
                              'AGCTCGGAG',
                              'AGTCCAGTG',
                              'GTACAAGAG',
                              'TACCCAGCG',
                              'GTAACAGAG',
                              'GCCCCCTAG',
                              'GGCCACAAG',
                              'CGGCCAGTG',
                              'GCCCCGCAG',
                              'CACCCCGAG',
                              'TGCCCAGGA',
                              'CGCCACGAG',
                              'GCTCCAGGG',
                              'GGCGCCGAG',
                              'GCCGAAGAG',
                              'AGCTCAGAG',
                              'GGCCCAATG',
                              'GAGCCAGAG',
                              'GGCCAGTAG',
                              'GGCTCTGGG',
                              'GCCTTAGAG',
                              'ATCCAAGAG',
                              'CGCCTAAAG',
                              'GGCCAAATG',
                              'GACCCACAC',
                              'GACCCACAA',
                              'CGCTCTGAG',
                              'TGCCCTGTG',
                              'GACTCAGTG',
                              'GGATCAGTG',
                              'GGACCCAAG',
                              'GGTCGATAG',
                              'GACCCACAT',
                              'GCCTCAGTG',
                              'GTCCATGAG',
                              'GGGCCTTAG',
                              'AGTCAAGAG',
                              'CACCAAGAG',
                              'GGCCCTCAG',
                              'GCTACAGAG',
                              'TTCCTAGAG',
                              'AGCACAAAG',
                              'GAACCAGGG',
                              'AGCATAGAG',
                              'GTTCCAAAG',
                              'GGCCGCGTG',
                              'GGCTCCGTG',
                              'GGTCAATAG',
                              'GGCGCATTG',
                              'GGCGGAGTG',
                              'GTCCGACAG',
                              'GTCCTAGAC',
                              'GTCCTAGAA',
                              'GTCCTAGAG',
                              'GGCGCGGGG',
                              'GGCCAATGG',
                              'GGCGCAGTA',
                              'GGCTCGCAG',
                              'GGAGCAAAG',
                              'CGCCCCAAG',
                              'AGCGCAAAG',
                              'GGTTCTGAG',
                              'GGGCCTGGG',
                              'GGCGCAGTT',
                              'CCCCGAGAG',
                              'GGCTTAGTG',
                              'GGGTCATAG',
                              'GGAACACAG',
                              'GGCACAACG',
                              'GCGCCACAG',
                              'TTCCCAGGG',
                              'TGCTCACAG',
                              'ACCACAGAG',
                              'AGTCCCGAG',
                              'GATGCAGAG',
                              'AGCGCTGAG',
                              'GGGTCACAG',
                              'GTCCCCCAG',
                              'TGCCAGGAG',
                              'GGCGAAGAG',
                              'GGCCGAGTC',
                              'GGCCGAGTA',
                              'GGCGAAGAC',
                              'GGCCGAGTG',
                              'GGTTCCGAG',
                              'GGCGCGGAT',
                              'GCCCGAGAT',
                              'GGCGCGGAC',
                              'GGCGCGGAA',
                              'GGCGAAGAT',
                              'GGCGCGGAG',
                              'GGCCGAGTT',
                              'GCCCGAGAA',
                              'GCCCGAGAC',
                              'GCCCGAGAG',
                              'GGATCCGAG',
                              'GCCCCTGGG',
                              'CTCCCAGCG',
                              'GGCCCTCAT',
                              'GCCCCTGAG',
                              'GGCCCTCAC',
                              'AGGACAGAG',
                              'GGCCCTCAA',
                              'AGTCCAGGG',
                              'AGCTCAGTG',
                              'GACCCTGGG',
                              'GGCCCGGGG',
                              'GGCCCGGGA',
                              'GGCCCGGGC',
                              'GGCACACAC',
                              'GGTCCTTAG',
                              'GGCACACAA',
                              'GGCACACAG',
                              'GTCCCAAAT',
                              'TGACCCGAG',
                              'GGCCCGGGT',
                              'GGCACACAT',
                              'GGCCCCAAA',
                              'GGCCCCAAC',
                              'TGTCCCGAG',
                              'GGCCCCAAG',
                              'TGCCCCGGG',
                              'GGTACACAG',
                              'ATCCTAGAG',
                              'GGCAAAAAG',
                              'GGCGGCGAG',
                              'AGCCTAGGG',
                              'TCCCCAGGG',
                              'GTTCCTGAG',
                              'GGCCCCTTG',
                              'GGCGCCCAG',
                              'GGCCTTGAC',
                              'GCAGCAGAG',
                              'GGCCAGGCG',
                              'CGCCCGGGG',
                              'GCCTCAGGG',
                              'GCCCTTGAG',
                              'GGATAAGAG',
                              'AGCCGAGAT',
                              'GATCCTGAG',
                              'GGCTGAGAT',
                              'AGCCGAGAC',
                              'CGCCCAGGG',
                              'AGCCGAGAA',
                              'GGCTCCGAA',
                              'AGCCGAGAG',
                              'CGTCCTGAG',
                              'GGCTGAGAA',
                              'GGCTGAGAC',
                              'GGCATAGTG',
                              'GGCCCCCCG',
                              'GGCTGAGAG',
                              'GACCCGGGG',
                              'GGTCAGGAG',
                              'GACCCCAAG',
                              'TGACTAGAG',
                              'GACGAAGAG',
                              'AGCCCATCG',
                              'TTCCCAGAA',
                              'GGGCCGGGG',
                              'AGGCCACAG',
                              'CTCTCAGAG',
                              'TGCACACAG',
                              'CGCCCACAT',
                              'GGCCTCGGG',
                              'GGGCCCAAG',
                              'GCCACAAAG',
                              'CGCCCACAA',
                              'CGCCCGCAG',
                              'CGCCCACAC',
                              'TGCCGAGAC',
                              'CGCCCACAG',
                              'TGACCAGCG',
                              'GAGCCATAG',
                              'GCCGCAGGG',
                              'GCACAAGAG',
                              'GTGCCACAG',
                              'GGCTAAGCG',
                              'GGCTCCGAT',
                              'GGCCAAGTA',
                              'GGTCCGTAG',
                              'GGCCAAGTC',
                              'GGGACTGAG',
                              'AGGCTAGAG',
                              'CGCGCAAAG',
                              'GGCCGCGGG',
                              'GGTCTTGAG',
                              'GAACAAGAG',
                              'GACGCAGGG',
                              'GGTCGCGAG',
                              'GACCAATAG',
                              'GGCCAAGTG',
                              'GGGCTAGTG',
                              'TACCTAGAG',
                              'GGTCCAGAC',
                              'GGTCCAGAA',
                              'GGTCCAGAG',
                              'GGCCAACGG',
                              'GGGCCATCG',
                              'GGCCCGATG',
                              'CGCCCAACG',
                              'GGACTAAAG',
                              'GGTCCAGAT',
                              'TGCCAATAG',
                              'GTCCCATGG',
                              'GTCCAAGTG',
                              'GGCCCAACG',
                              'GGCCCAACA',
                              'GTCACAGTG',
                              'GGCCCAACC',
                              'GGCCGGTAG',
                              'GGCCCAACT',
                              'CGACAAGAG',
                              'GGCCAAGTT',
                              'GGCTCTGTG',
                              'GGACGAAAG',
                              'GGCGTAGGG',
                              'GGCCCACAG',
                              'GGCCCACAA',
                              'GGCCCACAC',
                              'ACCCGAGAG',
                              'AGACCACAG',
                              'GACCCAGGC',
                              'GGCCTTGCG',
                              'GGCCCACAT',
                              'ACACCAGAG',
                              'GGCTCACAA',
                              'TGGCCATAG',
                              'AGCCCATAG',
                              'AGCCTATAG',
                              'AGCCCATAA',
                              'AGCCCATAC',
                              'GTCCGAGCG',
                              'AGCCCATAT',
                              'GGCCTACTG',
                              'GGCTCACAT',
                              'CGCTCAGCG',
                              'GGCTCCGCG',
                              'GGACACGAG',
                              'GTACCTGAG',
                              'GGCGCACGG',
                              'GCCCCAGCT',
                              'GTCGTAGAG',
                              'GGTGTAGAG',
                              'GCCCCAGCC',
                              'GCCCCAGCA',
                              'GGGGCAAAG',
                              'GCCCATGAG',
                              'GACCTAGCG',
                              'GTCCCGGAG',
                              'GTCCCGGAC',
                              'GTCCCGGAA',
                              'AGGCGAGAG',
                              'TGCCCAGCA',
                              'GACACACAG',
                              'TGCCCAGCC',
                              'TGCCCAGCG',
                              'GGCATAGGG',
                              'GTCCCGGAT',
                              'GGTCAAGTG',
                              'CACCCAGGG',
                              'CCCCCAGAA',
                              'CCCCCAGAC',
                              'TGCCCGGCG',
                              'CCCCCAGAG',
                              'GCCCTAGGG',
                              'GTCCCTCAG',
                              'GATCCAGGG',
                              'CGCCCTTAG',
                              'CCCCCAGAT',
                              'GGCATAGAT',
                              'GGCTCAAGG',
                              'TGCCTAGAT',
                              'TTCCCAGTG',
                              'GGTCCCGTG',
                              'GAGCCAGAA',
                              'GGTTTAGAG',
                              'GGCGCCGAC',
                              'GGCGCCGAA',
                              'TGCCTAGAG',
                              'TGCCTAGAA',
                              'TGCCTAGAC',
                              'GGCGCCGAT',
                              'GAGCCAGAT',
                              'CGTCCGGAG',
                              'GCCAGAGAG',
                              'CGCCCCGGG',
                              'GACTCGGAG',
                              'GGCTCTTAG',
                              'GGCCCAGCT',
                              'GTCCTGGAG',
                              'GGCCGAGAG',
                              'GCCTGAGAG',
                              'GGCCGAGAA',
                              'GGCCGAGAC',
                              'GGCGCTGCG',
                              'GGCCGTGTG',
                              'GGCCGAGAT',
                              'AGCCCAGCA',
                              'AGCCCAGCC',
                              'AGCCCAGCG',
                              'GGCAGATAG',
                              'AGCCCAGCT',
                              'GACCCGGTG',
                              'GGTTCAGGG',
                              'GGCTCGTAG',
                              'GGCTCAGGA',
                              'AGGCCTGAG',
                              'GGGCCACAA',
                              'GGGCCACAC',
                              'GCTCGAGAG',
                              'GGGCCACAG',
                              'GGAACATAG',
                              'GACCGACAG',
                              'GATCCCGAG',
                              'CTCCCTGAG',
                              'GGGCCACAT',
                              'GGGCCGGTG',
                              'GACATAGAG',
                              'ATCACAGAG',
                              'GCACCAGGG',
                              'GCCGCAGTG',
                              'GCCTCATAG',
                              'GGCAGAGCG',
                              'GGACCAGCT',
                              'GTCATAGAG',
                              'TGGCCCGAG',
                              'GGCGTTGAG',
                              'TACCGAGAG',
                              'CGCAAAGAG',
                              'GGACCAGCA',
                              'GGACCAGCC',
                              'TACGCAGAG',
                              'GGACCAGCG',
                              'TACCCAGAA',
                              'GGCCCTCTG',
                              'TACCCAGAC',
                              'GGGCCAACG',
                              'AGCTCAGGG',
                              'GAACCCGAG',
                              'GTCCGCGAG',
                              'TACCCAGAT',
                              'AGTCGAGAG',
                              'GGCTGAGTG',
                              'GGCCGTTAG',
                              'GCTCCAGAG',
                              'GGTGCACAG',
                              'CACCCAGCG',
                              'GACCCATTG',
                              'GGCCAAGAC',
                              'GTTCCACAG',
                              'GACCGAAAG',
                              'AACCGAGAG',
                              'GACCCCTAG',
                              'GGCCCAATT',
                              'TCCCTAGAG',
                              'TGCCAAGTG',
                              'AACCCAAAG',
                              'GACGCAGTG',
                              'GGCACATTG',
                              'GGCCCAATA',
                              'GGCCCAATC',
                              'TGCCCCCAG',
                              'GCCACGGAG',
                              'GGGCTAGGG',
                              'CGCAGAGAG',
                              'GGGCGAGAG',
                              'GGCCTGGAG',
                              'TGTCCAGGG',
                              'TGCCAAGCG',
                              'GGGGCAGAC',
                              'GGCCTGGAC',
                              'GTCGCAGTG',
                              'TGCCCCAAG',
                              'GGCCTGGAT',
                              'GGCCACGAA',
                              'TGACCTGAG',
                              'GGCCACGAC',
                              'GGCCACGAG',
                              'GGACAGGAG',
                              'GTCCCAATG',
                              'CGCCCGGTG',
                              'GGCCACGAT',
                              'GTCCCCTAG',
                              'GGGCCTAAG',
                              'GACCACGAG',
                              'GTCCCAGCT',
                              'TGACCGGAG',
                              'GGTCTATAG',
                              'AGCCAAGTG',
                              'GGCCGGCAG',
                              'CGACCCGAG',
                              'GCCGCCGAG',
                              'GGACAAGAT',
                              'GCCCTCGAG',
                              'GGCGCAGTC',
                              'GGCCTAAGG',
                              'TGACCACAG',
                              'GGCCTACGG',
                              'TCACCAGAG',
                              'GTCCCATTG',
                              'GGGCTATAG',
                              'GTCCAAGGG',
                              'TGCCCAGTA',
                              'GGTCCAGTG',
                              'TACCCAGAG',
                              'GGTCCAGTA',
                              'GGTCCAGTC',
                              'GGATCACAG',
                              'AGCCCGGGG',
                              'GGGTCAGCG',
                              'AGGGCAGAG',
                              'GGTCCAGTT',
                              'AGCGCGGAG',
                              'GGACCCTAG',
                              'TGGCAAGAG',
                              'AGCACAGAT',
                              'GTCTCTGAG',
                              'GGCCCGCGG',
                              'AGCACAGAA',
                              'GTTTCAGAG',
                              'AGCACAGAC',
                              'GGTCCAAAG',
                              'GGTCCCGGG',
                              'AGCACAGAG',
                              'GGCCTTGTG',
                              'AGTCCACAG',
                              'GACCGATAG',
                              'TGCCCCGAA',
                              'CGCCTACAG',
                              'TGCCCCTAG',
                              'GGCGCACTG',
                              'GAGCCAGAC',
                              'GGAGCACAG',
                              'GACCAGGAG',
                              'GGCAAATAG',
                              'GTCCGAAAG',
                              'TGCCCGGAG',
                              'GTCTCAGGG',
                              'GTCCTAAAG',
                              'GGACTAGAA',
                              'GCCCGAGCG',
                              'GGCGAAGGG',
                              'TGCGCGGAG',
                              'GGTCAAGCG',
                              'CTCCCAGAA',
                              'CTCCCAGAC',
                              'CGCCCCCAG',
                              'CTCCCAGAG',
                              'AGCGCACAG',
                              'GTCACAGAG',
                              'GTCACAGAC',
                              'GTCACAGAA',
                              'CTCCCAGAT',
                              'GGCCCTCCG',
                              'GTCACAGAT',
                              'GGAACAAAG',
                              'GGACGGGAG',
                              'GCGCCAAAG',
                              'CGACCAGGG',
                              'GCCCCCGTG',
                              'CGGCAAGAG',
                              'GATCCATAG',
                              'AGCAGAGAG',
                              'GGAGTAGAG',
                              'GCCACACAG',
                              'GACCCAGCT',
                              'GGACATGAG',
                              'CGCACCGAG',
                              'GTACCGGAG',
                              'GACCCAGCA',
                              'GACCCAGCC',
                              'GACCCAGCG',
                              'GGACTCGAG',
                              'CGTCCAGTG',
                              'TCCCCAGAG',
                              'GGCACAGCT',
                              'GGCGCGAAG',
                              'TCCCCAGAA',
                              'GCTCCAGAC',
                              'GGCACAGCA',
                              'TCCCCAGAT',
                              'AGCTCTGAG',
                              'TGAACAGAG',
                              'GCTCCAGAA',
                              'GGCAGAGTG',
                              'GGTCTAGCG',
                              'GGACCGGGG',
                              'GCACCAGTG',
                              'GACTCTGAG',
                              'CACCCAGTG',
                              'CGCGGAGAG',
                              'GGCCCCCAT',
                              'AGCCATGAG',
                              'GGCCCCCAG',
                              'GGCTGAGGG',
                              'GGCCCCCAC',
                              'GCTCCAGAT',
                              'GGCCCCCAA',
                              'GGGCAATAG',
                              'GACCCATGG',
                              'GGTTCAGTG',
                              'GAGACAGAG',
                              'CGGCCGGAG',
                              'CGGCCACAG',
                              'GGGCAAGCG',
                              'GGCGTAAAG',
                              'GGCCTAGTA',
                              'GGCCGATAT',
                              'GGCACATGG',
                              'GGCCTAGTG',
                              'AAACCAGAG',
                              'GCCCAAGTG',
                              'GGCCAAGGT',
                              'CGCCCACCG',
                              'GGCCGATAG',
                              'GGCCGATAC',
                              'GGCCGATAA',
                              'TGCCCTGGG',
                              'GGCCAAGGC',
                              'ATCCCGGAG',
                              'TGTCAAGAG',
                              'GGCCAAGGG',
                              'GCCGCAGAT',
                              'GGCACCTAG',
                              'GGTGCATAG',
                              'GACCATGAG',
                              'GGGACCGAG',
                              'GGAGGAGAG',
                              'GTCCGAGTG',
                              'TGCCGAGCG',
                              'GGCCAAAGG',
                              'GGCTTAGCG',
                              'GGTCCGGCG',
                              'TGCTCCGAG',
                              'AGAGCAGAG',
                              'GGCTGATAG',
                              'GACCCAGTT',
                              'AGCCGAAAG',
                              'GCCCCTAAG',
                              'TAACCAGAG',
                              'GGTCCAGCT',
                              'GGCCGTGGG',
                              'AACACAGAG',
                              'ACCCCGGAG',
                              'GGTCTACAG',
                              'GCCCTAAAG',
                              'GACCGCGAG',
                              'GATCCACAG',
                              'GGCCCCGTT',
                              'TGCCGAAAG',
                              'GGTCCCCAG',
                              'CGCCCAAAG',
                              'GGCGCAGGT',
                              'CGCCCAAAC',
                              'GGGAAAGAG',
                              'TGCGCATAG',
                              'GGCCAATTG',
                              'GGCGCAGGG',
                              'GGTCCAGCA',
                              'GCCCCCAAG',
                              'CGCCCAAAT',
                              'GGCGCAGGC',
                              'GGCGCAGGA',
                              'TGCGCAGGG',
                              'GCGACAGAG',
                              'GGACCACAT',
                              'GTGCCAAAG',
                              'GGCCCAAAT',
                              'GGACAAGGG',
                              'GGCGGACAG',
                              'GACCCAGTG',
                              'TGTCCGGAG',
                              'CGCCGAGTG',
                              'AGGCCGGAG',
                              'GGCCCAAAG',
                              'GACCCAGTA',
                              'GGCCCACCC',
                              'GGCCCACCA',
                              'AGTCCAGAG',
                              'GGCCCACCG',
                              'CGCGCAGAA',
                              'GTCCCACAT',
                              'CGCGCAGAC',
                              'GACCCAGTC',
                              'GCACCGGAG',
                              'CGCGCAGAG',
                              'AACCCAGAA',
                              'GGGCTCGAG',
                              'GGTCGTGAG',
                              'GGCCCACCT',
                              'GGTCTCGAG',
                              'CGCGCAGAT',
                              'GGGACACAG',
                              'GGCCCCTAT',
                              'GGAGCTGAG',
                              'GGACCTGTG',
                              'GTCCGTGAG',
                              'GCCCAGGAG',
                              'GGCCCCTAA',
                              'GGCCCCTAC',
                              'TGCACAAAG',
                              'GGCCCCTAG',
                              'GGGCCCCAG',
                              'GCCCCAGGG',
                              'GACCTAGAT',
                              'TGCCTAGCG',
                              'GGCCTGGAA',
                              'GGACGACAG',
                              'GACCTAGAA',
                              'GACCTAGAC',
                              'GACCTAGAG',
                              'GCCCCAGGC',
                              'GGCGTCGAG',
                              'GGCTCAGCT',
                              'CGCACAGGG',
                              'GGCGTAGAC',
                              'GGCTCAGCG',
                              'GGCTCAGCC',
                              'GGCTCAGCA',
                              'GTCCCGGCG',
                              'TGCCCATGG',
                              'GGCGTAGAG',
                              'CCCCCAGCG',
                              'AGACCAAAG',
                              'GGTACCGAG',
                              'GGCTCAAAT',
                              'TGCCCGGAT',
                              'GAACCATAG',
                              'GGACGATAG',
                              'GGCCTCCAG',
                              'TGCCCCGTG',
                              'TGCCCGGAA',
                              'GGTGCAGAA',
                              'GGCTCAAAC',
                              'GGCTCAAAA',
                              'GGCGCCGGG',
                              'CGCTCCGAG',
                              'GTTCGAGAG',
                              'GCCCCAGGT',
                              'GGGACAGGG',
                              'GGGGCAGAA',
                              'GACCGGGAG',
                              'GGTCTGGAG',
                              'GGGGCAGAG',
                              'CGTCCAGGG',
                              'GGCCCCATG',
                              'AGCCTAGTG',
                              'ATCCCATAG',
                              'GTCCCGAAG',
                              'GGGGCAGAT',
                              'GTCCCTGCG',
                              'GACACAAAG',
                              'ATCCCAGGG',
                              'GGCCGAGCC',
                              'GCCCGAGTG',
                              'GGCCGAGCA',
                              'GGCGAAGTG',
                              'GGCCGAGCG',
                              'GTCAGAGAG',
                              'GCCCAAAAG',
                              'GGGCTAAAG',
                              'GGATCGGAG',
                              'AGCCCAGAG',
                              'AGCCCAGAC',
                              'GGCCCTGTT',
                              'AGCCCAGAA',
                              'GGCCCTGTA',
                              'GGCCCTGTC',
                              'AGCCCAGAT',
                              'AGCCAGGAG',
                              'GGCCCTGTG',
                              'GACCCATCG',
                              'ACCCCTGAG',
                              'GATTCAGAG',
                              'ACCCAAGAG',
                              'GGGCCACCG',
                              'GCCCGTGAG',
                              'CGCCAAAAG',
                              'GTCTCCGAG',
                              'GGCCCTTTG',
                              'GGGGGAGAG',
                              'GTCGCAAAG',
                              'GCCCCACTG',
                              'ATACCAGAG',
                              'GGCATACAG',
                              'TGCTAAGAG',
                              'GGCAGAGAC',
                              'GGCAGAGAA',
                              'GGCAGAGAG',
                              'CGCCCTGAT',
                              'AGCAAAGAG',
                              'GGACCAGAT',
                              'CGCCCTGAC',
                              'CGCCCTGAA',
                              'CGCCCTGAG',
                              'GGCAGAGAT',
                              'GGACCAGAG',
                              'GGACCAGAC',
                              'CCTCCAGAG',
                              'GGACCAGAA',
                              'GGGCCAAAG',
                              'TACCCAGGG',
                              'CTCCCACAG',
                              'GGGCCAAAC',
                              'AGCTCAGAC',
                              'GGGCCAAAA',
                              'AGCTCAGAA',
                              'CCCGCAGAG',
                              'AGCTCAGAT',
                              'GGGCCAAAT',
                              'GGCACCGTG',
                              'AGGCCCGAG',
                              'TGCACCGAG',
                              'GCTCCAGCG',
                              'AGCCTACAG',
                              'GCCCCACAT',
                              'TGCCCATAT',
                              'GTGCCCGAG',
                              'GGCCTCGTG',
                              'GCGCCAGAG',
                              'GCCCCACAA',
                              'GCCCCACAC',
                              'GCCCCACAG',
                              'GGTCGGGAG',
                              'GTCTAAGAG',
                              'GGCACTGGG',
                              'GGCACATCG',
                              'CGCCGCGAG',
                              'AACCCACAG',
                              'TGCCAAGAA',
                              'CCGCCAGAG',
                              'TGCCAAGAC',
                              'TACCCATAG',
                              'GGACAAGTG',
                              'TGCCAAGAG',
                              'GGGCGGGAG',
                              'GGCCTAGGA',
                              'GGCCTAGGC',
                              'GGCCTAGGG',
                              'TGCCAAGAT',
                              'CGCGTAGAG',
                              'GATCCAGTG',
                              'GGCCTAGGT',
                              'GGCCACGCG',
                              'CGCCCAGCC',
                              'TGGACAGAG',
                              'GGCCCACTC',
                              'AGCCGGGAG',
                              'GCCCCAACG',
                              'GCCGCAAAG',
                              'GCTCCAAAG',
                              'GGACCTGGG',
                              'GGCCCGTTG',
                              'GGCCCGAGG',
                              'TGCCTATAG',
                              'TGGCCTGAG',
                              'GGGGTAGAG',
                              'TGCCCAGTT',
                              'AGCCCGTAG',
                              'GCCCCGGGG',
                              'TGCCCATTG',
                              'TGCCCAGTG',
                              'CGCTCAGTG',
                              'TGCCCAGTC',
                              'GTCTCGGAG',
                              'GTCCAAGAG',
                              'GCACTAGAG',
                              'GTCCAAGAC',
                              'GTCCAAGAA',
                              'GTCCAAGAT',
                              'GGCCGAGCT',
                              'GCCCCTCAG',
                              'GGGCCTCAG',
                              'TTCCCAGCG',
                              'CGCATAGAG',
                              'GGGCGTGAG',
                              'AGCACAGCG',
                              'GTTCCAGCG',
                              'GGACCAAAC',
                              'AACCCAGAT',
                              'TGCGCTGAG',
                              'GTCGCCGAG',
                              'GGCAACGAG',
                              'AACCCAGAC',
                              'AACCCAGAG',
                              'GAACCACAG',
                              'GGTCAAGAT',
                              'GGGACAGTG',
                              'GACCAAGGG',
                              'GGACGAGAG',
                              'GGTCAAGAA',
                              'GGTCAAGAC',
                              'TGATCAGAG',
                              'GGTCAAGAG',
                              'GGCGCTGTG',
                              'CTCCCAGGG',
                              'GTCACAGCG',
                              'GGAACAGAT',
                              'GCGCCAGAT',
                              'GGAACAGAA',
                              'AGTCCAGCG',
                              'GGAACAGAC',
                              'GGAACAGAG',
                              'GCGCCAGAA',
                              'TGTCTAGAG',
                              'GCGCCAGAC',
                              'GACCTAGTG',
                              'AGGCAAGAG',
                              'GGGTCAGAA',
                              'GGCGACGAG',
                              'GGGTCAGAC',
                              'GGCCATAAG',
                              'GGGTCAGAG',
                              'ACCCCAGGG',
                              'AGTCCAAAG',
                              'GGGTCAGAT',
                              'TGCCCCGCG',
                              'GGCTAATAG',
                              'GACCCAGAT',
                              'TGCCCGGTG',
                              'CCCCCAGTG',
                              'GGCGCAATG',
                              'ATCCCTGAG',
                              'GATCGAGAG',
                              'GACCCAGAG',
                              'GGGCGAGTG',
                              'GACCCGCAG',
                              'GACCCAGAC',
                              'GACCCAGAA',
                              'GGCCGTAAG',
                              'GGCACGGTG'
                              ])
        self.assertEqual(sorted(bioinfo1.neighbors(test_pattern, test_distance)), test_output)


class FrequentWordsWithMismatchesTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_frequent_words_mismatches_sample(self):
        """
        """
        test_text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        test_k = 4
        test_d = 1
        test_output = ['ATGC', 'ATGT', 'GATG']
        self.assertEqual(sorted(bioinfo1.frequent_words_with_mismatches(test_text, test_k, test_d)),
                         test_output)

    def test_frequent_words_mismatches_1(self):
        """
        This dataset checks that your code includes k-mers that do not actually appear in Text.
        Notice here that, although none of the output k-mers except for AA actually appear in Text,
        they are all valid because they appear in Text with up to 1 mismatch (i.e. 0 or 1 mismatch).
        """
        test_text = 'AAAAAAAAAA'
        test_k = 2
        test_d = 1
        test_output = sorted(['AA', 'AC', 'AG', 'CA', 'AT', 'GA', 'TA'])
        self.assertEqual(sorted(bioinfo1.frequent_words_with_mismatches(test_text, test_k, test_d)),
                         test_output)

    def test_frequent_words_mismatches_2(self):
        """
        This dataset makes sure that your code is not accidentally swapping k and d.
        """
        test_text = 'AGTCAGTC'
        test_k = 4
        test_d = 2
        test_output = sorted(['TCTC', 'CGGC', 'AAGC', 'TGTG', 'GGCC', 'AGGT', 'ATCC', 'ACTG',
                              'ACAC', 'AGAG', 'ATTA', 'TGAC', 'AATT', 'CGTT', 'GTTC', 'GGTA',
                              'AGCA', 'CATC'])
        self.assertEqual(sorted(bioinfo1.frequent_words_with_mismatches(test_text, test_k, test_d)),
                         test_output)

    def test_frequent_words_mismatches_3(self):
        """
        This dataset makes sure you are not finding patterns in the Reverse Complement of Text
        (that is the next problem, “Frequent Words with Mismatches and Reverse Complements
        Problem”).
        """
        test_text = 'AATTAATTGGTAGGTAGGTA'
        test_k = 4
        test_d = 0
        test_output = ['GGTA']
        self.assertEqual(sorted(bioinfo1.frequent_words_with_mismatches(test_text, test_k, test_d)),
                         test_output)

    def test_frequent_words_mismatches_4(self):
        """
        This dataset first checks that k-mers with exactly d mismatches are being found. Then, it
        checks that k-mers with less than d mismatches are being allowed (i.e. you are not only
        allowing k-mers with exactly d mismatches). Next, it checks that you are not returning too
        few k-mers. Last, it checks that you are not returning too many k-mers.
        """
        test_text = 'ATA'
        test_k = 3
        test_d = 1
        test_output = sorted(['GTA', 'ACA', 'AAA', 'ATC', 'ATA', 'AGA', 'ATT', 'CTA', 'TTA',
                              'ATG'])
        self.assertEqual(sorted(bioinfo1.frequent_words_with_mismatches(test_text, test_k, test_d)),
                         test_output)

    def test_frequent_words_mismatches_5(self):
        """
        This dataset checks that your code is not looking for k-mers in the Reverse Complement of
        Text.
        """
        test_text = 'AAT'
        test_k = 3
        test_d = 0
        test_output = ['AAT']
        self.assertEqual(sorted(bioinfo1.frequent_words_with_mismatches(test_text, test_k, test_d)),
                         test_output)

    def test_frequent_words_mismatches_6(self):
        """
        This dataset checks that your code correctly delimiting your output (i.e. using spaces) and
        verifies that our k-mers are actually of length k.
        """
        test_text = 'TAGCG'
        test_k = 2
        test_d = 1
        test_output = ['GG', 'TG']
        self.assertEqual(sorted(bioinfo1.frequent_words_with_mismatches(test_text, test_k, test_d)),
                         test_output)


if __name__ == "__main__":
        unittest.main()
