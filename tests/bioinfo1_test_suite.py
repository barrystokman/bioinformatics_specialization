import unittest

import bioinfo1_unit_test as bioinfo1


def suite():

    suite = unittest.TestSuite()

    suite.addTest(bioinfo1.PatternCountTests('test_pattern_count_1'))
    suite.addTest(bioinfo1.PatternCountTests('test_pattern_count_2'))
    suite.addTest(bioinfo1.PatternCountTests('test_pattern_count_3'))
    suite.addTest(bioinfo1.PatternCountTests('test_pattern_count_4'))
    suite.addTest(bioinfo1.PatternCountTests('test_pattern_count_5'))
    suite.addTest(bioinfo1.PatternCountTests('test_pattern_count_6'))

    suite.addTest(bioinfo1.FrequentWordsTests('test_frequent_words_1'))
    suite.addTest(bioinfo1.FrequentWordsTests('test_frequent_words_2'))
    suite.addTest(bioinfo1.FrequentWordsTests('test_frequent_words_3'))
    suite.addTest(bioinfo1.FrequentWordsTests('test_frequent_words_4'))

    suite.addTest(bioinfo1.ReverseComplementTests('test_reverse_complement_1'))

    suite.addTest(bioinfo1.PatternMatchingProblemTests('test_pattern_matching_problem_1'))
    suite.addTest(bioinfo1.PatternMatchingProblemTests('test_pattern_matching_problem_2'))
    suite.addTest(bioinfo1.PatternMatchingProblemTests('test_pattern_matching_problem_3'))
    suite.addTest(bioinfo1.PatternMatchingProblemTests('test_pattern_matching_problem_4'))

    suite.addTest(bioinfo1.ClumpFindingProblemTests('test_clump_finding_problem_1'))
    suite.addTest(bioinfo1.ClumpFindingProblemTests('test_clump_finding_problem_2'))
    suite.addTest(bioinfo1.ClumpFindingProblemTests('test_clump_finding_problem_3'))

    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
