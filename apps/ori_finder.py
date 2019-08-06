from bioinformatics_1.functions import minimum_skew
from bioinformatics_1.functions import frequent_words_with_mismatches_and_reverse_complement \
                                       as frequent_words
from bioinformatics_1.functions import approx_pattern_and_reverse_complement_count as kmers_count


class OriFinder:

    def __init__(self, genome):
        self.genome = genome.genome

    @property
    def ori_candidate(self):
        """
        finds starting position of ori candidate
        """
        return minimum_skew(self.genome)[0]

    def find_frequent_kmers(self, start=None, k_mer_length=9, window_size=500, mismatch=1):
        start = self.ori_candidate if start is None else start
        end = start + window_size
        text = self.genome[start:end]
        frequent_kmers = frequent_words(self.genome[start:end], k_mer_length, mismatch)
        kmer = frequent_kmers[0]
        kmers_count = kmers_count(kmer, text, mismatch)
