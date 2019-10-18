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
        frequent_kmers = frequent_words(text, k_mer_length, mismatch)
        count = kmers_count(frequent_kmers[0], text, mismatch)

        return frequent_kmers, count

    def count_kmers_in_window(self, kmer, start=None, window_size=500, mismatch=1):
        start = 0 if start is None else start
        end = start + window_size
        text = self.genome[start:end]
        return kmers_count(kmer, text, mismatch)

    @staticmethod
    def kmer_count(genome_segment, kmer):
        return
