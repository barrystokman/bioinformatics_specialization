"""
Class reads datasets specified in bioinfo CLI commands and returns an object that allows parsing of
variables.
Main class: ReadDataset: instantiating this class reads the file specified on the CLI
Subclasses: xxxDataset: returns variables <var1>, <var2>, ..., <varn> for function <func>
"""
from helper.constants import DATASET_PATH


class ReadDataset:

    def __init__(self, file_name):
        with open(DATASET_PATH + file_name) as f:
            self.lines = [line for line in f]

    def __str__(self):
        return_str = ''
        for line in self.lines:
            return_str += line
        return_str += '\nThis file contains ' + str(len(self.lines)) + ' lines.'

        return return_str


class PatternCountDataset(ReadDataset):

    def get_text(self):
        return str(self.lines[1]).rstrip()

    def get_pattern(self):
        return self.lines[2].rstrip()

    def get_expected_result(self):
        return int(self.lines[4])

    def get_text_challenge(self):
        return self.lines[0].rstrip()

    def get_pattern_challenge(self):
        return self.lines[1].rstrip()


class FrequentWordsDataset(ReadDataset):

    def get_text(self):
        return self.lines[1].rstrip()

    def get_k(self):
        return int(self.lines[2])

    def get_expected_result(self):
        return self.lines[4].rstrip()

    def get_text_challenge(self):
        return self.lines[0].rstrip()

    def get_k_challenge(self):
        return int(self.lines[1])


class ReverseComplementDataset(ReadDataset):

    def get_pattern(self):
        return self.lines[1].rstrip()

    def get_expected_result(self):
        return self.lines[3].rstrip()

    def get_pattern_challenge(self):
        return self.lines[0].rstrip()


class PatternMatchingDataset(ReadDataset):

    def get_pattern(self):
        return self.lines[1].rstrip()

    def get_genome(self):
        return self.lines[2].rstrip()

    def get_expected_result(self):
        return self.lines[4].rstrip()

    def get_pattern_challenge(self):
        return self.lines[0].rstrip()

    def get_genome_challenge(self):
        return self.lines[1].rstrip()


class ClumpFindingDataset(ReadDataset):

    def get_genome(self):
        return self.lines[1].rstrip()

    def get_variables(self):
        return [int(var) for var in self.lines[2].split(' ')]

    def get_k(self):
        return self.get_variables()[0]

    def get_l(self):
        return self.get_variables()[1]

    def get_t(self):
        return self.get_variables()[2]

    def get_expected_result(self):
        return self.lines[4].rstrip()

    def get_genome_challenge(self):
        return self.lines[0].rstrip()

    def get_variables_challenge(self):
        return [int(var) for var in self.lines[1].split(' ')]

    def get_k_challenge(self):
        return self.get_variables_challenge()[0]

    def get_l_challenge(self):
        return self.get_variables_challenge()[1]

    def get_t_challenge(self):
        return self.get_variables_challenge()[2]


class ComputingFrequenciesDataset(ReadDataset):

    def get_text(self):
        return self.lines[1].rstrip()

    def get_k(self):
        return int(self.lines[2])

    def get_expected_result(self):
        return self.lines[4].rstrip()

    def get_text_challenge(self):
        return self.lines[0].rstrip()

    def get_k_challenge(self):
        return int(self.lines[1])


class PatternToNumberDataset(ReadDataset):

    def get_pattern(self):
        return self.lines[1].rstrip()

    def get_expected_result(self):
        return int(self.lines[3])

    def get_pattern_challenge(self):
        return self.lines[0].rstrip()


class NumberToPatternDataset(ReadDataset):

    def get_number(self):
        return int(self.lines[1])

    def get_k(self):
        return int(self.lines[2])

    def get_expected_result(self):
        return self.lines[4]

    def get_number_challenge(self):
        return int(self.lines[0])

    def get_k_challenge(self):
        return int(self.lines[1])


class MinimumSkewDataset(ReadDataset):

    def get_genome(self):
        return self.lines[1].rstrip()

    def get_expected_result(self):
        return self.lines[3].rstrip()

    def get_genome_challenge(self):
        return self.lines[0].rstrip()


class HammingDistanceDataset(ReadDataset):

    def get_string1(self):
        return self.lines[1].rstrip()

    def get_string2(self):
        return self.lines[2].rstrip()

    def get_expected_result(self):
        return int(self.lines[4])

    def get_string1_challenge(self):
        return self.lines[0].rstrip()

    def get_string2_challenge(self):
        return self.lines[1].rstrip()


class ApproxMatchDataset(ReadDataset):

    def get_pattern(self):
        return self.lines[1].rstrip()

    def get_text(self):
        return self.lines[2].rstrip()

    def get_d(self):
        return int(self.lines[3])

    def get_expected_result(self):
        return self.lines[5].rstrip()

    def get_pattern_challenge(self):
        return self.lines[0].rstrip()

    def get_text_challenge(self):
        return self.lines[1].rstrip()

    def get_d_challenge(self):
        return int(self.lines[2])


def main():
    pass


if __name__ == '__main__':
    main()
