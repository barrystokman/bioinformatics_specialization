"""
Class reads datasets specified in bioinfo CLI commands and returns an object that allows parsing of
variables.
Main class: ReadDataset: instantiating this class reads the file specified on the CLI
Subclasses: xxxDataset: returns variables <var1>, <var2>, ..., <varn> for function <func>
"""
from apps.constants import DATASET_PATH, NUCLEOTIDES


class ReadDataset:

    def __init__(self, file_name, challenge=False):
        self.file_name = file_name
        self.challenge = challenge
        with open(DATASET_PATH + file_name) as f:
            self.lines = [line for line in f]

    def __str__(self):
        return_str = ''
        for line in self.lines:
            return_str += line
        return_str += '\nThis file contains ' + str(len(self.lines)) + ' lines.'
        return_str += '\nThe challenge flag is ' + str(self.challenge) + '.'

        return return_str


class PatternCountDataset(ReadDataset):

    @property
    def text(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return str(self.lines[1]).rstrip()

    @property
    def pattern(self):
        if self.challenge:
            return self.lines[1].rstrip()
        else:
            return self.lines[2].rstrip()

    @property
    def result(self):
        return int(self.lines[4])


class FrequentWordsDataset(ReadDataset):

    @property
    def text(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def k(self):
        if self.challenge:
            return int(self.lines[1])
        else:
            return int(self.lines[2])

    @property
    def result(self):
        return self.lines[4].rstrip()


class ReverseComplementDataset(ReadDataset):

    @property
    def pattern(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def result(self):
        return self.lines[3].rstrip()


class PatternMatchingDataset(ReadDataset):

    @property
    def pattern(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def genome(self):
        if self.challenge:
            return self.lines[1].rstrip()
        else:
            return self.lines[2].rstrip()

    @property
    def result(self):
        return self.lines[4].rstrip()


class ClumpFindingDataset(ReadDataset):

    @property
    def genome(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    def _get_variables(self):
        if self.challenge:
            return [int(var) for var in self.lines[1].split(' ')]
        else:
            return [int(var) for var in self.lines[2].split(' ')]

    @property
    def k(self):
        return self._get_variables()[0]

    @property
    def l(self):
        return self._get_variables()[1]

    @property
    def t(self):
        return self._get_variables()[2]

    @property
    def result(self):
        return self.lines[4].rstrip()


class ComputingFrequenciesDataset(ReadDataset):

    @property
    def text(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def k(self):
        if self.challenge:
            return int(self.lines[1])
        else:
            return int(self.lines[2])

    @property
    def result(self):
        return self.lines[4].rstrip()


class PatternToNumberDataset(ReadDataset):

    @property
    def pattern(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def result(self):
        return int(self.lines[3])


class NumberToPatternDataset(ReadDataset):

    @property
    def number(self):
        if self.challenge:
            return int(self.lines[0])
        else:
            return int(self.lines[1])

    @property
    def k(self):
        if self.challenge:
            return int(self.lines[1])
        else:
            return int(self.lines[2])

    @property
    def result(self):
        return self.lines[4]


class MinimumSkewDataset(ReadDataset):

    @property
    def genome(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def result(self):
        return self.lines[3].rstrip()


class HammingDistanceDataset(ReadDataset):

    @property
    def string1(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def string2(self):
        if self.challenge:
            return self.lines[1].rstrip()
        else:
            return self.lines[2].rstrip()

    @property
    def result(self):
        return int(self.lines[4])


class ApproxMatchDataset(ReadDataset):

    @property
    def pattern(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def text(self):
        if self.challenge:
            return self.lines[1].rstrip()
        else:
            return self.lines[2].rstrip()

    @property
    def d(self):
        if self.challenge:
            return int(self.lines[2])
        else:
            return int(self.lines[3])

    @property
    def result(self):
        return self.lines[5].rstrip()


class ApproxCountDataset(ReadDataset):

    @property
    def pattern(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def text(self):
        if self.challenge:
            return self.lines[1].rstrip()
        else:
            return self.lines[2].rstrip()

    @property
    def d(self):
        if self.challenge:
            return int(self.lines[2])
        else:
            return int(self.lines[3])

    @property
    def result(self):
        return int(self.lines[5])


class NeighborsDataset(ReadDataset):

    @property
    def pattern(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def d(self):
        if self.challenge:
            return int(self.lines[1])
        else:
            return int(self.lines[2])

    @property
    def result(self):
        result = []
        for neighbor in self.lines[4:len(self.lines)]:
            result.append(neighbor.rstrip())
        return '\n'.join(sorted(result))


class FrequentWordsMismatchesDataset(ReadDataset):

    @property
    def text(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    def _get_variables(self):
        if self.challenge:
            return [int(var) for var in self.lines[1].split(' ')]
        else:
            return [int(var) for var in self.lines[2].split(' ')]

    @property
    def k(self):
        return self._get_variables()[0]

    @property
    def d(self):
        return self._get_variables()[1]

    @property
    def result(self):
        return ' '.join(sorted([result.rstrip() for result in self.lines[4].split(' ')]))


class Genome(ReadDataset):

    @property
    def genome(self):
        if len(self.lines) == 1:
            return self.lines[0].rstrip()
        else:
            return ''.join([line.rstrip() for line in self.lines[1:]])


class MotifEnumerationDataset(ReadDataset):

    def __output_line_index(self):
        return self.lines.index('Output\n')

    @property
    def dna(self):

        dna = []

        if self.challenge:
            for line in self.lines[1:]:
                dna.append(line.rstrip())
        else:
            output_index = self.__output_line_index()
            for line in self.lines[2:output_index]:
                dna.append(line.rstrip())

        return dna

    @property
    def k(self):
        if self.challenge:
            return int(self.lines[0].split(' ')[0])
        else:
            return int(self.lines[1].split(' ')[0])

    @property
    def d(self):
        if self.challenge:
            return int(self.lines[0].split(' ')[1].rstrip())
        else:
            return int(self.lines[1].split(' ')[1].rstrip())

    @property
    def result(self):
        output_index = self.__output_line_index()
        result = []

        for line in self.lines[output_index + 1:]:
            result.append(line.rstrip())

        return ' '.join(result)


class DistanceBetweenPatternAndString(ReadDataset):

    def __output_line_index(self):
        return self.lines.index('Output\n')

    @property
    def pattern(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def dna(self):

        if self.challenge:
            return self.lines[1].rstrip().split(' ')
        else:
            return self.lines[2].rstrip().split(' ')

    @property
    def result(self):
        output_index = self.__output_line_index()

        return int(self.lines[output_index + 1])


class MedianString(ReadDataset):

    def __output_line_index(self):
        return self.lines.index('Output\n')

    @property
    def dna(self):

        dna = []

        if self.challenge:
            for line in self.lines[1:]:
                dna.append(line.rstrip())
        else:
            output_index = self.__output_line_index()
            for line in self.lines[2:output_index]:
                dna.append(line.rstrip())

        return dna

    @property
    def k(self):
        if self.challenge:
            return int(self.lines[0].rstrip())
        else:
            return int(self.lines[1].rstrip())

    @property
    def result(self):
        output_index = self.__output_line_index()

        return self.lines[output_index + 1].rstrip()


class ProfileMostProbableKmer(ReadDataset):

    @property
    def text(self):
        if self.challenge:
            return self.lines[0].rstrip()
        else:
            return self.lines[1].rstrip()

    @property
    def k(self):
        if self.challenge:
            return int(self.lines[1].rstrip())
        else:
            return int(self.lines[2].rstrip())

    @property
    def profile(self):
        if self.challenge:
            probabilities = [[float(p) for p in row.split(' ')] for row in self.lines[2:6]]
        else:
            probabilities = [[float(p) for p in row.split(' ')] for row in self.lines[3:7]]

        return dict(zip(NUCLEOTIDES, probabilities))

    @property
    def result(self):

        return self.lines[8].rstrip()


class GreedyMotifSearch(ReadDataset):

    def __output_line_index(self):
        return self.lines.index('Output\n')

    @property
    def dna(self):
        if self.challenge:
            return [value.rstrip() for value in self.lines[1:]]
        else:
            output_line_index = self.__output_line_index()
            return [value.rstrip() for value in self.lines[2:output_line_index]]

    @property
    def k(self):
        if self.challenge:
            return int(self.lines[0].rstrip().split(' ')[0])
        else:
            return int(self.lines[1].rstrip().split(' ')[0])

    @property
    def t(self):
        if self.challenge:
            return int(self.lines[0].rstrip().split(' ')[1])
        else:
            return int(self.lines[1].rstrip().split(' ')[1])

    @property
    def result(self):
        output_line_index = self.__output_line_index()

        return ' '.join([result.rstrip() for result in self.lines[output_line_index+1:]])


def main():
    pass


if __name__ == '__main__':
    main()
