import matplotlib.pyplot as plt
from bioinformatics_1.functions import skew


class SkewPlot:

    def __init__(self, genome):
        self.genome = genome
        self.skew = skew(genome)[1:]

    @property
    def position_axis(self):
        return [pos for pos in range(len(self.genome))]

    @property
    def skew_axis(self):
        return self.skew

    @property
    def nucleotides(self):
        return [nucleotide for nucleotide in self.genome]

    @property
    def minimum_skew(self):
        pass

    def generate_plot(self):
        plt.plot(self.position_axis, self.skew_axis)
        plt.xlim([0, len(self.position_axis)])
        plt.title('Skew Diagram')
        plt.xlabel('position')
        plt.ylabel('skew')
        plt.show()
