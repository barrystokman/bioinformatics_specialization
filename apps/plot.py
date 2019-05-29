import matplotlib.pyplot as plt
from bioinformatics_1.functions import skew, minimum_skew


class SkewPlot:

    def __init__(self, genome):
        self.genome = genome.genome
        self.organism_name = genome.file_name.split('.')[0].replace('_', ' ')

    @property
    def skew(self):
        return skew(self.genome)

    @property
    def position_axis(self):
        return [pos for pos in range(len(self.genome) + 1)]

    @property
    def skew_axis(self):
        return self.skew

    @property
    def nucleotides(self):
        return [nucleotide for nucleotide in self.genome]

    @property
    def minimum_skew_positions(self):
        return minimum_skew(self.genome)

    @property
    def minimum_skew_values(self):
        return [self.skew[position] for position in self.minimum_skew_positions]

    def generate_plot(self):
        position_axis = self.position_axis
        skew_axis = self.skew_axis
        plt.plot(position_axis, skew_axis, 'b')
        plt.xlim([0, len(position_axis)])
        plt.title(f'Skew Diagram for {self.organism_name}')
        plt.xlabel('position')
        plt.ylabel('skew')

    def show_minimum_skew(self):
        minimum_skew_positions = self.minimum_skew_positions
        minimum_skew_values = self.minimum_skew_values
        plt.plot(minimum_skew_positions, minimum_skew_values, 'ro')
        for x, y in zip(minimum_skew_positions, minimum_skew_values):
            plt.text(x, y, f"{x}, {y}")

    def show_plot(self):
        plt.show()
