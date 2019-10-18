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

    @staticmethod
    def skew_text(positions, values):
        for pos, val in enumerate(zip(positions, values)):
            pass

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

        x_fraction=0.6
        y_fraction = 0.8
        y_offset = -0.05

        for x, y in zip(minimum_skew_positions, minimum_skew_values):
            plt.annotate(
                f"({x}, {y})",
                xy=(x, y), xycoords='data',
                xytext=(x_fraction, y_fraction), textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', facecolor='black')
            )
            y_fraction += y_offset

    def show_ori_candidate(self):
        pass

    def show_plot(self):
        plt.show()


class CountBarPlot:

    def __init__(self, genome, count):
        self.genome = genome.genome
        self.organism_name = genome.file_name.split('.')[0].replace('_', ' ')
        self.count = count

    @property
    def position_axis(self):
        return [pos for pos in range(len(self.genome) - 500)]

    @property
    def count_axis(self):
        return self.count

    def generate_plot(self):
        position_axis = self.position_axis
        count_axis = self.count_axis
        plt.plot(position_axis, count_axis, 'b')
        plt.xlim([0, len(position_axis)])
        plt.title(f'Count Diagram for {self.organism_name}')
        plt.xlabel('position')
        plt.ylabel('count')

    def show_plot(self):
        plt.show()
