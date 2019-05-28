import click

from bioinformatics_1 import functions as bioinfo1
from apps.cli import result_color, cli_output, get_correct_result
from apps.plot import SkewPlot
import apps.dataset_reader as dsr


@click.group()
@click.option('--code-challenge', '-c', is_flag=True,
              help='Use this option if attempting a Code Challenge')
@click.pass_context
def biocli(context, code_challenge):
    """
    Collection of CLI callable bioinformatics functions to make Code Challenges quick and easy to
    perform.
    """
    if context.obj is None:
        context.obj = {}

    context.obj['CHALLENGE'] = code_challenge


@biocli.command('pattern-count')
@click.argument('dataset', required=True)
@click.pass_context
def pattern_count(context, dataset, sort_result=False, listing=False):
    """
    Runs pattern_count(text, pattern). The input variables 'text' and 'pattern' are read from the
    DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.PatternCountDataset(dataset, challenge)

    text = data.text
    pattern = data.pattern
    correct_result = get_correct_result(data, challenge)

    args = text, pattern
    func = bioinfo1.pattern_count

    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('frequent-words')
@click.argument('dataset', required=True)
@click.pass_context
def frequent_words(context, dataset, sort_result=False, listing=False):
    """
    Runs frequent_words_by_sorting(text, k). The input variables 'text' and 'k' are read from the
    DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.FrequentWordsDataset(dataset, challenge)

    text = data.text
    k = data.k
    correct_result = get_correct_result(data, challenge)

    args = text, k
    func = bioinfo1.frequent_words_by_sorting
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('reverse-complement')
@click.argument('dataset', required=True)
@click.pass_context
def reverse_complement(context, dataset, sort_result=False, listing=False):
    """
    Runs reverse_complement(pattern). The input variable 'pattern' is read from the DATASET
    argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.ReverseComplementDataset(dataset, challenge)

    pattern = data.pattern
    correct_result = get_correct_result(data, challenge)

    args = [pattern]
    func = bioinfo1.reverse_complement
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('pattern-matching')
@click.argument('dataset', required=True)
@click.pass_context
def pattern_matching(context, dataset, sort_result=False, listing=False):
    """
    Runs pattern_matching_problem(pattern, genome). The input variables 'pattern' and 'genome' are
    read from the DATASET argument, where DATASET is the text file containing the input data. This
    function is only available in Code Challenge mode.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.PatternMatchingDataset(dataset, challenge)

    pattern = data.pattern
    genome = data.genome
    correct_result = get_correct_result(data, challenge)

    args = pattern, genome
    func = bioinfo1.pattern_matching_problem
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('clump-finding')
@click.argument('dataset', required=True)
@click.pass_context
def clump_finding(context, dataset, sort_result=False, listing=False):
    """
    Runs better_clump_finding(genome, k, l, t). The input variables 'genome', 'k', 'l', and 't' are
    read from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.ClumpFindingDataset(dataset, challenge)

    genome = data.genome
    k = data.k
    l = data.l
    t = data.t
    correct_result = get_correct_result(data, challenge)

    args = genome, k, l, t
    func = bioinfo1.better_clump_finding
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('computing-frequencies')
@click.argument('dataset', required=True)
@click.pass_context
def computing_frequencies(context, dataset, sort_result=False, listing=False):
    """
    Runs computing_frequencies(text, k). The input variables 'text' and 'k' are read from the
    DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.ComputingFrequenciesDataset(dataset, challenge)

    text = data.text
    k = data.k
    correct_result = get_correct_result(data, challenge)

    args = text, k
    func = bioinfo1.computing_frequencies
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('pattern-to-number')
@click.argument('dataset', required=True)
@click.pass_context
def pattern_to_number(context, dataset, sort_result=False, listing=False):
    """
    Runs pattern_to_number(pattern). The input variable 'pattern' is read from the DATASET
    argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.PatternToNumberDataset(dataset, challenge)

    pattern = data.pattern
    correct_result = get_correct_result(data, challenge)

    args = [pattern]
    func = bioinfo1.pattern_to_number
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('number-to-pattern')
@click.argument('dataset', required=True)
@click.pass_context
def number_to_pattern(context, dataset, sort_result=False, listing=False):
    """
    Runs number_to_pattern(number, k). The input variables 'number' and 'k' are read from the
    DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.NumberToPatternDataset(dataset, challenge)

    number = data.number
    k = data.k
    correct_result = get_correct_result(data, challenge)

    args = number, k
    func = bioinfo1.number_to_pattern
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('minimum-skew')
@click.argument('dataset', required=True)
@click.pass_context
def minimum_skew(context, dataset, sort_result=False, listing=False):
    """
    Runs minimum_skew(genome). The input variable 'genome' is read from the DATASET argument, where
    DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.MinimumSkewDataset(dataset, challenge)

    genome = data.genome
    correct_result = get_correct_result(data, challenge)

    args = [genome]
    func = bioinfo1.minimum_skew
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('hamming-distance')
@click.argument('dataset', required=True)
@click.pass_context
def hamming_distance(context, dataset, sort_result=False, listing=False):
    """
    Runs hamming_distance(string1, string2). The input variables 'string1' and 'string2' are read
    from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.HammingDistanceDataset(dataset, challenge)

    string1 = data.string1
    string2 = data.string2
    correct_result = get_correct_result(data, challenge)

    args = string1, string2
    func = bioinfo1.hamming_distance
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('approx-matching')
@click.argument('dataset', required=True)
@click.pass_context
def approx_matching(context, dataset, sort_result=False, listing=False):
    """
    Runs approx_pattern_match(pattern, text, d). The input variables 'pattern' and 'text' are read from
    the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.ApproxMatchDataset(dataset, challenge)

    pattern = data.pattern
    text = data.text
    d = data.d
    correct_result = get_correct_result(data, challenge)

    args = pattern, text, d
    func = bioinfo1.approx_pattern_match
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('approx-count')
@click.argument('dataset', required=True)
@click.pass_context
def approx_count(context, dataset, sort_result=False, listing=False):
    """
    Runs approx_pattern_count(pattern, text, d). The input variables 'pattern', 'text' and 'd' are
    read from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.ApproxCountDataset(dataset, challenge)

    pattern = data.pattern
    text = data.text
    d = data.d
    correct_result = get_correct_result(data, challenge)

    args = pattern, text, d
    func = bioinfo1.approx_pattern_count
    cli_output(challenge, correct_result, sort_result, func, args)


@biocli.command('neighbors')
@click.argument('dataset', required=True)
@click.pass_context
def neighbors(context, dataset, sort_result=True, listing=True):
    """
    Runs neighbors(pattern, d). The input variables 'pattern' and 'd' are read from the DATASET
    argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.NeighborsDataset(dataset, challenge)

    pattern = data.pattern
    d = data.d
    correct_result = get_correct_result(data, challenge)

    args = pattern, d
    func = bioinfo1.neighbors
    cli_output(challenge, correct_result, sort_result, listing, func, args)


@biocli.command('frequent-words-mismatches')
@click.argument('dataset', required=True)
@click.pass_context
def frequent_words_mismatches(context, dataset, sort_result=False, listing=False):
    """
    Runs frequent_words_with_mismatches(pattern, k, d). The input variables 'pattern', 'k' and 'd'
    are read from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.FrequentWordsMismatchesDataset(dataset, challenge)

    text = data.text
    k = data.k
    d = data.d
    correct_result = get_correct_result(data, challenge)

    args = text, k, d
    func = bioinfo1.frequent_words_with_mismatches_sorting
    cli_output(challenge, correct_result, sort_result, listing, func, args)


@biocli.command('frequent-words-mismatches-and-reverse-complement')
@click.argument('dataset', required=True)
@click.pass_context
def frequent_words_mismatches_and_reverse_complement(context, dataset, sort_result=False, listing=False):
    """
    Runs frequent_words_with_mismatches_and_reverse_complement(pattern, k, d). The input variables
    'pattern', 'k' and 'd' are read from the DATASET argument, where DATASET is the text file
    containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.FrequentWordsMismatchesDataset(dataset, challenge)

    text = data.text
    k = data.k
    d = data.d

    correct_result = get_correct_result(data, challenge)

    args = text, k, d
    func = bioinfo1.frequent_words_with_mismatches_and_reverse_complement
    cli_output(challenge, correct_result, sort_result, listing, func, args)


@biocli.command('skew-plot')
@click.option('-m', '--minimum-skew', required=False, is_flag=True, help='adds minimum skew indicators to plot')
@click.argument('genome', required=True)
def skew_plot(minimum_skew, genome):
    """
    Generates a skew plot of the provided genome
    """
    skew_genome = dsr.Genome(genome).genome
    SkewPlot(skew_genome).generate_plot()
    if minimum_skew:
        SkewPlot(skew_genome).show_minimum_skew()
    SkewPlot(skew_genome).show_plot()


if __name__ == '__main__':
    biocli(obj={})
