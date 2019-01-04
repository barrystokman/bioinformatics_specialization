import click

from bioinformatics_1 import functions as bioinfo1
from helper.cli_helper import result_color
from helper.dataset_reader import (PatternCountDataset,
                                   FrequentWordsDataset,
                                   ReverseComplementDataset,
                                   PatternMatchingDataset,
                                   ClumpFindingDataset,
                                   ComputingFrequenciesDataset,
                                   PatternToNumberDataset,
                                   NumberToPatternDataset,
                                   MinimumSkewDataset,
                                   HammingDistanceDataset,
                                   ApproxMatchDataset,
                                   ApproxCountDataset,
                                   NeighborsDataset,
                                   FrequentWordsMismatchesDataset
                                   )


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
def pattern_count(context, dataset):
    """
    Runs pattern_count(text, pattern). The input variables 'text' and 'pattern' are read from the
    DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = PatternCountDataset(dataset, challenge)

    text = data.text
    pattern = data.pattern

    if challenge:
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.pattern_count(text, pattern)}",
                               fg="yellow", bold=True))
    else:
        correct_result = data.result
        result = bioinfo1.pattern_count(text, pattern)
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('frequent-words')
@click.argument('dataset', required=True)
@click.pass_context
def frequent_words(context, dataset):
    """
    Runs frequent_words_by_sorting(text, k). The input variables 'text' and 'k' are read from the
    DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = FrequentWordsDataset(dataset, challenge)

    text = data.text
    k = data.k

    if challenge:
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(sorted(bioinfo1.frequent_words_by_sorting(text, k)))}",
                               fg="yellow", bold=True))
    else:
        correct_result = data.result
        result = bioinfo1.frequent_words_by_sorting(text, k)
        # prepare result for output by sorting and joining
        result = ' '.join(sorted(result))
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('reverse-complement')
@click.argument('dataset', required=True)
@click.pass_context
def reverse_complement(context, dataset):
    """
    Runs reverse_complement(pattern). The input variable 'pattern' is read from the DATASET
    argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = ReverseComplementDataset(dataset, challenge)

    pattern = data.pattern

    if challenge:
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.reverse_complement(pattern)}",
                               fg="yellow", bold=True))
    else:
        correct_result = data.result
        result = bioinfo1.reverse_complement(pattern)
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('pattern-matching')
@click.argument('dataset', required=True)
@click.pass_context
def pattern_matching(context, dataset):
    """
    Runs pattern_matching_problem(pattern, genome). The input variables 'pattern' and 'genome' are
    read from the DATASET argument, where DATASET is the text file containing the input data. This
    function is only available in Code Challenge mode.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = PatternMatchingDataset(dataset, challenge)

    pattern = data.pattern
    genome = data.genome

    if challenge:
        result = bioinfo1.pattern_matching_problem(pattern, genome)
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, result))}", fg="yellow", bold=True))
    else:
        # no regular data set available for this function
        click.echo(click.style(f"THIS FUNCTION IS ONLY AVAILABLE IN CODE CHALLENGE MODE!",
                               fg='red', bold=True))


@biocli.command('clump-finding')
@click.argument('dataset', required=True)
@click.pass_context
def clump_finding(context, dataset):
    """
    Runs better_clump_finding(genome, k, l, t). The input variables 'genome', 'k', 'l', and 't' are
    read from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = ClumpFindingDataset(dataset, challenge)

    if context.obj['CHALLENGE']:
        genome = data.genome
        var_k = data.k
        var_l = data.l
        var_t = data.t
        result = bioinfo1.better_clump_finding(genome, var_k, var_l, var_t)
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, result))}", fg="yellow", bold=True))
    else:
        genome = data.genome
        var_k = data.k
        var_l = data.l
        var_t = data.t
        correct_result = data.result
        result = bioinfo1.better_clump_finding(genome, var_k, var_l, var_t)
        # prepare result for output by sorting and joining
        result = ' '.join(sorted(result))
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('computing-frequencies')
@click.argument('dataset', required=True)
@click.pass_context
def computing_frequencies(context, dataset):
    """
    Runs computing_frequencies(text, k). The input variables 'text' and 'k' are read from the
    DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = ComputingFrequenciesDataset(dataset, challenge)

    if context.obj['CHALLENGE']:
        text = data.text
        k = data.k
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, bioinfo1.computing_frequencies(text, k)))}",
                               fg="yellow", bold=True))
    else:
        text = data.text
        k = data.k
        correct_result = data.result
        result = bioinfo1.computing_frequencies(text, k)
        # prepare result for output by sorting and joining
        result = ' '.join(map(str, result))
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('pattern-to-number')
@click.argument('dataset', required=True)
@click.pass_context
def pattern_to_number(context, dataset):
    """
    Runs pattern_to_number(pattern). The input variable 'pattern' is read from the DATASET
    argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = PatternToNumberDataset(dataset, challenge)

    if context.obj['CHALLENGE']:
        pattern = data.pattern
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.pattern_to_number(pattern)}",
                               fg="yellow", bold=True))
    else:
        pattern = data.pattern
        correct_result = data.result
        result = bioinfo1.pattern_to_number(pattern)
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('number-to-pattern')
@click.argument('dataset', required=True)
@click.pass_context
def number_to_pattern(context, dataset):
    """
    Runs number_to_pattern(number, k). The input variables 'number' and 'k' are read from the
    DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = NumberToPatternDataset(dataset, challenge)

    if context.obj['CHALLENGE']:
        number = data.number
        k = data.k
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.number_to_pattern(number, k)}",
                               fg="yellow", bold=True))
    else:
        number = data.number
        k = data.k
        correct_result = data.result
        result = bioinfo1.number_to_pattern(number, k)
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('minimum-skew')
@click.argument('dataset', required=True)
@click.pass_context
def minimum_skew(context, dataset):
    """
    Runs minimum_skew(genome). The input variable 'genome' is read from the DATASET argument, where
    DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = MinimumSkewDataset(dataset, challenge)

    if context.obj['CHALLENGE']:
        genome = data.genome
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, bioinfo1.minimum_skew(genome)))}",
                               fg="yellow", bold=True))
    else:
        genome = data.genome
        correct_result = data.result
        result = bioinfo1.minimum_skew(genome)
        # prepare result for output by sorting and joining
        result = ' '.join(map(str, result))
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('hamming-distance')
@click.argument('dataset', required=True)
@click.pass_context
def hamming_distance(context, dataset):
    """
    Runs hamming_distance(string1, string2). The input variables 'string1' and 'string2' are read
    from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = HammingDistanceDataset(dataset, challenge)

    if context.obj['CHALLENGE']:
        string1 = data.string1
        string2 = data.string2
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.hamming_distance(string1, string2)}",
                               fg="yellow", bold=True))
    else:
        string1 = data.string1
        string2 = data.string2
        correct_result = data.result
        result = bioinfo1.hamming_distance(string1, string2)
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('approx-matching')
@click.argument('dataset', required=True)
@click.pass_context
def approx_matching(context, dataset):
    """
    Runs approx_matching(pattern, text, d). The input variables 'pattern' and 'text' are read from
    the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = ApproxMatchDataset(dataset, challenge)

    pattern = data.pattern
    text = data.text
    d = data.d

    if challenge:
        result = bioinfo1.approx_pattern_match(pattern, text, d)
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, result))}", fg="yellow", bold=True))
    else:
        correct_result = data.result
        result = bioinfo1.approx_pattern_match(pattern, text, d)
        result = ' '.join(map(str, (result)))
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('approx-count')
@click.argument('dataset', required=True)
@click.pass_context
def approx_count(context, dataset):
    """
    Runs approx_pattern_count(pattern, text, d). The input variables 'pattern', 'text' and 'd' are
    read from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = ApproxCountDataset(dataset, challenge)

    pattern = data.pattern
    text = data.text
    d = data.d

    if challenge:
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.approx_pattern_count(pattern, text, d)}",
                               fg="yellow", bold=True))
    else:
        correct_result = data.result
        result = bioinfo1.approx_pattern_count(pattern, text, d)
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('neighbors')
@click.argument('dataset', required=True)
@click.pass_context
def neighbors(context, dataset):
    """
    Runs neighbors(pattern, d). The input variables 'pattern' and 'd' are read from the DATASET
    argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = NeighborsDataset(dataset, challenge)

    pattern = data.pattern
    d = data.d

    if challenge:
        result = bioinfo1.neighbors(pattern, d)
        nl = '\n'
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{nl}{nl.join(sorted(result))}", fg="yellow", bold=True))
    else:
        correct_result = data.result
        result = bioinfo1.neighbors(pattern, d)
        correct_result = ' '.join(sorted(correct_result))
        result = ' '.join(sorted(result))
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('frequent-words-mismatches')
@click.argument('dataset', required=True)
@click.pass_context
def frequent_words_mismatches(context, dataset):
    """
    Runs frequent_words_with_mismatches(pattern, k, d). The input variables 'pattern', 'k' and 'd'
    are read from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = FrequentWordsMismatchesDataset(dataset, challenge)

    text = data.text
    k = data.k
    d = data.d

    if challenge:
        result = bioinfo1.frequent_words_with_mismatches_sorting(text, k, d)
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, result))}", fg="yellow", bold=True))
    else:
        correct_result = data.result
        result = bioinfo1.frequent_words_with_mismatches(text, k, d)
        correct_result = ' '.join(sorted(correct_result))
        result = ' '.join(sorted(result))
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('frequent-words-mismatches-and-reverse-complement')
@click.argument('dataset', required=True)
@click.pass_context
def frequent_words_mismatches_and_reverse_complement(context, dataset):
    """
    Runs frequent_words_with_mismatches_and_reverse_complement(pattern, k, d). The input variables
    'pattern', 'k' and 'd' are read from the DATASET argument, where DATASET is the text file
    containing the input data.
    """
    challenge = context.obj['CHALLENGE']
    click.clear()

    data = FrequentWordsMismatchesDataset(dataset, challenge)

    text = data.text
    k = data.k
    d = data.d

    if challenge:
        result = bioinfo1.frequent_words_with_mismatches_and_reverse_complement(text, k, d)
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, result))}", fg="yellow", bold=True))
    else:
        correct_result = data.result
        result = bioinfo1.frequent_words_with_mismatches_and_reverse_complement(text, k, d)
        correct_result = ' '.join(sorted(correct_result))
        result = ' '.join(sorted(result))
        text_color = result_color(result, correct_result)
        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


if __name__ == '__main__':
    biocli(obj={})
