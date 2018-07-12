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
    click.clear()

    data = PatternCountDataset(dataset)

    if context.obj['CHALLENGE']:
        text = data.get_text_challenge()
        pattern = data.get_pattern_challenge()

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.pattern_count(text, pattern)}",
                               fg="yellow", bold=True))
    else:
        text = data.get_text()
        pattern = data.get_pattern()
        correct_result = data.get_expected_result()

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
    click.clear()

    data = FrequentWordsDataset(dataset)

    if context.obj['CHALLENGE']:
        text = data.get_text_challenge()
        k = data.get_k_challenge()

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(sorted(bioinfo1.frequent_words_by_sorting(text, k)))}",
                               fg="yellow", bold=True))
    else:
        text = data.get_text()
        k = data.get_k()
        correct_result = data.get_expected_result()

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
    click.clear()

    data = ReverseComplementDataset(dataset)

    if context.obj['CHALLENGE']:
        pattern = data.get_pattern_challenge()

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.reverse_complement(pattern)}",
                               fg="yellow", bold=True))
    else:
        pattern = data.get_pattern()
        correct_result = data.get_expected_result()

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
    click.clear()

    data = PatternMatchingDataset(dataset)

    if context.obj['CHALLENGE']:
        pattern = data.get_pattern_challenge()
        genome = data.get_genome_challenge()
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
    click.clear()

    data = ClumpFindingDataset(dataset)

    if context.obj['CHALLENGE']:
        genome = data.get_genome_challenge()
        var_k = data.get_k_challenge()
        var_l = data.get_l_challenge()
        var_t = data.get_t_challenge()

        result = bioinfo1.better_clump_finding(genome, var_k, var_l, var_t)

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, result))}", fg="yellow", bold=True))
    else:
        genome = data.get_genome()
        var_k = data.get_k()
        var_l = data.get_l()
        var_t = data.get_t()

        correct_result = data.get_expected_result()

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
    click.clear()

    data = ComputingFrequenciesDataset(dataset)

    if context.obj['CHALLENGE']:
        text = data.get_text_challenge()
        k = data.get_k_challenge()

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, bioinfo1.computing_frequencies(text, k)))}",
                               fg="yellow", bold=True))
    else:
        text = data.get_text()
        k = data.get_k()
        correct_result = data.get_expected_result()

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
    click.clear()

    data = PatternToNumberDataset(dataset)

    if context.obj['CHALLENGE']:
        pattern = data.get_pattern_challenge()

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.pattern_to_number(pattern)}",
                               fg="yellow", bold=True))
    else:
        pattern = data.get_pattern()
        correct_result = data.get_expected_result()

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
    click.clear()

    data = NumberToPatternDataset(dataset)

    if context.obj['CHALLENGE']:
        number = data.get_number_challenge()
        k = data.get_k_challenge()

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.number_to_pattern(number, k)}",
                               fg="yellow", bold=True))
    else:
        number = data.get_number()
        k = data.get_k()
        correct_result = data.get_expected_result()

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
    click.clear()

    data = MinimumSkewDataset(dataset)

    if context.obj['CHALLENGE']:
        genome = data.get_genome_challenge()

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, bioinfo1.minimum_skew(genome)))}",
                               fg="yellow", bold=True))
    else:
        genome = data.get_genome()
        correct_result = data.get_expected_result()

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
    click.clear()

    data = HammingDistanceDataset(dataset)

    if context.obj['CHALLENGE']:
        string1 = data.get_string1_challenge()
        string2 = data.get_string2_challenge()

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.hamming_distance(string1, string2)}",
                               fg="yellow", bold=True))
    else:
        string1 = data.get_string1()
        string2 = data.get_string2()

        correct_result = data.get_expected_result()

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
    click.clear()

    data = ApproxMatchDataset(dataset)

    if context.obj['CHALLENGE']:
        pattern = data.get_pattern_challenge()
        text = data.get_text_challenge()
        d = data.get_d_challenge()

        result = bioinfo1.approx_pattern_match(pattern, text, d)

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, result))}", fg="yellow", bold=True))
    else:
        pattern = data.get_pattern()
        text = data.get_text()
        d = data.get_d()

        correct_result = data.get_expected_result()

        result = bioinfo1.approx_pattern_match(pattern, text, d)

        # prepare result for output by sorting and joining

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
    click.clear()

    data = ApproxCountDataset(dataset)

    if context.obj['CHALLENGE']:
        pattern = data.get_pattern_challenge()
        text = data.get_text_challenge()
        d = data.get_d_challenge()

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.approx_pattern_count(pattern, text, d)}",
                               fg="yellow", bold=True))
    else:
        pattern = data.get_pattern()
        text = data.get_text()
        d = data.get_d()

        correct_result = data.get_expected_result()

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
    click.clear()

    data = NeighborsDataset(dataset)

    if context.obj['CHALLENGE']:
        pattern = data.get_pattern_challenge()
        d = data.get_d_challenge()

        result = bioinfo1.neighbors(pattern, d)

        # backlashes in f-string expressions are not allowed. Assigning '\n' to a variable can
        # circumvent this:
        nl = '\n'

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{nl}{nl.join(sorted(result))}", fg="yellow", bold=True))
    else:
        pattern = data.get_pattern()
        d = data.get_d()

        correct_result = data.get_expected_result()

        result = bioinfo1.neighbors(pattern, d)

        # prepare result for output by sorting and joining

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
    click.clear()

    data = FrequentWordsMismatchesDataset(dataset)

    if context.obj['CHALLENGE']:
        text = data.get_text_challenge()
        k = data.get_k_challenge()
        d = data.get_d_challenge()

        # result = bioinfo1.frequent_words_with_mismatches(text, k, d)
        result = bioinfo1.frequent_words_with_mismatches_sorting(text, k, d)

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, result))}", fg="yellow", bold=True))
    else:
        text = data.get_text()
        k = data.get_k()
        d = data.get_d()

        correct_result = data.get_expected_result()

        result = bioinfo1.frequent_words_with_mismatches(text, k, d)

        # prepare result for output by sorting and joining

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
    click.clear()

    data = FrequentWordsMismatchesDataset(dataset)

    if context.obj['CHALLENGE']:
        text = data.get_text_challenge()
        k = data.get_k_challenge()
        d = data.get_d_challenge()

        result = bioinfo1.frequent_words_with_mismatches_and_reverse_complement(text, k, d)

        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, result))}", fg="yellow", bold=True))
    else:
        text = data.get_text()
        k = data.get_k()
        d = data.get_d()

        correct_result = data.get_expected_result()

        result = bioinfo1.frequent_words_with_mismatches_and_reverse_complement(text, k, d)

        # prepare result for output by sorting and joining

        correct_result = ' '.join(sorted(correct_result))
        result = ' '.join(sorted(result))

        text_color = result_color(result, correct_result)

        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


if __name__ == '__main__':
    biocli(obj={})
