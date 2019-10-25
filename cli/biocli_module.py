import click
from inspect import getmembers, isfunction

from bioinformatics_1 import functions as bioinfo1
from apps.cli import cli_output, get_correct_result, determine_padding
from apps.plot import SkewPlot, CountBarPlot
from apps.ori_finder import OriFinder
from apps.constants import WINDOW_SIZE
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
    import ipdb; ipdb.set_trace()
    challenge = context.obj['CHALLENGE']

    data = dsr.PatternCountDataset(dataset, challenge)

    text = data.text
    pattern = data.pattern
    correct_result = get_correct_result(data, challenge)

    args = text, pattern
    func = bioinfo1.pattern_count

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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

    cli_output(challenge, correct_result, sort_result, listing, func, args)


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


@biocli.command('motif-enumeration')
@click.argument('dataset', required=True)
@click.pass_context
def motif_enumeration(context, dataset, sort_result=True, listing=False):
    """
    Runs motif_enumeration(dna, k, d). The input variables 'dna', 'k' and 'd' are read from the
    DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.MotifEnumerationDataset(dataset, challenge)

    dna = data.dna
    k = data.k
    d = data.d

    correct_result = get_correct_result(data, challenge)

    args = dna, k, d
    func = bioinfo1.motif_enumeration

    cli_output(challenge, correct_result, sort_result, listing, func, args)


@biocli.command('distance-between-pattern-and-string')
@click.argument('dataset', required=True)
@click.pass_context
def distance_between_pattern_and_string(context, dataset, sort_result=True, listing=False):
    """
    Runs distance_between_pattern_and_string(pattern, dna). The input variables 'pattern' and 'dna'
    are read from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.DistanceBetweenPatternAndString(dataset, challenge)

    pattern = data.pattern
    dna = data.dna

    correct_result = get_correct_result(data, challenge)

    args = pattern, dna
    func = bioinfo1.distance_between_pattern_and_strings

    cli_output(challenge, correct_result, sort_result, listing, func, args)


@biocli.command('median-string')
@click.argument('dataset', required=True)
@click.pass_context
def median_string(context, dataset, sort_result=True, listing=False):
    """
    Runs median_string(dna, k). The input variables 'dna' and 'k'
    are read from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.MedianString(dataset, challenge)
    dna = data.dna
    k = data.k

    correct_result = get_correct_result(data, challenge)

    args = dna, k
    func = bioinfo1.median_string

    cli_output(challenge, correct_result, sort_result, listing, func, args)


@biocli.command('profile-most-probable-kmer')
@click.argument('dataset', required=True)
@click.pass_context
def profile_most_probable_kmer(context, dataset, sort_result=False, listing=False):
    """
    Runs profile_most_probable_kmer(text, k, profile). The input variables 'text', 'k' and 'profile'
    are read from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.ProfileMostProbableKmer(dataset, challenge)
    text = data.text
    k = data.k
    profile = data.profile

    correct_result = get_correct_result(data, challenge)

    args = text, k, profile
    func = bioinfo1.profile_most_probable_kmer

    cli_output(challenge, correct_result, sort_result, listing, func, args)


@biocli.command('greedy-motif-search')
@click.argument('dataset', required=True)
@click.pass_context
def greedy_motif_search(context, dataset, sort_result=False, listing=False):
    """
    Runs greedy_motif_search(dna, k, t). The input variables 'dna', 'k' and 't'
    are read from the DATASET argument, where DATASET is the text file containing the input data.
    """
    challenge = context.obj['CHALLENGE']

    data = dsr.GreedyMotifSearch(dataset, challenge)
    dna = data.dna
    k = data.k
    t = data.t

    correct_result = get_correct_result(data, challenge)

    args = dna, k, t
    func = bioinfo1.greedy_motif_search
    cli_output(challenge, correct_result, sort_result, listing, func, args)


@biocli.command('randomized-motif-search')
@click.argument('dataset', required=True)
@click.pass_context
def randomized_motif_search(context, dataset, sort_result=False, listing=True):
    """
    Loops over randomized_motif_search(dna, k, t) 1000 times. The input variables 'dna', 'k' and
    't' are read from the DATASET argument, where DATASET is the text file containing the input
    data.
    """
    challenge = context.obj['CHALLENGE']
    data = dsr.RandomizedMotifSearch(dataset, challenge)
    dna = data.dna
    k = data.k
    t = data.t

    correct_result = get_correct_result(data, challenge)

    args = dna, k, t
    func = bioinfo1.loop_randomized_motif_search
    cli_output(challenge, correct_result, sort_result, listing, func, args)


@biocli.command('subtle-motif-problem')
@click.argument('dataset', required=True)
@click.option('-n', '--number-of-iterations', required=False, default=10000, type=int, help='')
@click.pass_context
def subtle_motif_problem(context, dataset, sort_result=False, listing=True):
    """
    Tries to solve the subtle motif problem as stated in Week 4, Chapter 1.1, step 9:
    1) Run randomized_motif_search(dna, k, t) n (standard = 10 000) times by running loop_randomized_motif_search().
    2) Determine motif score (by row)
    3) Determine the concensus string
    """
    import ipdb;ipdb.set_trace() 
    challenge = True
    data = dsr.RandomizedMotifSearch(dataset, challenge)
    dna = data.dna
    k = data.k
    t = data.t

    args = dna, k, t, n
    # motifs = bioinfo1.loop_randomized_motif_search(*args)
    motifs = bioinfo1.loop_randomized_motif_search(dna, k, t, n)
    click.echo(click.style(f"Result of running randomized motif search:", fg="blue", bold=True))
    for motif in motifs:
        click.echo(click.style(f"{motif}", fg="yellow", bold=True))

    score = bioinfo1.motifs_score_by_rows(motifs)
    click.echo(click.style(f"Motif Score:", fg="blue", bold=True))
    click.echo(click.style(f"{score}", fg="yellow", bold=True))

    concensus = bioinfo1.motifs_concensus(motifs)
    click.echo(click.style(f"Concensus:", fg="blue", bold=True))
    click.echo(click.style(f"{concensus}", fg="yellow", bold=True))


@biocli.command('skew-plot')
@click.option('-m', '--minimum-skew', required=False, is_flag=True, help='adds minimum skew indicators to plot')
@click.argument('genome', required=True)
def skew_plot(minimum_skew, genome):
    """
    Generates a skew plot of the provided genome
    """
    skew_genome = dsr.Genome(genome)
    skewplot_obj = SkewPlot(skew_genome)
    skewplot_obj.generate_plot()
    if minimum_skew:
        skewplot_obj.show_minimum_skew()
    skewplot_obj.show_plot()


@biocli.command('find-ori')
@click.argument('genome', required=True)
def find_ori(genome):
    """
    Finds the origin of replication in a given microbial genome
    1. Read genome
    2. Plot skew diagram with minimum skew indicators
    3. Get ori candidate based on minimum skew
    4. Find most DnaA box candidates
    5. Plot bar graphs for DnaA box candidates
    """
    # 1 read genome
    ori_genome = dsr.Genome(genome)
    ori_obj = OriFinder(ori_genome)

    # 2 plot skew diagram
    plot_obj = SkewPlot(ori_genome)
    plot_obj.generate_plot()
    plot_obj.show_minimum_skew()
    plot_obj.show_plot()

    # 3 get ori candidate
    ori_candidate = ori_obj.ori_candidate
    print(f"ori candidate: {ori_candidate}")

    # 4.1 find starting point of first window
    first_position = ori_candidate - 2000

    # 4.2 find starting point of last window
    last_position = ori_candidate - 1999
    # last_position = ori_candidate + 2000 - WINDOW_SIZE

    print(f"{first_position}, {last_position}")

    # 4.3 find most frequent 9-mers with 1 mismatch
    for mismatch in [2]:
        max_kmer_freq = -999
        most_frequent_kmers_and_rc = []
        ori_position = [-999]
        for position in range(first_position, last_position):
            print('{:.2%}'.format((position - first_position) / (last_position - first_position)))
            result = ori_obj.find_frequent_kmers(start=position, mismatch=mismatch)
            ori_frequent_kmers = result[0]
            frequent_kmers_count = result[1]

            if frequent_kmers_count > max_kmer_freq:
                max_kmer_freq = frequent_kmers_count
                most_frequent_kmers_and_rc = [ori_frequent_kmers]
                ori_position = [position]
            elif frequent_kmers_count == max_kmer_freq and ori_frequent_kmers not in most_frequent_kmers_and_rc:
                most_frequent_kmers_and_rc.append(ori_frequent_kmers)
                if position - 500 > ori_position[-1]:
                    ori_position = [position]

        print(f"highest kmer frequency with {mismatch} mismatches: {most_frequent_kmers_and_rc}, occuring \
              {max_kmer_freq} times,on position {ori_position}")

    import ipdb; ipdb.set_trace()
    # most_frequent_kmers_and_rc = [['TCATGATCA', 'TGATCATGA'], ['ATGATCATG', 'CATGATCAT', 'TCATGATCA', 'TGATCATGA'], ['ATGATCATG', 'CATGATCAT']]
    most_frequent_kmers_and_rc = [['AATGATCAT', 'ATGATCATT'], ['AATGATCAT', 'ATCATGATC', 'ATGATCATT', 'GATCATGAT'], ['ATCATGATC', 'GATCATGAT']]
    most_frequent_kmers = list()

    for i in most_frequent_kmers_and_rc:
        for j in i:
            if j not in most_frequent_kmers and bioinfo1.reverse_complement(j) not in most_frequent_kmers:
                most_frequent_kmers.append(j)

    print(f"most frequent kmers with {mismatch} mismatch(es): {most_frequent_kmers}")

    # find count of all most frequent kmers per position in a size 500 window in whole genome!:

    candidates = ['TCATGATCA', 'ATGATCATG']
    # candidate = 'ATGATCATG'
    candidate = 'TCATGATCA'

    graph_count = list()
    for position in range(len(ori_genome.genome) - WINDOW_SIZE):
        number_of_kmers = ori_obj.count_kmers_in_window(candidate, start=position, mismatch=mismatch)
        graph_count.append(number_of_kmers)
        print('{:.2%}'.format(position / (len(ori_genome.genome) - WINDOW_SIZE)))


    import ipdb; ipdb.set_trace()
    # plot  diagram
    plot_obj = CountBarPlot(ori_genome, graph_count)
    plot_obj.generate_plot()
    plot_obj.show_plot()

    import csv
    with open('TCATGATCA-count.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(graph_count)


@biocli.command('quiz')
def quiz():
    """
    Easy quiz handling
    """
    func_list = getmembers(bioinfo1, isfunction)
    func_dict = dict(zip(range(len(func_list)), func_list))

    for key, value in func_dict.items():
        click.echo(f"{key}: {value[0]}")

    while True:
        try:
            choice = click.prompt("Please select a function from the list")
            function_name, function = func_dict[int(choice)]
        except (ValueError, KeyError):
            click.echo(f"Invalid input, try again")
            continue
        break

    click.echo(f"You have selected: {function_name}")
    arg_count = function.__code__.co_argcount
    args = function.__code__.co_varnames[:arg_count]
    annotations = function.__annotations__

    arg_values = {}

    for arg in args:
        arg_value = click.prompt(f"Enter a value for {arg}")
        if annotations[arg] == int:
            arg_values[arg] = int(arg_value)
        elif annotations[arg] == list:
            arg_values[arg] = arg_value.split(' ')
        elif annotations[arg] == dict:
            try:
                if function_name in ['profile_most_probable_kmer']:
                    split_arg_value = arg_value.split('|')
                    probabilities = [[float(p) for p in row.split(' ')] for row in split_arg_value]
                    arg_values[arg] = dict(zip(bioinfo1.NUCLEOTIDES, probabilities))
            except:
                click.echo(f"I'm not sure what to do with this dictionary")
        else:
            arg_values[arg] = arg_value

    result = function(**arg_values)
    if annotations['return'] == dict:
        padding = determine_padding(result)
        click.echo(f"Result:")
        for key, value in result.items():
            print(f"{key}:", end=" ")
            for frequency in value:
                frequency = str(frequency).rjust(padding, ' ')
                print(f"{frequency}", end=" ")
            print()
    if annotations['return'] == list:
        print(f"Result:")
        for item in result:
            print(f"{item}", sep=" ", end="\n")
    else:
        click.echo(f"Result: {result}")


if __name__ == '__main__':
    biocli(obj={})
