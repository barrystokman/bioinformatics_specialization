import click

from bioinformatics_1 import functions as bioinfo1

DATASET_PATH = '/home/barry/proj/bioinformatics/datasets/'


@click.group()
@click.option('--code-challenge', '-c', is_flag=True,
              help='Use this option if attempting a Code Challenge')
@click.pass_context
def biocli(context, code_challenge):
    """
    Collection of CLI callable bioinformatics functions to make Code Challenges quick and easy to
    perform
    """
    if context.obj is None:
        context.obj = {}

    context.obj['CHALLENGE'] = code_challenge


@biocli.command('pattern-count')
@click.argument('dataset', required=True)
@click.pass_context
def pattern_count(context, dataset):
    """
    Calls the function pattern_count(text, pattern).
    The input variables 'text' and 'pattern' are read from the DATASET argument, where DATASET is
    the text file containing the input data
    """
    click.clear()

    with open(DATASET_PATH + dataset) as f:
        read_data = f.read().splitlines()

    if context.obj['CHALLENGE']:
        # TODO: helper funtion read code challenge data
        text = read_data[0]
        pattern = read_data[1]
        # TODO: helper function print code challenge result
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.pattern_count(text, pattern)}",
                               fg="yellow", bold=True))
    else:
        # TODO: helper funtion read non-code challenge data
        text = read_data[1]
        pattern = read_data[2]
        correct_result = int(read_data[4])
        result = bioinfo1.pattern_count(text, pattern)

        # TODO: helper function print function result
        if result == correct_result:
            text_color = "green"
        else:
            text_color = "red"

        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('frequent-words')
@click.argument('dataset', required=True)
@click.pass_context
def frequent_words(context, dataset):
    """
    Calls the function faster_frequent_words(tekst, k).
    The input variables 'text' and 'k' are read from the DATASET argument, where DATASET is
    the text file containing the input data
    """
    click.clear()

    with open(DATASET_PATH + dataset) as f:
        read_data = f.read().splitlines()

    if context.obj['CHALLENGE']:
        # TODO: helper funtion read code challenge data
        text = read_data[0]
        k = int(read_data[1])
        # TODO: helper function print code challenge result
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(sorted(bioinfo1.faster_frequent_words(text, k)))}",
                               fg="yellow", bold=True))
    else:
        # TODO: helper funtion read non-code challenge data
        text = read_data[1]
        k = int(read_data[2])
        correct_result = read_data[4]
        result = bioinfo1.faster_frequent_words(text, k)

        # prepare result for output by sorting an joining

        result = ' '.join(sorted(result))

        # TODO: helper function print function result
        if result == correct_result:
            text_color = "green"
        else:
            text_color = "red"

        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('reverse-complement')
@click.argument('dataset', required=True)
@click.pass_context
def reverse_complement(context, dataset):
    """
    Calls the function reverse_complement(pattern).
    The input variable 'pattern' is read from the DATASET argument, where DATASET is
    the text file containing the input data
    """
    click.clear()

    with open(DATASET_PATH + dataset) as f:
        read_data = f.read().splitlines()

    if context.obj['CHALLENGE']:
        # TODO: helper funtion read code challenge data
        pattern = read_data[0]
        # TODO: helper function print code challenge result
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.reverse_complement(pattern)}",
                               fg="yellow", bold=True))
    else:
        # TODO: helper funtion read non-code challenge data
        pattern = read_data[1]
        correct_result = read_data[3]
        result = bioinfo1.reverse_complement(pattern)

        # TODO: helper function print function result
        if result == correct_result:
            text_color = "green"
        else:
            text_color = "red"

        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('pattern-matching')
@click.argument('dataset', required=True)
@click.pass_context
def pattern_matching(context, dataset):
    """
    Calls the function pattern_matching_problem(pattern, genome).
    The input variables 'pattern' and 'genome' are read from the DATASET argument, where DATASET is
    the text file containing the input data
    This function is only available in Code Challenge mode
    """
    click.clear()

    with open(DATASET_PATH + dataset) as f:
        read_data = f.read().splitlines()

    if context.obj['CHALLENGE']:
        # TODO: helper funtion read code challenge data
        pattern = read_data[0]
        genome = read_data[1]
        result = bioinfo1.pattern_matching_problem(pattern, genome)
        # TODO: helper function print code challenge result
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
    Calls the function better_clump_finding(genome, k, l, t).
    The input variables 'genome', 'k', 'l', and 't' are read from the DATASET argument, where
    DATASET is the text file containing the input data
    """
    click.clear()

    with open(DATASET_PATH + dataset) as f:
        read_data = f.read().splitlines()

    if context.obj['CHALLENGE']:
        # TODO: helper funtion read code challenge data
        genome = read_data[0]
        k_l_t = [int(i) for i in read_data[1].split(' ')]
        k = k_l_t[0]
        l = k_l_t[1]
        t = k_l_t[2]
        result = bioinfo1.better_clump_finding(genome, k, l, t)
        # TODO: helper function print code challenge result
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, result))}", fg="yellow", bold=True))
    else:
        # TODO: helper funtion read non-code challenge data
        genome = read_data[1]
        k_l_t = [int(i) for i in read_data[2].split(' ')]
        k = k_l_t[0]
        l = k_l_t[1]
        t = k_l_t[2]
        correct_result = read_data[4]
        result = bioinfo1.better_clump_finding(genome, k, l, t)

        # prepare result for output by sorting an joining

        result = ' '.join(sorted(result))

        # TODO: helper function print function result
        if result == correct_result:
            text_color = "green"
        else:
            text_color = "red"

        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('computing-frequencies')
@click.argument('dataset', required=True)
@click.pass_context
def computing_frequencies(context, dataset):
    """
    Calls the function computing_frequencies(tekst, k).
    The input variables 'text' and 'k' are read from the DATASET argument, where DATASET is
    the text file containing the input data
    """
    click.clear()

    with open(DATASET_PATH + dataset) as f:
        read_data = f.read().splitlines()

    if context.obj['CHALLENGE']:
        # TODO: helper funtion read code challenge data
        text = read_data[0]
        k = int(read_data[1])
        # TODO: helper function print code challenge result
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{' '.join(map(str, bioinfo1.computing_frequencies(text, k)))}",
                               fg="yellow", bold=True))
    else:
        # TODO: helper funtion read non-code challenge data
        text = read_data[1]
        k = int(read_data[2])
        correct_result = read_data[4]
        result = bioinfo1.computing_frequencies(text, k)

        # prepare result for output by sorting an joining

        result = ' '.join(map(str, result))

        # TODO: helper function print function result
        if result == correct_result:
            text_color = "green"
        else:
            text_color = "red"

        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('pattern-to-number')
@click.argument('dataset', required=True)
@click.pass_context
def pattern_to_number(context, dataset):
    """
    Calls the function pattern_to_number(pattern).
    The input variable 'pattern' is read from the DATASET argument, where DATASET is
    the text file containing the input data
    """
    click.clear()

    with open(DATASET_PATH + dataset) as f:
        read_data = f.read().splitlines()

    if context.obj['CHALLENGE']:
        # TODO: helper funtion read code challenge data
        pattern = read_data[0]
        # TODO: helper function print code challenge result
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.pattern_to_number(pattern)}",
                               fg="yellow", bold=True))
    else:
        # TODO: helper funtion read non-code challenge data
        pattern = read_data[1]
        correct_result = int(read_data[3])
        result = bioinfo1.pattern_to_number(pattern)

        # TODO: helper function print function result
        if result == correct_result:
            text_color = "green"
        else:
            text_color = "red"

        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


@biocli.command('number-to-pattern')
@click.argument('dataset', required=True)
@click.pass_context
def number_to_pattern(context, dataset):
    """
    Calls the function number_to_pattern(number, k).
    The input variables 'number' and 'k' are read from the DATASET argument, where DATASET is
    the text file containing the input data
    """
    click.clear()

    with open(DATASET_PATH + dataset) as f:
        read_data = f.read().splitlines()

    if context.obj['CHALLENGE']:
        # TODO: helper funtion read code challenge data
        number = int(read_data[0])
        k = int(read_data[1])
        # TODO: helper function print code challenge result
        click.echo(f"The result of the Coding Challenge is:")
        click.echo(click.style(f"{bioinfo1.number_to_pattern(number, k)}",
                               fg="yellow", bold=True))
    else:
        # TODO: helper funtion read non-code challenge data
        number = int(read_data[1])
        k = int(read_data[2])
        correct_result = read_data[4]
        result = bioinfo1.number_to_pattern(number, k)

        # TODO: helper function print function result
        if result == correct_result:
            text_color = "green"
        else:
            text_color = "red"

        click.echo(click.style(f"The result of this function is:"))
        click.echo(click.style(f"{result}", fg=text_color, bold=True))


if __name__ == '__main__':
    biocli(obj={})
