"""
A collection of helper functions to support the CLI module(s)
"""
import click
import pyperclip

import apps.decorators


def result_color(result, correct_result):
    """
    returns text_color for cli echo, 'green' if correct, else 'red'
    """
    if result == correct_result:
        return 'green'
    return 'red'

def result_to_clipboard(result):
    return pyperclip.copy(str(result))

def format_result(result, sort_result, listing):
    if listing:
        separator = '\n'
    else:
        separator = ' '

    if type(result) is list:
        if sort_result:
            result = separator.join(sorted(map(str, result)))
        else:
            result = separator.join(map(str, result))

    return result

def cli_output(challenge, correct_result, sort_result, listing, func, args):
    function_output = apps.decorators.cli_output_decorator(challenge, correct_result, sort_result,
                                                          listing)(func)
    return function_output(*args)

def get_correct_result(data, challenge):
    if not challenge:
        return data.result

def challenge_output(func, result):
    click.echo(click.style(f"The result of the Coding Challenge for function {func.__name__!r} is:", fg="white", bold=True))
    click.echo(click.style(f"{result}", fg="yellow", bold=True))

def no_challenge_output(func, result, correct_result):
    text_color = result_color(result, correct_result)
    click.echo(click.style(f"The result of function {func.__name__!r} is:", fg="white", bold=True))
    click.echo(click.style(f"{result}", fg=text_color, bold=True))
