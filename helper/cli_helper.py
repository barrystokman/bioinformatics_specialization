"""
A collection of helper functions to support the CLI module(s)
"""


def result_color(result, correct_result):
    """
    returns text_color for cli echo, 'green' if correct, else 'red'
    """
    if result == correct_result:
        return 'green'
    return 'red'
