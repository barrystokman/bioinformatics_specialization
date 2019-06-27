import functools

import apps.cli


def cli_output_decorator(challenge, correct_result, sort_result, listing):
    def cli_output(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            unformatted_result = func(*args, **kwargs)
            result = apps.cli.format_result(unformatted_result, sort_result, listing)
            apps.cli.result_to_clipboard(result)
            if challenge:
                apps.cli.challenge_output(func, result)
            else:
                apps.cli.no_challenge_output(func, result, correct_result)
            return result
        return wrapper
    return cli_output
