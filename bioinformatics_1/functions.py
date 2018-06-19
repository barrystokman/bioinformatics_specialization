"""
Collection of funtions developed for the course Bioinformatics I: Finding Hidden Messages in DNA.
Note that these functions return python objects. The code graded problems (code_graded_problems.py)
are formatted for stdout output.
"""

"""
Week 1:

pattern_count(text, pattern)
frequent_words_by_sorting(text, k)
reverse_complement(pattern)
pattern_matching_problem(pattern, genome)
better_clump_finding(genome, k, l, t)
computing_frequencies(text, k)
pattern_to_number(pattern)
number_to_pattern(number, k)

DEPRECATED:
frequent_words(text, k)
faster_frequent_words(text, k)
clump_finding_problem(genome, k, l, t)
clump_finding(genome, k, l, t)
number_to_pattern_non_recursive(number, k)
pattern_to_number_non_recursive(pattern)

Week 2:
"""


BASE_COMPLEMENT = {'A': 'T',
                   'C': 'G',
                   'G': 'C',
                   'T': 'A'
                   }

SYMBOL_TO_NUMBER = {'A': 0,
                    'C': 1,
                    'G': 2,
                    'T': 3
                    }

NUMBER_TO_SYMBOL = {0: 'A',
                    1: 'C',
                    2: 'G',
                    3: 'T'
                    }

NUMBER_OF_SYMBOLS = 4


def pattern_count(text: str, pattern: str) -> int:
    """
    count the number of times a string appears in a text
    input: text, pattern
    output: count
    """
    count = 0

    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count


def frequent_words(text: str, k: int) -> list:
    """
    find the most frequent k-mers in a string 'text'
    """
    frequent_patterns = []

    count = [0 for index in range(len(text) - k + 1)]

    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        count[i] = pattern_count(text, pattern)
    max_count = max(count)

    for i in range(len(text) - k + 1):
        if count[i] == max_count:
            frequent_patterns.append(text[i:i+k])

    return list(set(frequent_patterns))


def frequent_words_by_sorting(text: str, k: int) -> list:
    """
    Works as frequent_words, but uses a sorted index array to improve running time
    """
    frequent_patterns = []
    index = []
    count = [1 for index in range(len(text) - k + 1)]

    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        index.append(pattern_to_number(pattern))

    index.sort()

    for i in range(1, len(text) - k + 1):
        if index[i] == index[i - 1]:
            count[i] = count[i - 1] + 1

    max_count = max(count)

    for i in range(len(text) - k + 1):
        if count[i] == max_count:
            frequent_patterns.append(number_to_pattern(index[i], k))

    return frequent_patterns


def faster_frequent_words(text: str, k: int) -> list:
    """
    Works as frequent_words, but uses computing_frequencies to improve running time
    """
    frequent_patterns = []
    frequency_array = computing_frequencies(text, k)
    max_count = max(frequency_array)

    for i in range(NUMBER_OF_SYMBOLS**k):
        if frequency_array[i] == max_count:
            frequent_patterns.append(number_to_pattern(i, k))

    return frequent_patterns


def reverse_complement(pattern: str) -> str:
    """
    finds the reverse complement of a given string
    input: pattern
    output: reverse complement
    """
    reverse_complement = ""

    for i in range(len(pattern)):
        reverse_complement = BASE_COMPLEMENT[pattern[i]] + reverse_complement
    return reverse_complement


def pattern_matching_problem(pattern: str, genome: str) -> list:
    """
    finds the starting positions of where a pattern appears in a genome
    input: two strings, pattern and genome.
    output: a list of integers specifying all starting positions where pattern appears as a
    substring of genome.
    """
    start_index_list = []

    for i in range(len(genome) - len(pattern) + 1):
        if pattern == genome[i:i+len(pattern)]:
            start_index_list.append(i)

    return start_index_list


def clump_finding_problem(genome: str, k: int, l: int, t: int) -> list:
    """
    finds clumps of k-mers of length k in a region of the genome of size l, where the k-mer appears
    at least t times.
    input: genome: text representing genome
           k     : integer representing the length of the k-mer
           l     : integer representing the size of the region of the genome that is checked for
                   clumps
           t     : integer representing the minimum number of times the k-mer should appear in the
                   interval
    output: a list of all distinct k-mers forming (l, t)-clumps in a region of the genome
    """
    clumps = []

    for i in range(len(genome) - l + 1):
        for j in range(l - k + 1):
            kmer = genome[i+j:i+j + k]
            count = pattern_count(genome[i:i+l], kmer)
            if count >= t and kmer not in clumps:
                clumps.append(kmer)

    return clumps


def clump_finding(genome: str, k: int, l: int, t: int) -> list:
    """
    finds clumps of k-mers of length k in a region of the genome of size l, where the k-mer appears
    at least t times.
    input: genome: text representing genome
           k     : integer representing the length of the k-mer
           l     : integer representing the size of the region of the genome that is checked for
                   clumps
           t     : integer representing the minimum number of times the k-mer should appear in the
                   interval
    output: a list of all distinct k-mers forming (l, t)-clumps in a region of the genome
    """
    clump_patterns = []
    clumps = [0 for _ in range(NUMBER_OF_SYMBOLS**k)]

    for i in range(len(genome) - l + 1):
        text = genome[i:i + l]
        frequency_array = computing_frequencies(text, k)

        for index in range(NUMBER_OF_SYMBOLS**k):
            if frequency_array[index] >= t:
                clumps[index] = 1

    for i in range(NUMBER_OF_SYMBOLS**k):
        if clumps[i] == 1:
            clump_patterns.append(number_to_pattern(i, k))

    return clump_patterns


def better_clump_finding(genome: str, k: int, l: int, t: int) -> list:
    """
    finds clumps of k-mers of length k in a region of the genome of size l, where the k-mer appears
    at least t times.
    input: genome: text representing genome
           k     : integer representing the length of the k-mer
           l     : integer representing the size of the region of the genome that is checked for
                   clumps
           t     : integer representing the minimum number of times the k-mer should appear in the
                   interval
    output: a list of all distinct k-mers forming (l, t)-clumps in a region of the genome
    """
    clump_patterns = []
    clumps = [0 for _ in range(NUMBER_OF_SYMBOLS**k)]
    text = genome[0:l]
    frequency_array = computing_frequencies(text, k)

    for i in range(NUMBER_OF_SYMBOLS**k):
        if frequency_array[i] >= t:
            clumps[i] = 1

    for i in range(1, len(genome) - l + 1):
        first_pattern = genome[i - 1:i - 1 + k]
        index_first_pattern = pattern_to_number(first_pattern)
        frequency_array[index_first_pattern] -= 1

        last_pattern = genome[i + l - k:i + l]
        index_last_pattern = pattern_to_number(last_pattern)
        frequency_array[index_last_pattern] += 1

        if frequency_array[index_last_pattern] >= t:
            clumps[index_last_pattern] = 1

    for i in range(NUMBER_OF_SYMBOLS**k):
        if clumps[i] == 1:
            clump_patterns.append(number_to_pattern(i, k))

    return clump_patterns


def computing_frequencies(text: str, k: int) -> list:
    """
    generates a frequency array for all k-mers in a string
    input: text: represents a string of nucleotides
    output: a list representing a frequency array
    """
    frequency = [0 for i in range(NUMBER_OF_SYMBOLS**k)]

    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        j = pattern_to_number(pattern)
        frequency[j] += 1

    return frequency


def number_to_pattern_non_recursive(number: int, k: int) -> str:
    """
    Converts a number into a k-mer
    """
    pattern = ''

    for i in reversed(range(k)):
        pattern += NUMBER_TO_SYMBOL[number // NUMBER_OF_SYMBOLS**i]
        number = number % NUMBER_OF_SYMBOLS**i

    return pattern


def number_to_pattern(number: int, k: int) -> str:
    """
    Converts a number into a k-mer
    """
    if k == 1:
        return NUMBER_TO_SYMBOL[number]
    else:
        quotient = number // NUMBER_OF_SYMBOLS
        remainder = number % NUMBER_OF_SYMBOLS
        return number_to_pattern(quotient, k - 1) + NUMBER_TO_SYMBOL[remainder]


def pattern_to_number_non_recursive(pattern: str) -> int:
    """
    Converts a pattern of nucleotides to a number representing its lexicological index
    """
    number = 0

    for i, char in enumerate(pattern):
        number += SYMBOL_TO_NUMBER[char] * NUMBER_OF_SYMBOLS**(len(pattern) - i - 1)

    return number


def pattern_to_number(pattern: str) -> int:
    """
    Converts a pattern of nucleotides to a number representing its lexicological index
    This functions works recursively
    """
    if len(pattern) == 0:
        return 0
    else:
        symbol = pattern[-1:]
        prefix = pattern[:-1]
        return 4 * pattern_to_number(prefix) + SYMBOL_TO_NUMBER[symbol]


def main():
    pass


if __name__ == '__main__':
    main()
