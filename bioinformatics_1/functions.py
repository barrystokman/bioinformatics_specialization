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

skew(genome)
minimum_skew(genome)
hamming_distance(string1, string2)
approx_pattern_match(pattern, text, d)
approx_pattern_count(pattern, text, d)
neighbors(pattern, d)
computing_frequencies_with_mismatches(text, k, d)
frequent_words_with_mismatches(text, k, d)
immediate_neighbors(pattern)
frequent_words_with_mismatches_and_reverse_complement(text, k, d)
frequent_words_with_mismatches_sorting(text, k, d)

DEPRECATED:
count_d(pattern, text, d)
"""


NUCLEOTIDES = ['A', 'C', 'G', 'T']

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


def skew(genome: str) -> list:
    """
    We define skew_i(genome) as the difference between the total number of occurrences of G and the
    total number of occurrences of C in the first i nucleotides of Genome.

    Note that we can compute skew_i+1(genome) from skew_i(genome) according to the nucleotide in
    position i of genome. If this nucleotide is G, then skew_i+1(genome) = skew_i(genome) + 1; if
    this nucleotide is C, then skew_i+1(genome)= skew_i(genome) – 1; otherwise, skew_i+1(enome) =
    skew_i(genome).

    """
    skew_values = [0]

    for i in range(len(genome)):
        if genome[i] == 'G':
            skew_values.append(skew_values[i] + 1)
        elif genome[i] == 'C':
            skew_values.append(skew_values[i] - 1)
        else:
            skew_values.append(skew_values[i])

    return skew_values


def minimum_skew(genome: str) -> list:
    """
    Input: a DNA string genome.
    Output: All integer(s) i minimizing skew_i(genome) among all values of i (from 0 to len(genome))
    """
    _skew = skew(genome)
    minimum_skew = min(_skew)

    min_skew_positions = [i for i in range(len(genome) + 1) if _skew[i] == minimum_skew]

    return min_skew_positions


def hamming_distance(string1: str, string2: str) -> int:
    """
    Computes the Hamming distance between two strings.
    Input: two strings of equal length: string1 and string2
    Output: the Hamming distance between these strings
    """
    distance = 0

    for i in range(len(string1)):
        if string1[i] != string2[i]:
            distance += 1

    return distance


def approx_pattern_match(pattern: str, text: str, d: int) -> list:
    """
    Find all approximate occurences of a pattern in a string
    Input: string pattern and text, integer d
    Output: All starting positions where pattern appears as a substring of text with at most d
    mismatches.
    """
    start_index_list = []

    for i in range(len(text) - len(pattern) + 1):
        if hamming_distance(pattern, text[i:i+len(pattern)]) <= d:
            start_index_list.append(i)

    return start_index_list


def count_d(pattern: str, text: str, d: int) -> int:
    """
    Counts the total number of occurrences of Pattern in Text with at most d mismatches
    """
    return len(approx_pattern_match(pattern, text, d))


def approx_pattern_count(pattern: str, text: str, d: int) -> int:
    """
    Counts the total number of occurrences of Pattern in Text with at most d mismatches
    """
    count = 0

    for i in range(len(text) - len(pattern) + 1):
        sub_text = text[i: i + len(pattern) + 1]
        if hamming_distance(pattern, sub_text) <= d:
            count += 1

    return count


def neighbors(pattern: str, d: int) -> list:
    """
    """
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return ['A', 'C', 'G', 'T']

    neighborhood = []

    suffix_pattern = pattern[1:]
    first_symbol = pattern[0]

    suffix_neighbors = neighbors(suffix_pattern, d)

    for string in suffix_neighbors:
        if hamming_distance(suffix_pattern, string) < d:
            for nucleotide in NUCLEOTIDES:
                neighborhood.append(nucleotide + string)
        else:
            neighborhood.append(first_symbol + string)

    return neighborhood


def computing_frequencies_with_mismatches(text: str, k: int, d: int) -> list:
    """
    generates a frequency array for all k-mers in a string, with at most d mismatches
    input: text: represents a string of nucleotides, length of k-mer k, Hamming distance d
    output: a list representing a frequency array
    """
    frequency = [0 for i in range(NUMBER_OF_SYMBOLS**k)]

    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        neighborhood = neighbors(pattern, d)
        for approximate_pattern in neighborhood:
            j = pattern_to_number(approximate_pattern)
            frequency[j] += 1

    return frequency


def frequent_words_with_mismatches(text: str, k: int, d: int) -> list:
    """
    Works as frequent_words, but uses computing_frequencies, while allowing for mismatches to
    improve running time.
    """
    frequent_patterns = []
    frequency_array = computing_frequencies_with_mismatches(text, k, d)
    max_count = max(frequency_array)

    for i in range(NUMBER_OF_SYMBOLS**k):
        if frequency_array[i] == max_count:
            frequent_patterns.append(number_to_pattern(i, k))

    return frequent_patterns


def immediate_neighbors(pattern: str) -> list:
    """
    generates a list of strings consisting of pattern and all strings 1 Hamming distance away from
    pattern.
    """
    neighborhood = [pattern]

    for i in range(len(pattern)):
        symbol = pattern[i]
        for nucleotide in NUCLEOTIDES:
            if nucleotide != symbol:
                neighbor = pattern[:i] + nucleotide + pattern[i + 1:]
                neighborhood.append(neighbor)

    return neighborhood


def frequent_words_with_mismatches_and_reverse_complement(text: str, k: int, d: int) -> list:
    """
    Find the most frequent k-mers (with mismatches and reverse complements) in a string.
    Input: A DNA string Text as well as integers k and d.
    Output: All k-mers Pattern maximizing the sum Countd(Text, Pattern) + Countd(Text, Patternrc)
    over all possible k-mers.
    """
    frequent_patterns = []

    frequency_array_text = computing_frequencies_with_mismatches(text, k, d)
    frequency_array_rc = computing_frequencies_with_mismatches(reverse_complement(text), k, d)

    frequency_array = [frequency_array_text[i] +
                       frequency_array_rc[i] for i in range(len(frequency_array_text))]

    max_count = max(frequency_array)

    for i in range(NUMBER_OF_SYMBOLS**k):
        if frequency_array[i] == max_count:
            frequent_patterns.append(number_to_pattern(i, k))

    return frequent_patterns


def frequent_words_with_mismatches_sorting(text: str, k: int, d: int) -> list:
    """
    More efficient version of frequent_words_with_mismatches using sorting.
    """
    frequent_patterns = []
    index = []

    neighborhoods = [neighbors(text[i:i+k], d) for i in range(len(text) - k + 1)]

    neighborhood_array = [neighbor for neighborhood in neighborhoods for neighbor in
                          neighborhood]

    for i in range(len(neighborhood_array)):
        pattern = neighborhood_array[i]
        index.append(pattern_to_number(pattern))

    count = [1 for index in range(len(neighborhood_array))]
    sorted_index = sorted(index)

    for i in range(len(neighborhood_array) - 1):
        if sorted_index[i] == sorted_index[i + 1]:
            count[i + 1] = count[i] + 1

    max_count = max(count)

    for i in range(len(neighborhood_array)):
        if count[i] == max_count:
            frequent_patterns.append(number_to_pattern(sorted_index[i], k))

    return sorted(list(set(frequent_patterns)))


def main():
    pass


if __name__ == '__main__':
    main()
