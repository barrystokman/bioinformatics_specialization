import math
import random

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

Week 3:

motif_enumeration(dna, k, d)
distance_between_pattern_and_strings(pattern, dna)
median_string(dna, k)
motifs_score_by_columns(motifs)
motifs_score_by_rows(motifs)
motifs_count(motifs)
motifs_profile(motifs)
motifs_profile_laplace(motifs)
motifs_concensus(motifs)
motif_entropy(motifs)
profile_most_probable_kmer(text, k, profile)
greedy_motif_search(dna, t)

Week 4:
randomized_motif_search(dna, k, t)
gibbs_random(probability_distribution)

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


def approx_pattern_and_reverse_complement_count(pattern: str, text: str, d: int) -> int:
    """
    Finds the frequency of a pattern and its reverse complement in a text, and all patterns within
    a hamming distance d
    """
    pattern_frequency = approx_pattern_count(pattern, text, d)
    rc_pattern_frequency = approx_pattern_count(reverse_complement(pattern), text, d)

    return pattern_frequency + rc_pattern_frequency


def motif_enumeration(dna: list, k: int, d: int) -> set:
    """
    Brute force or exhaustive search algorithm to find an implanted motif.

    Given a collection of strings Dna and an integer d, a k-mer is a (k,d)-motif if it appears in
    every string from Dna with at most d mismatches

    Input: Integers k and d, followed by a collection of strings Dna.
    Output: All (k, d)-motifs in Dna.
    """
    patterns = set()

    first_string = dna[0]

    for i in range(len(first_string) - k + 1):
        pattern = first_string[i:i+k]
        for neighbor in neighbors(pattern, d):
            add_neighbor = []
            for string in dna:
                add_neighbor.append(any(_ in string for _ in neighbors(neighbor, d)))
            if all(add_neighbor):
                patterns.add(neighbor)

    return list(patterns)


def distance_between_pattern_and_strings(pattern: str, dna: list) -> int:
    """
    Input: A string pattern followed by a collection of strings dna.
    Output: d(pattern, dna).
    """
    k = len(pattern)
    distance = 0

    for string in dna:
        d = math.inf
        for i in range(len(string) - k + 1):
            pattern_in_string = string[i:i+k]
            if hamming_distance(pattern, pattern_in_string) < d:
                d = hamming_distance(pattern, pattern_in_string)
        distance += d

    return distance


def median_string(dna: list, k: int) -> str:
    """
    Input: A collection of strings Dna and an integer k.
    Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers.
    """
    distance = math.inf

    for i in range(NUMBER_OF_SYMBOLS ** k - 1):
        pattern = number_to_pattern(i, k)
        current_distance = distance_between_pattern_and_strings(pattern, dna)
        if distance > current_distance:
            distance = current_distance
            median = pattern

    return median


def motifs_score_by_columns(motifs: list) -> int:
    """
    Given a motif matrix motifs, motifs_score(motifs) determines the sum of the number of unpopular
    letters in each column of the matrix
    """
    count_matrix = motifs_count(motifs)
    score = 0
    columns = len(motifs[0])

    for i in range(columns):
        all_scores = [value[i] for value in count_matrix.values()]
        score += sum(all_scores) - max(all_scores)

    return score


def motifs_score_by_rows(motifs: list) -> int:
    """

    """
    concensus = motifs_concensus(motifs)
    score = 0

    for string in motifs:
        score += hamming_distance(concensus, string)

    return score


def motifs_count(motifs: list) -> dict:
    """

    """
    count_matrix = {nucleotide: [] for nucleotide in NUCLEOTIDES}
    columns = len(motifs[0])

    for i in range(columns):
        column = [string[i] for string in motifs]
        for nucleotide in NUCLEOTIDES:
            count_matrix[nucleotide].append(column.count(nucleotide))

    return count_matrix


def motifs_profile(motifs: list) -> dict:
    """

    """
    count_matrix = motifs_count(motifs)
    sum_scores = len(motifs)

    profile = {k: [score / sum_scores for score in v] for k, v in count_matrix.items()}

    return profile


def motifs_profile_laplace(motifs: list) -> dict:
    """

    """
    count_matrix = motifs_count(motifs)
    pseudo_count_matrix = {key: [v + 1 for v in values] for key, values in count_matrix.items()}
    sum_scores = sum([values[0] for values in list(pseudo_count_matrix.values())])

    profile = {k: [score / sum_scores for score in v] for k, v in pseudo_count_matrix.items()}

    return profile


def motifs_concensus(motifs: list) -> str:
    """

    """
    profile = motifs_profile(motifs)
    columns = len(motifs[0])
    nucleotides = list(profile.keys())
    concensus = str()

    for i in range(columns):
        profile_values = [value[i] for value in profile.values()]
        max_profile_value = max(profile_values)
        max_profile_index = profile_values.index(max_profile_value)
        concensus_nucleotide = nucleotides[max_profile_index]
        concensus += concensus_nucleotide

    return concensus


def motif_entropy(motifs: list) -> float:
    """

    """
    profile = motifs_profile(motifs)
    columns = len(motifs[0])
    entropy = 0

    for i in range(columns):
        profile_values = [value[i] for value in profile.values()]
        entropy += sum([value * math.log2(value) if value != 0 else 0 for value in profile_values])

    return entropy


def profile_most_probable_kmer(text: str, k: int, profile: dict) -> str:
    """
    Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.

    Input: A string Text, an integer k, and a 4 × k matrix Profile.
    Output: A Profile-most probable k-mer in Text.
    """
    text_length = len(text)
    highest_probability = -math.inf
    most_probable_kmer = str()

    for i in range(text_length - k + 1):
        pattern = text[i:i+k]
        p = 1
        for j in range(len(pattern)):
            p *= profile[pattern[j]][j]
        if p > highest_probability:
            highest_probability = p
            most_probable_kmer = pattern

    return most_probable_kmer


def greedy_motif_search(dna: list, k: int, t: int) -> list:
    """
    Input: Integers k and t, followed by a collection of strings dna.
    Output: A collection of strings best_motifs resulting from applying greedy_motif_search(dna, k,
    t). If at any step you find more than one profile-most probable k-mer in a given string, use
    the one occurring first.
    """
    first_dna_string = dna[0]
    best_motifs = [string[:k] for string in dna]

    for i in range(len(first_dna_string) - k + 1):
        motifs = [first_dna_string[i:i+k]]
        for j in range(1, t):
            profile = motifs_profile_laplace(motifs)
            best_motif = profile_most_probable_kmer(dna[j], k, profile)
            motifs.append(best_motif)
        score_motifs = motifs_score_by_rows(motifs)
        score_best_motifs = motifs_score_by_rows(best_motifs)
        if score_motifs < score_best_motifs:
            best_motifs = motifs

    return best_motifs


def randomized_motif_search(dna: list, k: int, t: int) -> list:
    """
    """
    motifs = []
    for i, string in enumerate(dna):
        random_index = random.randrange(len(string) - k)
        motifs.append(string[random_index:random_index+k])
    best_motifs = motifs.copy()
    while True:
        profile = motifs_profile_laplace(motifs)
        motifs = [profile_most_probable_kmer(string, k, profile) for string in dna]
        if motifs_score_by_rows(motifs) < motifs_score_by_rows(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


def loop_randomized_motif_search(dna: list, k: int, t: int, n: int=1000) -> list:
    """

    """
    current_motifs = randomized_motif_search(dna, k, t)

    for i in range(n):
        best_motifs = randomized_motif_search(dna, k, t)
        if motifs_score_by_rows(best_motifs) < motifs_score_by_rows(current_motifs):
            current_motifs = best_motifs
    return current_motifs


def gibbs_random(probability_distribution: list) -> int:
    """
    Given a probability distribution (p1, …, pn), this random number generator, denoted
    Random(p1, …, pn), models an n-sided biased die and returns integer i with probability pi. If
    the pi sum to some C > 0 instead, then Random(p1, …, pn) is defined as Random(p1/C, …, pn/C),
    where (p1/C, …, pn/C) is a probability distribution.
    """
    c = sum(probability_distribution)
    normalized_distribution = [p/c for p in probability_distribution]

    roll = random.random()

    for i, pc in enumerate(normalized_distribution):
        if roll < sum(normalized_distribution[:i+1]):
            return i


def gibbs_sampler(dna: list, k: int, t: int, n: int) -> list:
    """

    """

    motifs = []
    for i, string in enumerate(dna):
        random_index = random.randrange(len(string) - k)
        motifs.append(string[random_index:random_index+k])

    best_motifs = motifs.copy()

    for j in range(n):
        remove_index = random.randrange(t)
        motifs_one_removed = motifs[:remove_index] + motifs[remove_index+1:]
        profile = motifs_profile_laplace(motifs_one_removed)

        probabilities = []
        for kmer_index in range(len(dna[remove_index]) - k + 1):
            probability_index = 1
            kmer = dna[remove_index][kmer_index:kmer_index+k]
            for char_index, char in enumerate(kmer):
                probability_index *= profile[char][char_index]
            probabilities.append(probability_index)

        gibbs_random_index = gibbs_random(probabilities)

        insert_motif = dna[remove_index][gibbs_random_index:gibbs_random_index+k]
        motifs[remove_index] = insert_motif

        if motifs_score_by_rows(motifs) < motifs_score_by_rows(best_motifs):
            best_motifs = motifs.copy()

    return best_motifs


def gibbs_sampler_loop(dna: list, k: int, t: int, n: int) -> list:
    """

    """
    NUMBER_OF_STARTS = 20
    best_motifs = []
    best_motif_score = math.inf

    for i in range(NUMBER_OF_STARTS):
        current_motifs = gibbs_sampler(dna, k, t, n)
        if motifs_score_by_rows(current_motifs) < best_motif_score:
            best_motifs = current_motifs
            best_motif_score = motifs_score_by_rows(best_motifs)

    return best_motifs


if __name__ == '__main__':
    pass
