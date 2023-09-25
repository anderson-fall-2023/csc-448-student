import numpy as np

def count(text,pattern):
    count = 0
    # YOUR SOLUTION HERE
    return count

def frequent_words(text,k):
    frequent_patterns = []
    counts = []
    # YOUR SOLUTION HERE
    return list(np.unique(frequent_patterns)),max_count

def reverse_complement(text):
    text = text[::-1].lower()
    chars = []
    # YOUR SOLUTION HERE
    return "".join(chars)

def frequency_table(text,k):
    freq_map = {}
    n = len(text)
    for i in range(n-k+1):
        # YOUR SOLUTION HERE
        pass
    return freq_map

def better_frequent_words(text,k):
    frequent_patterns = []
    freq_map = frequency_table(text,k)
    # YOUR SOLUTION HERE
    return frequent_patterns,max_value

def skew(genome):
    skews = [0]
    # YOUR SOLUTION HERE
    return skews