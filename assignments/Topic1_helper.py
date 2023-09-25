import numpy as np

def count(text,pattern):
    count = 0
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    for i in range(len(text)-len(pattern)):
        if text[i:(i+len(pattern))] == pattern:
            count += 1
    ## END SOLUTION
    return count

def frequent_words(text,k):
    frequent_patterns = []
    counts = []
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    for i in range(len(text)-k):
        pattern = text[i:(i+k)]
        counts.append(count(text,pattern))
    max_count = max(counts)
    for i in range(len(text)-k):
        if counts[i] == max_count:
            frequent_patterns = frequent_patterns + [text[i:(i+k)]]
    ## END SOLUTION
    return list(np.unique(frequent_patterns)),max_count

def reverse_complement(text):
    text = text[::-1].lower()
    chars = []
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    for i,c in enumerate(text):
        if c == "a":
            chars.append("t")
        elif c == "t":
            chars.append("a")
        elif c == "c":
            chars.append("g")
        elif c == "g":
            chars.append("c")
        else:
            chars.append("?")
    ## END SOLUTION
    return "".join(chars)

def frequency_table(text,k):
    freq_map = {}
    n = len(text)
    for i in range(n-k+1):
        # YOUR SOLUTION HERE
        ## BEGIN SOLUTION
        pattern = text[i:(i+k)]
        if pattern not in freq_map:
            freq_map[pattern] = 0
        freq_map[pattern]+=1
        ## END SOLUTION
        pass
    return freq_map

def better_frequent_words(text,k):
    frequent_patterns = []
    freq_map = frequency_table(text,k)
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    values = freq_map.values()
    max_value = max(values)
    for key in freq_map.keys():
        if freq_map[key] == max_value:
            frequent_patterns.append(key)
    ## END SOLUTION
    return frequent_patterns,max_value

def skew(genome):
    skews = [0]
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    for i in range(0,len(genome)):
        if genome[i] == "G":
            skews.append(skews[i]+1)
        elif genome[i] == "C":
            skews.append(skews[i]-1)
        else:
            skews.append(skews[i])
    ## END SOLUTION
    return skews