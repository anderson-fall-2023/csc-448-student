---
jupyter:
  jupytext:
    formats: ipynb,md,py
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.8.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<!-- #region slideshow={"slide_type": "slide"} hideCode=false hidePrompt=false -->
# Comparing DNA Sequences

## Dynamic programming

In this assignment, we will:
1. Verify that your solutions for exercises in Topic 2 are implemented correctly
2. Extend these solutions
3. Apply them to scenarios of both success and failure
<!-- #endregion -->

## Verifying solutions to Topic 2 are correct
Please see Topic2 notebook for detailed information. You must copy Topic2_helper.py to this folder for these commands to work.

```python slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
%load_ext autoreload
%autoreload 2

## BEGIN SOLUTION
import joblib
answers = {}
## END SOLUTION

import pandas as pd
import numpy as np

import Assignment2_helper 
import Topic2_helper

from pathlib import Path
home = str(Path.home()) # all other paths are relative to this path. 
# This is not relevant to most people because I recommended you use my server, but
# change home to where you are storing everything. Again. Not recommended.
```

## Exercise 1

```python
## BEGIN SOLUTION
answers["answer_exercise_1"] = Topic2_helper.greedy_lcs("AACCTTGG","ACACTGTGA")
## END SOLUTION

print(Topic2_helper.greedy_lcs("AACCTTGG","ACACTGTGA"))

print(Topic2_helper.greedy_lcs("AACCTTGG","ACAC"))
```

## Exercise 2

```python
## BEGIN SOLUTION
answers["answer_exercise_2"] = Topic2_helper.greedy_alignment("AACCTTGG","ACACTGTGA")
## END SOLUTION

print(Topic2_helper.greedy_alignment("AACCTTGG","ACACTGTGA"))
print()
print(Topic2_helper.greedy_alignment("AACCTTGG","ACAC"))
```

## Exercise 3

```python
## BEGIN SOLUTION
answers["answer_exercise_3"] = Topic2_helper.min_num_coins(27,[6,5,1])
## END SOLUTION

Topic2_helper.min_num_coins(27,[6,5,1])
```

## Exercise 4

```python
## BEGIN SOLUTION
answers["answer_exercise_4"] = Topic2_helper.min_num_coins_dynamic(27,[6,5,1])
## END SOLUTION

Topic2_helper.min_num_coins_dynamic(27,[6,5,1])
```

## Exercise 5

```python
## BEGIN SOLUTION
answers["answer_exercise_5"] = Topic2_helper.align("AACCT","ACACTG")
## END SOLUTION

score, aligned_s1, aligned_s2 = Topic2_helper.align("AACCT","ACACTG")
print(score)
print(aligned_s1)
print(aligned_s2)
```

## Exercise 6

```python
## BEGIN SOLUTION
answers["answer_exercise_6"] = Topic2_helper.align_dynamic("AACCT","ACACTG")
## END SOLUTION

score = Topic2_helper.align_dynamic("AACCT","ACACTG")
score
```

## Exercise 7

```python
## BEGIN SOLUTION
answers["answer_exercise_7"] = Topic2_helper.align_dynamic2("AACCT","ACACTG")
## END SOLUTION

score,s1_aligned,s2_aligned = Topic2_helper.align_dynamic2("AACCT","ACACTG")
print(score)
print(s1_aligned)
print(s2_aligned)
```

## Exercise 8. Different scoring functions

For exercise 8, modify align_dynamic3 to correctly incorporate a match, mismatch, and gap penalty.

```python
s1="CGCAACCACAGCGCGCAGGGCAGGCGCGAGCTGTCTGAGCCCCGGCCTCGGACCGCCCACTGGACTCCCGGCACGCCCGGTGCCGCCTTCCGGCTCCAGTCCCCC"
s2="CGCAACGGCAGCGCGCAGGGCAGGCGCGAGCTGGCCTCTGAGCCCCGGCCTCGGACCGCCCACTCCACGCCCGGCAGGCCCGGTGCCGCCTTCCGGCTCCAGTCCCCCCGC"
score_1,aligned_s1_1,aligned_s2_1 = Assignment2_helper.align_dynamic3(s1,s2,match_score=1,mismatch_score=0,gap_score=0)
score_2,aligned_s1_2,aligned_s2_2 = Assignment2_helper.align_dynamic3(s1,s2,match_score=2,mismatch_score=-3,gap_score=-1)
## BEGIN SOLUTION
answers["answer_exercise_8"] = score_1,score_2
## END SOLUTION
score_1,score_2
```

**Problem 1:** Now that you have implemented the function and verified that your score works, please consider the resulting alignments. Discuss how the gap, match, and mismatch penalty have changed the alignment. Feel free to experiment with other parameter values and sequences. Speculate on how a biologist would decide on the best match, mismatch, and penalties.

```python
Assignment2_helper.print_alignment(aligned_s1_1,aligned_s2_1)
```

```python
Assignment2_helper.print_alignment(aligned_s1_2,aligned_s2_2)
```

**Your answer here**

```python slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
## BEGIN SOLUTION
joblib.dump(answers,"../tests/answers_Assignment2.joblib");
## END SOLUTION
# Don't forget to push!
```

```python hideCode=false hidePrompt=false

```
