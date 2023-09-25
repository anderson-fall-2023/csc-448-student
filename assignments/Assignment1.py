# ---
# jupyter:
#   jupytext:
#     formats: ipynb,md,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.8.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# + [markdown] slideshow={"slide_type": "slide"} hideCode=false hidePrompt=false
# # Origin of Genome Replication
#
# In this assignment, we will:
# 1. Verify that your solutions for exercises in Topic 1 are implemented correctly
# 2. Extend these solutions
# 3. Apply them to scenarios of both success and failure
# -

# ## Verifying solutions to Topic 1 are correct
# Please see Topic1 notebook for detailed information. You must copy Topic1_helper.py to this folder for these commands to work.

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
# %load_ext autoreload
# %autoreload 2


import pandas as pd
import numpy as np

# Put all your solutions into Assignment1_helper.py as this is the script which is autograded
import Assignment1_helper 
import Topic1_helper

from pathlib import Path
home = str(Path.home()) # all other paths are relative to this path. 
# This is not relevant to most people because I recommended you use my server, but
# change home to where you are storing everything. Again. Not recommended.
# -

Topic1_helper.count("ACAACTATGCATACTATCGGGAACTATCCT","ACTAT")

print(Topic1_helper.frequent_words("ACAACTATGCATACTATCGGGAACTATCCT",5))
print(Topic1_helper.frequent_words("ACAACTATGCATACTATCGGGAACTATCCT",4))

Topic1_helper.reverse_complement("cagt")

text = "atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"
freq_map = Topic1_helper.frequency_table(text,3)
pd.Series(freq_map)

Topic1_helper.better_frequent_words(text,9)

data = pd.read_table("http://bioinformaticsalgorithms.com/data/realdatasets/Rearrangements/E_coli.txt",header=None)
genome = data.values[0,0]
skews = Topic1_helper.skew(genome)
skews = pd.Series(Topic1_helper.skew(genome))
skews

# %matplotlib inline
skews.plot.line();

# ## Moving past prokaryotic organisms
#
# We developed a bioinformatics algorithm that at a glance seems to detect origin sequences in a genome, but a closer inspection reveals that we have only tested this algorithm for a single prokaryotic organism. This *skew* algorithm will perform well on organisms with cirucular chromosomes, but what happens when we try it on an organism with linear chromosomes? Let's take a look at Yeast. 
#
# ### Yeast Origin Database
# It is a difficult biological question when you consider more complicated genomes. Even genomes that are well studied such as Yeast. Here is a link to a database with origin sequences: 
#
# http://cerevisiae.oridb.org/search.php?chr=all&confirmed=true&likely=true&dubious=true&name= 
#
# Notice anything? I notice the number of confirmed and likely origin locations on each chromosome. Oh yeah. We are now dealing with multiple chromosomes as well of course. Let's see how our algorithm does on this data. 
#
# This is a good time to talk about file formats. Bioinformatics loves file formats and they often are not everything we would desire. A common sequence format is called a FASTA file. I have downloaded the Yeast genome for us, and we can take a look at the top:

# !head {home}/csc-448-student/data/GCF_000146045.2_R64_genomic.fna

# The first line showing the name of the chromosome (but could be anything) begins with a ``>`` sign. If we look through the file, we will find a few of these lines, which means that this file actually contains multiple sequences (in this case multiple chromosomes).
#
# The other thing you will notice is that the sequence is actually broken up over multiple lines. The destinction between lines is purely formatting based. It has nothing to do biologically, so when we read in the file, we must remove newlines.

# !grep chromosome {home}/csc-448-student/data/GCF_000146045.2_R64_genomic.fna

# **Exercise 7.** Write a function that reads in the sequences in a FASTA file. You must handle the newlines and the headers.

file = f"{home}/csc-448-student/data/GCF_000146045.2_R64_genomic.fna"
headers,sequences = Assignment1_helper.read_fasta(file)
headers

pd.Series(sequences)

# **Problem 1:** What is the average length of a chromosome in Yeast?

avg = None
#### Your solution here
print(f"{avg} bp")

# **Problem 2:** Apply your skew function to the first 2000 bp of the first chromosome. Graph the skew as previously done. Are there any clear origin of replication locations? If not, do you have any guesses as to why? We will discuss in class with Dr. Davidson. Write your guesses down here, so we can impress her or alternatively give her a good laugh :)
#
# Don't forget that you can also look into the database of known sites here: http://cerevisiae.oridb.org/search.php?chr=all&confirmed=true&likely=true&dubious=true&name= 

import pandas as pd
skews = pd.Series(Topic1_helper.skew(sequences[0][:2000]))
skews.plot.line();

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
# Don't forget to push!

# + hideCode=false hidePrompt=false

