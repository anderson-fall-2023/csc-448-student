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

# + [markdown] slideshow={"slide_type": "slide"}
# # Assignment 3 - Genome Assembly
# ## Graph Algorithms

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Learning Outcomes
# * Understand genome assembly bioinformatics motivation
# * Apply, analyze, and evaluate genome sequence alignment graph based algorithms
#
# In this assignment, we will:
# 1. Verify that your solutions for exercises in Topic 3 are implemented correctly
# 2. Extend these solutions
# 3. Apply them to scenarios of both success and failure

# + slideshow={"slide_type": "skip"}
# %matplotlib inline
# %load_ext autoreload
# %autoreload 2


import pandas as pd
import numpy as np

import Assignment3_helper 

from pathlib import Path
home = str(Path.home()) # all other paths are relative to this path. 
# This is not relevant to most people because I recommended you use my server, but
# change home to where you are storing everything. Again. Not recommended.

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 1** String composition problem: Generate the k-mer composition of a string.
#     
# Input: An integer $k$ and a string $text$.
#
# Output: The lexicographic ordered collection of all k-mer substrings in text (including repeated k-mers) ``composition(k,text)``, where the k-mers are arranged in lexicographic order.
#
# Example: ``composition(3,TATGGGGTGC)`` => [ATG, GGG, GGG, GGT, GTG, TAT, TGC, TGG]. Notice how there are two entries for GGG, so we keep duplicates.
#
# Learning objective: Build up the skills and vocabulary we need to solve the reconstruction problem.

# + slideshow={"slide_type": "subslide"}
kmers = Assignment3_helper.composition(3,"TATGGGGTGC")
kmers

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 2** Construct the de Bruijn graph from a set of k-mers.
#
# Input: a collection of $k$-mers ``patterns``.
#
# Output: NetworkX directed graph that corresponds to the adjacency list of the de Bruijn graph. 

# + slideshow={"slide_type": "subslide"}
dB = Assignment3_helper.de_bruijn(["AAT","ATG","ATG","ATG","CAT","CCA","GAT","GCC","GGA","GGG","GTT","TAA","TGC","TGG","TGT"])
Assignment3_helper.show(dB)

# + slideshow={"slide_type": "subslide"}
Assignment3_helper.to_adj(dB)

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 3** Given a Eulerian directed graph, find the Eulerian cycle.
#
# Given: An Eulerian directed graph, in the form of a MultiDiGraph networkx object.
#
# Return: An Eulerian cycle in this graph returned as a list.

# + slideshow={"slide_type": "skip"}
# Print example with random start
cycle = Assignment3_helper.eulerian_cycle(Assignment3_helper.G)
Assignment3_helper.show(Assignment3_helper.G)
print("Random start")
print(cycle)

# Example with start node specified
cycle = Assignment3_helper.eulerian_cycle(Assignment3_helper.G,start=6)

print("Start at specific node")
print(cycle) # Should also result in a shifted but equivalent answer. 

# If your code isn't working, please carefully read the section on the ant 
# moving around the graph and getting stuck :)

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Eulerian path
# Consider the following graph that does not contain a Eulerian cycle.

# + slideshow={"slide_type": "fragment"}
Assignment3_helper.show(Assignment3_helper.G2)

# + [markdown] slideshow={"slide_type": "subslide"}
# Consider the following function that finds the in and out degree of every node.

# + slideshow={"slide_type": "fragment"}
Assignment3_helper.calc_in_out(Assignment3_helper.G2)

# + [markdown] slideshow={"slide_type": "subslide"}
# Is there "nearly" a cycle? If so, how would you find it? What two nodes would you connect?

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 4** Given a Eulerian directed graph, find the Eulerian path.
#
# Given: An Eulerian directed graph, in the form of a MultiDiGraph networkx object.
#
# Return: An Eulerian cycle in this graph returned as a list.
#
# Hint: Can you find the correct starting node from the output of ``in_out``?

# + slideshow={"slide_type": "subslide"}
path = Assignment3_helper.eulerian_path(Assignment3_helper.G2)
path

# + [markdown] slideshow={"slide_type": "subslide"}
# We can finally put everything together!
#
# **Exercise 5:** Reconstruct a string from its k-mer composition.
#
# Given: Given a list of $k$-mers patterns.
#
# Return: A string ``text`` with $k$-mer composition equal to patterns.

# + slideshow={"slide_type": "subslide"}
assembly = Assignment3_helper.reconstruct(Assignment3_helper.kmers)
assembly
# -

# Don't forget to push!

# ## Important notes
# Problems below are not autograded.
#
# They do require you to copy your align_dynamic2 function into the Assignment3_helper.py in order to get the print outs. It also requires you to copy your function for reading a FASTA file.

# We've written the bones of a modern assembly algorithm.

seq="TATGGGGTGC"
kmers = Assignment3_helper.composition(3,seq)
assembled_seq=Assignment3_helper.reconstruct(kmers)
assembled_seq

# We have an alignment algorithm so we can check it out:

sc,s1,s2=Assignment3_helper.align_dynamic2(seq,assembled_seq)
Assignment3_helper.print_alignment(s1,s2)

# Looks like we solved it! Time to go home...
#
# Let's take a look at Yeast's genome. We will read in all the chromosomes, but for simplicity we will focus on the first 1,000 nucleotides of the first chromosome.

headers,seqs=Assignment3_helper.read_fasta(f'{home}/csc-448-student/data/GCF_000146045.2_R64_genomic.fna')
headers[0]

ch1_first_1000 = seqs[0][:1000].upper()

ch1_first_1000

# **Problem 1:** Do you think we will get the correct assembly using the same parameters and functions that were so successful above? Explain your answer and your reasoning.
#
# **Your answer here**

kmers = Assignment3_helper.composition(3,ch1_first_1000)
assembled_seq=Assignment3_helper.reconstruct(kmers)
assembled_seq

Assignment3_helper.print_alignment(ch1_first_1000,assembled_seq,num_to_print=50)

# How does that look to you? Take a look at all C's in a row. And all the T's and the altnerating AC's. Ahhh....

# The real question is not whether this assembly is wrong (it clearly is), but the question is what caused this problem? There are a number of things to consider:
# 1. It could be a bug in our program (it isn't :)
# 2. It could be that our program returned one of many possible assemblies (ding ding ding!). 
#
# It is the second complication. With a small $k$-mer value we run the risk of too much variability in the assemblies that are possible. 

# **Problem 2:** See how the results change as you increase the $k$ value:

kmers=Assignment3_helper.composition(11,ch1_first_1000)
assembled_seq=Assignment3_helper.reconstruct(kmers)
assembled_seq

Assignment3_helper.print_alignment(ch1_first_1000,assembled_seq,num_to_print=50)

# Well that's better! But there are still a couple of errors. Try increasing the k-value and see if we can correct the problem.

kmers=Assignment3_helper.composition(13,ch1_first_1000)
assembled_seq=Assignment3_helper.reconstruct(kmers)
assembled_seq

Assignment3_helper.print_alignment(ch1_first_1000,assembled_seq,num_to_print=50)

# So now our question becomes, is this it? Can we call it a day? Well don't forget that Yeast has more than one chromosome. If you construct kmers from two chromosomes, then our algorithm falls apart and needs significant improvements.

ch2_first_30 = seqs[1][:30].upper()
ch2_first_30

# **Warning:** The following cell is supposed to produce an error as the graph produced by combining the kmers from the two chromosomes is not structured correctly for our assembler.

kmers=Assignment3_helper.composition(13,ch1_first_1000)+Assignment3_helper.composition(13,ch2_first_30)
assembled_seq=Assignment3_helper.reconstruct(kmers)
assembled_seq
