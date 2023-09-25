# -*- coding: utf-8 -*-
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
# # Topic 1
# # Where in the Genome Does DNA Replication Begin?
# ## Secondary Title: Algorithm Warmup
#
# Motivation and some exercises are variations on those available in Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.

# + [markdown] slideshow={"slide_type": "subslide"}
# ### Learning objectives for the week and lab:
# 1. Build our mental model of biology (two ways)
# 2. Warm up our algorithm abilities after summer break
# 3. Gain experience translating a biological problem into a problem we can solve via code

# + [markdown] slideshow={"slide_type": "slide"} hideCode=false hidePrompt=false
# ## Molecular biology for bioinformatics
# Short introduction into some of the concepts we must familiarize ourselves with in this class. We will add to our biological knowledge a little bit at a time.
#
# Credit for some of this content comes from: http://web.stanford.edu/class/cs173/papers/primer.pdf
#
# <img src="https://www.thoughtco.com/thmb/pTDTZI1GH12gJr-DEyozAcVOGvA=/1500x1039/filters:fill(auto,1)/genes-DNA-573388755f9b58723d5eb4fd.jpg" width=800>

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Goals of Molecular Biology
# * Sequencing and comparing full genomes of organisms.
# * Identifying the genes and determining the foundations of the proteins they encode.
# * Understanding gene expression.
# * Understanding genetic diseases.
# * Understanding evolution and evolutionary history.
# * Understanding proteins, which means predicting the folding of the amino acid
# sequence, and characterizing the function of the protein based on this folding.
# * Constructing synthetic proteins, which means creating amino acid sequences,
# such that the protein produced from these have a desired function.

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Polymers
# Three types of polymers will play a roll in this class and in biology: DNA, RNA, and proteins
#
# DNA sequences are the information-containing molecules and are composed of
# nucleotides from an alphabet of four letters: A, C, G and T. The DNA of an organism plays a central role in its existence. It is arranged in
# the form of chromosomes. These strings may be millions of nucleotides long,
# measured in base pairs (bp).
#
# The entire set of genetic information of an organism is called its genome. Genome sizes vary for different species.
#
# <img src="https://ib.bioninja.com.au/_Media/genome-size-table_med.jpeg" widt=300>

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# Proteins, which are the operational molecules, are composed of chains of amino
# acids, called polypeptides, each from an alphabet of 20 letters:
#     
# <img src="https://qph.fs.quoracdn.net/main-qimg-aaa1523b121f7374bb1bdc294d0107ba">

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# RNA sequences, which stand between DNA and protein, are composed of nucleotides from an alphabet of four letters: A, C, G, U.
#
# The Central Dogma of Molecular Biology describes the interaction of these polymers:
# - DNA acts as a template to replicate itself;
# - DNA is also transcribed into RNA; and
# - RNA is translated into protein.

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# <img src="https://ib.bioninja.com.au/_Media/central-dogma_med.jpeg">

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# 1. Replication of DNA.
# Each strand in a DNA is a chemical ”mirror image” of the other. If there is
# an a on one strand, there will always be a t in the same position on the other
# strand, and vice versa; if there is a c on the one strand, its ”partner” on the
# other strand will always be a g, and vice versa.
# When a cell divides to form daughter cells, DNA is replicated by untwisting the
# two strands and using each strand as a template to produce its chemical mirror
# image.
# 2. Transcription of DNA.
# DNA also act as a blueprint for RNA, more exactly three main types of RNA:
# messenger RNA (mRNA), transfer RNA (tRNA), and ribosomal RNA (rRNA).
# They carry information from the genome to the ribosomes, the protein synthesis
# apparatus in a cell.
# 3. Translation of mRNA.
# The information in an mRNA will be translated into a sequence of amino acids,
# creating a polypeptide molecule.3

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# #### More on proteins
#
# "Organic chemistry is the chemistry of carbon compounds. Biochemistry is the study of carbon compounds that crawl." - Mike Adam
#
# Example, human insulin is composed by two words (chains of amino acids):
#
# A: gly ile val glu gln cys cys thr ser ile cys ser leu tyr glu leu glu asn tyr cys asn.
#
# B: phe val asn gln his leu cys gly ser his leu val glu ala leu tyr leu val cys gly glu arg
# gly phe phe tyr thr pro lys thr.
#
# The function of a protein is a direct consequence of its three-dimensional structure

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### More on genes
# Historically, the heritable factors which determine much of the physical make up of
# organisms are called genes.
#
# Usually there are several different forms one gene can have. These forms are called
# alleles.
#
# A combination of alleles describes the make-up of an individual, more exactly:
# * The genetic make-up of an individual is its genotype.
# * The expression of the genes of an individual is its phenotype.

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# **Problem 1.** Find as many different pictures of human insulin protein and insert them here.
#
# **Problem 2.** Why are there different pictures?
#
# **Problem 3.** Now find the gene sequence of insulin and paste it below.

# + [markdown] slideshow={"slide_type": "slide"} hideCode=false hidePrompt=false
# ## Genome Replication
# * One of most important tasks carried out in the cell. 
# * Must be carried out before cell division
# * In 1953, James Watson and Francis Crick ended their paper on DNA double helix with:
#
# "It has not escaped our notice that the specific pairing we have postulated immediately suggests a possible copying mechanism for the genetic material."
#
# <img src="https://images.fineartamerica.com/images/artworkimages/mediumlarge/1/watson-and-crick-a-barrington-brown-and-photo-researchers.jpg" alt="drawing" width="400"/>

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ## Let's talk some biology
# <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/37/Difference_DNA_RNA-EN.svg/1200px-Difference_DNA_RNA-EN.svg.png" width="600"/>
#
#

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ## DNA and RNA code
# * Not a binary alphabet
# * DNA alphabet: AGCT
# * RNA alphabet: AGCU
# * Nucleotides are complementary (A binds to T and G binds to C)
# * Replication begins at replication origin (*ori*)
# * **Binary is a base-2 system, what is DNA/RNA?**

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ## Why should I care?
# There are molecular copy machines known as DNA polymerases that start by locating a *ori*. Some gene therapy methods use genetically engineered mini-genomes, which are called **viral vectors** because they are able to penetrate cell walls. Viral vectors carry artificial genes that have been used to engineer frost-resistant tomatoes and pesticide-resistant corn. In 1990, gene therapy was successfully performed on humans when it saved the life of a four year old girl suffered from Severe Combined Immunodeficiency Disorder. To ensure the treatment works, scientists must know the location of *ori* and avoid disrupting this site.

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ## Looking for *ori*
# Verified *ori* of Vibrio cholerae, the bacterium that causes cholera (~500 nucleotides):
# <pre>
# atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaac
# ctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgacca
# cggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgactt
# gtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggatt
# acgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttagga
# tagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaat
# tgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaag
# atcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtt
# tccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc
# </pre>

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ## DnaA box
# * There is a hidden message in *ori* that orders the cell to begin replication here.
# * We know that the initiation of replication is mediated by a protein called **DnaA** that looks for a short segment within *ori*.
# * This short segment is known as a *DnaA box*
# * Biologists want to find this hidden message, but is that clearly defined enough for us CS/STAT/MATH/EGR folks?

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ## Counting words
# * Turns out that the patterns in our DNA are not random. 
# * Some patterns are more common than others. 
# * Biologically speaking this helps because certain protins can only bind to DNA if a specific string of nucleotides is present and if that string is more prevelant then we have a greater chance of success (and less likely a mutation will cause problems). 
# * We are going to refer to a *k*-mer as a string of length *k*.

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Why? Why? Why?
# "Nothing in biology makes sense except in the light of evolution." - Theodosius Dobzhansky

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# We are looking for surprisingly frequent substrings (contiguous strings appearing within) this *ori*.
# <pre>
# atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaac
# ctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgacca
# cggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgactt
# gtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggatt
# acgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttagga
# tagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaat
# tgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaag
# atcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtt
# tccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc
# </pre>
# Are there any substrings that occur more frequent than others?
#
# Before we go about searching for unknown substrings, we'll write a function that counts the number of occurances of a specific substring.

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
# %load_ext autoreload
# %autoreload 2


# Put all your solutions into Lab1_helper.py as this script which is autograded
import Topic1_helper 

from pathlib import Path
home = str(Path.home()) # all other paths are relative to this path. 
# This is not relevant to most people because I recommended you use my server, but
# change home to where you are storing everything. Again. Not recommended.

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# **Exercise 1.**
# A *k*-mer is a string of length ``k``. For this exercise, define a function ``count(text, pattern)`` as the number of times that a k-mer ``pattern`` appears as a substring of ``text``. For example,
#
# For example:
# <pre>
# count("ACAACTATGCATACTATCGGGAACTATCCT","ACTAT")=3.
# </pre>
# Please note that count("CGATATATCCATAG", "ATA") is equal to 3 (not 2) since we should account for overlapping occurrences of ``pattern`` in ``text``.

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
Topic1_helper.count("ACAACTATGCATACTATCGGGAACTATCCT","ACTAT")

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# **Exercise 2.** Find the most frequent *k*-mers in a string.
# * Input: A string ``text`` and an integer ``k``.
# * Output: All most frequent *k*-mers in ``text`` and their count.
# * Requirements: Do not use a dictionary/map

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
print(Topic1_helper.frequent_words("ACAACTATGCATACTATCGGGAACTATCCT",5))
print(Topic1_helper.frequent_words("ACAACTATGCATACTATCGGGAACTATCCT",4))
# -

# **Problem 4:** What is the Big-O of frequent words? Define |text| as the length of text. Assume the unit of measurement is comparing a single charater (i.e., comparing ABC to DEF costs 3 units).
#
# A. '|text|^2'<br>
# B. '|text|^2*k'<br>
# C. 'k^2'<br>
#
#
# **Your answer here**

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Now let's look at the *ori* and see what 9-mers appear

# + slideshow={"slide_type": "fragment"} hideCode=false hidePrompt=false
text = "atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"
Topic1_helper.frequent_words(text,9)

# + [markdown] slideshow={"slide_type": "fragment"} hideCode=false hidePrompt=false
# Notice anything interesting about the sequences?

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# As previously stated, nucleotides only bind to their complement, so A and T bind and G and C bind. It is also true that DNA is read in specific direction. Very much in the same way we read left to right. DNA is read from what is called the 5' end to the 3' end.
#
# <img src="https://image.slidesharecdn.com/dna-replication-lin-140210083429-phpapp02/95/dna-replicationlin-4-638.jpg?cb=1392021295" width=400/>
#
# So we can now understand and look for something very important called a reverse complement. The definition of which is right there in the name. ACTG is the reverse complement of CAGT. Let's now write a simple funciton to find the reverse complement.

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# **Exercise 3.** Write a function that find the reverse complement of a DNA sequence.
# * Input: A string ``text`` representing DNA.
# * Output: The reverse complement of ``text``.

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
Topic1_helper.reverse_complement("cagt")

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Back to our 9-mers

# + slideshow={"slide_type": "fragment"} hideCode=false hidePrompt=false
solutions = Topic1_helper.frequent_words(text,9)
print(solutions)
print("Reverse complement of first 9-mer:",Topic1_helper.reverse_complement(solutions[0][0]))

# + [markdown] slideshow={"slide_type": "fragment"} hideCode=false hidePrompt=false
# **Problem 4.** What is interesting about the reverse complement of the first 9-mer?
#
# **Your answer here**

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Writing faster code
# **Exercise 4.** Let's now write faster code that produces a frequency map. 
# * Input: A string ``text`` representing DNA and integer ``k``.
# * Output: a frequency map (Python dictionary) that maps every pattern of size ``k`` to the number of times that pattern occurs.

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
freq_map = Topic1_helper.frequency_table(text,3)

# + [markdown] slideshow={"slide_type": "subslide"}
# ## A word about packages
# I try to limit the number of Python packages that you need for this class. They are roughly pandas+numpy and networkx.
#
# In your project, you are welcome to use bioinformatics Python and non-Python packages. You are encouraged to do so.

# + slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# I'm only using pandas here so the output is reasonable, you can remove it of course and see the full dictionary
import pandas as pd
pd.Series(Topic1_helper.frequency_table(text,3))

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Write better frequent words
# **Exercise 5.** Write a function that finds the frequent patterns using a dictionary/map. 
# * Input: A string ``text`` representing DNA and integer ``k``.
# * Output: All most frequent *k*-mers in ``text`` and their count.
# * Requirements: Use your frequency_table function (i.e., use the dictionary).

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
Topic1_helper.better_frequent_words(text,9)

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Clump Finding Problem
# * Imagine you are trying to find *ori* in a newly sequenced genome
# * Old frequent hiddent messages won't be useful
# * One solution is to use a sliding window and look for a region where a $k$-mer appears several times in short succession
# * For example if TGCA forms a (25,3)-clump then it appears at least 3 times in a window of length 25

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# Even if we solve the clump finding problem, we still have an issue
# * Specifically, for the *E. coli* genome we find hundreds of different 9-mers forming (500,3)-clumps
# * This makes it absolutely unclear which of these 9-mers might represent a DnaA box in the bacterium’s *ori* region.
# * Please read the next sections entitled "The Simplest Way to Replicate DNA" and "Asymmetry of Replication". Dig into the biology, but the abstract model/representation we are using in this class does not require you to understand that biology in detail. Chat with me in Slack about what you find confusing and interesting. 

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Statistics of the Foward and Reverse Half-Strands
# The most important consequence for us from the discussion of DNA replication is that we now have four pieces
#     1. Forward half-strand x 2
#     2. Reverse half-strand x 2
#
# <img src="http://bioinformaticsalgorithms.com/images/Replication/half_strands.png" width=400>
#
#

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Why does this matter?
# Consider the genome of *Thermotoga petrophila*. If we count the nucleotides in the forward and reverse half strands, then we get the following:
#
# <img src="http://bioinformaticsalgorithms.com/images/Replication/forward_reverse_nucleotide_counts.png" width=400>
#
# Notice that the number of C's and G's is different in the reverse and forward half-strand. Why is this?

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# Take a minute to read this and then we will discuss together and then we will discuss:
#
# "It turns out that we observe these discrepancies because cytosine (C) has a tendency to mutate into thymine (T) through a process called deamination. Deamination rates rise 100-fold when DNA is single-stranded, which leads to a decrease in cytosine on the forward half-strand. Also, since C-G base pairs eventually change into T-A base pairs, deamination results in the observed decrease in guanine (G) on the reverse half-strand (recall that a forward parent half-strand synthesizes a reverse daughter half-strand, and vice-versa)." - Bioinformatics Algorithms 3rd Edition

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# ### Minimum skew problem
# We can use this statistic to find the *ori*.
#
# Our idea is to traverse the genome, keeping a running total of the difference between the counts of G and C. If this difference starts increasing, then we guess we are on the forward half-strand. If this difference starts decreasing, then we guess that we are on the reverse half-strand.
#
#
# <img src="http://bioinformaticsalgorithms.com/images/Replication/increasing_decreasing_skew.png" width=600>
#

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# We define $Skew_i(Genome)$ as the difference between the total number of occurrences of G and the total number of occurrences of C in the first $i$ nucleotides of Genome. 
#
# Note that we can compute $Skew_i(Genome)$ incrementally.  
#
# If the next nucleotide is G, then $Skew_{i+1}(Genome)$ = $Skew_i(Genome)$ + 1
#
# if this nucleotide is C, then $Skew_{i+1}(Genome)$ = $Skew_i(Genome)$ – 1
#
# otherwise, $Skew_{i+1}(Genome)$ = $Skew_i(Genome)$.

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# **Exercise 6:** Compute the skew at every position of a Genome
#
# Input: A DNA string Genome.
#
# Output: An array that computes the $Skew_i(Genome)$. You can assume $Skew_0(Genome)$=0

# + [markdown] slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# **Read in the *E coli* genome**

# + slideshow={"slide_type": "fragment"} hideCode=false hidePrompt=false
import pandas as pd
data = pd.read_table("http://bioinformaticsalgorithms.com/data/realdatasets/Rearrangements/E_coli.txt",header=None)
genome = data.values[0,0]

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
skews = Topic1_helper.skew(genome)

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
# Again I'm using pandas because the display would be horrible otherwise.
# I will either give you the pandas code or teach you that pandas/numpy code
skews = pd.Series(Topic1_helper.skew(genome))
skews

# + slideshow={"slide_type": "subslide"} hideCode=false hidePrompt=false
# %matplotlib inline
skews.plot.line();

# + [markdown] slideshow={"slide_type": "fragment"} hideCode=false hidePrompt=false
# Where do you think the *ori* is located?

# + slideshow={"slide_type": "fragment"} hideCode=false hidePrompt=false
print('Position:',skews.idxmin()+1)

# + slideshow={"slide_type": "skip"} hideCode=false hidePrompt=false
# Don't forget to push!

# + hideCode=false hidePrompt=false

