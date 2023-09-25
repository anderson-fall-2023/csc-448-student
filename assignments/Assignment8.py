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

# + [markdown] slideshow={"slide_type": "slide"}
# # Assignment 8 - Was T-Rex a giant chicken?
#
# ## Proteomics

# + slideshow={"slide_type": "skip"}
# %matplotlib inline
# %load_ext autoreload
# %autoreload 2

    
import pandas as pd
import numpy as np

import Assignment8_helper 

from pathlib import Path
home = str(Path.home()) # all other paths are relative to this path. 
# This is not relevant to most people because I recommended you use my server, but
# change home to where you are storing everything. Again. Not recommended.

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 1**: Spectrum graph construction
#
# We represent the masses in a spectrum as a sequence Spectrum of integers in increasing order. We define a labeled graph Graph(Spectrum) by forming a node for each element of ``spectrum``, then connecting nodes $s_i$ and $s_j$ by a directed edge labeled by an amino acid $a$ if $s_j−s_i$ is equal to the mass of $a$. We do not distinguish between amino acids having the same integer masses (i.e., the pairs K/Q and I/L). You'll need to add in a zero mass node at the beginning and the sum of all the masses as the last element.
#
# Input: A list of integers ``spectrum``.
#
# Output: A networkx graph that represents the graph described above.

# + slideshow={"slide_type": "subslide"}
spectrum1 = [57,71,154,185,301,332,415,429,486]
graph1 = Assignment8_helper.spectrum_graph_construction(spectrum1)

Assignment8_helper.show(graph1)

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 2**: Ideal spectrum
#
# Input: A ``peptide``
#
# Output: A list that represents the masses of the ideal spectrum

# + slideshow={"slide_type": "subslide"}
peptide1 = "GPG"
spectrum2 = Assignment8_helper.ideal_spectrum(peptide1)
fragments = []
spectrum3 = Assignment8_helper.ideal_spectrum("REDCA",fragments=fragments)

print(f"Spectrum for {peptide1}")
print(spectrum2)
print(f"Fragments for REDCA")
print(fragments)
print(f"Spectrum for REDCA")
print(spectrum3)

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 3**: Decoding the ideal spectrum
#
# Input: A ``spectrum``
#
# Output: All matching strings representing the peptide that corresponds to a path from *source* to *sink* in Graph(spectrum) whose ideal spectrum is equal to ``spectrum``.

# + slideshow={"slide_type": "subslide"}
import numpy as np

spectrum5 = [57,114,128,215,229,316,330,387,444]
peptides5 = Assignment8_helper.decoding_ideal_spectrum(spectrum5)

print(peptides5)

# + [markdown] slideshow={"slide_type": "slide"}
# ## From Ideal to Real Spectra

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 4**: Convert a peptide into a binary peptide vector.
#
# Input: A peptide P.
#
# Output: The peptide vector of P as a numpy array.

# + slideshow={"slide_type": "subslide"}
import numpy as np

peptide_v1 = Assignment8_helper.construct_peptide_vector("XZZXX")

print(peptide_v1)
display(Assignment8_helper.construct_peptide_vector("XZZXX",verbose=True))

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 5**: Convert a binary vector into a peptide (do the reverse).
#
# Input: A numpy binary vector ``p``.
#
# Output: A peptide whose binary peptide vector matches ``p``. For masses with more than one amino acid, any choice may be used.

# + slideshow={"slide_type": "subslide"}
p = np.array([0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1])
peptides6 = Assignment8_helper.construct_peptide_from_vector(p)

peptides6

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 6:** Peptide Sequencing Problem
#
# Input: A spectral vector ``s``.
#
# Output: A peptide with maximum score against ``s``. For masses with more than one amino acid, any choice may be used.

# + slideshow={"slide_type": "subslide"}
p2 = [0,0,0,4,-2,-3,-1,-7,6,5,3,2,1,9,3,-8,0,3,1,2,1,0]
peptide7 = Assignment8_helper.max_peptide(p2,debug=True)

peptide7

# + slideshow={"slide_type": "skip"}
# Don't forget to push!
