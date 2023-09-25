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
# # Topic 8 - Was T. Rex just a big chicken?
# ## Computational Proteomics

# + slideshow={"slide_type": "skip"}
# %matplotlib inline
# %load_ext autoreload
# %autoreload 2

    
import pandas as pd
import numpy as np
from IPython.display import Image

from pathlib import Path
home = str(Path.home()) # all other paths are relative to this path. 
# This is not relevant to most people because I recommended you use my server, but
# change home to where you are storing everything. Again. Not recommended.

# + slideshow={"slide_type": "skip"}
import networkx as nx
import pandas as pd
import numpy as np
import copy

import matplotlib.pyplot as plt

def draw(A):
    return Image(A.draw(format='png', prog='dot'))

def to_adj(T):
    df = pd.DataFrame(nx.adjacency_matrix(T).todense(),index=T.nodes(),columns=T.nodes())
    for i in range(len(df)):
        for j in range(len(df)):
            if df.iloc[i,j] == 1:
                data = T.get_edge_data(df.index[i],df.columns[j])
                df.iloc[i,j] = data['label']
            else:
                df.iloc[i,j] = ""
    return df

def show(G):
    # same layout using matplotlib with no labels
    pos = nx.drawing.nx_agraph.graphviz_layout(G, prog='neato')
    #print(edge_labels)
    # Modify node fillcolor and edge color.
    #D.node_attr.update(color='blue', style='filled', fillcolor='yellow')
    #D.edge_attr.update(color='blue', arrowsize=1)
    A = nx.nx_agraph.to_agraph(G)
    A.graph_attr["rankdir"] = "LR"
    # draw it in the notebook
    display(draw(A))


# + [markdown] slideshow={"slide_type": "subslide"}
# ## Learning Outcomes
# * Understand proteomics applications and data
# * Apply, analyze, and evaluate proteomics algorithms

# + [markdown] slideshow={"slide_type": "slide"}
# # Background
#
# ## Paleontology meets computing (and statistics)
# Controversy ensues (make sure you note the names in the coming articles)...
#
#

# + [markdown] slideshow={"slide_type": "subslide"}
# * Protein sequencing was very difficult in 1950s but DNA sequencing was impossible
# * Today DNA sequencing is essentially trivial while protein sequencing remains difficult
# * Dinosaur drama (https://www.wired.com/2009/06/ff-originofspecies/)
#
# <img src="https://images-na.ssl-images-amazon.com/images/I/81ZfeGANs9L._AC_SL1500_.jpg" width=400>

# + [markdown] slideshow={"slide_type": "subslide"}
# <img src="https://www.researchgate.net/profile/Jeovanis_Gil_Valdes2/publication/321347984/figure/fig2/AS:581671432593414@1515692684716/General-workflow-for-the-mass-spectrometry-based-proteomics-analysis-of-acetylated.png">

# + [markdown] slideshow={"slide_type": "slide"}
# ## Matching Ideal Spectra
# This first section of exercises will allow you to implement algorithms under ideal circumstances where masses are not lost and there is no "noise" mass.
#
# In reality... we have to deal with false and missing masses. One source of this misinformation is that when a mass spectrometer breaks a peptide, small parts of the resulting fragments may be lost, thus lowering their masses (location on the x axis). When breaking REDCA into RE, DCA, RE might lose a water molecule and DCA might lost an ammonia.

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 1**: Spectrum graph construction
#
# We represent the masses in a spectrum as a sequence Spectrum of integers in increasing order. We define a labeled graph Graph(Spectrum) by forming a node for each element of ``spectrum``, then connecting nodes $s_i$ and $s_j$ by a directed edge labeled by an amino acid $a$ if $s_j−s_i$ is equal to the mass of $a$. We do not distinguish between amino acids having the same integer masses (i.e., the pairs K/Q and I/L). You'll need to add in a zero mass node at the beginning and the sum of all the masses as the last element.
#
# Input: A list of integers ``spectrum``.
#
# Output: A networkx graph that represents the graph described above.

# + [markdown] slideshow={"slide_type": "subslide"}
# ### OK... What?
#
# Let's break down what we know. We know the weights of amino acids

# + slideshow={"slide_type": "subslide"}
a_mass = {
    "G": 57,
    "A": 71,
    "S": 87,
    "P": 97,
    "V": 99,
    "T": 101,
    "C": 103,
    "I": 113,
    "L": 113,
    "N": 114,
    "D": 115,
    "K": 128,
    "Q": 128,
    "E": 129,
    "M": 131,
    "H": 137,
    "F": 147,
    "R": 156,
    "Y": 163,
    "W": 186
}
import pandas as pd
pd.Series(a_mass).plot.bar();

# + [markdown] slideshow={"slide_type": "subslide"}
# ### We can put that information in a dictionaries that map amino acids to masses and masses to amino acids

# + slideshow={"slide_type": "fragment"}
import networkx as nx

mass_a = {}
for key in a_mass.keys():
    mass = a_mass[key]
    if mass not in mass_a:
        mass_a[mass] = []
    mass_a[mass].append(key)
    
pd.Series(mass_a) # Just for printing

# + [markdown] slideshow={"slide_type": "subslide"}
# ### What about a real spectrum?
# <img src="http://bioinformaticsalgorithms.com/images/Proteomics/dinosaur_spectrum_unannotated.png">

# + [markdown] slideshow={"slide_type": "subslide"}
# ### You could take a look at a spectrum in a very simplified manner

# + slideshow={"slide_type": "subslide"}
spectrum1 = [57,71,154,185,301,332,415,429,486]
pd.Series(1+0*np.array(spectrum1),index=spectrum1).plot.bar();


# + slideshow={"slide_type": "subslide"}
def spectrum_graph_construction(spectrum,mass_a=mass_a):
    spectrum = copy.copy(spectrum)
    spectrum.insert(0,0)
    G = nx.DiGraph()
    for i,s in enumerate(spectrum):
        G.add_node(s)
    # Your solution here
    return G

spectrum1 = [57,71,154,185,301,332,415,429,486]
graph1 = spectrum_graph_construction(spectrum1)

show(graph1)


# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 2**: Ideal spectrum
#
# Input: A ``peptide``
#
# Output: A list that represents the masses of the ideal spectrum

# + slideshow={"slide_type": "subslide"}
# fragments is for debugging purposes
def ideal_spectrum(peptide,a_mass=a_mass,prefix=True,suffix=True,fragments=None):
    if fragments is None:
        fragments = []
    ideal = [0]
    # Your solution here
    ideal.sort()
    return ideal

peptide1 = "GPG"
spectrum2 = ideal_spectrum(peptide1)
fragments = []
spectrum3 = ideal_spectrum("REDCA",fragments=fragments)

print(f"Spectrum for {peptide1}")
print(spectrum2)
print(f"Fragments for REDCA")
print(fragments)
print(f"Spectrum for REDCA")
print(spectrum3)
# -

fragments=[]
print(ideal_spectrum('GASDGG',fragments=fragments))
print(fragments)

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 3**: Decoding the ideal spectrum
#
# Input: A ``spectrum``
#
# Output: All matching strings representing the peptide that corresponds to a path from *source* to *sink* in Graph(spectrum) whose ideal spectrum is equal to ``spectrum``.
# -

ideal_spectrum('NTDN')

# + slideshow={"slide_type": "subslide"}
import numpy as np

def decoding_ideal_spectrum(spectrum,a_mass=a_mass,debug=False):
    mass_a = {}
    for key in a_mass.keys():
        mass = a_mass[key]
        if mass not in mass_a:
            mass_a[mass] = []
        mass_a[mass].append(key)
    G = spectrum_graph_construction(spectrum,mass_a=mass_a)
    if debug:
        show(G)
    # Your solution here
    matches = []
    return matches

spectrum5 = [57,114,128,215,229,316,330,387,444]
peptides5 = decoding_ideal_spectrum(spectrum5,debug=True)

print(peptides5)


# + [markdown] slideshow={"slide_type": "subslide"}
# ### Again... What?
#
# Take one of the solutions and take a look at how much it weights.

# + slideshow={"slide_type": "fragment"}
s=0
for c in 'GGTTQ':
    s += a_mass[c]
s

# + [markdown] slideshow={"slide_type": "subslide"}
# ### But does it match the spectrum?
#
# <pre>
# spectrum5 = [57,114,128,215,229,316,330,387,444]</pre>

# + slideshow={"slide_type": "fragment"}
ideal_spectrum("GGTTQ")
# -

# ## From Ideal to Real Spectra

# **Exercise 4**: Convert a peptide into a binary peptide vector.
#
# Input: A peptide P.
#
# Output: The peptide vector of P as a numpy array.

# +
import numpy as np

def construct_peptide_vector(peptide,a_mass={"X":4,"Z":5},verbose=False):
    total_mass = sum([a_mass[c] for c in peptide])
    vector = np.zeros((total_mass),dtype=int)
    # Your solution here
    return vector

peptide_v1 = construct_peptide_vector("XZZXX")

print(peptide_v1)
display(construct_peptide_vector("XZZXX",verbose=True))


# -

# **Exercise 5**: Convert a binary vector into a peptide (do the reverse).
#
# Input: A numpy binary vector ``p``.
#
# Output: A peptide whose binary peptide vector matches ``p``. For masses with more than one amino acid, any choice may be used.

# +
def construct_peptide_from_vector(p,a_mass={"X":4,"Z":5}):
    peptides = []
    # Your solution here
    return peptides

p = np.array([0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1])
peptides6 = construct_peptide_from_vector(p)

peptides6


# -

# **Exercise 6:** Peptide Sequencing Problem
#
# Input: A spectral vector ``s``.
#
# Output: A peptide with maximum score against ``s``. For masses with more than one amino acid, any choice may be used.

# +
def max_peptide(s,a_mass={"X":4,"Z":5},debug=False):
    peptide = ""
    mass_a = {}
    for key in a_mass.keys():
        mass = a_mass[key]
        if mass not in mass_a:
            mass_a[mass] = []
        mass_a[mass].append(key)
    # Your solution here
    return peptide

p2 = [0,0,0,4,-2,-3,-1,-7,6,5,3,2,1,9,3,-8,0,3,1,2,1,0]
peptide7 = max_peptide(p2,debug=True)

peptide7

# + [markdown] slideshow={"slide_type": "subslide"}
# ### So again... what?
#
# The node labels mean something. The first integer is the mass. There is a transition between 9:6 and 13:1 because X is equal to 4 and so is 13-9. Why do we have a graph? Well. We want to score different options that are paths from this graph. How do we score? Once we have a path we can sum it up using the second integers from the nodes. So a path that went from 9:6 to 13:1 would get a score of ...+6+1+...
# -


