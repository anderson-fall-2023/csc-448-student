import networkx as nx
import pandas as pd
import numpy as np
import copy
from IPython.display import Image

import matplotlib.pyplot as plt

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

mass_a = {}
for key in a_mass.keys():
    mass = a_mass[key]
    if mass not in mass_a:
        mass_a[mass] = []
    mass_a[mass].append(key)

def spectrum_graph_construction(spectrum,mass_a=mass_a):
    spectrum = copy.copy(spectrum)
    spectrum.insert(0,0)
    G = nx.DiGraph()
    for i,s in enumerate(spectrum):
        G.add_node(s)
    # Your solution here
    ## BEGIN SOLUTION
    masses = list(mass_a.keys())
    for i,s1 in enumerate(spectrum):
        for j,s2 in enumerate(spectrum):
            if s2-s1 in masses:
                G.add_edge(s1,s2,label="/".join(mass_a[s2-s1]))
            
    ## END SOLUTION
    return G

# fragments is for debugging purposes
def ideal_spectrum(peptide,a_mass=a_mass,prefix=True,suffix=True,fragments=None):
    if fragments is None:
        fragments = []
    ideal = [0]
    # Your solution here
    ## BEGIN SOLUTION
    i = 0
    if prefix:
        for j in range(i,len(peptide)):
            s = 0
            k = i
            fragment = ""
            while k < (i+j+1):
                c = peptide[k]
                s += a_mass[c]
                fragment += c
                k += 1
            ideal.append(s)
            fragments.append(fragment)
    else:
        j = len(peptide)-1
    if suffix:
        for i in range(1,len(peptide)):
            s = 0
            k = i
            fragment = ""
            while k < len(peptide):#(i+j+1):
                c = peptide[k]
                s += a_mass[c]
                fragment+=c
                k += 1
            #for c in peptide[i:(i+j+1)]:
            #    s += a_mass[c]
            fragments.append(fragment)
            ideal.append(s)        
        
    #print(peptide[i:(len(peptide)-j)])
    ## END SOLUTION
    ideal.sort()
    return ideal

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
    ## BEGIN SOLUTION
    for path in nx.all_simple_paths(G,0,spectrum[-1]):
        peptides = [""]
        for i in range(len(path)-1):
            a = G.get_edge_data(path[i],path[i+1])['label']
            if "/" in a:
                a1,a2 = a.split("/")
                for i in range(len(peptides)):
                    peptides.append(peptides[i]+a2)
                    peptides[i] += a1
            else:
                for i in range(len(peptides)):
                    peptides[i] += a
        for peptide in peptides:
            if np.all(ideal_spectrum(peptide,a_mass=a_mass)==[0]+spectrum):
                matches.append(peptide)
    #print(peptide[i:(len(peptide)-j)])
    ## END SOLUTION
    return matches

def construct_peptide_vector(peptide,a_mass={"X":4,"Z":5},verbose=False):
    total_mass = sum([a_mass[c] for c in peptide])
    vector = np.zeros((total_mass),dtype=int)
    # Your solution here
    ## BEGIN SOLUTION
    fragments=[]
    index = ["" for _ in range(len(vector))]
    spectrum = ideal_spectrum(peptide,a_mass=a_mass,suffix=False,fragments=fragments)
    for i,m in enumerate(spectrum):
        vector[m-1] = 1
        index[m-1] = fragments[i-1]
    if verbose:
        return pd.Series(vector,index=index)
    ## END SOLUTION
    return vector

def construct_peptide_from_vector(p,a_mass={"X":4,"Z":5}):
    peptides = []
    # Your solution here
    ## BEGIN SOLUTION
    prefix_spectrum = list(np.where(p == 1)[0]+1)
    suffix_spectrum = []
    for m in prefix_spectrum:
        suffix_spectrum.append(prefix_spectrum[-1]-m)
    spectrum = prefix_spectrum + suffix_spectrum
    spectrum = list(np.sort(spectrum))[1:]
    peptides = decoding_ideal_spectrum(spectrum,a_mass=a_mass)
    ## END SOLUTION
    return peptides

def max_peptide(s,a_mass={"X":4,"Z":5},debug=False):
    peptide = ""
    mass_a = {}
    for key in a_mass.keys():
        mass = a_mass[key]
        if mass not in mass_a:
            mass_a[mass] = []
        mass_a[mass].append(key)
    # Your solution here
    ## BEGIN SOLUTION
    G = nx.DiGraph()
    s.insert(0,0)
    for i,si in enumerate(s):
        G.add_node("%d:%d"%(i,si))
    masses = list(mass_a.keys())
    for i,s1 in enumerate(s):
        for j,s2 in enumerate(s):
            if j-i in masses:
                G.add_edge("%d:%d"%(i,s1),"%d:%d"%(j,s2),label="/".join(mass_a[j-i]))
    if debug:
        show(G)
     
    peptides = []
    scores = []
    for path in nx.all_simple_paths(G,"%d:%d"%(0,0),"%d:%d"%(len(s)-1,0)):
        peptide = ""
        pscore = 0
        for i in range(len(path)-1):
            ix,score = path[i].split(":")
            pscore += int(score)
            a = G.get_edge_data(path[i],path[i+1])['label']
            peptide += a
        peptides.append(peptide)
        scores.append(pscore)

    ix = np.argmax(scores)
    peptide = peptides[ix]
    ## END SOLUTION
    return peptide





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