import numpy as np
import pandas as pd

import networkx as nx
import pandas as pd
import copy

import matplotlib.pyplot as plt

display_available = True
try:
    from IPython.display import Image
except:
    display_available = False
try:
    import pygraphviz
    graphviz_installed = True # Set this to False if you don't have graphviz
except:
    graphviz_installed = False

def draw(A):
    return Image(A.draw(format='png', prog='dot'))

patterns1 = ['ATAGA','ATC','GAT']
patterns2 = ['ananas','and','antenna','banana','bandana','nab','nana','pan']

# Inputs: G - networkx graph, current - node name, c - character on edge
# Output: a neighbor of current that is reached by an edge that has label==c; otherwise None
def find_edge(G,current,c): 
    for n in G.neighbors(current):
        if n is None:
            return None
        data = G.get_edge_data(current,n)
        if data['label'] == c:
            return n
    return None

def trie_construction(patterns):
    G = nx.DiGraph()
    G.add_node('root')
    # Your solution here
    return G

def trie_matching(text,trie):
    positions = []
    # YOUR SOLUTION HERE
    return positions



def prefix_trie_matching(text,trie):
    # Your solution here
    return None

def suffix_trie(text):
    G = nx.DiGraph()
    G.add_node('root')
    # Your solution here
    return G

# Inputs: G - networkx graph, current - node name, c - character on edge
# Output: a neighbor of current that is reached by an edge that has label==c; otherwise None
def modified_find_edge(G,current,c):
    cv,j = c.split(",")
    j = int(j)
    for n in G.neighbors(current):
        if n is None:
            return None
        data = G.get_edge_data(current,n)
        cw,i = data['label'].split(",")
        i = int(i)
        if cw == cv and j > i:
            return n
    return None

def modified_suffix_trie(text):
    G = nx.DiGraph()
    G.add_node('root')
    leaf_nodes = []
    # Your solution here
    return G,leaf_nodes

def suffix_tree_construction(text):
    trie,leaf_nodes = modified_suffix_trie(text)
    return trie

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
    if graphviz_installed:
        # same layout using matplotlib with no labels
        pos = nx.drawing.nx_agraph.graphviz_layout(G, prog='dot')
        #print(edge_labels)
        # Modify node fillcolor and edge color.
        #D.node_attr.update(color='blue', style='filled', fillcolor='yellow')
        #D.edge_attr.update(color='blue', arrowsize=1)
        A = nx.nx_agraph.to_agraph(G)
        # draw it in the notebook
        if display_available:
            display(draw(A))
        else:
            print(A)
    else:
        if display_available:
            display(to_adj(G))
        else:
            print(to_adj(G))
            
            