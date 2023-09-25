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
    ## BEGIN SOLUTION
    cnt = 1
    for pattern in patterns:
        current = "root"
        for i in range(len(pattern)):
            c = pattern[i]
            n = find_edge(G,current,c)
            if n is None:
                G.add_edge(current,str(cnt),label=c)
                current = str(cnt)
                cnt += 1
            else:
                current = n
            
    ## END SOLUTION
    return G

def trie_matching(text,trie):
    positions = []
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    i = 0
    while i < len(text):
        r = prefix_trie_matching(text[i:],trie)
        if r is not None:
            positions.append(i)
        i += 1
    ## END SOLUTION
    return positions



def prefix_trie_matching(text,trie):
    # Your solution here
    ## BEGIN SOLUTION
    symbol = ""
    v = "root"
    i = 0
    while i < len(text):            
        if len(list(trie.neighbors(v))) == 0:
            return symbol
        else:
            w = find_edge(trie,v,text[i])
            if w is None:
                return None
            symbol += text[i]
            i += 1
            v = w
    ## END SOLUTION
    return None

def suffix_trie(text):
    G = nx.DiGraph()
    G.add_node('root')
    # Your solution here
    ## BEGIN SOLUTION
    i = 0
    cnt = 1
    while i < len(text):
        current = "root"
        suffix = text[i:]
        for c in suffix:
            if c == "$":
                G.add_edge(current,"[%d]"%i,label=c)
            else:
                n = find_edge(G,current,c)
                if n is None:
                    G.add_edge(current,str(cnt),label=c)
                    current = str(cnt)
                    cnt += 1
                else:
                    current = n  
        i+=1
    ## END SOLUTION
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
    ## BEGIN SOLUTION
    i = 0
    cnt = 1
    while i < len(text):
        current = "root"
        suffix = text[i:]
        for j,c in enumerate(suffix):
            if c == "$":
                G.add_edge(current,"[%d]"%i,label="%s,%d"%(c,i+j))
                leaf_nodes.append("[%d]"%i)
            else:
                n = modified_find_edge(G,current,"%s,%d"%(c,i+j))
                if n is None:
                    G.add_edge(current,str(cnt),label="%s,%d"%(c,i+j))
                    current = str(cnt)
                    cnt += 1
                else:
                    current = n  
        i+=1
    ## END SOLUTION
    return G,leaf_nodes

def suffix_tree_construction(text):
    trie,leaf_nodes = modified_suffix_trie(text)
    ## BEGIN SOLUTION
    to_remove = []
    for leaf_node in leaf_nodes:
        current = "root"
        full_paths = list(nx.all_simple_paths(trie,current,leaf_node))
        for full_path in full_paths: # this loop isn't necessary as there is only a single simple path
            if len(full_path) == 2:
                data = trie.get_edge_data(full_path[-2],full_path[-1])
                cw,i = data['label'].split(",")
                trie.remove_node(full_path[-1])
                trie.add_edge(full_path[-2],full_path[-1],label=cw)
            else:
                while True:
                    path = [full_path[0]]
                    for k,n in enumerate(full_path[1:]):
                        path += [n]
                        if len(list(trie.neighbors(n))) > 1:
                            break
                    if len(path) > 2:
                        data = trie.get_edge_data(path[0],path[1])
                        cw,i = data['label'].split(",")
                        i = int(i)
                        data = trie.get_edge_data(path[-2],path[-1])
                        cv,j = data['label'].split(",")
                        j = int(j)
                        trie.add_edge(path[0],path[-1],label=text[i:j+1])
                        for n in path[1:-1]:
                            to_remove.append(n)
                    else:
                        if path[-1] not in leaf_nodes:
                            data = trie.get_edge_data(path[-2],path[-1])
                            if "," in data["label"]:
                                cw,i = data['label'].split(",")
                                nx.set_edge_attributes(trie, {(path[-2],path[-1]): {'label':cw}})
                    if k+1 < len(full_path) and full_path[k+1] not in leaf_nodes:
                        current = full_path[k+1]
                        full_path = full_path[k+1:]
                    else:
                        break
    for n in set(to_remove):
        trie.remove_node(n)
    ## END SOLUTION
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
            
            