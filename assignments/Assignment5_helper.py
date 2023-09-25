import copy
import pandas as pd
import numpy as np

import networkx as nx

G = nx.Graph()

G.add_edge('v1', 'v5', weight=11)
G.add_edge('v2', 'v5', weight=2)
G.add_edge('v5', 'v6', weight=4)
G.add_edge('v6', 'v3', weight=6)
G.add_edge('v6', 'v4', weight=7)

names = ["v1","v2","v3","v4"]
D = pd.DataFrame([[0,13,21,22],[13,0,12,13],[21,12,0,13],[22,13,13,0]],index=names,columns=names)

def show(T):
    T = copy.deepcopy(T)
    labels = nx.get_edge_attributes(T,'weight')
    max_value = 0
    for n1,n2 in T.edges():
        if T[n1][n2]['weight'] > max_value:
            max_value = T[n1][n2]['weight']
    for n1,n2 in T.edges():
        T[n1][n2]['weight']=max_value - T[n1][n2]['weight'] + 3
    pos=nx.spring_layout(T)
    nx.draw(T,pos,with_labels=True)
    nx.draw_networkx_edge_labels(T,pos,edge_labels=labels);
    
def show_adj(T):
    if len(T.nodes()) == 0:
        return pd.DataFrame()
    return pd.DataFrame(nx.adjacency_matrix(T).todense(),index=T.nodes(),columns=T.nodes())

def compute_d(G):
    d = {}
    for nodei in G.nodes():
        for nodej in G.nodes():
            d[nodei,nodej] = 0
    # Fill in all adjacent values
    for nodei,nodej,data in G.edges(data=True):
        d[nodei,nodej] = data['weight']
        d[nodej,nodei] = d[nodei,nodej]
    for nodei in G.nodes():
        for nodej in G.nodes():
            if d[nodei,nodej] == 0 and nodei != nodej:
                dij = 0
                ## YOUR SOLUTION HERE
                # networkx has a function to compute simple paths. You can uncomment out the line
                # below in order to see this function working. I'll be reviewing your solutions though
                # and if you don't write this from scratch, then I will consider this an attempt to
                # circumvent the autograder
                #path = list(nx.all_simple_paths(G,nodei,nodej))[0] # get the first simple path
                path = ["v1","v2"] #remove this when you are ready
                a = path[0]
                for b in path[1:]:
                    dij += d[a,b]
                    a = b
                d[nodei,nodej] = dij
                d[nodej,nodei] = dij
    d = pd.DataFrame(d.values(),index=d.keys(),columns=['d']).unstack()
    d.columns = [n for l,n in d.columns]
    return d

def limb(D,j):
    min_length = np.Inf
    nodes = D.drop(j).index
    for ix,i in enumerate(nodes):
        for kx in range(ix+1,len(nodes)):
            k = nodes[kx]
    return min_length

def find(D,n):
    D = copy.copy(D)
    nodes = D.drop(n).index
    for ix,i in enumerate(nodes):
        for kx in range(ix+1,len(nodes)):
            # Your solution here
            pass
    return None,None

def base_case(D):
    T = nx.Graph()
    ## YOUR SOLUTION HERE
    return T

def additive_phylogeny(D,new_number):
    D = copy.deepcopy(D)
    if len(D) == 2:
        return base_case(D) # Implemented correctly above
    n = D.index[-1]
    limbLength = limb(D,n) # our algorithm will choose the last node
    print(limbLength)
    Dtrimmed = D.drop(n).drop(n,axis=1)
    for j in Dtrimmed.index:
        D.loc[j,n] = D.loc[j,n] - limbLength
        D.loc[n,j] = D.loc[j,n]

    Dtrimmed = D.drop(n).drop(n,axis=1)
    T = additive_phylogeny(Dtrimmed,new_number+1)

    i,k = find(D,n) # Implemented correctly above
    #if D.loc[j,n] < D.loc[i,n]:
    #    i,k = k,i
    v = "v%s"%new_number
    ## Your solution here
    # This is definitely the most complicated thing conceptually
    # You'll need to add edges and remove edges (T.add_edge and T.remove_edge)
    return T

def compute_path_cost(T,i,k):
    cost = 0
    try:
        path = list(nx.all_simple_paths(T,i,k))[0]
    except:
        return -1
    a = path[0]
    cost = 0
    A = show_adj(T)
    for b in path[1:]:
        cost += A.loc[a,b]
        a = b
    return cost