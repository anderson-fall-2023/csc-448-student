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

<!-- #region slideshow={"slide_type": "slide"} -->
# Topic 5 - Which animal gave us SARS?
## Secondary Title: Evolutionary Tree Construction

Motivation and some exercises are variations on those available in Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Learning Outcomes
* Understand evolutionary trees and their uses in biology
* Apply, analyze, and evaluate evolutionary tree algorithms
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Ice breaker

What's your favorite hobby picked up during the pandemic?
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
# Fastest Outbreak?
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## SARS crosses Pacific Ocean within a week
* Feb 21, 2003 - Chinese doctor Liu Jianlu flys to Hong Kong to attend a wedding
* Two weeks later Dr. Liu Jianlu dies but tells doctors that he recently treated sick patients in Guangdong Province
* Feb 23, 2003 - Man staying across the hall from Liu Jianlu travels to Hanoi and dies after infecting 80 people
* Feb 26, 2003 - Woman traveled home to Toronto bringing the disease and initiating an outbreak there.

* It took 5 years for the Black Death to travel from Constantinople to Kiev
* It took HIV two decades to circle the global
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Reaction
The reaction was immense with Chinese officials threatening to execute infected patients who violated quarantine.

If international travel spread the disease then international collaboration would contain it.

Within weeks bioligists identified a virus that caused the epdemic and sequences the genome. The new disease earned the name Severe Acute Respiratory Syndrome or SARS

<a href="https://theconversation.com/the-mysterious-disappearance-of-the-first-sars-virus-and-why-we-need-a-vaccine-for-the-current-one-but-didnt-for-the-other-137583">Article on SARS and COVID19</a>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Evolution of SARS
* SARS belongs to a family of viruses called conoraviruses (named after latin word *corona* meaning crown)
* Before SARS the coronaviruses were thought to only cause minor problems like the common cold

<table>
    <tr>
        <td>
            <img src="http://bioinformaticsalgorithms.com/images/Evolution/coronavirus.png" width=200>
        </td>
        <td>
            <img src="http://bioinformaticsalgorithms.com/images/Evolution/eclipse.png" width=200>
        </td>
    </tr>
    </table>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Sequencing
* By fall 2003, researchers sequenced many strains from patients in various contries
* OK ... so we keep talking about sequencing...
<center>
<img src="https://www.genengnews.com/wp-content/uploads/2018/10/Aug23_GEN_Baker_NGSChallenges_GettyImages537618080_ktsimage_DNATestSangerSeq2004017011.jpg" width=400>
</center>
* <a href="https://www.youtube.com/watch?v=zBPKj0mMcDg&ab_channel=ThermoFisherScientific">But what are some modern ways to do sequencing?</a>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Questions
* How did SARS-CoV cross the species barrier to humans? 
* When and where did it happen? How did SARS spread around the world, and who infected whom?

Each of these questions about SARS is ultimately related to the problem of constructing evolutionary trees (also known as phylogenies). Here is an evolutionary tree for HIV, but what algorithms do we use for this?

<img src="http://bioinformaticsalgorithms.com/images/Evolution/HIV_phylogeny.png" width=500>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
Our weekly injection of biology from a biologist :)

Let's watch this video and try to answer:
1. What is cladistics?
2. What is parsimony?
3. Do we compare the whole genome or part of it? Why or why not?
4. What do we use ancestoral genes?
5. Is it one size algorithm fits all?

<a href="https://calpoly.zoom.us/rec/share/hY2IiOUQZlOwwt54utQtXUuwdIcWsERw_ylxc48HcNT-BZ4VhxPwWCueNA1hLAzB.1Fhte5u1Cs-Kn4-2">Dr. Jean Davidson Video</a> Passcode: q$20al&M
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
# Distance Matrices to Evolutionary Trees
* Scientists started sequencing cornonavirus from various species to determine which is most similar to SARS
* Comparing (multiple alignment) of entire viral genomes is tricky because viral genes are often rearranged, inserted, and deleted
* Scientists focused on one of six genes in SARS-CoV
* Gene that encodes Spike protein
    * Identifies and binds to receptor site on host's cell membrane
    * Spike protein is 1,255 amino acids long and rather weak similarity with Spike proteins in other coronaviruses
    * Even subtle similarities turned out to be ssufficient for constructing a multiple alignment (comparison) across coronaviruses
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Distance matrix
Consider the DNA sequences shown below for 4 different species. 
<img src="http://bioinformaticsalgorithms.com/images/Evolution/mammal_alignment_distance_matrix.png" width=400>
The above example defined a distance matrix of +1 for every mismatched position. In general, $D$ must satisfy three properties. It must be symmetric (for all $i$ and $j$, $D_{i,j}$ = $D_{j,i}$), non-negative (for all $i$ and $j$, $D_{i,j}$ $\ge$ 0) and satisfy the triangle inequality (for all $i$, $j$, and $k$, $D_{i,j} + D_{j,k} \ge D_{i,k}$ ).
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
<img src="http://bioinformaticsalgorithms.com/images/Evolution/tree_of_life.png" width=700>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Rooted trees
<img src="http://bioinformaticsalgorithms.com/images/Evolution/rooted_tree_time.png" width=700>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## What are we aiming for?
We say that a weighted unooted tree $T$ fits a distance matrix $D$ if $d_{i,j}=D_{i,j}$ for every pair of leaves $i$ and $j$. Example:

<img src="http://bioinformaticsalgorithms.com/images/Evolution/additive_distance_matrix.png" width=300>

<img src="http://bioinformaticsalgorithms.com/images/Evolution/simple_tree_fitting_additive_matrix.png" width=300>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
## Technical Detour (networkx and pandas): 

``networkx`` is a very useful Python package. It provides a great visualization for many different applications. In your projects you can use other evolutionary tree tools.
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
%matplotlib inline 

import networkx as nx

G = nx.Graph()

G.add_edge('v1', 'v5', weight=11)
G.add_edge('v2', 'v5', weight=2)
G.add_edge('v5', 'v6', weight=4)
G.add_edge('v6', 'v3', weight=6)
G.add_edge('v6', 'v4', weight=7)
```

```python slideshow={"slide_type": "subslide"}
import copy
import pandas as pd

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
```

```python slideshow={"slide_type": "subslide"}
show(G)
```

```python slideshow={"slide_type": "subslide"}
show_adj(G)
```

<!-- #region slideshow={"slide_type": "slide"} -->
## Into the algorithms
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Let's talk about a given graph and neighboring leaves
<br>
<center>
${\color{red}{d_{k,m}}} = \dfrac{({\color{purple}{d_{i,m}}} + {\color{red}{d_{k,m}}}) + ({\color{blue}{d_{j,m}}} + {\color{red}{d_{k,m}}}) - ({\color{purple}{d_{i,m}}} + {\color{blue}{d_{j,m}}})}{2} = \dfrac{d_{i,k} +d_{j,k} - d_{i,j}}{2}$

<img src="http://bioinformaticsalgorithms.com/images/Evolution/neighboring_leaves_equality.png" width=500>
</center>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 1** Compute the distances between leaves in a weighted tree

Input: A weighted tree defined by the package networkx

Output: $n \times n$ matrix ($d_{i,j}$), where $d_{i,j}$ is the length of the path between leaves $i$ and $j$.

Learning objectives:
1. Refresh memory of graph traversal (path finding)
2. Understand the difference between $d_{i,j}$ and $D_{i,j}$.
3. Gain exposure and work with networkx python package.
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
import pandas as pd

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
                ## BEGIN SOLUTION
                path = list(nx.all_simple_paths(G,nodei,nodej))[0] # get the first simple path
                ## END SOLUTION
                a = path[0]
                for b in path[1:]:
                    dij += d[a,b]
                    a = b
                d[nodei,nodej] = dij
                d[nodej,nodei] = dij
    d = pd.DataFrame(d.values(),index=d.keys(),columns=['d']).unstack()
    d.columns = [n for l,n in d.columns]
    return d

compute_d(G)
```

<!-- #region slideshow={"slide_type": "subslide"} -->
## What if you need to find the graph?
<img src="http://bioinformaticsalgorithms.com/images/Evolution/additive_phylogeny.png" width=450>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 2** Implement limb length algorightm.

Input: An addititve distance matrix $D$ and a node $j$

Output: The length of the limb connect leaf $j$ to its parent in $Tree(D)$.

Learning outcomes:
1. Understanding why this function is needed when we just computed the paths weights previously.
2. Understanding the Limb Length Theorem in Chapter 7.
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
import pandas as pd
import numpy as np

def limb(D,j):
    min_length = np.Inf
    nodes = D.drop(j).index
    for ix,i in enumerate(nodes):
        for kx in range(ix+1,len(nodes)):
            k = nodes[kx]
            ## BEGIN SOLUTION
            length = (D.loc[i,j] + D.loc[j,k] - D.loc[i,k])/2
            if length < min_length:
                min_length = length
            ## END SOLUTION
    return min_length

names = ["v1","v2","v3","v4"]
D = pd.DataFrame([[0,13,21,22],[13,0,12,13],[21,12,0,13],[22,13,13,0]],index=names,columns=names)
print(D)
limb(D,"v4")

```

<!-- #region slideshow={"slide_type": "subslide"} -->
## Building up additive phylogeny
One piece at a time... First thing will be to find out where we could insert a noode.

Learning outcomes:
1. Using recursion and by extension debugging recursion
2. Understand additive versus non-additive $D$
3. Constructing our first evolutionary tree
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 3a** Implement a portion of ``AdditivePhylogeny`` algorithm from Chapter 7.

Input: Distance matrix $D$ and node name $n$.

Output: Return the node names $i,k$ that satisfy $D_{i,k} = D_{i,n} + D_{n,k}$. In other words, where can you insert $n$ back in.
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
Dorig = copy.copy(D)
```

```python slideshow={"slide_type": "subslide"}
def find(D,n):
    nodes = D.drop(n).index
    for ix,i in enumerate(nodes):
        for kx in range(ix+1,len(nodes)):
            ## BEGIN SOLUTION
            k = nodes[kx]
            if D.loc[i,k] == D.loc[i,n] + D.loc[n,k]:
                return i,k
            ## END SOLUTION
            # Your solution here
            pass
    return None,None

D = copy.copy(Dorig)
print("Starting D")
print(D)
limbLength = limb(D,D.index[-1]) # our algorithm will choose the last node
n = D.index[-1]
print("Node to remove:",n)
Dtrimmed = D.drop(n).drop(n,axis=1)
for j in Dtrimmed.index:
    D.loc[j,n] = D.loc[j,n] - limbLength
    D.loc[n,j] = D.loc[j,n]
print("New D")
print(D)
print(n)
find(D,n)
```

**Exercise 3b** Implement a portion of ``AdditivePhylogeny`` algorithm from Chapter 7.

Input: Distance matrix $D$ of size $2 \times 2$.

Output: Return a networkx graph with the correct weight.

```python
def base_case(D):
    T = nx.Graph()
    ## YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    T.add_edge(D.index[0], D.index[1], weight=D.loc[D.index[0],D.index[1]])
    ## END SOLUTION
    return T

base_G = base_case(D.iloc[:2,:].iloc[:,:2])
show(base_G)
```

```python slideshow={"slide_type": "subslide"}
def additive_phylogeny(D,new_number):
    D = D.copy()
    if len(D) == 2:
        return base_case(D) # Implemented correctly above
    n = D.index[-1]
    limbLength = limb(D,n) # our algorithm will choose the last node
    Dtrimmed = D.drop(n).drop(n,axis=1)
    for j in Dtrimmed.index:
        D.loc[j,n] = D.loc[j,n] - limbLength
        D.loc[n,j] = D.loc[j,n]

    Dtrimmed = D.drop(n).drop(n,axis=1)
    T = additive_phylogeny(Dtrimmed,new_number+1)

    i,k = find(D,n) # Implemented correctly above
    print(n,i,k)
    #if D.loc[j,n] < D.loc[i,n]:
    #    i,k = k,i
    v = "v%s"%new_number
    ## Your solution here
    # This is definitely the most complicated thing conceptually
    # You'll need to add edges and remove edges (T.add_edge and T.remove_edge)
    ## BEGIN SOLUTION
    path = list(nx.all_simple_paths(T,i,k))[0]
    a = path[0]
    cost = 0
    A = show_adj(T)
    for b in path[1:]:
        segment_cost = A.loc[a,b]
            
        cost += segment_cost
        
        if D.loc[i,n] == cost:
            raise Exception("Not implemented")
            break
        elif D.loc[i,n] < cost: # along an edge
            cost = cost - segment_cost
            break
                        
        a = b

    T.add_edge(v,n,weight=limbLength)
    T.add_edge(a,v,weight=D.loc[i,n] - cost)
    T.add_edge(v,b,weight=segment_cost - (D.loc[i,n] - cost))
    T.remove_edge(a,b)
    ## END SOLUTION
    return T

D = copy.copy(Dorig)
G2 = additive_phylogeny(D,len(D)+1)
show(G2)
```

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 4.** Run your new algorithm on SARS data derived from multiple alignment of Spike proteins.
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
import os.path

file = None
locations = ['../data/coronavirus_distance_matrix_additive.txt','../csc-448-student/data/coronavirus_distance_matrix_additive.txt']
for f in locations:
    if os.path.isfile(f):
        file = f
        break
print('Opening',file)
D_sars = pd.read_csv(file,index_col=0)
D_sars
```

```python slideshow={"slide_type": "subslide"}
from pylab import rcParams
rcParams['figure.figsize'] = 10, 5

G3 = additive_phylogeny(D_sars,len(D_sars)+1)
show(G3)
```

```python slideshow={"slide_type": "skip"}
show_adj(G3)
```

<!-- #region slideshow={"slide_type": "skip"} -->
**Helper function to check things out**
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
def compute_path_cost(T,i,k):
    cost = 0
    try:
        path = list(nx.all_simple_paths(T,i,k))[0]
    except:
        return -1
    a = path[0]
    cost = 0
    A = show_adj(G3)
    for b in path[1:]:
        cost += A.loc[a,b]
        a = b
    return cost
```

```python slideshow={"slide_type": "skip"}
compute_path_cost(G3,'Human','Turkey')
```

```python slideshow={"slide_type": "fragment"}
compute_path_cost(G3,'Human','Turkey') == D_sars.loc['Human','Turkey']
```

```python slideshow={"slide_type": "skip"}
# Don't forget to push!
```

```python

```

```python

```
