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
# Assignment 5 - Which animal gave us SARS?
## Secondary Title: Evolutionary Tree Construction
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
%matplotlib inline
%load_ext autoreload
%autoreload 2

## BEGIN SOLUTION
import joblib
answers = {}
## END SOLUTION

import pandas as pd
import numpy as np

import Assignment5_helper 

from pathlib import Path
home = str(Path.home()) # all other paths are relative to this path. 
# This is not relevant to most people because I recommended you use my server, but
# change home to where you are storing everything. Again. Not recommended.
```

<!-- #region slideshow={"slide_type": "subslide"} -->
## Learning Outcomes
* Understand evolutionary trees and their uses in biology
* Apply, analyze, and evaluate evolutionary tree algorithms
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
d=Assignment5_helper.compute_d(Assignment5_helper.G)
## BEGIN SOLUTION
answers["exercise_1"] = d
## END SOLUTION
d
```

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 2** Implement limb length algorightm described in Chapter 7.

Input: An addititve distance matrix $D$ and a node $j$

Output: The length of the limb connect leaf $j$ to its parent in $Tree(D)$.

Learning outcomes:
1. Understanding why this function is needed when we just computed the paths weights previously.
2. Understanding the Limb Length Theorem in Chapter 7.
<!-- #endregion -->

```python
Assignment5_helper.D
```

```python slideshow={"slide_type": "subslide"}
length = Assignment5_helper.limb(Assignment5_helper.D,"v4")
## BEGIN SOLUTION
answers["exercise_2"] = length
## END SOLUTION
length
```

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 3** Implement a portion of ``AdditivePhylogeny`` algorithm from Chapter 7.

Input: Distance matrix $D$ and node name $n$.

Output: Return the node names $i,k$ that satisfy $D_{i,k} = D_{i,n} + D_{n,k}$. In other words, where can you insert $n$ back in.
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
D=Assignment5_helper.D.copy()
print("Starting D")
print(D)
limbLength = Assignment5_helper.limb(D,D.index[-1]) # our algorithm will choose the last node
n = D.index[-1]
print("Node to remove:",n)
Dtrimmed = D.drop(n).drop(n,axis=1)
for j in Dtrimmed.index:
    D.loc[j,n] = D.loc[j,n] - limbLength
    D.loc[n,j] = D.loc[j,n]
print("New D")
print(D)
i,k = Assignment5_helper.find(D,"v4")
## BEGIN SOLUTION
answers["exercise_3"] = i,k
## END SOLUTION
i,k
```

**Exercise 4** Implement a portion of ``AdditivePhylogeny`` algorithm from Chapter 7. The base case when there are only two nodes.

Input: Distance matrix $D$ of size $2 \times 2$.

Output: Return a networkx graph with the correct weight.

```python
base_G = Assignment5_helper.base_case(Assignment5_helper.D.iloc[:2,:].iloc[:,:2])
## BEGIN SOLUTION
import networkx as nx
answers["exercise_4"] = nx.adjacency_matrix(base_G).todense()
## END SOLUTION
Assignment5_helper.show(base_G)
```

**Exercise 5:** We are ready to put everything together! Implement the full additive phylogeny algorithm from Chapter 7. 

```python slideshow={"slide_type": "subslide"}
G2 = Assignment5_helper.additive_phylogeny(Assignment5_helper.D,len(D)+1)
## BEGIN SOLUTION
answers["exercise_5"] = nx.adjacency_matrix(G2).todense()
## END SOLUTION
Assignment5_helper.show(G2)
```

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 6** Run your new algorithm on SARS data derived from multiple alignment of Spike proteins.
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
import os.path

file = f'{home}/csc-448-student/data/coronavirus_distance_matrix_additive.txt'
print('Opening',file)
D_sars = pd.read_csv(file,index_col=0)
D_sars
```

```python slideshow={"slide_type": "subslide"}
from pylab import rcParams
rcParams['figure.figsize'] = 10, 5

G3 = Assignment5_helper.additive_phylogeny(D_sars,len(D_sars)+1)
## BEGIN SOLUTION
answers["exercise_6"] = nx.adjacency_matrix(G3).todense()
## END SOLUTION
Assignment5_helper.show(G3)
```

```python slideshow={"slide_type": "skip"}
Assignment5_helper.show_adj(G3)
```

<!-- #region slideshow={"slide_type": "skip"} -->
**Helper function to help you check your algorithms**
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
Assignment5_helper.compute_path_cost(G3,'Human','Turkey')
```

```python slideshow={"slide_type": "fragment"}
Assignment5_helper.compute_path_cost(G3,'Human','Turkey') == D_sars.loc['Human','Turkey']
```

```python slideshow={"slide_type": "skip"}
## BEGIN SOLUTION
joblib.dump(answers,"../tests/answers_Assignment5.joblib");
## END SOLUTION
# Don't forget to push!
```

```python

```
