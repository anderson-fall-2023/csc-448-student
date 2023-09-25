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
# Topic 3 - Genome Assembly
## Graph Algorithms
Material and embedded exercises.

Motivation and some exercises are variations on those available in Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Learning Outcomes
* Understand genome assembly bioinformatics motivation
* Apply, analyze, and evaluate genome sequence alignment graph based algorithms
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
## BEGIN SOLUTION
import joblib
answers = {}
## END SOLUTION
%load_ext autoreload
%autoreload 2
```

<!-- #region slideshow={"slide_type": "slide"} -->
# History and motivation
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Why do we sequence different species?

<table>
    <tr>
                <td><img src="https://www.akcchf.org/assets/images/Dog-Genome-Nature.jpg"></td>
                <td><img src="http://www.scienceagainstevolution.org/images/v10i1g1.jpg" width=300></td>
                <td><img src="https://www.ncbi.nlm.nih.gov/genome/guide/bee/cover_nature.jpg"></td>
    </tr>
    </table>
Many applications in medicine, agriculture, biotechnology, etc begin with a sequenced and annotated genome.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Brief History of Genome Sequencing
* 1977 - Walter Gilbert and Frederick Sanger develop independent DNA sequencing methods.
* 1980 - They share the Nobel Prize.
* Sequencing methods were too expensive ($3 billion to sequence the human genome)
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## The Race to Sequence the Human Genome
* 1990 - The public Human Genome Project, headed by Francis Collins, aims to sequence the human genome by 2005.
* 1997 - Craig Venter founds Celera Genomics, a private firm with the same goal :)
* Genome become available 5 years early in 2000
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
<img src="https://thumbs-prod.si-cdn.com/CZYSL1uFzN85JMfOk1XBKFSA9uI=/fit-in/1600x0/https://public-media.si-cdn.com/filer/f5/45/f545637b-13e4-42a8-927b-a1f59d326e0b/biogenome.jpg">
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
# Exploding Newspapers
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Consider a silly example
<img src="http://bioinformaticsalgorithms.com/images/Assembly/newspaper_blowup.png">
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "notes"} -->
"Imagine that we stack a hundred copies of the June 27, 2000 edition of the New York Times on a pile of dynamite, and then we light the fuse. We ask you to further suspend your disbelief and assume that the newspapers are not all incinerated but instead explode cartoonishly into smoldering pieces of confetti. How could we use the tiny snippets of newspaper to figure out what the news was on June 27, 2000? We will call this crazy conundrum the Newspaper Problem." - Bioinformatics Chapter 3
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
<img src="http://bioinformaticsalgorithms.com/images/Assembly/overlapping_newspaper.png" width=400>

Use overlapping shreds of paper to figure out the news.

"Fine, you ask, but what do exploding newspapers have to do with biology? Determining the order of nucleotides in a genome, or genome sequencing, presents a fundamental task in bioinformatics. Genomes vary in length; your own genome is roughly 3 billion nucleotides long, whereas the genome of Amoeba dubia, an amorphous unicellular organism, is approximately 200 times longer! This unicellular organism competes with the rare Japanese flower Paris japonica for the title of species with the longest genome." - Bioinformatics Algorithms, Chapter 3.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Traditional Method for Genome Sequencing

<img src="http://bioinformaticsalgorithms.com/images/Assembly/sequencing_overview.png">
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## What makes genome sequencing difficult?
* Modern sequencing machines cannot read an entire genome one nucleotide at a time from beginning to end (like we read a book)
* They can only shred the genome and generate short  reads (though we can now mix longer reads with these shorter reads)
* The genome assembly is not the same as a jigsaw puzzle: we must use overlapping reads to reconstruct the genome, a  giant overlap puzzle!
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
## Into the language of computer scientists: String Reconstruction Problem
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 1** String composition problem: Generate the k-mer composition of a string.
    
Input: An integer $k$ and a string $text$.

Output: The lexicographic ordered collection of all k-mer substrings in text (including repeated k-mers) ``composition(k,text)``, where the k-mers are arranged in lexicographic order.

Example: ``composition(3,TATGGGGTGC)`` => [ATG, GGG, GGG, GGT, GTG, TAT, TGC, TGG]. Notice how there are two entries for GGG, so we keep duplicates.

Learning objective: Build up the skills and vocabulary we need to solve the reconstruction problem.
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
def composition(k,text):
    patterns = []
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    for i in range(0,len(text)-k+1):
        patterns.append(text[i:(i+k)])
    patterns.sort()
    ## END SOLUTION
    return patterns

## BEGIN SOLUTION
answers["answer_exercise_1"] = tuple(composition(3,"TATGGGGTGC"))
## END SOLUTION
composition(3,"TATGGGGTGC")
```

<!-- #region slideshow={"slide_type": "subslide"} -->
Well that's great and pretty straightfoward, but we need to solve the inverse problem...

**String Reconstruction Problem:** Reconstruct a string from its k-mer composition.

Input: An integer $k$ and a collection ``patterns`` of $k$-mers

Output: A string ``text`` with $k$-mer composition equal to ``patterns`` if such a string exists.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
### Easy string reconstruction problem first
Consider ``patterns`` = [AAT ATG GTT TAA TGT]

Most straightfoward way to solve this is to "connect" a pair of $k$-mers if they overlap in $k$-1 symbols. Where should we start? Let's start with ``TAA`` because there is no 3-mer ending in ``TA``. The next $k$-mer has to be ``AAT`` because it is the only one that satisfies. We can continue with this pattern easily for this example:
<pre>
TAA    
 AAT   
  ATG  
   TGT 
    GTT
TAATGTT
</pre>
Nice! We did it, so time to move on to another chapter...
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
### A harder string reconstruction problem

AAT  ATG  ATG  ATG  CAT  CCA  GAT  GCC  GGA  GGG  GTT  TAA  TGC  TGG  TGT

Let's start with ``TAA`` again:
<pre>
TAA  
 AAT 
  ATG
</pre>
So far so good! Let's try to continue. ``ATG`` can be extended by either ``TGC`` or ``TGG`` or ``TGT``. How do we choose? In other words, how good are you at chess because we need to be able to look ahead. For now, let's choose ``TGT``.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
After we choose ``TGT`` our only choose is ``GTT`` so we get:
<pre>
TAA   
 AAT  
  ATG 
   TGT
    GTT
</pre>
But we are now stuck! There is nothing that matches :(

If you are good at looking ahead you could have extended ``ATG`` by ``TGC`` instead.

<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
If we make that change, then we obtain the following assembly:
<pre>
TAA             
 AAT            
  ATG           
   TGC          
    GCC         
     CCA        
      CAT       
       ATG      
        TGG     
         GGA    
          GAT   
           ATG  
            TGT 
             GTT
TAATGCCATGGATGTT
</pre>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## So what was the problem?
Repeats complicate genome assembly!

Our previous problem came from the fact that ``ATG`` is repeated three times which causes us to have the three choices by which to extend ``ATG``. This doesn't pose too much of a problem for small examples, but consider having this happen when we have millions of reads!
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
# Bioinformatics from a biologist
## Genome Assembly

https://calpoly.zoom.us/rec/share/AMLFJgx84KoyoXpm4FXy8mViguXg21gvFp1QBf1_J0DXS3ownJriANDFZ_jqzxiI.DMc-nPbPVRe2u1E8
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
# Bioinformatics from a computer scientist
## Genome Assembly

https://calpoly.zoom.us/rec/share/PBxQoxqNsGMdR0l3jccQlSopZpuEyM_4hpU19FjbRJHO6VNmoI8BLgQKum1Z_HXi.fp9FrRd2-4fakeRV 

Passcode: Sg7.^4XU
<!-- #endregion -->

## Seven Bridge of Konigsberg

<img src="https://bitnine.net/wp-content/uploads/2016/07/20160713_3.jpg">

<!-- #region slideshow={"slide_type": "subslide"} -->
**A resident asks:** Can I set out from my home and cross each bridge once and return home at the end?
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
Since you know a graph algorithm is coming, how would you represent the above as a graph?
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "fragment"} -->
Answer: colors as nodes or vertices and bridges as edges.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
In comes the father of graph theory: Leonhard Euler. 

Graphs have edges and vertices. 

This is not a graph today:

<img src="https://www.betterevaluation.org/sites/default/files/scatterplot2.gif" width=200>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
<img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fnbt.2023/MediaObjects/41587_2011_Article_BFnbt2023_Fig1_HTML.gif">
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**Hamiltonian Paths and Cycles:**

> Hamiltonian path is a path in an undirected or directed graph that visits each vertex exactly once. A Hamiltonian cycle is a Hamiltonian path that is a cycle. Determining whether such paths and cycles exist in graphs is the Hamiltonian path problem, which is NP-complete. - https://en.wikipedia.org/wiki/Hamiltonian_path
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**How does this relate to bioinformatics?**

For years bioinformaticians came up with heuristic/approximate algorithms for finding Hamiltonian paths and cycles in order to perform genome assembly. But it is a NP-complete problem.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**After moderate but hard won success the field realized that there was another way...**

And that is where are Bridges of Konigsberg comes into play.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
A **Eulerian cycle** is a cycle that visits each edge once but does not matter how many times a vertices or edge is visited.

Big difference is that you can visit a vertices more than once. This seemingly small difference has a dramatic affect on the computability of genome assembly.

The way Euler came up with an algorithm to determine if a Eulerian cycle exists, he came up with a method for finding them. 
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
## String Reconstruction: A reverse approach/order from the book
The book builds up from walking overlap graphs to an *Algorithm for Finding Eulerian Cycles* (wait what are those words!).

We are going to work our way backwards because it is a good thing to have different complementary approaches. Please please read the book and listen to the textbook online lectures.

Our big goal is the following...
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**Solve the string reconstruction problem** using the algorithm for finding Eulerian Cycles. Our approach is as follows:
<pre>
def string_reconstruction(patterns):
    dB = de_bruijn(patterns)
    path = eulerian_path(dB)
    text = path_to_genome(path)
    return text
</pre>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
## Constructing de Bruijn graphs from $k$-mer composition
Given a collection of $k$-mers Patterns, the nodes of de_bruijn(k,patterns) are simply all unique (k−1)-mers occurring as a prefix or suffix in Patterns. For example, say we are given the following collection of 3-mers:
<pre>
AAT   ATG   ATG   ATG    CAT   CCA   GAT   GCC   GGA   GGG   GTT   TAA   TGC   TGG   TGT
</pre>
Then the set of eleven unique 2-mers is:
<pre>
AA   AT   CA   CC   GA   GC   GG   GT   TA   TG   TT
</pre>
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
For every $k$-mer in ``patterns`` we will connect its prefix node to its suffix node by a directed edge in order to produce our final graph.

<img src="http://bioinformaticsalgorithms.com/images/Assembly/debruijn_graph_alternate_rendering.png">
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 2** Construct the de Bruijn graph from a set of k-mers.

Input: a collection of $k$-mers ``patterns``.

Output: NetworkX directed graph that corresponds to the adjacency list of the de Bruijn graph. 

Now that you know how to use ``networkx``... use it to make your life easier for representing your graph.
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
# Just code to help you out here
%matplotlib inline

import networkx as nx
import pandas as pd
import copy
from collections import Counter

def to_adj(T):
    try:
        return pd.DataFrame(nx.adjacency_matrix(T).todense(),index=T.nodes(),columns=T.nodes())
    except:
        print("Cannot convert to adjacency matrix")
    return None

def show(T):
    T = copy.deepcopy(T)
    width_dict = Counter(T.edges())
    edge_width = [ (u, v, {'width': value}) 
                  for ((u, v), value) in width_dict.items()]
    
    G_new = nx.DiGraph()
    G_new.add_edges_from(edge_width)
    pos=nx.kamada_kawai_layout(G_new)
    #pos=nx.spring_layout(G_new)
    nx.draw(G_new, pos)
    edge_labels=dict([((u,v,),d['width'])
                 for u,v,d in G_new.edges(data=True)])
    
    nx.draw(G_new,pos,with_labels=True)
    nx.draw_networkx_edges(G_new, pos=pos)
    nx.draw_networkx_edge_labels(G_new, pos, edge_labels=edge_labels,
                                 label_pos=0.55, font_size=10)
```

```python slideshow={"slide_type": "subslide"}
def de_bruijn(patterns):
    dB = nx.MultiDiGraph()
    # dB.add_edge("AA","AT") # sample edge in case you want to run the code without implementing your solution
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    k = len(patterns[0])
    for kmer in patterns:
        prefix = kmer[:(k-1)]
        suffix = kmer[1:]
        dB.add_edge(prefix,suffix)
    ## END SOLUTION
    return dB

dB = de_bruijn(["AAT","ATG","ATG","ATG","CAT","CCA","GAT","GCC","GGA","GGG","GTT","TAA","TGC","TGG","TGT"])
## BEGIN SOLUTION
answers["answer_exercise_2"] = nx.adjacency_matrix(dB).todense()
## END SOLUTION
show(dB)
```

```python slideshow={"slide_type": "subslide"}
to_adj(dB)
```

<!-- #region slideshow={"slide_type": "subslide"} -->
We now have step 1 completed. The next step is to find the Eulerian path, but we will start by talking about and finding Eulerian cycles. More undefined words! Don't forget. For a different bottom up approach (we are taking top down in this lab), then see the textbook. 

A Eulerian path is a path in a graph traversing each edge of a graph exactly once. A cycle that traverses each edge of a graph exactly once is called an Eulerian cycle.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
The book works in a lot more detail, but I will boil it down to what I view are the essentials.

If we assume an arbitrary directed graph is balanced (in-degree of each node is equal to out-degree) and it is strongly connected
* Then we can find the Eulerian cycle by just walking around! We just need to make sure we always leave out of a different edge.

<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 3** Given a Eulerian directed graph, find the Eulerian cycle.

Given: An Eulerian directed graph, in the form of a MultiDiGraph networkx object.

Return: An Eulerian cycle in this graph returned as a list.
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
def eulerian_cycle(G,start=None):
    # YOUR SOLUTION HERE
    cycle = None
    ## BEGIN SOLUTION
    edges = {}
    cnt = 0
    for u,v in G.edges():
        if u not in edges:
            edges[u] = []
        edges[u].append(v)
        cnt += 1
    if start is None:
        start = list(edges.keys())[0]
    cycle = [start]
    current = start
    ecnt = 0
    while ecnt < cnt:
        while True:
            neighbors = edges[current]
            if len(neighbors) > 0:
                break
            # we need to keep shifting until we get to a point that has unused portions
            cycle = [cycle[-2]] + cycle[0:-2] + [cycle[-2]]
            current = cycle[-1]
            
        current = neighbors.pop()
        cycle.append(current)
        ecnt += 1
    ## END SOLUTION
    return cycle
    
G = nx.MultiDiGraph()
G.add_edge(0,3)
G.add_edge(1,0)
G.add_edge(2,1)
G.add_edge(2,6)
G.add_edge(3,2)
G.add_edge(4,2)
G.add_edge(5,4)
G.add_edge(6,5)
G.add_edge(6,8)
G.add_edge(7,9)
G.add_edge(8,7)
G.add_edge(9,6)

## BEGIN SOLUTION
answers["answer_exercise_3"] = tuple(eulerian_cycle(G))
## END SOLUTION

show(G)
print(eulerian_cycle(G))

## BEGIN SOLUTION
answers["answer_exercise_3"] = tuple(eulerian_cycle(G,start=6))
## END SOLUTION

print(eulerian_cycle(G,start=6)) # Should also result in a shifted but equivalent answer. 
# If your code isn't working, please carefully read the section on the ant moving around the graph and getting stuck :)
```

<!-- #region slideshow={"slide_type": "subslide"} -->
## Eulerian path
Consider the following graph that does not contain a Eulerian cycle.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"}
G2 = nx.MultiDiGraph()
G2.add_edge(0,2);G2.add_edge(1,3);G2.add_edge(2,1);G2.add_edge(3,0);G2.add_edge(3,4);G2.add_edge(6,3);G2.add_edge(6,7);G2.add_edge(7,8);G2.add_edge(8,9);G2.add_edge(9,6)
show(G2)
```

<!-- #region slideshow={"slide_type": "subslide"} -->
Consider the following function that finds the in and out degree of every node.
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
def calc_in_out(G):
    in_deg = {}
    out_deg = {}
    for u,v in G.edges():
        if v not in in_deg:
            in_deg[v] = 0
        if u not in out_deg:
            out_deg[u] = 0
        in_deg[v] += 1
        out_deg[u] += 1
    in_out = pd.Series(in_deg,name="in").to_frame().join(pd.Series(out_deg,name="out").to_frame(),how='outer')
    return in_out.fillna(0).astype(int)
```

```python slideshow={"slide_type": "fragment"}
calc_in_out(G2)
```

<!-- #region slideshow={"slide_type": "subslide"} -->
Is there "nearly" a cycle? If so, how would you find it? What two nodes would you connect?
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "subslide"} -->
**Exercise 4** Given a Eulerian directed graph, find the Eulerian path.

Given: An Eulerian directed graph, in the form of a MultiDiGraph networkx object.

Return: An Eulerian cycle in this graph returned as a list.

Hint: Can you find the correct starting node from the output of ``in_out``?
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
def eulerian_path(G):
    # YOUR SOLUTION HERE
    path = []
    ## BEGIN SOLUTION
    in_out = calc_in_out(G)
    diff = in_out["out"] - in_out["in"]
    
    end = list(in_out.index[diff < 0])[0]
    diff = in_out["in"] - in_out["out"]
    start = list(in_out.index[diff < 0])[0]
    G2 = copy.deepcopy(G)
    G2.add_edge(end,start)
    cycle = eulerian_cycle(G2)
    while not (cycle[0] == start and cycle[-2] == end):
        cycle = [cycle[-2]] + cycle[0:-2] + [cycle[-2]]
    path = cycle[:-1]
    ## END SOLUTION
    return path

## BEGIN SOLUTION
answers["answer_exercise_4"] = tuple(eulerian_path(G2))
## END SOLUTION
eulerian_path(G2)
```

<!-- #region slideshow={"slide_type": "subslide"} -->
We can finally put everything together!

**Exercise 5:** Reconstruct a string from its k-mer composition.

Given: Given a list of $k$-mers patterns.

Return: A string ``text`` with $k$-mer composition equal to patterns.
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
def reconstruct(kmers):
    dB = de_bruijn(kmers)
    path = eulerian_path(dB)
    text = ""
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    prefix = path.pop(0)
    text = [prefix]
    while len(path) > 0:
        suffix = path.pop(0)
        text.append(suffix[-1])
    text = "".join(text)
    ## END SOLUTION
    return text
    
kmers = ["CTTA","ACCA","TACC","GGCT","GCTT","TTAC"]
## BEGIN SOLUTION
answers["answer_exercise_5"] = reconstruct(kmers)
## END SOLUTION
reconstruct(kmers)
```

```python
kmers = ['AAT','ATG','ATG','ATG','CAT','CCA','GAT','GCC','GGA','GGG','GTT','TAA','TGC','TGG','TGT']
de_bruijn(kmers)
show(dB)
```

```python slideshow={"slide_type": "skip"}
## BEGIN SOLUTION
joblib.dump(answers,"../tests/answers_Topic3.joblib");
## END SOLUTION
# Don't forget to push!
```
