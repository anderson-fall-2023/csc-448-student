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
# # Topic 4 - Are there fragile regions in the human genome?
#
# ## Clustering Algorithms - Chapter 8
#
# Motivation and some exercises are variations on those available in Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Assignments for week
# * Labs and keep on making progress on the project (keep up the good work)!

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Slack ice breaker
# Best meal you've ever eaten?

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Prelude
#
# There are so many ways we could slice and dice this chapter. Fundamentally, we've got a lot of different angles. We could approach this chapter from a biological and biochemical perspective and focus on the chemsitry and biology necessary to perform transcriptomics. Or we could dive into the statistical approaches necessary to accurately quantify gene expression values. Or we could focus more on alignment algorithms that power the heart of this analysis. We could also focus on the problem beginning at a gene expression values and then focus on algorithms that analyze data. We are going to try to strike a balance in the following order:
# 1. Discuss some of the biochemistry that makes modern sequencing possible
# 2. Discuss some of the different ways scientists investigate what is going on inside a cell
# 3. Discuss the clustering algorithms that are the first lines of the analysis

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Illumina Sequencing
#
# <a href="https://www.youtube.com/watch?v=womKfikWlxM&ab_channel=Illumina">Click here for video</a>

# + [markdown] slideshow={"slide_type": "subslide"}
# <img src="https://edu.t-bio.info/wp-content/uploads/2020/01/Molecular-Data-Cascade.jpg">

# + [markdown] slideshow={"slide_type": "slide"}
# ## Biology from a biologist
#
# <a href="https://calpoly.zoom.us/rec/share/l5hfMH_OtdbAo4-ow76eOgR2F4Lh92mG9YHkEPh6CSElixfS2awWOKcEW34XxrbT.mZKbQHX4hni3U4IO?startTime=1604435224000">Transcriptomics perspective</a>
#
# Watch the first 3 minutes and then answer in your group: "What is the transcriptome?" We'll discuss together about 5 minutes is up.

# + [markdown] slideshow={"slide_type": "slide"}
# ## Clustering
# Clustering or partitioning data into sets is not specific to bioinformatics. Let's first talk about clustering in a generic sense.
#
# <a href="http://anderson-data-science.com/csc_448_2020_fall/clustering.pptx">Slides available here</a>

# + [markdown] slideshow={"slide_type": "slide"}
# ## What's all this about yeast and wine?

# + [markdown] slideshow={"slide_type": "subslide"}
# The species of yeast that we will consider in this chapter is Saccharomyces cerevisiae. Why?

# + [markdown] slideshow={"slide_type": "fragment"}
# **It can brew wine because it converts the glucose found in fruit into ethanol**

# + [markdown] slideshow={"slide_type": "subslide"}
# **Our question:**<br>
#
# If S. cerevisiae often lives on grapevines, **why must crushed grapes be stored in tightly sealed barrels in order to make wine?

# + [markdown] slideshow={"slide_type": "subslide"}
# * If the supply of glucose runs out, S. cerevisiae must do something to survive
# * It will then invert its metabolism, with the ethanol (alcohol) that it just produced becoming its new food supply. 
# * This metabolic inversion, called the diauxic shift, can only occur in the presence of oxygen. 
# * Without oxygen, S. cerevisiae hibernates until either glucose or oxygen becomes available. 
#
# In conclusion, if winemakers don’t seal their barrels, then the yeast in the barrel will metabolize the ethanol that it just produced, ruining the wine.

# + [markdown] slideshow={"slide_type": "subslide"}
# The diauxic shift is a complex process that affects the expression of many genes. 

# + [markdown] slideshow={"slide_type": "slide"}
# ## Our data
# >In 1997, Joseph DeRisi conducted the first massive gene expression experiment by sampling an S. cerevisiae culture every two hours for the six hours before and after the diauxic shift. Since there are approximately 6,400 genes in S. cerevisiae, and there were seven time points, this experiment resulted in a 6,400 × 7 gene expression matrix. 

# + slideshow={"slide_type": "subslide"}
import pandas as pd
df=pd.read_csv('http://bioinformaticsalgorithms.com/data/realdatasets/Clustering/diauxic_raw_ratios_RG.txt',sep='\t')
df

# + [markdown] slideshow={"slide_type": "subslide"}
# Values above 1 in expression vectors correspond to increased expression, while values below 1 correspond to decreased expression.

# + slideshow={"slide_type": "subslide"}
import altair as alt
plot_df = df.set_index('ORF').drop('Name',axis=1).loc[['YPR055W','YLR258W','YPL012W']]
plot_df.columns.name = 'Sample Point'
plot_df = plot_df.stack().to_frame()
plot_df.columns=["Ratio"]
plot_df = plot_df.reset_index()
alt.Chart(plot_df).mark_line().encode(
    x='Sample Point:N',
    y='Ratio',
    color='ORF'
)

# + [markdown] slideshow={"slide_type": "subslide"}
# **Stop and think:** What is the interpretation of this plot?

# + [markdown] slideshow={"slide_type": "subslide"}
# ### YOUR SOLUTION HERE

# + [markdown] slideshow={"slide_type": "subslide"}
# Consider what to do about the other genes?

# + slideshow={"slide_type": "fragment"}
df.shape

# + [markdown] slideshow={"slide_type": "subslide"}
# **Stop and think:** Considering the dataset above and what you now know about clustering, what questions could you ask?

# + [markdown] slideshow={"slide_type": "subslide"}
# ### YOUR SOLUTION HERE
# ### YOUR SOLUTION HERE

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Sample of genes

# + slideshow={"slide_type": "fragment"}
import altair as alt
plot_df = df.set_index('ORF').drop('Name',axis=1).sample(n=100)
plot_df.columns.name = 'Sample Point'
plot_df = plot_df.stack().to_frame()
plot_df.columns=["Ratio"]
plot_df = plot_df.reset_index()
alt.Chart(plot_df).mark_line().encode(
    x='Sample Point:N',
    y='Ratio',
    color='ORF'
)

# + [markdown] slideshow={"slide_type": "subslide"}
# **Stop and think:** What are your observations about this graph?

# + [markdown] slideshow={"slide_type": "subslide"}
# ### YOUR SOLUTION HERE
# ### YOUR SOLUTION HERE

# + [markdown] slideshow={"slide_type": "subslide"}
# Let's remove genes that are not of interest. This is done in the textbook by removing genes that don't go up or down by a significant amount. 

# + slideshow={"slide_type": "fragment"}
df_subset = pd.read_csv('http://bioinformaticsalgorithms.com/data/realdatasets/Clustering/230genes_log_expression.txt',sep='\t')
df_subset

# + [markdown] slideshow={"slide_type": "subslide"}
# ### Redo the plot

# + slideshow={"slide_type": "fragment"}
plot_df = df_subset.set_index('ORF').drop('Name',axis=1).sample(n=100)
plot_df.columns.name = 'Sample Point'
plot_df = plot_df.stack().to_frame()
plot_df.columns=["Ratio"]
plot_df = plot_df.reset_index()
alt.Chart(plot_df).mark_line().encode(
    x='Sample Point:N',
    y='Ratio',
    color='ORF'
)

# + [markdown] slideshow={"slide_type": "subslide"}
# **Stop and think:** Now that you know about k-means clustering, what is a good $k$ value?

# + slideshow={"slide_type": "subslide"}
# our standard imports
import numpy as np
import pandas as pd

from sklearn.cluster import KMeans

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 1:** Using your k value, cluster the genes using k-means. You may use sklearn's version of kmeans. Color the plot above using your clusters.

# + slideshow={"slide_type": "subslide"}
clusterer = KMeans(n_clusters=2, random_state=10)
df_subset["Cluster"] = clusterer.predict(df_subset.drop(['ORF','Name'],axis=1))
df_subset

# + [markdown] slideshow={"slide_type": "subslide"}
# **Problem 2:** Plot all of the genes with color according to their cluster. 

# + slideshow={"slide_type": "subslide"}
### YOUR SOLUTION HERE
### YOUR SOLUTION HERE

# + [markdown] slideshow={"slide_type": "subslide"}
# **Stop and think:** How can you now if you selected the right number of clusters?

# + [markdown] slideshow={"slide_type": "subslide"}
# ### YOUR SOLUTION HERE
# ### YOUR SOLUTION HERE

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 2:** Analyze two clusterings (k=2 and k=3) by calculating the silhouette score.

# + slideshow={"slide_type": "subslide"}
from sklearn.metrics import silhouette_score

### YOUR SOLUTION HERE
cluster2 = clusterer2.predict(df_subset.drop(['ORF','Name','Cluster'],axis=1))
cluster3 = clusterer3.predict(df_subset.drop(['ORF','Name','Cluster'],axis=1))
print('Score for k=2',silhouette_score(df_subset.drop(['ORF','Name','Cluster'],axis=1), cluster2))
print('Score for k=3',silhouette_score(df_subset.drop(['ORF','Name','Cluster'],axis=1), cluster3))
# -


