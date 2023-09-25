def greedy_lcs(s1,s2):
    lcs = ""
    # YOUR SOLUTION HERE
    return lcs

def greedy_alignment(s1,s2):
    s1_new = ""
    s2_new = ""
    # YOUR SOLUTION HERE
    return "\n".join([s1_new,s2_new])

import numpy as np

def min_num_coins(money,coins):
    min_coins = np.Inf
    # YOUR SOLUTION HERE
    return min_coins

def min_num_coins_dynamic(money,coins):
    min_coins = {0:0} # Base case, no coins needed for no money
    # YOUR SOLUTION HERE
    return min_coins[m]

import pandas as pd
import numpy as np
def align(s1,s2):
    # Below are the exact base cases that I want you to use
    if len(s1) == 0:
        aligned_s1 = "".join(["-" for i in range(len(s2))])
        return 0,aligned_s1,s2
    if len(s2) == 0: # no way to match
        aligned_s2 = "".join(["-" for i in range(len(s1))])
        return 0,s1,aligned_s2
    
    # You don't have to use my dataframe that helps with the choices, but ... I recommend it
    choices_df = pd.DataFrame({
        "remainder(s1)":[s1[1:],s1[1:],s1],
        "remainder(s2)":[s2[1:],s2,s2[1:]],
        "s1_part":[s1[0],s1[0],"-"],
        "s2_part":[s2[0],"-",s2[0]],
        "score(s1_part,s2_part)":[int(s1[0]==s2[0]),0,0]})
    max_score = -np.Inf
    aligned_s1 = None
    aligned_s2 = None
    for i,choice in choices_df.iterrows():
        # here is how to get these values into base Python
        rem_s1,rem_s2,s1_part,s2_part,score = choice.values
        # YOUR SOLUTION HERE
        # print(rem_s1,rem_s2)
    return max_score,aligned_s1,aligned_s2

def align_dynamic(s1,s2):
    scores = pd.DataFrame(index=["-"]+[s1[:i+1] for i in range(len(s1))],columns=["-"]+[s2[:i+1] for i in range(len(s2))])
    for s2_part in scores.columns:
        scores.loc["-",s2_part] = 0
    for s1_part in scores.index:
        scores.loc[s1_part,"-"] = 0
    
    nrows,ncols = scores.shape
    for i in range(1,nrows):
        for j in range(1,ncols):
            # What are our three options
            opt1_s1 = scores.index[i-1] # remember the rows are representative of s1
            opt1_s2 = scores.columns[j-1] # remember the columns are representative of s2
            score_opt1 = -np.Inf # FIX THIS!
            
            opt2_s1 = scores.index[i-1]
            opt2_s2 = scores.columns[j]
            score_opt2 = -np.Inf # FIT THIS!
            
            opt3_s1 = scores.index[i]
            opt3_s2 = scores.columns[j-1]
            score_opt3 = -np.Inf # FIT THIS!
            
            scores.loc[scores.index[i],scores.columns[j]] = max(score_opt1,score_opt2,score_opt3)
            
    return scores.loc[s1,s2]

def align_dynamic2(s1,s2,verbose=False):
    scores = pd.DataFrame(index=["-"]+[s1[:i+1] for i in range(len(s1))],columns=["-"]+[s2[:i+1] for i in range(len(s2))])
    aligned = pd.DataFrame(index=["-"]+[s1[:i+1] for i in range(len(s1))],columns=["-"]+[s2[:i+1] for i in range(len(s2))])
    for s2_part in scores.columns:
        scores.loc["-",s2_part] = 0
        if s2_part == "-":
            aligned.loc["-","-"] = ("","")
        else:
            aligned.loc["-",s2_part] = ("".join(["-" for i in range(len(s2_part))]),s2_part)
    for s1_part in scores.index:
        scores.loc[s1_part,"-"] = 0
        if s1_part == "-":
            aligned.loc["-","-"] = ("","")
        else:
            aligned.loc[s1_part,"-"] = (s1_part,"".join(["-" for i in range(len(s1_part))]))
    if verbose:
        display(aligned)
    
    nrows,ncols = scores.shape
    for i in range(1,nrows):
        for j in range(1,ncols):
            # What are our three options
            opt1_s1 = scores.index[i-1] # remember the rows are representative of s1
            opt1_s2 = scores.columns[j-1] # remember the columns are representative of s2
            score_opt1 = -np.Inf # FIX THIS!
            s1_aligned_opt1 = "" # FIX THIS!
            s2_aligned_opt1 = "" # FIX THIS!
            
            opt2_s1 = scores.index[i-1]
            opt2_s2 = scores.columns[j]
            score_opt2 = -np.Inf # FIT THIS!
            s1_aligned_opt2 = "" # FIX THIS!
            s2_aligned_opt2 = "" # FIX THIS!
            
            opt3_s1 = scores.index[i]
            opt3_s2 = scores.columns[j-1]
            score_opt3 = -np.Inf # FIT THIS!
            s1_aligned_opt3 = "" # FIX THIS!
            s2_aligned_opt3 = "" # FIX THIS!
            
            scores.loc[scores.index[i],scores.columns[j]] = max(score_opt1,score_opt2,score_opt3)
            if max(score_opt1,score_opt2,score_opt3) == score_opt1:
                aligned.loc[scores.index[i],scores.columns[j]] = (s1_aligned_opt1,s2_aligned_opt1)
            elif max(score_opt1,score_opt2,score_opt3) == score_opt2:
                aligned.loc[scores.index[i],scores.columns[j]] = (s1_aligned_opt2,s2_aligned_opt2)
            else:
                aligned.loc[scores.index[i],scores.columns[j]] = (s1_aligned_opt3,s2_aligned_opt3)
    if verbose:
        display(scores)
        display(aligned)
    return scores.loc[s1,s2],aligned.loc[s1,s2][0],aligned.loc[s1,s2][1]