import sys
sys.path.append(".")

# Import the student solutions
import Assignment5_helper

import pandas as pd

import pathlib
DIR=pathlib.Path(__file__).parent.absolute()

import joblib
answers = joblib.load(str(DIR)+"/answers_Assignment5.joblib")

import networkx as nx
import numpy as np
import copy

file = f'{DIR}/../data/coronavirus_distance_matrix_additive.txt'
D_sars = pd.read_csv(file,index_col=0)

def test_exercise_1():
    assert (Assignment5_helper.compute_d(Assignment5_helper.G) == answers['exercise_1']).all().all()

def test_exercise_2():
    assert Assignment5_helper.limb(Assignment5_helper.D,"v4") == answers['exercise_2']

def test_exercise_3():
    D=Assignment5_helper.D.copy()
    limbLength = Assignment5_helper.limb(D,D.index[-1]) # our algorithm will choose the last node
    n = D.index[-1]
    Dtrimmed = D.drop(n).drop(n,axis=1)
    for j in Dtrimmed.index:
        D.loc[j,n] = D.loc[j,n] - limbLength
        D.loc[n,j] = D.loc[j,n]
    assert Assignment5_helper.find(D,"v4") == answers['exercise_3']

def test_exercise_4():
    assert np.all(nx.adjacency_matrix(Assignment5_helper.base_case(Assignment5_helper.D.iloc[:2,:].iloc[:,:2])).todense() == answers['exercise_4'])

def test_exercise_5():
    assert np.all(nx.adjacency_matrix(Assignment5_helper.additive_phylogeny(Assignment5_helper.D,len(Assignment5_helper.D)+1)).todense() == answers['exercise_5'])

def test_exercise_6():
    assert np.all(nx.adjacency_matrix(Assignment5_helper.additive_phylogeny(D_sars,len(D_sars)+1)).todense() == answers['exercise_6'])

# git clone https://github.com/anderson-github-classroom/csc-448-student ../csc-448-student && sudo -H pip3 install -r ../csc-448-student/requirements.txt

# pytest ../csc-448-student/tests/test_Lab1.py::test_exercise_1
