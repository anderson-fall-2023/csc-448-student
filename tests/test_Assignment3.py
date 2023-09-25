import sys
sys.path.append(".")

# Import the student solutions
import Assignment3_helper

import pathlib
DIR=pathlib.Path(__file__).parent.absolute()

import joblib
answers = joblib.load(str(DIR)+"/answers_Assignment3.joblib")
print("Keys",answers.keys())

import numpy as np
import networkx as nx

def cycles_equal(cycle1,cycle2):
    assert len(cycle1) == len(cycle2)
    assert cycle1[0] == cycle1[-1]
    assert cycle2[0] == cycle2[-1]
    if np.all(cycle1==cycle2):
        return True
    for i in range(len(cycle2)):
        cycle1 = [cycle1[-2]] + cycle1[0:-2] + [cycle1[-2]]
        print(cycle1,cycle2)
        if np.all(np.array(cycle1)==np.array(cycle2)):
            return True
    return False

def test_exercise_1():
    assert tuple(answers['exercise_1']) == tuple(Assignment3_helper.composition(3,"TATGGGGTGC"))

def test_exercise_2():
    assert np.all(answers['exercise_2'] == nx.adjacency_matrix(Assignment3_helper.de_bruijn(["AAT","ATG","ATG","ATG","CAT","CCA","GAT","GCC","GGA","GGG","GTT","TAA","TGC","TGG","TGT"])))

def test_exercise_3():
    assert cycles_equal(answers['exercise_3'],Assignment3_helper.eulerian_cycle(Assignment3_helper.G,start=6))

def test_exercise_4():
    assert tuple(answers['exercise_4']) == tuple(Assignment3_helper.eulerian_path(Assignment3_helper.G2))

def test_exercise_5():
    assert answers['exercise_5'] == Assignment3_helper.reconstruct(Assignment3_helper.kmers)