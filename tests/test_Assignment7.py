import sys
sys.path.append(".")

# Import the student solutions
import Assignment7_helper

import pathlib
DIR=pathlib.Path(__file__).parent.absolute()

import joblib
answers = joblib.load(str(DIR)+"/answers_Assignment7.joblib")
print("Keys",answers.keys())

import numpy as np
answer_tol = 1e-10

def test_exercise_1():
    answer = Assignment7_helper.trie_construction(Assignment7_helper.patterns2)
    assert np.all(Assignment7_helper.to_adj(answer).values == answers['answer_exercise_1'].values)

def test_exercise_2():
    trie2 = Assignment7_helper.trie_construction(Assignment7_helper.patterns2)
    answer = Assignment7_helper.trie_matching("bananablahblahantennanabnablkjdf",trie2)
    assert np.all(tuple(answer) == tuple(answers['answer_exercise_2']))

def test_exercise_3():
    answer = Assignment7_helper.suffix_trie("panamabananas$")
    assert np.all(Assignment7_helper.to_adj(answer).values == answers['answer_exercise_3'].values)

def test_exercise_4():
    answer,discard = Assignment7_helper.modified_suffix_trie("panamabananas$")
    assert np.all(Assignment7_helper.to_adj(answer).values == answers['answer_exercise_4'].values)

def test_exercise_5():
    answer = Assignment7_helper.to_adj(Assignment7_helper.suffix_tree_construction("panamabananas$"))
    answer.index = [str(c) for c in answer.index]
    answer.columns = [str(c) for c in answer.columns]
    instructor_answer = answers['answer_exercise_5']
    answer_order = list(instructor_answer.index)
    answer = answer.loc[answer_order,answer_order]
    assert np.all(answer.values == instructor_answer.values)
