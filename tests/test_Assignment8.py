import sys
sys.path.append(".")

# Import the student solutions
import Assignment8_helper

import pathlib
DIR=pathlib.Path(__file__).parent.absolute()

import joblib
answers = joblib.load(str(DIR)+"/answers_Assignment8.joblib")
print("Keys",answers.keys())

import numpy as np
answer_tol = 1e-10

spectrum1 = [57,71,154,185,301,332,415,429,486]
spectrum5 = [57,114,128,215,229,316,330,387,444]


def test_exercise_1():
    answer = Assignment8_helper.spectrum_graph_construction(spectrum1)
    assert np.all(Assignment8_helper.to_adj(answer).values == answers['answer_exercise_1'].values)

def test_exercise_2():
    answer = Assignment8_helper.ideal_spectrum("REDCA")
    assert np.all(answer == answers['answer_exercise_2'])

def test_exercise_3():
    answer = Assignment8_helper.decoding_ideal_spectrum(spectrum5,debug=False)
    assert len(set(answer).difference(set(answers['answer_exercise_3']))) == 0

def test_exercise_4():
    answer = Assignment8_helper.construct_peptide_vector("XZZXX")
    assert np.all(answer == answers['answer_exercise_4'])

def test_exercise_5():
    p = np.array([0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1])
    answer = Assignment8_helper.construct_peptide_from_vector(p)
    assert np.all(answer == answers['answer_exercise_5'])

def test_exercise_6():
    p2 = [0,0,0,4,-2,-3,-1,-7,6,5,3,2,1,9,3,-8,0,3,1,2,1,0]
    answer = Assignment8_helper.max_peptide(p2,debug=False)
    assert answer == answers['answer_exercise_6']
