import sys
sys.path.append(".")

# Import the student solutions
import Topic2_helper
import Assignment2_helper

import pathlib
DIR=pathlib.Path(__file__).parent.absolute()

import joblib
answers = joblib.load(str(DIR)+"/answers_Assignment2.joblib")
print("Keys",answers.keys())

def test_exercise_1():
    assert Topic2_helper.greedy_lcs("AACCTTGG","ACACTGTGA") == answers['answer_exercise_1']

def test_exercise_2():
    assert Topic2_helper.greedy_alignment("AACCTTGG","ACACTGTGA") == answers['answer_exercise_2']

def test_exercise_3():
    assert Topic2_helper.min_num_coins(27,[6,5,1]) == answers['answer_exercise_3']

def test_exercise_4():
    assert Topic2_helper.min_num_coins_dynamic(27,[6,5,1]) == answers['answer_exercise_4']

def test_exercise_5():
    assert Topic2_helper.align("AACCT","ACACTG") == answers['answer_exercise_5']

def test_exercise_6():
    assert Topic2_helper.align_dynamic("AACCT","ACACTG") == answers['answer_exercise_6']

def test_exercise_7():
    assert Topic2_helper.align_dynamic2("AACCT","ACACTG") == answers['answer_exercise_7']

def test_exercise_8():
    s1="CGCAACCACAGCGCGCAGGGCAGGCGCGAGCTGTCTGAGCCCCGGCCTCGGACCGCCCACTGGACTCCCGGCACGCCCGGTGCCGCCTTCCGGCTCCAGTCCCCC"
    s2="CGCAACGGCAGCGCGCAGGGCAGGCGCGAGCTGGCCTCTGAGCCCCGGCCTCGGACCGCCCACTCCACGCCCGGCAGGCCCGGTGCCGCCTTCCGGCTCCAGTCCCCCCGC"
    score_1,aligned_s1_1,aligned_s2_1 = Assignment2_helper.align_dynamic3(s1,s2,match_score=1,mismatch_score=0,gap_score=0)
    score_2,aligned_s1_2,aligned_s2_2 = Assignment2_helper.align_dynamic3(s1,s2,match_score=2,mismatch_score=-3,gap_score=-1)

    assert (score_1 == answers['answer_exercise_8'][0]) and (score_2 == answers['answer_exercise_8'][1])


