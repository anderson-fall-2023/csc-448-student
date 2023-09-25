import sys
sys.path.append(".")

# Import the student solutions
import Topic1_helper

import pathlib
DIR=pathlib.Path(__file__).parent.absolute()

import joblib
answers = joblib.load(str(DIR)+"/answers_Topic1.joblib")

import pandas as pd
data = pd.read_table("http://bioinformaticsalgorithms.com/data/realdatasets/Rearrangements/E_coli.txt",header=None)
genome = data.values[0,0]

text = "atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"

def test_exercise_1():
    assert Topic1_helper.count("ACAACTATGCATACTATCGGGAACTATCCT","ACTAT") == answers['answer_exercise_1a']
    
def test_exercise_2():
    assert Topic1_helper.frequent_words("ACAACTATGCATACTATCGGGAACTATCCT",4) == answers['answer_exercise_2a']

def test_exercise_3():
    assert Topic1_helper.reverse_complement("cagt") == answers['answer_exercise_3']
    
def test_exercise_4():
    assert Topic1_helper.frequency_table(text,3) == answers["answer_exercise_4"]

def test_exercise_5():
    assert Topic1_helper.better_frequent_words(text,9) == answers["answer_exercise_5"]
    
def test_exercise_6():
    assert Topic1_helper.skew(genome) == answers["answer_exercise_6"]
