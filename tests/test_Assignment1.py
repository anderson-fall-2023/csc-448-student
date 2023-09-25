import sys
sys.path.append(".")

# Import the student solutions
import Assignment1_helper

import pathlib
DIR=pathlib.Path(__file__).parent.absolute()

import joblib
answers = joblib.load(str(DIR)+"/answers_Assignment1.joblib")

import pandas as pd
data = pd.read_table("http://bioinformaticsalgorithms.com/data/realdatasets/Rearrangements/E_coli.txt",header=None)
genome = data.values[0,0]

text = "atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"

file = f"{DIR}/../data/GCF_000146045.2_R64_genomic.fna"
    
def test_exercise_1():
    headers,sequences = Assignment1_helper.read_fasta(file)
    c = 0
    for seq in sequences:
        c += len(seq)
    avg1 = c/len(sequences)
    headers,sequences = answers["answer_exercise_1"]
    c = 0
    for seq in sequences:
        c += len(seq)
    avg2 = c/len(sequences)
    assert abs(avg1 - avg2) < 0.00001
    
def test_exercise_2():
    sequences = ["ccacacca","cacccacacacccacacaccacaccacacaccacacca","cacacacaccacacccacacca","caccacaccacacccacacccaca"]
    avg = Assignment1_helper.avg_length(sequences)
    assert abs(avg - answers["answer_exercise_2"]) < 0.00001
