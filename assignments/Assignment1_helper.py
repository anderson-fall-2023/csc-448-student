def read_fasta(file):
    seqs = []
    headers = []
    # implement a read fasta file
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    f = open(file)
    header = None
    sequence = []
    for line in f:
        line = line.strip()
        if line[0] == ">":
            if header is not None:
                headers.append(header)
                seq = "".join(sequence)
                seqs.append(seq.strip())
                sequence = []
            header = line
        else:
            sequence.append(line.strip())
    if header is not None:
        headers.append(header)
        seq = "".join(sequence)
        seqs.append(seq.strip())
        sequence = []
    ## END SOLUTION    
    return headers,seqs


def avg_length(sequences):
    avg = None
    ### BEGIN SOLUTION
    c = 0
    for seq in sequences:
        c += len(seq)
    avg = c/len(sequences)
    ## END SOLUTION
    return avg