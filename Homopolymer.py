from itertools import groupby

# Error rates Sources from MESA
# https://github.com/umr-ds/mesa_dna_sim/blob/master/simulators/error_sources/homopolymers.py

def error_func(homopolymer_length, base= None):
    if homopolymer_length < 3:
        return 0.0
    elif homopolymer_length < 4:
        return 0.3
    elif homopolymer_length < 5:
        return 0.6
    elif homopolymer_length < 6:
        return 0.9
    else:
        return 1.0

# take a sequence and error_fun
# return list of tuple of (Base,homopolymer_error probability)
def homopolymer(sequence):
    result = []

    # k (The Key): This is the character that is currently repeating (e.g., 'A')
    # g (The Group): This is an iterator containing all the items in that specific streak.
    # Because it is an iterator,you usually have to convert it to a list(g) to see the contents.

    max_homopolymers = [list(g) for k, g in groupby(sequence)] # just collecting all homopolymers here

    for seq in max_homopolymers:
        length = len(seq)
        error = error_func(length)
        if error > 0.0:
            result.append((seq,error))

    return result

if __name__ == "__main__":
    print(homopolymer("TGGCTCATTTCACAATCGGTCAACCAATACCTTCACCGGAGTGTCTACTCAAGATGAGAGATATATCGGCAGAATCTTACATAGCGTCGTTGCAGGGCGGACGGCGGCCGAGTACTGCCGGATCATAAATGACCTGCCGTGCAA"))






