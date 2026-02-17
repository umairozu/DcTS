
# Error rate --> SOURCED FROM MESA sim
# https://github.com/umr-ds/mesa_dna_sim/blob/master/simulators/error_sources/kmer.py
def error_func(kmer_amount):
    return kmer_amount ** 2 * 0.000002

# Accumulating all the Kmer > 1 errors in a sequence
# Right now, not removing any sequence solely on kmer error sum (error not too big and serious)
# moving on we can do the following:
# total = w1*k_mer + w2*gc + w3*homopolymers
# and remove if total > T

def k_mer(sequence, k):
    k_mer_dict = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i: i + k]
        if kmer not in k_mer_dict:
            k_mer_dict[kmer] = 1
        else:
            k_mer_dict[kmer] += 1

    sum_kmer_error = 0
    for index, value in k_mer_dict.items():
        if value > 1: # meaning more than 1 kmer found
            error = error_func(value)
            sum_kmer_error += error

    return sum_kmer_error
if __name__ == "__main__":
    print(k_mer("TGGCTCATTTCACAATCGGTCCAAACCATGCCTCCGTAAGCTACTATTGCTATGGCCCGCGTGCCACTCGCGCAAAGTTTATAATTGGGCATGATGAGCCGAGAAGGTCGTCCAGAATCAGATCATAAATGACCTGCCGTGCAA",8))




