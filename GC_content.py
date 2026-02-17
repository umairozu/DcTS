
"""
function to calculate error probabilities based on the gc_content.
"""
import os


# MESA GC contents values taken
# goes well with wet-lab steps (PCR, sequencing, synthesis) -->  behaves best at intermediate GC
def error_func(gc_percentage):
    gc_percentage = 1.0 * gc_percentage / 100.0
    if 0.5 <= gc_percentage <= 0.6:
        return 0
    elif gc_percentage > 0.6:
        return (gc_percentage - 0.6) * 2.5
    else:
        return (0.5 - gc_percentage) * 2


def global_gc_content(sequence):
    g_count = 0
    c_count = 0

    for i, _ in enumerate(sequence):
        if sequence[i].upper() == 'G':
            g_count += 1
        elif sequence[i].upper() == 'C':
            c_count += 1

    gc_sum = g_count + c_count

    gc_content_percentage = (gc_sum / len(sequence)) * 100

    return gc_content_percentage



# Currently output is like:
# [(gc_sequence, start_pos, end_pos, gc_content), ( , , , ) , ...]
# [('GCGCATATGCGCGCGCATATGCGCGCGCATATGCGCGCGC', 0, 40, 0.7), ('CGCATATGCGCGCGCATATGCGCGCGCATATGCGCGCGCA', 1, 41, 0.7), ('GCATATGCGCGCGCATATGCGCGCGCATATGCGCGCGCAT', 2, 42, 0.7), ...]

def local_gc_content(sequence, k):
    violating_sequences = []
    for i in range(len(sequence) - k + 1):
        g_count = 0
        c_count = 0
        gc_seq = sequence[i: i + k]
        for j, _ in enumerate(gc_seq):
            if sequence[j].upper() == 'G':
                g_count += 1
            elif sequence[j].upper() == 'C':
                c_count += 1
        gc_sum = g_count + c_count
        current_window_gc = gc_sum / len(gc_seq)
        if 0.70 > current_window_gc > 0.30:
            pass
        else:
            #print(f"local gc window [{gc_seq}] exceeded limit")
            violating_sequences.append((gc_seq,i, i+k,current_window_gc))
    return violating_sequences


def gc_error_probability(sequence):
    return error_func(global_gc_content(sequence))

if __name__ == "__main__":
    #print(gc_error_probability("GCGCATATGCGC"))
    #print(overall_gc_content("TGGCTCATTTCACAATCGGTTAAAGGATTAATGAGGTAGCGCTCCGGAGGGATCGCCCTCCGAACTGAAATGCTCTCAGCTCCCTGGTATCTCCTTATTACCTTCGCGACTGCGGGATCCGATCATAAATGACCTGCCGTGCAA"))
    #print(local_gc_content("GCGCATATGCGCGCGCATATGCGCGCGCATATGCGCGCGCATATGCGCGCGCATATGCGCGCGCATATGCGCGCGCATATGCGCGCGCATATGCGCGCGCATATGCGCGCGCATATGCGC",40))

    # This sequence is from Turkish_anthem.tar.gz.dna_order file, sequence is good, so no local gc outlier, empty list is returned
    # print(local_gc_content("TGGCTCATTTCACAATCGGTCAACCAATACCTTCACCGGAGTGTCTACTCAAGATGAGAGATATATCGGCAGAATCTTACATAGCGTCGTTGCAGGGCGGACGGCGGCCGAGTACTGCCGGATCATAAATGACCTGCCGTGCAA",40))

    input_file = "Turkish_anthem.tar.gz.dna_order.txt"
    max_range = 0.55
    min_range = 0.45
    path = fr'{os.getcwd()}\dna-fountain\{input_file}'
    with open(path, "r") as file:
        lines = []
        for line in file:
            # after primers and cutting site are added, ensuring that GC is still in stable range of 45% - 55%
            # by removing sequence that lie outside this range
            if  max_range > global_gc_content(line)/ 100 > min_range:
                lines.append(line)

        with open(path,"w+") as f:
            f.writelines(lines)




