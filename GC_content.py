
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


def overall_gc_content(sequence):
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

def gc_error_probability(sequence):
    return error_func(overall_gc_content(sequence))

if __name__ == "__main__":
    #print(gc_error_probability("GCGCATATGCGC"))
    #print(overall_gc_content("TGGCTCATTTCACAATCGGTTAAAGGATTAATGAGGTAGCGCTCCGGAGGGATCGCCCTCCGAACTGAAATGCTCTCAGCTCCCTGGTATCTCCTTATTACCTTCGCGACTGCGGGATCCGATCATAAATGACCTGCCGTGCAA"))

    input_file = "Turkish_anthem.tar.gz.dna_order.txt"
    path = fr'{os.getcwd()}\dna-fountain\{input_file}'
    with open(path, "r") as file:
        lines = []
        for line in file:
            # after primers and cutting site are added, ensuring that GC is still in stable range of 45% - 55%
            # by removing sequence that lie outside this range
            if  0.55 > overall_gc_content(line)/ 100 > 0.45:
                lines.append(line)

        with open(path,"w+") as f:
            f.writelines(lines)



