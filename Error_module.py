import random
from itertools import chain
import numpy as np
from ordered_set import OrderedSet

from Homopolymer import homopolymer

# sequencing_error.py

# Error_rates taken from mesa
# https://github.com/umr-ds/mesa_dna_sim/blob/master/simulators/sequencing/sequencing_error.py

err_rates = {"1": {"raw_rate": 0.0021, "mismatch": 0.81, "deletion": 0.0024, "insertion": 0.0013},
             "2": {"raw_rate": 0.0032, "mismatch": 0.79, "deletion": 0.0018, "insertion": 0.0011},
             "3": {"raw_rate": 0.02, "mismatch": 0.75, "deletion": 0.20, "insertion": 0.05},
             "4": {"raw_rate": 0.14, "mismatch": 0.37, "deletion": 0.21, "insertion": 0.42},
             "5": {"raw_rate": 0.2, "mismatch": 0.48, "deletion": 0.37, "insertion": 0.15},
             "6": {"raw_rate": 0.13, "mismatch": 0.41, "deletion": 0.36, "insertion": 0.23}}

mutation_attributes = {"1": {"deletion": {"position": {"random": 1},
                                          "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}},
                             "insertion": {"position": {"random": 1},
                                           "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}},
                             "mismatch": {"pattern": {"A": {"G": 0.50, "T": 0.25, "C": 0.25},
                                                      "T": {"G": 0.50, "A": 0.25, "C": 0.25},
                                                      "C": {"G": 0.50, "A": 0.25, "T": 0.25},
                                                      "G": {"T": 0.50, "A": 0.25, "C": 0.25}}}},
                       "2": {"deletion": {"position": {"random": 1},
                                          "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}},
                             "insertion": {"position": {"random": 1},
                                           "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}},
                             "mismatch": {"pattern": {"A": {"G": 0.50, "T": 0.25, "C": 0.25},
                                                      "T": {"G": 0.50, "A": 0.25, "C": 0.25},
                                                      "C": {"G": 0.50, "A": 0.25, "T": 0.25},
                                                      "G": {"T": 0.50, "A": 0.25, "C": 0.25}}}},
                       "3": {"deletion": {"position": {"homopolymer": 0.85, "random": 0.15},
                                          "pattern": {"G": 0.35, "C": 0.35, "A": 0.15, "T": 0.15}},
                             "insertion": {"position": {"homopolymer": 0.85, "random": 0.15},
                                           "pattern": {"A": 0.35, "T": 0.35, "C": 0.15, "G": 0.15}},
                             "mismatch": {"pattern": {"CG": {"CA": 0.5, "TG": 0.5}}}},
                       "4": {"deletion": {"position": {"homopolymer": 0.85, "random": 0.15},
                                          "pattern": {"G": 0.35, "C": 0.35, "A": 0.15, "T": 0.15}},
                             "insertion": {"position": {"homopolymer": 0.85, "random": 0.15},
                                           "pattern": {"A": 0.35, "T": 0.35, "C": 0.15, "G": 0.15}},
                             "mismatch": {"pattern": {"CG": {"CA": 0.5, "TG": 0.5}}}},
                       "5": {"deletion": {"position": {"homopolymer": 0.46, "random": 0.54},
                                          "pattern": {"G": 0.35, "C": 0.35, "A": 0.15, "T": 0.15}},
                             "insertion": {"position": {"homopolymer": 0.46, "random": 0.54},
                                           "pattern": {"A": 0.35, "T": 0.35, "C": 0.15, "G": 0.15}},
                             "mismatch": {"pattern": {"TAG": "TGG", "TAC": "TGC"}}},
                       "6": {"deletion": {"position": {"homopolymer": 0.46, "random": 0.54},
                                          "pattern": {"G": 0.35, "C": 0.35, "A": 0.15, "T": 0.15}},
                             "insertion": {"position": {"homopolymer": 0.46, "random": 0.54},
                                           "pattern": {"A": 0.35, "T": 0.35, "C": 0.15, "G": 0.15}},
                             "mismatch": {"pattern": {"TAG": "TGG", "TAC": "TGC"}}}}

class sequencingError:

    # user can provide their own mutation_attributes and error rates but should be of the same format
    def __init__(self, seq, process,attribute = None,error_rate = None, seed= None):
        self.bases = ['A', 'T', 'C', 'G']
        self.seq = seq
        self.process = process
        self.attributes = attribute
        self.error_rates = error_rate
        self.seed = seed if seed else np.random.seed()
        self.visited_bases = [{"base": char, "visited": False} for char in self.seq]


    @staticmethod
    def get_attributes(indels_type):
        # position is "Random" or "Homopolymer" location
        try:
            position = np.random.choice(list(indels_type["position"].keys()), p = list(indels_type["position"].values()))
        except KeyError:
            position = None

        try:
            pattern = indels_type["pattern"]
        except KeyError:
            pattern = None

        position_range = [20,100] # this is just an example initialization
        return position, pattern, position_range

    def insertion(self, ins_attr_dict = None):
        # optional insertion attribute dictionary
        if ins_attr_dict:
            position, pattern, position_range = self.get_attributes(ins_attr_dict)
        # default dictionary
        else:
            ins_dict = self.attributes['insertion']
            position, pattern, position_range = self.attributes(ins_dict)

        if not position or position == 'random':
            return self.indel(pattern,position_range,mode = 'insertion')

        if position == 'homopolymer':
            poly = homopolymer(self.seq)
            if poly:
                return self.indel_homopolymer(poly, pattern, mode = 'insertion')
        # if not position or position != random or position == homopolymer but poly is empty --> then random insertion
        return self.indel(pattern, position_range, mode='insertion')


    def indel_homopolymer(self, poly, pattern, mode):
        print(f"provided polymer: {poly} ")
        print(f"provided pattern: {pattern} ")

        poly_bases = list(OrderedSet(chain.from_iterable(poly)))
        print(f"Poly bases: {poly_bases}" )
        new_pattern = []
        new_pattern_weights = {}

        if not poly: # if no homopolymer available then go back to normal indel method
            return self.indel(pattern = None, position_range = None, mode = mode)

        for base in poly_bases:
            if base in pattern:
                new_pattern.append(base)  # this for loop is removing those bases from
                                          # pattern that are not in any of our homopolymers

            """        
            if not pattern: # if no pattern probabilities provided, then give equal weights to each base
            for base in new_pattern:
                new_pattern_weights[base] = 1 / len(new_pattern)
            else:
            """
            #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            # instead i could make a validate pattern/ mutation_attributes method or class
            #~~~~~~~~~~~~~~~IMPORTANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        #sum of bases in new_pattern
        sum_bases = 0
        for base in new_pattern:
            sum_bases += pattern[base]

        """normalizing bases weight to 1"""
        #new_pattern_weights = {}
        for base in new_pattern:
            new_pattern_weights[base] = pattern[base]/sum_bases

        print(f"new pattern: {new_pattern}")
        print(f"new pattern weights: {new_pattern_weights}")


        chosen_base = np.random.choice(list(new_pattern_weights.keys()), p= list(new_pattern_weights.values()))
        print(f"chosen base from new pattern: {chosen_base}")

        #poly is like [['T', 'T', 'T', 'T', 'T'], ['T', 'T', 'T'], ['G', 'G', 'G'], ['A', 'A', 'A']]
        # check homopolymer.py for details
        possible_mutables = ["".join(item) for item in poly if chosen_base in item]
        print(f"possible mutables: {possible_mutables}")

        chosen_mutable = random.choice(possible_mutables) # a list here right now e.g ['T', 'T', 'T']
        chosen_mutable = "".join(chosen_mutable) # converted to a string
        print(f"chosen mutable: {chosen_mutable}")

        index = random.randrange(len(chosen_mutable))
        print(f"chosen index:  {index} ")

        base = random.choice(self.bases)
        print(f"chosen base:  {base} ")

        # homopolymer is being mutated here based on mode
        if mode == 'insertion':
            chosen_mutable = chosen_mutable[:index] + base + chosen_mutable[index:]
        elif mode == 'deletion':
            chosen_mutable = chosen_mutable[:index] + " " + chosen_mutable[index + 1:]


        """Output for insertion as an example:"""
        """
        provided polymer: [['T', 'T', 'T', 'T', 'T'], ['T', 'T', 'T'], ['G', 'G', 'G'], ['A', 'A', 'A']] 
        provided pattern: {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25} 
        Poly bases: ['T', 'G', 'A']
        new pattern: ['T', 'G', 'A']
        new pattern weights: {'T': 0.3333333333333333, 'G': 0.3333333333333333, 'A': 0.3333333333333333}
        chosen base from new pattern: T
        possible mutables: ['TTTTT', 'TTT']
        chosen mutable: TTT
        chosen index:  0 
        chosen base:  G 
        mutated homopolymer: GTTT
        """

        return chosen_mutable


    """
    If a pattern is provided e.g; {"pattern": {"G": 0.35, "C": 0.35, "A": 0.15, "T": 0.15}} :
    pick a target base type and delegate it to def random_indel()
    Else:
        If pattern not provided:
            check If positon range is provided:
                if so, generate a random position base ensuring a empty " " position is not retrieved,
            Otherwise: pick a random position from the entire sequence
            -delegate position to def indel_sub_base()
    """

    def indel(self, pattern, position_range, mode):
        if not pattern:
            if position_range:
                pos = random.randrange(position_range[0], position_range[1] + 1)
                while self.seq[pos] == " ": # if the chosen index is " ", chose a different index then
                    pos = random.randrange(position_range[0], position_range[1] + 1)
            else:
                pos = random.randrange(len(self.seq))
            return self.indel_sub_base(pos, mode)
        else:
            target_base = np.random.choice(list(pattern.keys()), p = list(pattern.values()))
            return self.random_indel(target_base,position_range,mode)

    """
    This method is trying to randomly find a index(with some rules ofc) that matches the target base
    - then pass that index/ pos to def indel_sub_base for mutation based on the mode ('insertion', deletion') specified
    """
    def random_indel(self,target_base, position_range, mode, count = 0):
        if position_range:
            sequence_indices = range(position_range[0], position_range[1] + 1)
        else:
            sequence_indices = range(len(self.seq))

        valid_indices = [i for i in sequence_indices if self.seq[i] != " "]

        if valid_indices is None:
            return None ####

        indices = [i for i in valid_indices if self.seq[i] == target_base]

        pos = random.choice(valid_indices)

        if indices:
            pick_index = random.choice(indices)
            return self.indel_sub_base(pick_index, mode)
        else:
            if count < 12: # recursively try finding a random valid index using other bases (bases chosen randomly here also)
                target_base = random.choice(self.bases)
                return self.random_indel(target_base,position_range,mode, count = count + 1)

        # if unable to find the target base from amongst the valid sequence_indices, then mutate the
        # pos = random.choice(sequence_indices) anyways
        return self.indel_sub_base(pos,mode)



    """Insertion, deletion, substitution implementing methods"""
    def indel_sub_base(self, pos, mode):

        """2 improvements could be done later(if needed):
        - currently If position is already visited, nothing happens (silent skip),
            we can close this condition for increase mutation rate if needed
        - currently if base equals the original base, no substitution will happen, we can improve this
            by allowing any of the other 3 basses to go at that substitution place
        """

        assert mode in ("insertion", "deletion", "substitution")
        assert  0 <= pos <= len(self.seq)

        base = random.choice(self.bases)
        new_mutation = {"base": base, "visited": True}

        if self.visited_bases[pos]["visited"] == False:
            if mode == 'insertion':  # pre-insertion to be specific!
                self.seq = self.seq[:pos] + base + self.seq[pos:]
                self.visited_bases.insert(pos, new_mutation)  # not touching original base (at pos) again
            elif mode == 'deletion':
                self.seq = self.seq[:pos] + " " + self.seq[pos + 1:]
                self.visited_bases[pos] = {"base": " ", "visited": True}
            else:  # substitution
                self.seq = self.seq[:pos] + base + self.seq[pos + 1:]
                self.visited_bases[pos] = new_mutation


            #Working on HOMOPOLYMER SIDE OF INSERTION


if __name__ == "__main__":
    #my_poly = [['T', 'T', 'T', 'T', 'T'], ['T', 'T', 'T'], ['G', 'G', 'G'], ['A', 'A', 'A']]
    my_pattern = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    ins_mode = "insertion"
    del_mode = "deletion"

    sequence = "TGGCTCATTTCACAATCGGTAAAGAAAGGGAAGGAATAGGTTACTAGGCCCAACCGCAAGCCCTTTGGTCAACCGCAGTGGAAAGAAGGGCTAATAGGTCCTGGTAGATTTACCACTGAAGATCATAAATGACCTGCCGTGCAA"
    my_poly2 = homopolymer(sequence)

    sE = sequencingError(sequence, "sequencing", mutation_attributes,err_rates)

    indel_in_homopolymer = sE.indel_homopolymer(my_poly2,my_pattern,ins_mode)
    print(indel_in_homopolymer)

