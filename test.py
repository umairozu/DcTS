from TapeFS import TapeFS
from Partition_Module.Partition import Partition
from DNA_Payload import DNA_Payload
from OligoSequence import OligoSequence

if __name__ == "__main__":

                                        # Degradation module

        """#fig 5E & G
        xlsx = "RawData.xlsx"

        sim = CassetteTapeDecay.from_xlsx(xlsx)
        plot_class = PlotClass()
        plot_class.model = sim


        print("  D-DNA:", sim.k_dDNA)
        print("  E-DNA:", sim.k_eDNA)
        print("Ea (J/mol):", "D:", sim.ea_d, "E:", sim.ea_e)


        plot_class.plot_fig5G(out_png="Fig5G_sim.png", temp_min=-20.0, temp_max=80.0)
        plot_class.plot_fig5E(out_png="Fig5E_sim.png")"""

        #<----------------------------------------------------------------------------------->#



                                        #partition module

        #Oligos
        oligo1 = OligoSequence("A" * 20, "C" * 20, "ACGT" * 16, "G" * 16, "TTAA", "T" * 20).sequence()
        oligo2 = OligoSequence("C" * 20, "G" * 20, "TGCA" * 16, "A" * 16, "CCGG", "A" * 20).sequence()
        oligo3 = OligoSequence("G" * 20, "T" * 20, "AGCT" * 16, "C" * 16, "GGCC", "C" * 20).sequence()
        oligo4 = OligoSequence("T" * 20, "A" * 20, "TCGA" * 16, "T" * 16, "ATGC", "G" * 20).sequence()
        oligo5 = OligoSequence("ACGT" * 5, "TGCA" * 5, "GATT" * 16, "ACGT" * 4, "CTAG", "CGTA" * 5).sequence()
        oligo6 = OligoSequence("GTAC" * 5, "CAGT" * 5, "TTAA" * 16, "GTAC" * 4, "AATT", "TACG" * 5).sequence()
        oligo7 = OligoSequence("AAAAACCCCCGGGGGTTTTT", "CCCCCGGGGGAAAAATTTTT", "ATCG" * 16, "GCGCGCGCGCGCGCGC", "CGCG","TTTTTGGGGGCCCCCAAAAA").sequence()
        oligo8 = OligoSequence("AGCT" * 5, "TCGA" * 5, "CGAT" * 16, "TATATATATATATATA", "TATA", "GCTA" * 5).sequence()
        oligo9 = OligoSequence("CTAG" * 5, "GATC" * 5, "AAGG" * 16, "CCCCAAAAGGGGTTTT", "AGCT", "TAGC" * 5).sequence()
        oligo10 = OligoSequence("GCGC" * 5, "ATAT" * 5, "CCGG" * 16, "ATGCATGCATGCATGC", "GCAT", "CGCG" * 5).sequence()

        """invalid Oligo"""
        #invalid_oligo = OligoSequence("M"*20,"C" * 20, "ACGT" * 16, "G" * 16, "TTAA", "T" * 20)

        # DNA Payloads
        dna_payload_1 = DNA_Payload([(oligo1,1,True)])
        dna_payload_2 = DNA_Payload([(oligo2, 1, True)])
        dna_payload_3 = DNA_Payload([(oligo3, 1, True)])
        dna_payload_4 = DNA_Payload([(oligo4, 1, True)])
        dna_payload_5 = DNA_Payload([(oligo5, 1, True)])
        dna_payload_6 = DNA_Payload([(oligo6, 1, True)])
        dna_payload_7 = DNA_Payload([(oligo7, 4, True)])
        dna_payload_8 = DNA_Payload([(oligo8, 0, True)])
        dna_payload_9 = DNA_Payload([(oligo9, 0, True)])
        dna_payload_10 = DNA_Payload([(oligo10, 0, True)])

        """invalid Oligo Payload"""
        #invalid_dna_payload = DNA_Payload([(invalid_oligo, 2, False)])

        dna_payload_bulk = [dna_payload_1,dna_payload_2,dna_payload_3,dna_payload_4,dna_payload_5,dna_payload_6,dna_payload_7,dna_payload_8,dna_payload_9,dna_payload_10]

        fs = TapeFS()
        #fs.createFolder("d")
        #for payload in dna_payload_bulk:
        #        fs.deposit("a_1", payload)
        #fs.deposit("a_0", dna_payload_7)
        #fs.empty_partition("a_0")
        #fs.empty_partition("a_1")

        #fs.retrieval("a_1","AAAAACCCCCGGGGGTTTTT","TTTTTGGGGGCCCCCAAAAA")
        #fs.retrieval("a_1","AAAAAAAAAAAAAAAAAAAA","TTTTTTTTTTTTTTTTTTTT")

        #fs.removal("a_1","GGGGGGGGGGGGGGGGGGGG","CCCCCCCCCCCCCCCCCCCC")

        #fs.list_folder("a1")

        #fs.edit_Oligo("a_1","AAAAACCCCCGGGGGTTTTT","TTTTTGGGGGCCCCCAAAAA","TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")

        print(fs.barcode_IDs("d"))
        fs.scan_barcode("l_191")


                        #<--------------------------------------------------------->#

