from CassetteTapeDecay import CassetteTapeDecay
from PlotClass import PlotClass

from Partition_Module.TapeFS import TapeFS
from Partition_Module.Partition import Partition
from DNA_Payload import DNA_Payload
from OligoSequence import OligoSequence

if __name__ == "__main__":

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
        oligo1 = OligoSequence("ADL_","RS_CODE_","ACGTACTC_","SEED_","ENZYME_","ADR_")
        dna = DNA_Payload([(oligo1,2,True)])
        # OR
        #dna = DNA_Payload(True, [(oligo1.sequence(), 2)])


        fs = TapeFS()
        """assert Partition.partitions_for_label("a") == 12
        assert Partition.partitions_for_label("a1") == 15
        assert Partition.partitions_for_label("a1a") == 18"""

        myfolder = fs.createFolder("a")
        #depo_1 = fs.deposit("a_0", oligo1, 1, True)
        depo_6 = fs.deposit("a_1", oligo1, 1, True)
        print("<-------->")
        depo_6_retrieval = fs.retrieve("a_1","ADL_","ADR_")
        print("<-------->")
        #fs.removal("a_1","ADL_", "ADR_")
        print("<-------->")
        #depo_6_retrieval_02 = fs.retrieve("a_11","ADL_","ADR_")
        fs.insertion("a_1","ADL_", "ADR_")

        print("<---------------------------------------->")

        #label2 = fs.deposit("a_0", oligo1, 1, True)
        #label3 = fs.deposit("a_0", oligo1, 1, True)
        #label4 = fs.deposit("a_0", oligo1, 1, True)
        #should not be deposited
        #label5 = fs.deposit("a_0", oligo1, 1, True)



        #label7 = fs.deposit("a_1", oligo1, 1, True)
        #label8 = fs.deposit("a_1", oligo1, 1, True)
        #label9 = fs.deposit("a_1", oligo1, 1, True)
        # should not be deposited
        #label10 = fs.deposit("a_1", oligo1, 1, True)


        #recovered_label1 = fs.retrieve(label1)
        #print(f"recovered DNA data from {label1}: {recovered_label1}")

        #recovered_label2 = fs.retrieve(label2)
        #print(f"recovered DNA data from {label2}: {recovered_label2}")

        print("Folder occupancy:", fs.list_folder("a")[:5])

         #<--------------------------------------------------------->#

