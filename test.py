#from CassetteTapeDecay import CassetteTapeDecay
#from PlotClass import PlotClass

from Partition_Module.TapeFS import TapeFS
from Partition_Module.Partition import Partition

if __name__ == "__main__":

        #fig 5E & G
        """xlsx = "RawData.xlsx"

        sim = CassetteTapeDecay.from_xlsx(xlsx)
        plot_class = PlotClass()
        plot_class.model = sim


        print("  D-DNA:", sim.k_dDNA)
        print("  E-DNA:", sim.k_eDNA)
        print("Ea (J/mol):", "D:", sim.ea_d, "E:", sim.ea_e)


        plot_class.plot_fig5G(out_png="Fig5G_sim.png", temp_min=-20.0, temp_max=80.0)
        plot_class.plot_fig5E(xlsx_path=xlsx, out_png="Fig5E_sim.png")"""

        #partition module
        fs = TapeFS()
        assert Partition.partitions_for_label("a") == 12
        assert Partition.partitions_for_label("a1") == 15
        assert Partition.partitions_for_label("a1a") == 18

        label1 = fs.deposit("JK Li","ACCGTC")
        label2 = fs.deposit("JK Li","AGTCATGC")

        recovered_label1 = fs.retrieve(label1)
        print(f"recovered DNA data from {label1}: {recovered_label1}")
        recovered_label2 = fs.retrieve(label2)
        print(f"recovered DNA data from {label2}: {recovered_label2}")

        print("Folder occupancy:", fs.list_folder("JK Li")[:10], "...")

                #<--------------------------------------------------------->


