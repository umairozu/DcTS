from CassetteTapeDecay import CassetteTapeDecay
from PlotClass import PlotClass



if __name__ == "__main__":
        xlsx = "RawData.xlsx"

        sim = CassetteTapeDecay.from_xlsx(xlsx)
        plot_class = PlotClass()
        plot_class.model = sim


        print("  D-DNA:", sim.k_dDNA)
        print("  E-DNA:", sim.k_eDNA)
        print("Ea (J/mol):", "D:", sim.ea_d, "E:", sim.ea_e)


        plot_class.plot_fig5G(out_png="Fig5G_sim.png", temp_min=-20.0, temp_max=80.0)
        plot_class.plot_fig5E(xlsx_path=xlsx, out_png="Fig5E_sim.png")

                #<--------------------------------------------------------->


