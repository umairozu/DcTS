import matplotlib.pyplot as plt
from typing import Optional
import numpy as np
import math
from new_DcTS_cli import CassetteTapeDecay


class PlotClass:
    model: CassetteTapeDecay

    def plot_fig5E(self, xlsx_path: str, out_png: Optional[str] = None) -> None:
        temps = [60, 65, 70]
        weeks = [0, 1, 2, 3]

        plt.figure(figsize=(8, 6))

        for temp in temps:
            if temp >= 70:
                weeks_D = [0, 1, 2]
            else:
                weeks_D = weeks

            sim_D = [self.model.remaining_dna_frac(temp, False, w) for w in weeks_D]
            plt.plot(weeks_D, sim_D, linestyle=":", marker=None, label=f"D-DNA sim {temp}°C")

            sim_E = [self.model.remaining_dna_frac(temp, True, w) for w in weeks]
            plt.plot(weeks, sim_E, linestyle="-", marker=None, label=f"E-DNA sim {temp}°C")

        plt.yscale("linear")
        plt.xlabel("Time (weeks)")
        plt.ylabel("Remaining fraction (C/C0)")
        plt.title("Fig 5E)")
        plt.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1.3, 1))
        plt.tight_layout()

        if out_png:
            plt.savefig(out_png, dpi=200)
        plt.show()








