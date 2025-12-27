
from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
from typing import Optional

from CassetteTapeDecay import CassetteTapeDecay

def isDna(seq: str)-> bool:
    if not seq:
        return False
    for char in seq.upper():
        if char not in "ACGT":
            return False
    return True

def fig_5E(
        model: CassetteTapeDecay, temp_C: float, encapsulated: bool,
        weeks: float, graph_points: int, out_png: str
            ):

    t_grid = np.linspace(0.0,weeks,graph_points)
    fractions = np.array([model.remaining_dna_frac(temp_C,encapsulated,w) for w in t_grid], dtype= float)
    k_values = model.k(temp_C,encapsulated)
    half_values = model.half_life(temp_C, encapsulated)

    mode = "E-DNA" if encapsulated else "D-DNA"

    plt.figure(figsize=(7,5))
    plt.plot(t_grid,fractions)
    plt.ylim(1e-6, 1.2)
    plt.yticks([1, 1e-2, 1e-4, 1e-6], ["1", "0.01", "0.0001", "0.000001"])
    plt.yscale("linear")
    plt.xlabel("Weeks")
    plt.ylabel("DNA content C / C0")
    plt.title(f"{mode} at {temp_C:g}°C | k={k_values:.3e} s⁻¹ | t½={half_values:.5g} years")
    plt.tight_layout()
    plt.savefig(out_png, dpi = 200)



def fig5G(
        model: CassetteTapeDecay, temp_min: float, temp_max: float,
        step_c: float, mark_temp: Optional[float],out_png: str
            ):

    temps = np.arange(float(temp_min), float(temp_max), float(step_c))

    half_D = np.array([model.half_life(float(T), encapsulated=False) for T in temps], dtype=float)
    half_E = np.array([model.half_life(float(T), encapsulated=True)  for T in temps], dtype=float)

    plt.figure(figsize=(7, 5))
    plt.plot(temps, half_D, linestyle="--", label="D-DNA")
    plt.plot(temps, half_E, linestyle="-",  label="E-DNA")

    if mark_temp is not None:
        half_D_pt = model.half_life(mark_temp, encapsulated=False)
        half_E_pt = model.half_life(mark_temp, encapsulated=True)
        plt.scatter([mark_temp], [half_D_pt], marker="s", label=f"D @ {mark_temp:g}°C")
        plt.scatter([mark_temp], [half_E_pt], marker="o", label=f"E @ {mark_temp:g}°C")

    plt.yscale("log")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Half-life (years)")
    plt.title(" Half-life vs Temperature (°C)")
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)

def main():
    p = argparse.ArgumentParser(description="DNA cassette tape simulator CLI")

    #p.add_argument("--xlsx", required= True, help="Path to RawData.xlsx")
    p.add_argument("--seq", required=True, help="DNA sequence (A/C/G/T)")
    p.add_argument("--temp", required= True, type= float, help="temperature in celsius for Fig5E and 5G")
    p.add_argument("--encapsulated", action= "store_true", help= "if encapsulated use E-DNA, else D-DNA")

    p.add_argument("--weeks", type= float, default= 3.0, help= "Simulating upto mentioned number of weeks (default = 3)")

    p.add_argument("--temp_min", type= float, default= -20.0, help="Min temp for fig5G (default = -20.0)")
    p.add_argument("--temp_max", type= float, default= 80.0, help= "Max temp for fig5G (default = 80.0)")

    p.add_argument("--save", default= "output", help= "prefix for output file (e.g run1_fig5G.png)")

    args = p.parse_args()

    if not isDna(args.seq):
        sys.exit("Input sequence must only contain A/C/G/T characters")

    model = CassetteTapeDecay.from_xlsx("RawData.xlsx")

    out_fig5E = f"{args.save}_fig5E.png"
    out_fig5G = f"{args.save}_fig5G.png"

    fig_5E(
        model = model,
        temp_C = args.temp,
        encapsulated= args.encapsulated,
        weeks = args.weeks if args.weeks > 0.0 else sys.exit("Invalid num weeks"),
        graph_points= 200,
        out_png= out_fig5E
    )

    fig5G(
        model = model,
        temp_min = args.temp_min,
        temp_max = args.temp_max,
        step_c = 1.0,
        mark_temp = args.temp,
        out_png= out_fig5G
    )

    """Terminal outputs"""
    mode = "E-DNA (encapsulated DNA)" if args.encapsulated else "D-DNA (decapsulated DNA)"
    print("=== DcTS --DNA cassette tape simulator ===")
    print(f"Sequence length: {len(args.seq)} nt")
    print(f"Fig 5E mode: {mode}")
    print(f"Temp for fig 5E: {args.temp} °C" )
    print(f"Saved fig 5E as: {out_fig5E}")
    print(f"Sved fig 5G as: {out_fig5G}")

if __name__ == "__main__":
        main()





