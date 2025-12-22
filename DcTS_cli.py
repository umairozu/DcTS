
from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
import math
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
    plt.yscale("log")
    plt.xlabel("Weeks")
    plt.ylabel("DNA content C / C0")
    plt.title(f"{mode} @ {temp_C:g}°C | k={k_values:.3e} s⁻¹ | t½={half_values:.3g} years")
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











