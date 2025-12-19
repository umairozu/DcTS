
#import argparse
import math
import random
from cmath import exp
from pathlib import Path
from typing import Optional, List

import numpy as np
import matplotlib.pyplot as plt


R = 8.31446261815324 # gas const in J/(mol * k)
eDNA_Ea =  133880.342  # E-DNA activation energy value from data
dDNA_Ea = 90813.822     # D-DNA activation energy value from data

"""assuming activation energy of D-DNA <= E-DNA 
(need to refine this with probably stats from MESA
 or DeSP)"""
#dDna_Ea = eDNA_Ea - 1

Temp_ref_C = 20.0
Temp_ref_K = 273.15 + Temp_ref_C

#k = -(math.log(C/Co,math.e)/ t)
#half = math.log(2,math.e)/k/60/60/24/365

"""half-life stats from DNA cassette paper @20 Celsius"""
Half_eDNA = 345.0
Half_dDNA = Half_eDNA / 40.5 # pg 7. eDNa 10x retention vs dDNa

num_weeks_yearly = 52.177456

"""user enters a temp value in Celsius, plus whether the tape is
Encapsulated or not"""
"""

-> Step 1: we get the right activation energy & half life value
-> Step 2: we calculate decay rate const (K_ref)
-> Step 3: using the formula in Supp text 6, log K(t) is calculated
-> Step 4: remove log and get k_T
-> Step 5: re-calculate half_life 
"""
def half_life_years(Temp_C: float, encapsulated: bool) -> float:
    Temp_K = 273.15 + Temp_C
    if encapsulated:
        Ea = eDNA_Ea
        half_life = Half_eDNA
    else:
        Ea = dDNA_Ea
        half_life = Half_dDNA

    k_ref = math.log(2.0) / half_life

    log_k_T = math.log(k_ref) - Ea/R * (1 / Temp_K - 1 / Temp_ref_K)

    k_T = math.exp(log_k_T)

    half_life = math.log(2.0) / k_T
    return  half_life

"""
Using Eq2 from Supp text 6
f(t) = C(t) / C_o = e power -k(T)t
"""
def remaining_dna_fraction(Temp_C: float, encapsulated: bool, num_weeks: float) -> float:
    half_life = half_life_years(Temp_C,encapsulated)
    k_T = math.log(2.0) / half_life
    time_years  = num_weeks / num_weeks_yearly
    return math.exp(-k_T * time_years)

#<---------------------------------------------------------->#

# ---------- Half life vs temperature (degree Celsius)(Fig. 5G style) ---------- #
