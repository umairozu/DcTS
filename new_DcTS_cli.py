import math
from dataclasses import dataclass
from typing import Dict

import numpy as np
from openpyxl import load_workbook

R = 8.31446261815324
SEC_PER_WEEK = 7 * 24 * 3600
SEC_PER_YEAR = 365 * 24 * 3600

@dataclass
class CassetteTapeDecay:

    k_d: Dict[int, float]  # D-DNA tape: temp_C -> k   (key, value)
    k_e: Dict[int, float]  # E-DNA tape: temp_C -> k   (key, value)
    ea_d: float = 90813.822     # Activation Energy decapsulated DNA
    ea_e: float = 133880.342    # Activation Energy encapsulated DNA
    lnA_d: float = None
    lnA_e: float = None

    @staticmethod
    def from_rawdata_xlsx(xlsx_path: str) -> "CassetteTapeDecay":
        wb = load_workbook(xlsx_path, data_only=True)
        ws = wb["Arrhenius Calculate"]

        temps = [60, 65, 70]

        # D-DNA tape table data fetch from xlsx
        t_d = np.array([ws["H66"].value, ws["H67"].value, ws["H68"].value], dtype=float)
        cols_d = ["I", "J", "K"]
        k_d = {}
        for temp, col in zip(temps, cols_d):
            y = np.array([ws[f"{col}66"].value, ws[f"{col}67"].value, ws[f"{col}68"].value], dtype=float)
            slope, _ = np.polyfit(t_d, y, 1)
            k_d[temp] = float(-slope)

        # E-DNA tape table data fetch from xlsx
        t_e = np.array([ws["N66"].value, ws["N67"].value, ws["N68"].value, ws["N69"].value], dtype=float)
        cols_e = ["O", "P", "Q"]
        k_e = {}
        for temp, col in zip(temps, cols_e):
            y = np.array([ws[f"{col}66"].value, ws[f"{col}67"].value, ws[f"{col}68"].value, ws[f"{col}69"].value], dtype=float)
            slope, _ = np.polyfit(t_e, y, 1)
            k_e[temp] = float(-slope)

        model = CassetteTapeDecay(k_d=k_d, k_e=k_e)

        model.lnA_d = model.fit_lnA(model.k_d, model.ea_d)
        model.lnA_e = model.fit_lnA(model.k_e, model.ea_e)

        return model


    """    # ln k = ln A - Ea/(R T)  => ln A = ln k + Ea/(R T)  """
    @staticmethod
    def fit_lnA(k_table: Dict[int, float], ea: float) -> float:
        vals = []
        for temp_C, k in k_table.items():
            T = 273.15 + float(temp_C)
            vals.append(math.log(k) + ea / (R * T))
        return float(sum(vals) / len(vals))

    def k(self, temp_C: float, encapsulated: bool) -> float:
        temp_C_int = int(round(temp_C))
        if encapsulated and temp_C_int in self.k_e:
            return self.k_e[temp_C_int]
        if (not encapsulated) and temp_C_int in self.k_d:
            return self.k_d[temp_C_int]

        T = 273.15 + float(temp_C)
        if encapsulated:
            return math.exp(self.lnA_e - self.ea_e / (R * T))
        return math.exp(self.lnA_d - self.ea_d / (R * T))

    def remaining_fraction(self, temp_C: float, encapsulated: bool, weeks: float) -> float:
        t_sec = float(weeks) * SEC_PER_WEEK
        return math.exp(-self.k(temp_C, encapsulated) * t_sec)

    def half_life_years(self, temp_C: float, encapsulated: bool) -> float:
        k_val = self.k(temp_C, encapsulated)
        return (math.log(2) / k_val) / SEC_PER_YEAR

