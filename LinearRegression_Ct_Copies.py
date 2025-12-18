import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


df = pd.read_excel("RawData.xlsx", sheet_name="Sheet1")
df_std = df[["Ct", "Copies"]].apply(pd.to_numeric, errors="coerce").dropna()

df_std = df_std[df_std["Copies"] > 0]

ct = df_std["Ct"].to_numpy()
log_copies = np.log10(df_std["Copies"].to_numpy())

slope, intercept, r_value, p_value, std_err = linregress(log_copies, ct)

print(f"Slope: {slope:.3f}")
print(f"Intercept: {intercept:.3f}")
print(f"R²: {r_value**2:.4f}")

efficiency = (10 ** (-1 / slope) - 1) * 100
print(f"PCR Efficiency: {efficiency:.1f}%")

x_fit = log_copies
y_fit = slope * x_fit + intercept

plt.figure(figsize=(6, 4))
plt.scatter(log_copies, ct, label="Standards")
plt.plot(x_fit, y_fit, label="Linear fit")
plt.gca().invert_yaxis()
plt.xlabel("log10(Copies)")
plt.ylabel("Ct")
plt.title("qPCR Standard Curve")
plt.legend()
plt.tight_layout()
plt.show()
