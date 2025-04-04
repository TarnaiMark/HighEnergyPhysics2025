import numpy as np
import matplotlib
matplotlib.use('Agg')  # headless environment fix
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Load simulation output file
data = np.loadtxt("output/ising_sim.txt", skiprows=1)


T = data[:, 0]
R42_simple = data[:, 1]         # C4/C2
R42_subtracted = data[:, 2]     # C4/C2 - 3*C2

# Define critical temperature (approximate) 2/ln(1+sqrt)
T_CEP = 2.26

# Interpolation
interp_simple = interp1d(T, R42_simple, kind='cubic')
interp_sub = interp1d(T, R42_subtracted, kind='cubic')

R_simple_CEP = interp_simple(T_CEP)
R_sub_CEP = interp_sub(T_CEP)

# Plot side-by-side
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# --- First plot: C4/C2 ---
axs[0].plot(T, R42_simple, 'bo-', label=r'$C_4/C_2$')
axs[0].axvline(T_CEP, color='red', linestyle='--', linewidth=1)
axs[0].plot(T_CEP, R_simple_CEP, 'ro', label='CEP')
axs[0].set_xlabel(r'Temperature $T$')
axs[0].set_ylabel(r'$C_4/C_2$')
axs[0].set_title("Cumulant Ratio: $C_4/C_2$")
axs[0].grid(True)
axs[0].legend()

# --- Second plot: C4/C2 - 3*C2 ---
axs[1].plot(T, R42_subtracted, 'go-', label=r'$C_4/C_2 - 3C_2$')
axs[1].axvline(T_CEP, color='red', linestyle='--', linewidth=1)
axs[1].plot(T_CEP, R_sub_CEP, 'ro', label='CEP')
axs[1].set_xlabel(r'Temperature $T$')
axs[1].set_ylabel(r'$C_4/C_2 - 3C_2$')
axs[1].set_title("Critical Ratio: $C_4/C_2 - 3C_2$")
axs[1].grid(True)
axs[1].legend()

plt.tight_layout()
plt.savefig("output/ising_sim.png")
