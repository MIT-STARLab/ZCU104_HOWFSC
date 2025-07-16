### Script to plot timing results for the QR givens decompositions

import matplotlib.pyplot as plt
import numpy as np

# Data
arr_size        = [128, 256, 512, 1024, 2048]
fpga_hw         = [1.40E-02, 1.18E-01, 1.02E+00, 8.13E+00, 6.55E+01]
zcu_cortex_a53  = [3.92E-01, 3.14E+00, 2.53E+01, 2.03E+02, 1.64E+03]
t1040           = [1.24E-02, 1.86E-01, 5.66E+00, 1.17E+02, 1.58E+03]


# Convert to log base 2
log_n = np.log2(arr_size)
log_time = np.log(fpga_hw)  # natural log, or use np.log2(fpga_time_sec) if you want base-2

# Linear fit: log(time) = a * log2(n) + b
a, b = np.polyfit(log_n, log_time, 1)


print("Angular Spectrum: log(time) = a * log2(n) + b")
print(f"a (slope): {a:.3f}")
print(f"b (intercept): {b:.3f}")

print(np.array(zcu_cortex_a53)/np.array(fpga_hw))


# Plot
plt.figure(figsize=(12, 8))  # Increase figure size for better readability

# plt.plot(arr_size, t1040, 'o-', label="T1040 — Numerical Recipies in C",color='blue', markersize=10, linewidth=2)
plt.plot(arr_size, fpga_hw, 's--', label="ZCU104 FPGA — Givens Rotations", color='orange', markersize=10, linewidth=2)
plt.plot(arr_size, zcu_cortex_a53, 'x-', label="ZCU104 MPCPU Cortex A53 — Givens Rotations", color='red', markersize=8, linewidth=2)

# Plot settings
plt.xscale('log', base=2)  # Log scale for matrix size
plt.yscale('log')          # Log scale for time
plt.xlabel("Array Size", fontsize=20)
plt.ylabel("Time [ms]", fontsize=20)
plt.title("QR Timing Analysis", fontsize=24, fontweight='bold')

plt.tick_params(axis='both', which='major', labelsize=16, width=2, length=6, color='black', labelcolor='black')
plt.tick_params(axis='both', which='minor', labelsize=16, width=2, length=4, color='black', labelcolor='black')

plt.legend(loc='best', fontsize=20)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()  # Adjust plot to fit labels and title

# Save and show
plt.savefig('./qr_timing_log.png', format='png', dpi=300)
plt.show()
