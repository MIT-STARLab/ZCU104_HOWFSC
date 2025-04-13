import matplotlib.pyplot as plt
import numpy as np

def read_data(filename):
    with open(filename, 'r') as file:
        # Read the first line to get the value of n
        n = int(file.readline().strip())
        real = np.zeros(n*n)
        imag = np.zeros(n*n)
        
        # Read the subsequent lines for real and imaginary values
        for idx, line in enumerate(file):
            r, i = map(float, line.strip().split())
            real[idx] = r
            imag[idx] = i
            
    return n, real, imag


# Read the data
n, software_real, software_imag = read_data('angular_spectrum_krnl/hls/csim/build/software_output.txt')
_, hardware_real, hardware_imag = read_data('angular_spectrum_krnl/hls/csim/build/hardware_output.txt')

# Reshape the data
software_real = software_real.reshape((n, n))
software_imag = software_imag.reshape((n, n))
hardware_real = hardware_real.reshape((n, n))
hardware_imag = hardware_imag.reshape((n, n))

# Compute differences
real_diff = np.abs(software_real - hardware_real)
imag_diff = np.abs(software_imag - hardware_imag)

# Calculate global min and max for real and imaginary components
real_min = min(np.min(software_real), np.min(hardware_real))
real_max = max(np.max(software_real), np.max(hardware_real))

imag_min = min(np.min(software_imag), np.min(hardware_imag))
imag_max = max(np.max(software_imag), np.max(hardware_imag))

# Plot the results
fig, axs = plt.subplots(3, 2, figsize=(15, 15))
titles = [
    "Real Component Hardware",     "Imag Component Hardware",
    "Real Component Software",     "Imag Component Software",
    "Difference (Real)",           "Difference (Imag)"
]

data = [
    hardware_real, hardware_imag, 
    software_real, software_imag,
    real_diff, imag_diff
]

# Create consistent colormap scaling
for ax, title, d in zip(axs.flat, titles, data):
    if "Difference" in title:
        cmap = "coolwarm"  # Different colormap for differences
        vmin, vmax = np.min(d), np.max(d)  # Use min/max for the difference
    else:
        cmap = "viridis"  # Consistent colormap for real and imaginary components
        if "Real" in title:
            vmin, vmax = real_min, real_max  # Same scale for real components
        else:
            vmin, vmax = imag_min, imag_max  # Same scale for imaginary components

    im = ax.imshow(d, extent=(0, n, 0, n), origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title)
    fig.colorbar(im, ax=ax)

plt.tight_layout()
plt.show()
