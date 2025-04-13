import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
data = pd.read_csv('fft_timing_results_mac.csv')

# Separate data for 1D and 2D FFT
data_1d = data[data['Type'] == '1D']
data_2d = data[data['Type'] == '2D']

# Plot the results
plt.figure(figsize=(10, 6))

plt.plot(data_1d['Data Size'], data_1d['Avg Time (s)'], label='1D FFT', marker='o')
plt.plot(data_2d['Data Size'], data_2d['Avg Time (s)'], label='2D FFT', marker='x')

plt.xlabel('Data Size')
plt.ylabel('Average Time (s)')
plt.title('FFT Timing Results')
plt.legend()
plt.grid(True)
plt.yscale('log')
plt.xscale('log')

# Save the plot as a PNG file
plt.savefig('fft_timing_results.png')

# Show the plot
plt.show()
