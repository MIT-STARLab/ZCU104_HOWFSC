import numpy as np
import time

def fft2d_timer(input_data, m, n, invert):
    input_data = np.random.uniform(100.0, 1000.0, (m, n)) + 1j * np.random.uniform(100.0, 1000.0, (m, n))

    start_time = time.time()
    if invert:
        output_data = np.fft.ifft2(input_data)
    else:
        output_data = np.fft.fft2(input_data)
    end_time = time.time()

    time_elapsed = end_time - start_time
    return time_elapsed

if __name__ == "__main__":
    n = 1024
    m = 1024
    trials = 10

    total_trials_time = 0
    for _ in range(trials):
        input_data = np.empty((m, n), dtype=np.complex128)
        total_trials_time += fft2d_timer(input_data, m, n, invert=False)

    avg = total_trials_time / trials

    print(f"Timing 2D FFT on matrix with size ({m}, {n}) for number of trials = {trials}")
    print(f"Total time = {total_trials_time}")
    print(f"Average time over one run = {avg}")