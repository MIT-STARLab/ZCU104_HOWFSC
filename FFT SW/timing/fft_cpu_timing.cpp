#include <time.h>
#include <iostream>
#include <fstream>
#include "../src/fft.h"

/////// Timing Helpers ///////
#define TIMER(label)  timespec label; clock_gettime(CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0

using namespace std;

int main() {
    ofstream outputFile;
    outputFile.open("fft_timing_results_mac.csv");

    outputFile << "Type,Data Size,Avg Time (s)\n";

    cout << "*******      TIMING 1D FFT       *******\n" << endl;
    cout << "DATA SIZE     |     avg time over 10 rounds" << endl;
    cout << "-------------------------------------------" << endl;

    for (int nfft = 1; nfft < 12; nfft++) {
        int fft_length = 1 << nfft;

        cmpx_data_t* software_generated_input_data = new cmpx_data_t[fft_length];
        cmpx_data_t* software_generated_ouput_data = new cmpx_data_t[fft_length];
        bool invert = false;

        fft_random_data_generator(software_generated_input_data, software_generated_ouput_data, fft_length, invert); //no scale

        int rounds = 10;
        TIMER(TIME_START);
        for (int i = 0; i < rounds; i++) {
            fft(software_generated_input_data, fft_length, invert); //no scale
        }
        TIMER(TIME_END);
        
        double elapsedTime = ELAPSED(TIME_END, TIME_START) / rounds;
        cout << fft_length << "           |     " << scientific << elapsedTime << endl;
        outputFile << "1D," << fft_length << "," << elapsedTime << "\n";

        delete[] software_generated_input_data;
        delete[] software_generated_ouput_data;
    }

    cout << "\n*******      TIMING 2D FFT       *******" << endl;
    cout << "DATA SIZE     |     avg time over 10 rounds" << endl;
    cout << "-------------------------------------------" << endl;

    for (int nfft = 1; nfft < 12; nfft++) {
        int fft_length_n = 1 << nfft;
        int fft_length_m = 1 << nfft;

        cmpx_data_t* software_generated_input_data = new cmpx_data_t[fft_length_n * fft_length_m];
        cmpx_data_t* software_generated_ouput_data = new cmpx_data_t[fft_length_n * fft_length_m];
        bool invert = false;

        fft2d_random_data_generator(software_generated_input_data, software_generated_ouput_data, fft_length_m, fft_length_n, invert); //no scale

        int rounds = 10;
        TIMER(TIME_START);
        for (int i = 0; i < rounds; i++) {
            fft2d(software_generated_input_data, fft_length_m, fft_length_n, invert); //no scale
        }
        TIMER(TIME_END);
        
        double elapsedTime = ELAPSED(TIME_END, TIME_START) / rounds;
        cout << fft_length_n << "           |     " << scientific << elapsedTime << endl;
        outputFile << "2D," << fft_length_n << "," << elapsedTime << "\n";

        delete[] software_generated_input_data;
        delete[] software_generated_ouput_data;
    }

    outputFile.close();

    return 0;
}
