

/*
 * MIT STAR Lab
 * Written by Subhi in Jun 24
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

#include "fft_top.h"


using namespace std;






/**
 * Function to split a string by delimiter and store in vector
 */
vector<string> split(const string &s, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);

    while (getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}


/**
 * Checks if two complex numbers are equal within a tolarance
 */
bool approx_equal(const cmpx_data_t &a, const cmpx_data_t &b, double tolerance) {
    return (std::abs(a.real() - b.real()) <= tolerance) &&
           (std::abs(a.imag() - b.imag()) <= tolerance);
}




int main() {

    //////////   LOAD DATA   ///////////

    cmpx_data_t input_array[FFT_LENGTH];
    cmpx_data_t expected_output_array[FFT_LENGTH];
    cmpx_data_t actual_output_array[FFT_LENGTH];


    ifstream file("ftt_tb_data.csv");                             //File generated with a python script

    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        return 1;
    }

    string line;

    getline(file, line);                                          //Read header line; skip
    for (int i = 0; i < FFT_LENGTH; ++i) {

        if (!getline(file, line)) {
            cerr << "Error reading line " << i + 2 << " from file!" << endl;
            return 1;
        }

        // Split line into tokens
        // The first two cols contain the real and imaginary parts of the x inputs
        // The last two have the expected outpur
        vector<string> tokens = split(line, ',');

        input_array[i] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
        expected_output_array[i] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
    }

    file.close();
    cout << "Loaded data:" << endl;



    /////////////////  Call FFT function   //////////////
    bool ovflo;
    const bool FORWARD_FFT = true;                                   // Use true for forward FFT, false for inverse FFT
    fft_top(FORWARD_FFT, input_array, actual_output_array, &ovflo);



    /////////////////  Compare Results   //////////////    
    bool passed = true;
    double tolerance = 1e-3;

    for (int j = 0; j < FFT_LENGTH; j++) {
        if (!approx_equal(expected_output_array[j], actual_output_array[j], tolerance)) {
            cerr << "Error at index " << j 
                 << ": Expected (" << expected_output_array[j].real() << ", " << expected_output_array[j].imag() 
                 << "), Actual ("  <<   actual_output_array[j].real() << ", " <<   actual_output_array[j].imag() 
                 << ")" << endl;

            passed = false;
        }
    }

    if (passed) {
        cout << "All tests passed!" << endl;
        return 0;
    } else if (ovflo){
        cout << " (OVERFLOW!!!)" << endl;
        return 0;
    } else{
        cout << "Tests failed!" << endl;
        return 1;
    }
}
