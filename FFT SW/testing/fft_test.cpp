
#include <iostream>
#include "../fft.hxx"
#include <fstream>



const string dataFileName = "fft_testing_data.csv";
const int FFT_LENGTH = 2048;

///////  Helpers ///////

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



///////  Tests ///////
/**
 * Testing FFT
 */
void test_fft(){
    cout << "Testing FFT" << endl;

    vector<cmpx_data_t> input_data(FFT_LENGTH);
    vector<cmpx_data_t> expected_output_data(FFT_LENGTH);


    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int i = 0; i < FFT_LENGTH; ++i) {
        if (!getline(file, line)) {
            cerr << "Error reading line " << i + 2 << " from file!" << endl;
        }

        vector<string> tokens = split(line, ',');

        input_data[i] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
        expected_output_data[i] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
    }

    file.close();
    cout << "Loaded data successfully" << endl;


    //////////////////////   Run FFT    //////////////////////
    fft(input_data, 0); // 0 for forward


    ////////////////////// Compare Results ///////////////////
    bool passed = true;
    double tolerance = 1e-2;

    for (int j = 0; j < FFT_LENGTH; ++j) {
        if (!approx_equal(expected_output_data[j], input_data[j], tolerance)) {
            cerr << "Error at index " << j << ": Expected (" << expected_output_data[j].real() << ", "
                 << expected_output_data[j].imag() << "), Actual (" << input_data[j].real() << ", "
                 << input_data[j].imag() << ")" << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "All FFT tests passed!" << endl;
    } else {
        cout << "FFT Tests failed!" << endl;
    }

}


/**
 * Testing InverseFFT
 */
void test_Inverse_fft(){
    cout << "Testing InverseFFT" << endl;

    vector<cmpx_data_t> input_data(FFT_LENGTH);
    vector<cmpx_data_t> expected_output_data(FFT_LENGTH);

    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int i = 0; i < FFT_LENGTH; ++i) {
        if (!getline(file, line)) {
            cerr << "Error reading line " << i + 2 << " from file!" << endl;
        }

        vector<string> tokens = split(line, ',');

        input_data[i] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
        expected_output_data[i] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
    }

    file.close();
    cout << "Loaded data successfully" << endl;


    //////////////////////   Run FFT    //////////////////////
    fft(input_data, 1); // 1 for inverse


    ////////////////////// Compare Results ///////////////////
    bool passed = true;
    double tolerance = 1e-2;

    for (int j = 0; j < FFT_LENGTH; ++j) {
        if (!approx_equal(expected_output_data[j], input_data[j], tolerance)) {
            cerr << "Error at index " << j << ": Expected (" << expected_output_data[j].real() << ", "
                 << expected_output_data[j].imag() << "), Actual (" << input_data[j].real() << ", "
                 << input_data[j].imag() << ")" << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "All FFT tests passed!" << endl;
    } else {
        cout << "FFT Tests failed!" << endl;
    }

}


/**
 * Testing Data Generator
 */
void test_data_generator(){

    cout << "Testing Data Generator" << endl;

    vector<cmpx_data_t> expected_input_data(FFT_LENGTH);
    vector<cmpx_data_t> expected_output_data(FFT_LENGTH);

    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int i = 0; i < FFT_LENGTH; ++i) {
        if (!getline(file, line)) {
            cerr << "Error reading line " << i + 2 << " from file!" << endl;
        }

        vector<string> tokens = split(line, ',');

        expected_input_data[i] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
        expected_output_data[i] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
    }

    file.close();
    cout << "Loaded data successfully" << endl;

    //////////////    Run Generator  ////////////////////
    vector<float> frequencies = {1, 4, 7};
    vector<float> amplitudes = {3, 1, 0.5};
    vector<float> time_range = {0, 0.01* FFT_LENGTH};

    vector<cmpx_data_t> input_data_generator(FFT_LENGTH, 0.);
    vector<cmpx_data_t> output_data_generator(FFT_LENGTH, 0.);

    fft_data_generator(input_data_generator,output_data_generator, frequencies, amplitudes, time_range, 0);
    
    //////////////    Compare Results  ////////////////////

    bool passed = true;
    double tolerance = 1e-2;

    cout << "Comparing input" << endl;
    for (int j = 0; j < FFT_LENGTH; ++j) {
        if (!approx_equal(expected_input_data[j], input_data_generator[j], tolerance)) {
            cerr << "Error at index " << j 
                 <<  ": Expected (" << expected_input_data[j].real() << ", " << expected_input_data[j].imag() 
                 << "), Actual (" << input_data_generator[j].real() << ", " << input_data_generator[j].imag() 
                 << ")" << endl;
            passed = false;
        }
    }

    cout << "Comparing Output" << endl;
    for (int j = 0; j < FFT_LENGTH; ++j) {
        if (!approx_equal(expected_output_data[j], output_data_generator[j], tolerance)) {
            cerr << "Error at index " << j 
                 <<  ": Expected (" << expected_output_data[j].real() << ", " << expected_output_data[j].imag() 
                 << "), Actual (" << output_data_generator[j].real() << ", " << output_data_generator[j].imag() 
                 << ")" << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "All FFT tests passed!" << endl;
    } else {
        cout << "FFT Tests failed!" << endl;
    }

}






int main(){

    test_fft();
    test_Inverse_fft();
    test_data_generator();

}

