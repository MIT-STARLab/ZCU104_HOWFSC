
#include <iostream>
#include <fstream>
#include "../src/fft.h"




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

    bool real_compare = (a.real()<tolerance && b.real()<tolerance) || (abs((a.real() - b.real())/max(abs(a.real()), abs(b.real()))) < tolerance);
    bool imag_compare = (a.imag()<tolerance && b.imag()<tolerance) || (abs((a.imag() - b.imag())/max(abs(a.imag()), abs(b.imag()))) < tolerance);

    if (!(real_compare && imag_compare)){
        cout << "relative error in real = " << abs((a.real() - b.real())/max(abs(a.real()), abs(b.real()))) <<endl;
        cout << "relative error in imag = " << abs((a.imag() - b.imag())/max(abs(a.imag()), abs(b.imag()))) <<endl;
    }


    return  real_compare && imag_compare;
}




///////  Tests 1D ///////
/**
 * Target fft_vector Forward
 */
void test_fft_vector(int data_length, string dataFileName, double tolerance, bool scale){

    vector<cmpx_data_t> input_data(data_length);
    vector<cmpx_data_t> expected_output_data(data_length);


    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int i = 0; i < data_length; ++i) {
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
    fft_vector(input_data, 0, scale); // 0 for forward


    ////////////////////// Compare Results ///////////////////
    bool passed = true;

    for (int j = 0; j < data_length; ++j) {
        if (!approx_equal(expected_output_data[j], input_data[j], tolerance)) {
            cerr << "Error at index " << j << ": Expected (" << expected_output_data[j].real() << ", "
                 << expected_output_data[j].imag() << "), Actual (" << input_data[j].real() << ", "
                 << input_data[j].imag() << ")" << endl;
        
            cout << "relative real error is " << std::abs(( input_data[j].real() - expected_output_data[j].real())/ expected_output_data[j].real()) << endl;
            cout << "relative imag error is " << std::abs(( input_data[j].imag() - expected_output_data[j].imag())/ expected_output_data[j].imag()) << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "PASSED:    fft_vector" << endl;
    } else {
        cout << "FAILED:    fft_vector" << endl;
    }

}


/**
 * Testing fft_vector Inverse
 */
void test_fft_vector_inverse(int data_length, string dataFileName, double tolerance, bool scale){

    vector<cmpx_data_t> input_data(data_length);
    vector<cmpx_data_t> expected_output_data(data_length);

    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int i = 0; i < data_length; ++i) {
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
    fft_vector(input_data, 1, scale); // 1 for inverse


    ////////////////////// Compare Results ///////////////////
    bool passed = true;

    for (int j = 0; j < data_length; ++j) {
        if (!approx_equal(expected_output_data[j], input_data[j], tolerance)) {
            cerr << "Error at index " << j << ": Expected (" << expected_output_data[j].real() << ", "
                 << expected_output_data[j].imag() << "), Actual (" << input_data[j].real() << ", "
                 << input_data[j].imag() << ")" << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "PASSED:    fft_vector inverse" << endl;
    } else {
        cout << "FAILED:    fft_vector inverse" << endl;
    }

}


/**
 * Target fft Forward
 */
void test_fft(int data_length, string dataFileName, double tolerance, bool scale){

    cmpx_data_t input_data[data_length];
    cmpx_data_t expected_output_data[data_length];


    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int i = 0; i < data_length; ++i) {
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
    fft(input_data,data_length, 0, scale); // 0 for forward


    ////////////////////// Compare Results ///////////////////
    bool passed = true;

    for (int j = 0; j < data_length; ++j) {
        if (!approx_equal(expected_output_data[j], input_data[j], tolerance)) {
            cerr << "Error at index " << j << ": Expected (" << expected_output_data[j].real() << ", "
                 << expected_output_data[j].imag() << "), Actual (" << input_data[j].real() << ", "
                 << input_data[j].imag() << ")" << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "PASSED:    fft" << endl;
    } else {
        cout << "FAILED:    fft" << endl;
    }

}


/**
 * Testing fft_vector Inverse
 */
void test_fft_inverse(int data_length, string dataFileName, double tolerance, bool scale){

    cmpx_data_t input_data[data_length];
    cmpx_data_t expected_output_data[data_length];

    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int i = 0; i < data_length; ++i) {
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
    fft(input_data, data_length, 1, scale); // 1 for inverse


    ////////////////////// Compare Results ///////////////////
    bool passed = true;

    for (int j = 0; j < data_length; ++j) {
        if (!approx_equal(expected_output_data[j], input_data[j], tolerance)) {
            cerr << "Error at index " << j << ": Expected (" << expected_output_data[j].real() << ", "
                 << expected_output_data[j].imag() << "), Actual (" << input_data[j].real() << ", "
                 << input_data[j].imag() << ")" << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "PASSED:    fft" << endl;
    } else {
        cout << "FAILED:    fft inverse" << endl;
    }

}




/**
 * Testing Data Generator vector
 */
void test_data_generator_vector(int data_length, string dataFileName, double tolerance, bool scale){

    cout << "Testing Data Generator Vector" << endl;

    vector<cmpx_data_t> expected_input_data(data_length);
    vector<cmpx_data_t> expected_output_data(data_length);

    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int i = 0; i < data_length; ++i) {
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
    float sampling_rate = 0.01;
    std::vector<float> frequencies = {0.97, 10.44 , 50.44}; // {20 * fund_freq, 250 * fund_freq,  1000 * fund_freq};
    std::vector<float> amplitudes = {203,124,34};
    std::vector<float> time_range = {0, sampling_rate * data_length};

    vector<cmpx_data_t> input_data_generator(data_length, 0.);
    vector<cmpx_data_t> output_data_generator(data_length, 0.);

    fft_data_generator_vector(input_data_generator,output_data_generator, frequencies, amplitudes, time_range, 0, scale);
    
    //////////////    Compare Results  ////////////////////

    bool passed = true;

    cout << "Comparing input" << endl;
    for (int j = 0; j < data_length; ++j) {
        if (!approx_equal(expected_input_data[j], input_data_generator[j], tolerance)) {
            cerr << "Error at index " << j 
                 <<  ": Expected (" << expected_input_data[j].real() << ", " << expected_input_data[j].imag() 
                 << "), Actual (" << input_data_generator[j].real() << ", " << input_data_generator[j].imag() 
                 << ")" << endl;
            passed = false;
        }
    }

    cout << "Comparing Output" << endl;
    for (int j = 0; j < data_length; ++j) {
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



/**
 * Testing Data Generator
 */
void test_data_generator(int data_length, string dataFileName, double tolerance, bool scale){

    cout << "Testing Data Generator" << endl;

    cmpx_data_t expected_input_data[data_length];
    cmpx_data_t expected_output_data[data_length];

    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int i = 0; i < data_length; ++i) {
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
    float sampling_rate = 0.01;
    std::vector<float> frequencies = {0.97, 10.44 , 50.44}; // {20 * fund_freq, 250 * fund_freq,  1000 * fund_freq};
    std::vector<float> amplitudes = {203,124,34};
    std::vector<float> time_range = {0, sampling_rate * data_length};

    cmpx_data_t input_data_generator[data_length];
    cmpx_data_t output_data_generator[data_length];

    fft_data_generator(input_data_generator,output_data_generator,data_length , frequencies, amplitudes, time_range, 0, scale);
    
    //////////////    Compare Results  ////////////////////

    bool passed = true;

    cout << "Comparing input" << endl;
    for (int j = 0; j < data_length; ++j) {
        if (!approx_equal(expected_input_data[j], input_data_generator[j], tolerance)) {
            cerr << "Error at index " << j 
                 <<  ": Expected (" << expected_input_data[j].real() << ", " << expected_input_data[j].imag() 
                 << "), Actual (" << input_data_generator[j].real() << ", " << input_data_generator[j].imag() 
                 << ")" << endl;
            passed = false;
        }
    }

    cout << "Comparing Output" << endl;
    for (int j = 0; j < data_length; ++j) {
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



/**
 * Testing fft2d_test_generator_vector
 */
void fft_test_random_generator(int n, double tolerance, bool scale){
    cout << "Testing 2d fft data generator vector" << endl;

    cmpx_data_t input_array[n];
    cmpx_data_t expected_output_array[n];

    fft_random_data_generator(input_array, expected_output_array,n, 0, scale); //0 for forward

    fft(input_array,n, 0, scale);

    bool passed = true;

    for (int i =0; i<n; i++){
        if (!approx_equal(expected_output_array[i], input_array[i], tolerance)){
            cerr << "Error at location (" << i
            << ": Expected (" << expected_output_array[i].real() << ", " << expected_output_array[i].imag() 
            << "), Actual (" << input_array[i].real() << ", " << input_array[i].imag() 
            << ")" << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "All FFT data generator tests passed!" << endl;
    } else {
        cout << "FFT data generator Tests failed!" << endl;
    }

}




///////  Tests 2D ///////

/**
 * Testing FFT2D vector
 */
void test_fft2d_vector(int n, int m, string dataFileName, double tolerance, bool scale){
    cout << "Testing forward fft 2d vector" << endl;
    vector<vector<cmpx_data_t>> input_array;
    vector<vector<cmpx_data_t>> expected_output_array;

    for (int i = 0; i< m; i++){
        vector<cmpx_data_t> v1;
        input_array.push_back(v1);
        vector<cmpx_data_t> v2;
        expected_output_array.push_back(v2);

        for (int j= 0; j< n; j++){
            input_array[i].push_back(0);
            expected_output_array[i].push_back(0);
        }

    }
    cout << "Array allocated." << endl;
    //load testing data

    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
        return;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int line_index = 0; line_index < m*n; ++line_index) {
        if (!getline(file, line)) {
            cerr << "Error reading line " << line_index + 2 << " from file!" << endl;
        }

        vector<string> tokens = split(line, ',');

        int row_index = line_index / n ;
        int col_index = line_index % n ;

        input_array[row_index][col_index] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
        expected_output_array[row_index][col_index] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
    }

    file.close();
    cout << "Loaded data successfully" << endl;

    //////////////////////   Run FFT    //////////////////////
    fft2d_vector(input_array, 0, scale); // 0 for forward


    ////////////////////// Compare Results ///////////////////
    bool passed = true;
    for (int i =0; i<m; i++){
        for (int j =0; j<n; j++){
            if (!approx_equal(expected_output_array[i][j], input_array[i][j], tolerance)){
                cerr << "Error at location (" << i << ", " << j << ")"
                << ": Expected (" << expected_output_array[i][j].real() << ", " << expected_output_array[i][j].imag() 
                << "), Actual (" << input_array[i][j].real() << ", " << input_array[i][j].imag() 
                << ")" << endl;
                passed = false;
            }
        }
    }

    if (passed) {
        cout << "All 2D FFT tests passed!" << endl;
    } else {
        cout << "2D FFT Tests failed!" << endl;
    } 

}


/**
 * Testing InverseFFT2D vector
 */
void test_fft2d_inverse_vector(int n, int m, string dataFileName,  double tolerance, bool scale){
    cout << "Testing inverse fft 2d vector" << endl;

    vector<vector<cmpx_data_t>> input_array;
    vector<vector<cmpx_data_t>> expected_output_array;

    for (int i = 0; i< m; i++){
        vector<cmpx_data_t> v1;
        input_array.push_back(v1);
        vector<cmpx_data_t> v2;
        expected_output_array.push_back(v2);

        for (int j= 0; j< n; j++){
            input_array[i].push_back(0);
            expected_output_array[i].push_back(0);
        }

    }


    //load testing data

    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int line_index = 0; line_index < m*n; ++line_index) {
        if (!getline(file, line)) {
            cerr << "Error reading line " << line_index + 2 << " from file!" << endl;
        }

        vector<string> tokens = split(line, ',');

        int row_index = line_index / n ;
        int col_index = line_index % n ;


        expected_output_array[row_index][col_index] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
        input_array[row_index][col_index] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
    }

    file.close();
    cout << "Loaded data successfully" << endl;

    //////////////////////   Run FFT    //////////////////////
    fft2d_vector(input_array, 1, scale); // 1 for inverse


    ////////////////////// Compare Results ///////////////////
    bool passed = true;

    for (int i =0; i<m; i++){
        for (int j =0; j<n; j++){
            if (!approx_equal(expected_output_array[i][j], input_array[i][j], tolerance)){
                cerr << "Error at location (" << i << ", " << j << ")"
                << ": Expected (" << expected_output_array[i][j].real() << ", " << expected_output_array[i][j].imag() 
                << "), Actual (" << input_array[i][j].real() << ", " << input_array[i][j].imag() 
                << ")" << endl;
                passed = false;
            }
        }
    }

    if (passed) {
        cout << "All inverse 2D FFT tests passed!" << endl;
    } else {
        cout << "FFT inverse 2D Tests failed!" << endl;
    }

}



/**
 * Testing InverseFFT2D vector
 */
void fft2d_test_generator_vector(int n, int m, string dataFileName, double tolerance, bool scale){
    cout << "Testing 2d fft data generator vector" << endl;

    vector<vector<cmpx_data_t>> input_array;
    vector<vector<cmpx_data_t>> expected_output_array;
    for (int i = 0; i< m; i++){
        vector<cmpx_data_t> v1;
        input_array.push_back(v1);
        vector<cmpx_data_t> v2;
        expected_output_array.push_back(v2);

        for (int j= 0; j< n; j++){
            input_array[i].push_back(0);
            expected_output_array[i].push_back(0);
        }
    }

    fft2d_random_data_generator_vector(input_array, expected_output_array, 0, scale); //0 for forward

    fft2d_vector(input_array, 0, scale);

    bool passed = true;

    for (int i =0; i<m; i++){
        for (int j =0; j<n; j++){
            if (!approx_equal(expected_output_array[i][j], input_array[i][j], tolerance)){
                cerr << "Error at location (" << i << ", " << j << ")"
                << ": Expected (" << expected_output_array[i][j].real() << ", " << expected_output_array[i][j].imag() 
                << "), Actual (" << input_array[i][j].real() << ", " << input_array[i][j].imag() 
                << ")" << endl;
                passed = false;
            }
        }
    }

    if (passed) {
        cout << "All 2D FFT data generator tests passed!" << endl;
    } else {
        cout << "2D FFT data generator Tests failed!" << endl;
    }

}


/**
 * Testing FFT2D
 */
void test_fft2d(int n, int m, string dataFileName, double tolerance, bool scale){
    cout << "Testing forward fft 2d vector" << endl;
    cmpx_data_t input_array[n*m];
    cmpx_data_t expected_output_array[n*m];

    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
        return;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int line_index = 0; line_index < m*n; ++line_index) {
        if (!getline(file, line)) {
            cerr << "Error reading line " << line_index + 2 << " from file!" << endl;
        }

        vector<string> tokens = split(line, ',');

        input_array[line_index] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
        expected_output_array[line_index] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
    }

    file.close();
    cout << "Loaded data successfully" << endl;

    //////////////////////   Run FFT    //////////////////////
    fft2d(input_array,m,n, 0, scale); // 0 for forward


    ////////////////////// Compare Results ///////////////////
    bool passed = true;
    for (int i =0; i<m*n; i++){
        if (!approx_equal(expected_output_array[i], input_array[i], tolerance)){
            cerr << "Error at location (" << i 
            << ": Expected (" << expected_output_array[i].real() << ", " << expected_output_array[i].imag() 
            << "), Actual (" << input_array[i].real() << ", " << input_array[i].imag() 
            << ")" << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "All 2D FFT tests passed!" << endl;
    } else {
        cout << "2D FFT Tests failed!" << endl;
    } 

}


/**
 * Testing FFT2D
 */
void test_fft2d_inverse(int n, int m, string dataFileName, double tolerance, bool scale){
    cout << "Testing inverse fft 2d vector" << endl;
    cmpx_data_t input_array[n*m];
    cmpx_data_t expected_output_array[n*m];

    ///////////////////// LOADING TESTING DATA /////////////////////
    ifstream file(dataFileName);   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
        return;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int line_index = 0; line_index < m*n; ++line_index) {
        if (!getline(file, line)) {
            cerr << "Error reading line " << line_index + 2 << " from file!" << endl;
        }

        vector<string> tokens = split(line, ',');

        expected_output_array[line_index] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
        input_array[line_index] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
    }

    file.close();
    cout << "Loaded data successfully" << endl;

    //////////////////////   Run FFT    //////////////////////
    fft2d(input_array,m,n, 1, scale); // 0 for forward


    ////////////////////// Compare Results ///////////////////
    bool passed = true;
    for (int i =0; i<m*n; i++){
        if (!approx_equal(expected_output_array[i], input_array[i], tolerance)){
            cerr << "Error at location (" << i 
            << ": Expected (" << expected_output_array[i].real() << ", " << expected_output_array[i].imag() 
            << "), Actual (" << input_array[i].real() << ", " << input_array[i].imag() 
            << ")" << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "All 2D FFT tests passed!" << endl;
    } else {
        cout << "2D FFT Tests failed!" << endl;
    } 

}




/**
 * Testing fft2d_test_generator_vector
 */
void fft2d_test_generator_vector(int n, int m, double tolerance, bool scale){
    cout << "Testing 2d fft data generator vector" << endl;

    vector<vector<cmpx_data_t>> input_array;
    vector<vector<cmpx_data_t>> expected_output_array;
    for (int i = 0; i< m; i++){
        vector<cmpx_data_t> v1;
        input_array.push_back(v1);
        vector<cmpx_data_t> v2;
        expected_output_array.push_back(v2);

        for (int j= 0; j< n; j++){
            input_array[i].push_back(0);
            expected_output_array[i].push_back(0);
        }
    }

    fft2d_random_data_generator_vector(input_array, expected_output_array, 0, scale); //0 for forward

    fft2d_vector(input_array, 0, scale);

    bool passed = true;

    for (int i =0; i<m; i++){
        for (int j =0; j<n; j++){
            if (!approx_equal(expected_output_array[i][j], input_array[i][j], tolerance)){
                cerr << "Error at location (" << i << ", " << j << ")"
                << ": Expected (" << expected_output_array[i][j].real() << ", " << expected_output_array[i][j].imag() 
                << "), Actual (" << input_array[i][j].real() << ", " << input_array[i][j].imag() 
                << ")" << endl;
                passed = false;
            }
        }
    }

    if (passed) {
        cout << "All 2D FFT data generator tests passed!" << endl;
    } else {
        cout << "2D FFT data generator Tests failed!" << endl;
    }

}



/**
 * Testing fft2d_test_generator_vector
 */
void fft2d_test_generator(int n, int m, double tolerance, bool scale){
    cout << "Testing 2d fft data generator vector" << endl;

    cmpx_data_t input_array[n*m];
    cmpx_data_t expected_output_array[n*m];

    fft2d_random_data_generator(input_array, expected_output_array,m,n, 0, scale); //0 for forward

    fft2d(input_array,m,n, 0, scale);

    bool passed = true;

    for (int i =0; i<m*n; i++){
        if (!approx_equal(expected_output_array[i], input_array[i], tolerance)){
            cerr << "Error at location (" << i
            << ": Expected (" << expected_output_array[i].real() << ", " << expected_output_array[i].imag() 
            << "), Actual (" << input_array[i].real() << ", " << input_array[i].imag() 
            << ")" << endl;
            passed = false;
        }
    }

    if (passed) {
        cout << "All 2D FFT data generator tests passed!" << endl;
    } else {
        cout << "2D FFT data generator Tests failed!" << endl;
    }

}




///////  RUN TESTS ///////

using namespace std;


int n = 512;
int m = 128;
string file1 = "file1.csv";
string file2 = "file2.csv";
double tolerance = 5e-2;
bool scale = true; // we are  comparing against numpy and numpy scales by dimentions when inverse fft

/////////// Main ////////////
int main(){

    cout << "/nTESTING 1D FFT RELATED FUNCTIONS" << endl;
    cout << "********************************" << endl;
    test_fft_vector(n,file1, tolerance, scale);
    test_fft_vector_inverse(n,file1, tolerance, scale);
    test_fft(n,file1, tolerance, scale);
    test_fft_inverse(n,file1, tolerance, scale);
    test_data_generator_vector(n, file1, tolerance, scale);
    test_data_generator(n, file1, tolerance, scale);
    fft_test_random_generator(n, tolerance, scale);


    cout << "\nTESTING 2D FFT RELATED FUNCTIONS" << endl;
    cout << "********************************" << endl;
    test_fft2d_vector(n,m, file2, tolerance, scale);
    test_fft2d_inverse_vector(n,m, file2, tolerance, scale);
    test_fft2d(n,m, file2, tolerance, scale);
    test_fft2d_inverse(n,m, file2, tolerance, scale);
    fft2d_test_generator_vector(n,m, tolerance, scale);
    fft2d_test_generator(n,m, tolerance, scale);

}


