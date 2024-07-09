#include <iostream>
#include "../fft.hxx"
#include <fstream>


using namespace std;

const int m = 1024;
const int n = m;
const string dataFileName = "fft2d_tb_data.csv"; 
double tolerance = 1e-4;


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




//////////// Tests ////////////

/**
 * Testing FFT2D
 */
void test_fft2d(){
    cout << "Testing forward fft 2d" << endl;
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
        return;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int line_index = 0; line_index < m*n; ++line_index) {
        if (!getline(file, line)) {
            cerr << "Error reading line " << line_index + 2 << " from file!" << endl;
        }

        vector<string> tokens = split(line, ',');

        int row_index = line_index / m ;
        int col_index = line_index % m ;


        input_array[row_index][col_index] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
        expected_output_array[row_index][col_index] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
    }

    file.close();
    cout << "Loaded data successfully" << endl;

    //////////////////////   Run FFT    //////////////////////
    fft2d(input_array, 0); // 0 for forward


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
 * Testing InverseFFT2D
 */
void test_fft2d_inverse(){
    cout << "Testing inverse fft 2d" << endl;

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

        int row_index = line_index / m ;
        int col_index = line_index % m ;


        expected_output_array[row_index][col_index] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
        input_array[row_index][col_index] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
    }

    file.close();
    cout << "Loaded data successfully" << endl;

    //////////////////////   Run FFT    //////////////////////
    fft2d(input_array, 1); // 1 for inverse


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



void test_generator(){
    cout << "Testing 2d fft data generator" << endl;

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

    fft2d_random_data_generator(input_array, expected_output_array, 0); //0 for forward

    fft2d(input_array, 0);

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


/////////// Main ////////////
int main(){
    test_fft2d();
    test_fft2d_inverse();
    test_generator();
}

