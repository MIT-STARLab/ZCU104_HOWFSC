#include "vadd_new.h"
#include <iostream>

int main() {
    Test1:
    unsigned int in1[] = {1, 2, 3, 4, 5, 6, 7};
    unsigned int in2[] = {1, 2, 3, 4, 5, 6, 7};
    unsigned int expected[] = {2, 4, 6, 8, 10, 12, 14};

    unsigned int result[7] = {0}; // Initialize result array

    // Call the vadd function
    vadd_new(in1, in2, result);

    // Check the result
    bool success = true;
    for (int i = 0; i < 7; ++i) {
        if (result[i] != expected[i]) {
            std::cout << "Test failed at index " << i << ": got " << result[i] << ", expected " << expected[i] << std::endl;
            success = false;
        }
    }

    if (success) {
        std::cout << "Test passed!" << std::endl;
    }

    return 0;
}
