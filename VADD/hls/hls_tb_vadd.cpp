/*
# Copyright © 2023 Advanced Micro Devices, Inc. All rights reserved.
# SPDX-License-Identifier: MIT
*/

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "vadd_include.h"


int main(int argc, char* argv[])
{

  data_t *ptr_inp_a = (data_t *) malloc(DATA_SIZE * sizeof(data_t));
  data_t *ptr_inp_b = (data_t *) malloc(DATA_SIZE * sizeof(data_t));
  data_t *ptr_res   = (data_t *) malloc(DATA_SIZE * sizeof(data_t));

  //setting input data
  FILE *fin1=fopen("vector_inputs.txt", "w");

  // Initialization, should only be called once.

  srand(41);
  for(int i = 0 ; i< DATA_SIZE; i++)
  {
    // to avoid overflows let us use short ints
    data_t r1 = (data_t) rand();
    data_t r2 = (data_t) rand();
    ptr_inp_a[i] = r1;
    ptr_inp_b[i] = r2;
    fprintf(fin1, "%10f\t%10f\n", r1, r2);
  }
  fclose(fin1);

  //call DUT
  krnl_vadd(ptr_inp_a, ptr_inp_b, ptr_res, DATA_SIZE);


  // Verify the result
  FILE *fin2=fopen("vector_out.txt", "w");
  int match = 0;
  for (int i = 0; i < DATA_SIZE; i++)
  {
    data_t host_result = ptr_inp_a[i] + ptr_inp_b[i];
    fprintf(fin2, "%10f\t%10f\n", host_result, ptr_res[i]);
    if (ptr_res[i] != host_result) {
        printf("MISMATCH ERROR at %d: expected %f got %f\n", i, host_result, ptr_res[i]);
        match = 1;
        //break;
    }
  }

  fclose(fin2);
  free(ptr_inp_a);
  free(ptr_inp_b);
  free(ptr_res);

  std::cout << "TEST " << (match ? "FAILED" : "PASSED") << std::endl;
  return match;

}