/*
 Copyright © 2023 Advanced Micro Devices, Inc. All rights reserved.
 SPDX-License-Identifier: MIT
*/

/*
  Modified by Subhi, Jul 2024
*/

#ifndef _H_VADD_INCLUDE_H_
#define _H_VADD_INCLUDE_H_

// Includes
#include <hls_stream.h>

#define DATA_SIZE 4096

typedef double data_t;

//avoid memory contention
#define MEMORY_DWIDTH 64
#define MEMORY_DWIDTH_BYTES (64/8)

extern "C" void krnl_vadd(const data_t* in1, const data_t* in2, data_t* out, int size);

#endif
