/**
* Copyright (C) 2019-2021 Xilinx, Inc
*
* Licensed under the Apache License, Version 2.0 (the "License"). You may
* not use this file except in compliance with the License. A copy of the
* License is located at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
* License for the specific language governing permissions and limitations
* under the License.
*/
#include <sys/syscall.h>      /* Definition of SYS_* constants */
#include <unistd.h>
#include <stdio.h>
#include <time.h>
// #include <algorithm>
// #include <vector>
// #include <chrono>
#include "appaccel.h"x
#include "matrixvector.h"

// Define the matrix dimensions and vector size
#define MATRIX_ROWS 256 // Should be divisible by PARALLELROWS
#define MATRIX_COLS 256*20 //

#define TIMER(label)  timespec label; syscall(SYS_clock_gettime, CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0

void timestamp(){
  time_t x;
  time ( &x );
  printf ( "Timestamp: %s\n", ctime (&x));
}


const char* kernelName = "krnl_inner_product";
//char* kernelName;

inline data_t at(data_t *A,unsigned M, unsigned N, unsigned i,unsigned j, bool rowmajor = true)
{
  return rowmajor ? A[i*N+j] : A[j*M+i] ;
}

inline void set(data_t *A,unsigned M, unsigned N, unsigned i,unsigned j,data_t x, bool rowmajor = true)
{
  A[rowmajor ? i*N+j : j*M+i] = x; 
}


void printMatrix(data_t* A, unsigned M, unsigned N)
{
    for (unsigned i=0; i<M*N; i++)
    {
        printf("%.0f ",A[i]);
        if (i % N == N-1) printf("\n");
    }
}

void printVector(data_t* C, unsigned M)
{
    for (unsigned j=0; j<M; j++)
        printf("%.0f ",C[j]);
    printf("\n");
}

void matrixMultiply(data_t C[MATRIX_ROWS], data_t A[MATRIX_ROWS*MATRIX_COLS], data_t B[MATRIX_COLS], bool rowmajor = true)
{
	unsigned i,j,k;
	data_t acc;
	LOOP_M: for(i=0, k=0; i < MATRIX_ROWS; i++)
	{
		acc = 0;
		LOOP_N: for(j=0; j < MATRIX_COLS; j++, k++)
		{
			acc += at(A,MATRIX_ROWS,MATRIX_COLS,i,j,rowmajor) * B[j]; // A[k++] for a row-major matrix
		}
		C[i] = acc;
	}
}


int krnl_eval(cl::CommandQueue& queue, cl::Buffer& buffer_B, cl::Kernel& matrixvector, cl::Buffer& buffer_output){
    cl_int err;
    TIMER(startTime);

    // Copy input data to device global memory. Should be no copy if allignedAllocate was used
    err = queue.enqueueMigrateMemObjects({buffer_B}, 0 /* 0 means from host*/);
    if (err) {printf("Error %d migrating input\n",err); return 10;}
    queue.finish();
    
    TIMER(loadTime);

    err = queue.enqueueTask(matrixvector);
    if (err) {printf("Error %d enqueueTask\n",err); return 11;}
    queue.finish();

    TIMER(runTime);

    // Copy Result from Device Global Memory to Host Local Memory
    err = queue.enqueueMigrateMemObjects({buffer_output}, CL_MIGRATE_MEM_OBJECT_HOST);
    if (err) {printf("Error %d migrating output\n",err); return 12;}
	queue.finish();

    TIMER(doneTime);
  
    printf("Time total %f load %f run %f readout %f\n", ELAPSED(doneTime,startTime),ELAPSED(loadTime,startTime),ELAPSED(runTime,loadTime),ELAPSED(doneTime,runTime));
    return 0;
}


int main(int argc, char** argv) 
{

    printf("Timing M x V with size %d ROWS and %d COLS \n",MATRIX_ROWS,MATRIX_COLS);
    printf("%s %s\n",__DATE__,__TIME__);
	if (argc < 2) 
    {
    	printf("Usage: %s <xclbin file> [num iterations] [kernel name]\n",argv[0]);
        return 1;
    }

    const char* programName = argv[1];
    
    unsigned iterations = 1;

    if (argc > 2)
    {
    	iterations = atoi(argv[2]);
    	if (is_emulation() && iterations > 2) iterations = 2; //don't take all day
    }

    if (argc > 3)
    {
    	kernelName = argv[3];
    }

    unsigned i,j;
    cl_int err;
    cl::Context context;
    cl::CommandQueue queue;
    cl::Program program;

    printf("Programming with %s\n",programName);

    if (!programXCL(programName, program, context, queue)) return 2;

    //it's possible to load multiple kernels and run them in parallel
    cl::Kernel matrixvector(program, kernelName, &err);
    if (err) {printf("Error %d kernel %s not found\n",err,kernelName); return 3;}

    // Allocate matrices in Host Memory
    unsigned sizeA = sizeof(data_t)*MATRIX_ROWS*MATRIX_COLS;
    unsigned sizeB = sizeof(data_t)*MATRIX_COLS;
    unsigned sizeC = sizeof(data_t)*MATRIX_ROWS;

    data_t* A;
   
    data_t* C2 = (data_t*)alignedAllocate(sizeC);
    data_t y=1;


    cl::Buffer buffer_output(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_WRITE_ONLY, sizeC,NULL, &err);
    if (err) {printf("Error %d assigning buffer C size %d\n",err,sizeC); return 4;}
    cl::Buffer buffer_in1(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, sizeA,NULL, &err);
    if (err) {printf("Error %d assigning buffer A size %d\n",err,sizeA); return 5;}
    cl::Buffer buffer_in2(context, CL_MEM_ALLOC_HOST_PTR |  CL_MEM_READ_ONLY, sizeB,NULL, &err);
    if (err) {printf("Error %d assigning buffer B size %d\n",err,sizeB); return 6;}

    err = matrixvector.setArg(0, buffer_output); if (err) {printf("Error %d setting argument 0\n",err); return 7;}
    err = matrixvector.setArg(1, buffer_in1); if (err) {printf("Error %d setting argument 1\n",err); return 8;}
    err = matrixvector.setArg(2, buffer_in2); if (err) {printf("Error %d setting argument 2\n",err); return 9;}
    err = matrixvector.setArg(3, MATRIX_ROWS); if (err) {printf("Error %d setting argument 3\n",err); return 10;}
    err = matrixvector.setArg(4, MATRIX_COLS); if (err) {printf("Error %d setting argument 4\n",err); return 11;}

    data_t *C = (data_t*) queue.enqueueMapBuffer(buffer_output, CL_TRUE, CL_MAP_READ, 0, sizeC, NULL, NULL, &err);
    if (err) {printf("Error %d mapping buffer C size %d\n",err,sizeC); return 14;}

    A = (data_t*) queue.enqueueMapBuffer(buffer_in1, CL_TRUE, CL_MAP_WRITE, 0, sizeA, NULL, NULL, &err);

    if (err) {printf("Error %d mapping buffer A size %d\n",err,sizeA); return 15;}
    data_t *B = (data_t*) queue.enqueueMapBuffer(buffer_in2, CL_TRUE, CL_MAP_WRITE, 0, sizeB, NULL, NULL, &err);
    if (err) {printf("Error %d mapping buffer B size %d\n",err,sizeB); return 16;}
    queue.finish();
    // Create the test data
    for (i=0; i<MATRIX_ROWS; i++)
      for (j=0; j<MATRIX_COLS; j++)
        set(A,MATRIX_ROWS,MATRIX_COLS,i,j,y++);
    for (j=0; j<MATRIX_COLS; j++)
      B[j]=y++;

    // comparison
    matrixMultiply(C2, A, B);


    //printf("Migrate input memory\n");
    TIMER(startTime);
    //printf("Time of start\n");
    //timestamp();

    // Copy input data to device global memory. Should be no copy if allignedAllocate was used
    err = queue.enqueueMigrateMemObjects({buffer_in1, buffer_in2}, 0 /* 0 means from host*/);
    if (err) {printf("Error %d migrating input\n",err); return 10;}
    queue.finish();
    
    //printf("Start kernel\n");
    TIMER(loadTime);
    //printf("Time of migrate out finish\n");
    //timestamp();

    // Launch the Kernel
    // For HLS kernels global and local size is always (1,1,1). Always use enqueueTask()
    err = queue.enqueueTask(matrixvector);
    if (err) {printf("Error %d enqueueTask\n",err); return 11;}
    queue.finish();

    TIMER(runTime);
    //printf("Time of task finish\n");
    //timestamp();
    //printf("Migrate output memory\n");

    // Copy Result from Device Global Memory to Host Local Memory
    err = queue.enqueueMigrateMemObjects({buffer_output}, CL_MIGRATE_MEM_OBJECT_HOST);
    if (err) {printf("Error %d migrating output\n",err); return 12;}
	queue.finish();

    TIMER(doneTime);
    //printf("Time of all done\n");
    //timestamp();
  
    printf("Time total %f load %f run %f readout %f\n", ELAPSED(doneTime,startTime),ELAPSED(loadTime,startTime),ELAPSED(runTime,loadTime),ELAPSED(doneTime,runTime));
    printf("First result %f\n",C[0]);


    //printVector(C, MATRIX_ROWS);
    //printVector(C2, MATRIX_ROWS);
    printf("First two entries hardware and application: %f, %f\n",C[0],C2[0]);

    // compare
    for (unsigned i=0; i<MATRIX_ROWS; i++)
    {
        if (((C2[i] - C[i])/(C2[i])) > 1e-10)
        {
            printf("Application and kernel %f, %f\n",C2[i],C[i]);
            printf("Test failed at %d %f ~= %f\n", i ,C2[i], C[i]); 
            return 20;
        }
    }

    printf("TEST PASSED\n");

    printf("Going to run everything one more time without matrix movement! \n");

    for (j=0; j<MATRIX_COLS; j++)
      B[j]++;

    matrixMultiply(C2, A, B);

    err = krnl_eval(queue, buffer_in2, matrixvector, buffer_output);
    if (err) {printf("Error %d in krnl_eval\n",err); return 99;}

    // compare
    for (unsigned i=0; i<MATRIX_ROWS; i++)
    {
        if (((C2[i] - C[i])/(C2[i])) > 1e-10)
        {
            printf("Second test failed at %d %f ~= %f\n", i ,C2[i], C[i]); 
            return 20;
        }
    }

    printf("SECOND TEST PASSED\n");

    printf("Going to run 10 times\n");
    TIMER(multiStartTime);
    for (int i = 0; i<10; i++){
        // Copy input data to device global memory. Should be no copy if allignedAllocate was used
        err = queue.enqueueMigrateMemObjects({buffer_in2}, 0 /* 0 means from host*/);
        if (err) {printf("Error %d migrating input\n",err); return 10;}
        
        // Launch the Kernel
        // For HLS kernels global and local size is always (1,1,1). Always use enqueueTask()
        err = queue.enqueueTask(matrixvector);
        if (err) {printf("Error %d enqueueTask\n",err); return 11;}


        // Copy Result from Device Global Memory to Host Local Memory
        err = queue.enqueueMigrateMemObjects({buffer_output}, CL_MIGRATE_MEM_OBJECT_HOST);
        if (err) {printf("Error %d migrating output\n",err); return 12;}
        queue.finish();
    }
    TIMER(multiDoneTime);

    printf("Time for 10 is %f\n",ELAPSED(multiDoneTime,multiStartTime));
    
    // alignedFree(B);
    // alignedFree(C);
    alignedFree(C2);
    // alignedFree(R);


err = queue.enqueueUnmapMemObject(buffer_in1, A); if (err) {printf("Error %d unmapping buffer A\n",err);return 25;}



err = queue.enqueueUnmapMemObject(buffer_output, C); if (err) {printf("Error %d unmapping buffer C\n",err);return 24;}
err = queue.enqueueUnmapMemObject(buffer_in2, B); if (err) {printf("Error %d unmapping buffer B\n",err);return 26;}
queue.finish();
    return 0;
}