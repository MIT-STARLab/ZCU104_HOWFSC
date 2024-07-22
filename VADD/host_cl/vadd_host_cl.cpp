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
#include "appaccel.h"
#include "vadd_include.h"



#define DATA_SIZE 4096
const char* kernelName = "krnl_vadd";
typedef double data_t;


#define TIMER(label)  timespec label; syscall(SYS_clock_gettime, CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0


void timestamp(){
  time_t x;
  time ( &x );
  printf ( "Timestamp: %s\n", ctime (&x));
}






int main(int argc, char** argv) 
{

    printf("Timing VADD with size %d \n",DATA_SIZE);
    printf("%s %s\n",__DATE__,__TIME__);
	if (argc != 2) 
    {
    	printf("Usage: %s <xclbin file>\n",argv[0]);
        return 1;
    }

    const char* programName = argv[1];
    

    unsigned i,j;
    cl_int err;
    cl::Context context;
    cl::CommandQueue queue;
    cl::Program program;

    printf("Programming with %s\n",programName);

    if (!programXCL(programName, program, context, queue)) return 2;

    //it's possible to load multiple kernels and run them in parallel
    cl::Kernel vadd(program, kernelName, &err);
    if (err) {printf("Error %d kernel %s not found\n",err,kernelName); return 3;}

    // Allocate matrices in Host Memory
    unsigned sizeA = sizeof(data_t)*DATA_SIZE;
    unsigned sizeB = sizeof(data_t)*DATA_SIZE;
    unsigned sizeC = sizeof(data_t)*DATA_SIZE;


    cl::Buffer buffer_output(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_WRITE_ONLY, sizeC,NULL, &err);
    if (err) {printf("Error %d assigning buffer C size %d\n",err,sizeC); return 4;}
    cl::Buffer buffer_in1(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, sizeA,NULL, &err);
    if (err) {printf("Error %d assigning buffer A size %d\n",err,sizeA); return 5;}
    cl::Buffer buffer_in2(context, CL_MEM_ALLOC_HOST_PTR |  CL_MEM_READ_ONLY, sizeB,NULL, &err);
    if (err) {printf("Error %d assigning buffer B size %d\n",err,sizeB); return 6;}

    err = vadd.setArg(0, buffer_in1); if (err) {printf("Error %d setting argument 0\n",err); return 7;}
    err = vadd.setArg(1, buffer_in2); if (err) {printf("Error %d setting argument 1\n",err); return 8;}
    err = vadd.setArg(2, buffer_output); if (err) {printf("Error %d setting argument 2\n",err); return 9;}
    err = vadd.setArg(3, DATA_SIZE); if (err) {printf("Error %d setting argument 3\n",err); return 10;}

    data_t *C = (data_t*) queue.enqueueMapBuffer(buffer_output, CL_TRUE, CL_MAP_READ, 0, sizeC, NULL, NULL, &err);
    if (err) {printf("Error %d mapping buffer C size %d\n",err,sizeC); return 14;}
    data_t* A = (data_t*) queue.enqueueMapBuffer(buffer_in1, CL_TRUE, CL_MAP_WRITE, 0, sizeA, NULL, NULL, &err);
    if (err) {printf("Error %d mapping buffer A size %d\n",err,sizeA); return 15;}
    data_t *B = (data_t*) queue.enqueueMapBuffer(buffer_in2, CL_TRUE, CL_MAP_WRITE, 0, sizeB, NULL, NULL, &err);
    if (err) {printf("Error %d mapping buffer B size %d\n",err,sizeB); return 16;}
    queue.finish();


    // Create the test data
    data_t* C2 = (data_t*)alignedAllocate(sizeC);
    data_t y=1;
    for (i=0; i<DATA_SIZE; i++){
        B[i]=y;
        A[i]=y;
        C2[i] = A[i] + B[i];
    }

    TIMER(startTime);

    printf("Migrate input memory\n"); 
    err = queue.enqueueMigrateMemObjects({buffer_in1, buffer_in2}, 0 /* 0 means from host*/);       // Copy input data to device global memory. Should be no copy if allignedAllocate was used
    if (err) {printf("Error %d migrating input\n",err); return 10;}
    queue.finish();
    

    TIMER(loadTime);

    // Launch the Kernel
    // For HLS kernels global and local size is always (1,1,1). Always use enqueueTask()
    printf("Start kernel\n");
    err = queue.enqueueTask(vadd);
    if (err) {printf("Error %d enqueueTask\n",err); return 11;}
    queue.finish();

    TIMER(runTime);

    printf("Migrate output memory\n");
    err = queue.enqueueMigrateMemObjects({buffer_output}, CL_MIGRATE_MEM_OBJECT_HOST);      // Copy Result from Device Global Memory to Host Local Memory
    if (err) {printf("Error %d migrating output\n",err); return 12;}
	queue.finish();

    TIMER(doneTime);


    printf("Time total %f load %f run %f readout %f\n", ELAPSED(doneTime,startTime),ELAPSED(loadTime,startTime),ELAPSED(runTime,loadTime),ELAPSED(doneTime,runTime));
    printf("First result %f\n",C[0]);
    printf("First two entries hardware and application: %f, %f\n",C[0],C2[0]);

    // compare
    for (unsigned i=0; i<DATA_SIZE; i++)
    {
        if (((C2[i] - C[i])/(C2[i])) > 1e-10)
        {
            printf("Application and kernel %f, %f\n",C2[i],C[i]);
            printf("Test failed at %d %f ~= %f\n", i ,C2[i], C[i]); 
            return 20;
        }
    }

    printf("TEST PASSED\n");

    alignedFree(C2);
    err = queue.enqueueUnmapMemObject(buffer_in1, A); if (err) {printf("Error %d unmapping buffer A\n",err);return 25;}
    err = queue.enqueueUnmapMemObject(buffer_output, C); if (err) {printf("Error %d unmapping buffer C\n",err);return 24;}
    err = queue.enqueueUnmapMemObject(buffer_in2, B); if (err) {printf("Error %d unmapping buffer B\n",err);return 26;}
    queue.finish();
    return 0;
}