
#include <CL/cl.h>
#include <cstddef>
#include <sys/syscall.h>      /* Definition of SYS_* constants */
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include "appaccel.h"
#include "vadd_new.h"
#include <cmath>


using namespace std;


#define TIMER(label)  timespec label; syscall(SYS_clock_gettime, CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0

void timestamp(){
  time_t x;
  time ( &x );
  printf ( "Timestamp: %s\n", ctime (&x));
}


#define DATA_SIZE 7


const char* kernelName = "vadd_new";



int main(int argc, char** argv){

    printf("Timing Vector adder kernel");

    if (argc < 2) 
    {
    	printf("Usage: %s <xclbin file> [num iterations] [kernel name]\n",argv[0]);
        return 1;
    }

    const char* programName = argv[1];
    

    cl_int err;
    cl::Context context;
    cl::CommandQueue queue;
    cl::Program program;


    printf("Programming with %s\n",programName);

    if (!programXCL(programName, program, context, queue)) return 2;

    cl::Kernel vadd_new(program, kernelName, &err);
    if (err){printf("Error %d kernel %s\n", err, kernelName);}
    
    unsigned data_vector_size = sizeof(int) * DATA_SIZE;


    // Memory allocation
    cl::Buffer buffer_input1(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, data_vector_size, NULL, &err);
    if (err) {printf("Error %d assigning input buffer 1 size %d\n",err,data_vector_size); return 5;}

    cl::Buffer buffer_input2(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, data_vector_size, NULL, &err);
    if (err) {printf("Error %d assigning input buffer 2 size %d\n",err,data_vector_size); return 6;}

    cl::Buffer buffer_output(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_WRITE_ONLY, data_vector_size, NULL, &err);
    if (err) {printf("Error %d assigning output buffer size %d\n",err,data_vector_size); return 4;}


    err = vadd_new.setArg(0, buffer_input1); if (err) {printf("Error %d setting argument 0\n",err); return 7;}
    err = vadd_new.setArg(1, buffer_input2); if (err) {printf("Error %d setting argument 1\n",err); return 8;}
    err = vadd_new.setArg(2, buffer_output); if (err) {printf("Error %d setting argument 2\n",err); return 9;}

    int* x = (int*) queue.enqueueMapBuffer(buffer_input1, CL_TRUE, CL_MAP_WRITE, 0, data_vector_size, NULL, NULL, &err);
    if (err) {printf("Error %d mapping buffer input 1 size %d\n",err,data_vector_size); return 10;}

    int* y = (int*) queue.enqueueMapBuffer(buffer_input2, CL_TRUE, CL_MAP_WRITE, 0, data_vector_size, NULL, NULL, &err);
    if (err) {printf("Error %d mapping buffer input 2 size %d\n",err,data_vector_size); return 11;}

    int* z = (int*) queue.enqueueMapBuffer(buffer_output, CL_TRUE, CL_MAP_READ, 0, data_vector_size, NULL, NULL, &err);
    if (err) {printf("Error %d mapping buffer output size %d\n",err,data_vector_size); return 12;}

    queue.finish();

    int * z2 = (int*) alignedAllocate(data_vector_size);

    // create testing data
    for (int i = 0; i< DATA_SIZE; i++){
        x[i] = i+1;
        y[i] = i+1;
        z2[i] = 2 * (i+1);
    }


    TIMER(startTime);

    err = queue.enqueueMigrateMemObjects({buffer_input1, buffer_input2}, 0);    // 0 means from host  
    if (err) {printf("Error %d migrating input\n",err); return 13;}
    queue.finish();


    TIMER(loadTime);

    err = queue.enqueueTask(vadd_new);
    if (err) {printf("Error %d enqueueTask\n",err); return 14;}
    queue.finish();


    TIMER(runTime);

    err = queue.enqueueMigrateMemObjects({buffer_output}, CL_MIGRATE_MEM_OBJECT_HOST);
    if (err) {printf("Error %d migrating output\n",err); return 15;}
	queue.finish();

    TIMER(doneTime);

    printf("Time total %f load %f run %f readout %f\n", ELAPSED(doneTime,startTime),ELAPSED(loadTime,startTime),ELAPSED(runTime,loadTime),ELAPSED(doneTime,runTime));


    const float tolerance = 1e-3;
    bool passed = true;
    // Validating Results
    for (int i =0 ; i< DATA_SIZE; i++){
        if (abs(float(z[i]) - float(z2[i])) > tolerance){
            printf("Error at index %d, expected %d but got %d", i, z2[i], z[i]);
            passed = false;
        }

    }

    if (passed) {
        printf("All passed");
    } else {
        printf("All failed");
    }

    return 0;
}





