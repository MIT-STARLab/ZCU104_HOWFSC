#include "mtimesv.h"

void matrixTimesVector1(data_t C[M], data_t A[M*N], data_t B[N])
{

// #pragma HLS interface mode=s_axilite
	unsigned i,j,k;
	// static data_t cache[N];
	data_t acc, x;
#pragma HLS array_partition factor=16 type=cyclic variable=A
#pragma HLS array_partition factor=16 type=cyclic variable=B
#pragma HLS array_partition factor=16 type=cyclic variable=C
  // LOOP_Cache: for (i=0; i<N; i++) cache[i] = B[i];
	
	LOOP_M: for(i=0, k=0; i < M; i++)
	{
		acc = 0;
		LOOP_N: for(j=0; j < N; j++, k++)
		{
			acc += A[k] * B[j];
		}
		C[i] = acc;
	}
}
