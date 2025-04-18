#include "mtimesv.h"

extern "C" void matrixTimesVector1(data_t C[M], data_t A[M*N], data_t B[N])
{

//#pragma HLS INTERFACE mode=m_axi port=A bundle=aximm1
//#pragma HLS INTERFACE mode=m_axi port=B bundle=aximm2
//#pragma HLS INTERFACE mode=m_axi port=C bundle=aximm3

//#pragma HLS array_partition factor=16 type=cyclic variable=A
//#pragma HLS array_partition factor=16 type=cyclic variable=B
//#pragma HLS array_partition factor=16 type=cyclic variable=C
//#pragma HLS bind_op variable=acc op=dadd impl=fabric latency=2
//#pragma HLS bind_op variable=acc op=dmul impl=fabric latency=2
//#pragma HLS pipeline II=1

  // LOOP_Cache: for (i=0; i<N; i++) cache[i] = B[i];
	unsigned i,j,k;
	// static data_t cache[N];
	data_t acc, x;
	
	LOOP_M: for(i=0, k=0; i < M; i++)
	{
		acc = 0;
		LOOP_N: for(j=0; j < N; j++, k++)
		{
#pragma HLS unroll factor=16
			acc += A[k] * B[j];
		}
		C[i] = acc;
	}
}
