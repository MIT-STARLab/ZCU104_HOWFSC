#include "mtimesv.h"
 
#define VERIFY 1
#define ROWMAJOR 0

inline data_t at(data_t *A,unsigned i,unsigned j)
{
#if ROWMAJOR
  return A[j*M+i];
#else
  return A[i*N+j];
#endif
}
inline void set(data_t *A,unsigned i,unsigned j,data_t x)
{
#if ROWMAJOR
  A[j*M+i] = x; 
#else
  A[i*N+j] = x;
#endif
}
void matrxToPipelined(data_t dest[M*N], data_t source[M*N])
{
  unsigned i,j,k;
}

void gold(data_t C[M], data_t A[M*N], data_t B[N])
{
	unsigned i,j,k;
	data_t acc, x;
	LOOP_M: for(i=0, k=0; i < M; i++)
	{
		acc = 0;
		LOOP_N: for(j=0; j < N; j++, k++)
		{
			acc += A[k] * B[j];
		}
		C[i] = acc;
	}
 return;
//  unsigned i,j,k;
//#if ROWMAJOR
//  for(i=0, k=0; i < M; i++)
//  {
//    C[i] = 0;
//    for(j=0; j < N; j++, k++)
//      C[i] += A[k] * B[j];
//  }
//#else
//  for(j=0; j < N; j++)
//  {
//    C[i] = 0;
//    for(i=0; i < M; i++)
//      C[i] += A[j*M+i] * B[j];
//  }
//#endif
}

int main () // C = A * B
{
  data_t A[N*M];
  data_t B[N];
  data_t C1[M];
  data_t C[M];
  data_t accum;
  unsigned i,j,k;  

#if VERIFY
  // Create input data
  // for (k=0; k<N*M; k++) A[k] = (data_t)k;
  k=0;
  for (i=0; i<M; i++)
    for (j=0; j<N; j++) 
#if ROWMAJOR
      A[i*N+j] = (data_t)k++;
#else
      A[j*M+i] = (data_t)k++;
#endif      
  for (i=0; i<N; i++) B[i] = (data_t)i;
  gold(C1,A,B);

  ofstream result;
  result.open ("result.dat");
  result << "B" << endl;
  for (j=0; j<N; j++) result << B[j] << " ";
#if ROWMAJOR  
  result << endl << "A [" << M << " x " << N << "]" << endl;
  for (i=0,k=0; i<M; i++)
  {
    for (j=0; j<N; j++, k++) result << A[k] << " ";
    result << endl;
  }
#else
  result << endl << "A' [" << N << " x " << M << "]" <<  endl;
  for (i=0,k=0; i<N; i++)
  {
    for (j=0; j<M; j++, k++) result << A[k] << " ";
    result << endl;
  }
#endif

#endif

  // Call the function
  for(i=0; i<NUM_ITERATIONS; i++) 
  {
    matrixTimesVector1(C,A,B);
  }

#if VERIFY
  result << "C" << endl;
  for (i=0; i<M; i++) result << C[i] << " ";
  result << endl;
  result << "gold" << endl;
  for (i=0; i<M; i++) result << C1[i] << " ";
  result << endl;
  result.close();

  for (j=0; j<M; j++)
  {
    if (C1[j] != C[j]) 
     { 
      cout << "Test failed at " << j << " " << C1[j] << "~=" << C[j] << endl; 
      return 1;
    }
  }
#endif

  cout << "Test passed" << endl;
  return 0;
}

