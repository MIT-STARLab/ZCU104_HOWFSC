#include "matrixvector.h"
#include <iostream>
#include <fstream>
using namespace std;
// #include "hls_stream.h"
// #include "ap_int.h"

inline data_t at(data_t *A,unsigned M, unsigned N, unsigned i,unsigned j, bool rowmajor = true)
{
  return rowmajor ? A[i*N+j] : A[j*M+i] ;
}

inline void set(data_t *A,unsigned M, unsigned N, unsigned i,unsigned j,data_t x, bool rowmajor = true)
{
  A[rowmajor ? i*N+j : j*M+i] = x; 
}

void normalMatrixToRake(data_t *A, data_t *R, unsigned M, unsigned N, unsigned K, bool rowmajor = true)
{
    for (unsigned i=0; i<M; i++)
	for (unsigned j=0; j<N; j++)
	R[i%K + j*K + (i/K)*(K*N)] = at(A,M,N,i,j,rowmajor);
}

void rakeToNormalMatrix(data_t *A, data_t *R, unsigned M, unsigned N, unsigned K, bool rowmajor = true)
{
    for (unsigned x=0; x<N*M; x++)
    {
        unsigned i = x % K + (x / (K*N))*K;
        unsigned j = (x % (K*N)) / K;
        set(A,M,N,i,j,R[x],rowmajor);
    }
}

void matrixMultiply(data_t C[MATRIX_ROWS], data_t A[MATRIX_ROWS*MATRIX_COLS], data_t B[MATRIX_COLS], bool rowmajor = true)
{
	unsigned i,j,k;
	data_t acc, x;
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

void printMatrix(data_t* A, unsigned M, unsigned N)
{
    // for (unsigned i=0; i<MATRIX_ROWS; i++)
    // {        
    //     for (unsigned j=0; j<MATRIX_COLS; j++)
    //         cout << at(A,MATRIX_ROWS,MATRIX_COLS,i,j) << " ";
    //     cout << endl;
    // }

    for (unsigned i=0; i<M*N; i++)
    {
        cout << A[i] << " ";
        if (i % N == N-1) cout << endl;
    }
}

void printVector(data_t* C, unsigned M)
{
    for (unsigned j=0; j<M; j++)
        cout << C[j] << " ";
    cout << endl;
}
int main() 
{
    unsigned i,j;
    bool printDebug = MATRIX_ROWS*MATRIX_COLS < 100;
    // C = A * B
    data_t A[MATRIX_ROWS*MATRIX_COLS];
    data_t R[MATRIX_ROWS*MATRIX_COLS];
    data_t B[MATRIX_COLS];
    data_t C[MATRIX_ROWS];
    data_t C2[MATRIX_ROWS];

    data_t y = 1;

    // initialize A as a row-major matrix
    for (i=0; i<MATRIX_ROWS; i++)
        for (j=0; j<MATRIX_COLS; j++)
            set(A,MATRIX_ROWS,MATRIX_COLS,i,j,y++);

    // display
    if (printDebug) printMatrix(A, MATRIX_ROWS, MATRIX_COLS);
    
    // initialize B 
    for (j=0; j<MATRIX_COLS; j++)
        B[j]=y++;

    // comparison
    matrixMultiply(C2, A, B);

    normalMatrixToRake(A, R, MATRIX_ROWS, MATRIX_COLS, MATRIX_RAKE);

    // display rake
    if (printDebug) printMatrix(R, MATRIX_COLS * MATRIX_ROWS/MATRIX_RAKE, MATRIX_RAKE);
    
    // hardware version
    for(i=0; i<NUM_ITERATIONS; i++) 
        matrixvector(C,R,B);
        
    if (printDebug) printVector(C2,MATRIX_COLS);
    if (printDebug) printVector(C,MATRIX_COLS);

    // compare
    for (i=0; i<MATRIX_ROWS; i++)
    {
        if (C2[i] != C[i])
        {
            cout << "Test failed at " << i << " " << C2[i] << "~=" << C[i] << endl; 
            return 1;
            
        }

    }
    return 0;
}