#include "realChunkedArray.h"
#include "params.h"
#include <stdio.h>
#include <stdlib.h>

void RealChunkedArray::moveHostToDevice()
{
#if ( PITTPACKACC )
#pragma acc update device( P [0:arraySize] )
#endif
}

PittPackResult RealChunkedArray::allocate( int *n, int nbl )
{
    Nx = n[0];
    Ny = n[1];
    Nz = n[2];

    arraySize =  Nx * Ny * Nz;

    // cout << " allocated P with size " << arraySize << endl;
    P = new ( std::nothrow ) double[arraySize];

    if ( nbl <= 0 )
    {
        cout << " Exit Code " << BLOCK_NUMBER_FAIL << endl;
        cout << BLUE << PittPackGetErrorEnum( BLOCK_NUMBER_FAIL ) << RESET << endl;
        exit( 1 );
    }

    cout << " allocated P with size " << arraySize << endl;

    for ( int i = 0; i < arraySize; i++ )
    {
        P[i] = 0.0;
    }
    nChunk    = nbl;
    chunkSize = arraySize / nbl; /*! chunk size refers to the size of the chunks, that is the whole
                                    array divided by nchunks, this is to get rid of multiplication by two */

#if ( PITTPACKACC )
#pragma acc enter data create( this [0:1] ) async(2)
#pragma acc update device( this ) 
#pragma acc enter data create( P [0:arraySize] ) async(3)
#endif

#if ( PITTPACKACC )
acc_async_wait_all();
#endif

    if ( P == NULL )
    {
        return ( ALLOCATION_FAIL );
    }
    else
    {
        return ( SUCCESS );
    }


}

int  RealChunkedArray::getChunkSize() { return ( chunkSize ); }
void RealChunkedArray::setCoords( double *X )
{
    Xa = X[0];
    Xb = X[1];
    Ya = X[2];
    Yb = X[3];
    Za = X[4];
    Zb = X[5];
}

void RealChunkedArray::setDirection( int dir )
{
    orientation = dir;

    double coeff[3]    = {1.0, 1.0, 1.0};
    coeff[orientation] = 1. / nChunk;

    nx = coeff[0] * Nx;
    ny = coeff[1] * Ny;
    nz = coeff[2] * Nz;
#if ( PITTPACKACC )
#pragma acc update device( nx, ny, nz )
#endif
    cout << "each chunk dims=" << nx << " " << ny << " " << nz << endl;
}

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
int                 RealChunkedArray::size() { return ( arraySize ); }

void RealChunkedArray::moveDeviceToHost()
{
#if ( PITTPACKACC )
#pragma acc update self( P [0:arraySize] )
#endif
}

#if ( PITTPACKACC )
#pragma acc routine seq
// index 0 retrieves real and indez zero retrieves
inline double &RealChunkedArray::operator()( int i, int j, int k) { return P[ ( i + nx * j + nx * ny * k ) ]; }
#else
double &RealChunkedArray::operator()( int i, int j, int k) { return P[( i + nx * j + nx * ny * k ) ]; }
#endif

#if ( PITTPACKACC )
#pragma acc routine seq
inline double &RealChunkedArray::operator()( int i, int j, int k, int dir )
{
    if ( dir == 0 )
    {
        return ( P[ ( i + nx * j + nx * ny * k ) ] );
    }
    else if ( dir == 1 )
    {
        return ( P[ ( j + ny * i + nx * ny * k ) ] );
    }
    else
    {
        return ( P[index] );
    }
 
}
#else

double &RealChunkedArray::operator()( int i, int j, int k, int dir )
{
    if ( dir == 0 )
    {
        return ( P[ ( i + nx * j + nx * ny * k ) ] );
    }
    else if ( dir == 1 )
    {
        return ( P[ ( j + ny * i + nx * ny * k ) ] );
    }
    else
    {
        cout<<"incorrect index in RealChunkedArray operator " <<endl;
        exit(1);  
        return ( P[index] );
    }
 
}

#endif

#if ( PITTPACKACC )
#pragma acc routine seq
inline double &RealChunkedArray::operator()( int chunkId, int dir, int i, int j, int k )
{
    /*! brief: rearrangement is consistent with the planes where (1,0,0) (0,1,0) and (0,0,1) around which rotation occurs */

    if ( dir == 0 )
    {
        return ( P[ ( chunkId * chunkSize  + i + nx * j + nx * ny * k ) ] );
    }
    else if ( dir == 1 )
    {
        // return( P[2 * ( chunkId*chunkSize/2 + j + ny * k + nx * ny * i ) ]);
        return ( P[  ( chunkId * chunkSize  + j + ny * i + nx * ny * k ) ] );
    }
    else if ( dir == 2 )
    {
        return ( P[( chunkId * chunkSize  + k + nz * j + nz * ny * i ) ] );
    }
    /* this one is added to help solve the Thomas in place and while in y- major direction, y , x and i and j replace each other*/

    else if ( dir == 3 )
    {
        return ( P[ ( chunkId * chunkSize + i + ny * j + nx * ny * k ) ] );
    }
    // note that at rotated to y-dir, ny corresponds to i, nx corresponds to j
    else if ( dir == 4 )
    {
        return ( P[ ( chunkId * chunkSize  + k + nz * j + nz * ny * i ) ] );
    }
 

}
#else
double &RealChunkedArray::operator()( int chunkId, int dir, int i, int j, int k )
{
    /*! brief: rearrangement is consistent with the planes where (1,0,0) (0,1,0) and (0,0,1) around which rotation occurs */

    if ( dir == 0 )
    {
        return ( P[( chunkId * chunkSize  + i + nx * j + nx * ny * k ) ] );
    }
    else if ( dir == 1 )
    {
        return ( P[ ( chunkId * chunkSize  + j + ny * i + nx * ny * k ) ] );
    }
    else if ( dir == 2 )
    {
        return ( P[ ( chunkId * chunkSize + k + nz * j + nz * ny * i ) ] );
    }
    /* this one is added to help solve the Thomas in place and while in y- major direction, y , x and i and j replace each other*/

    else if ( dir == 3 )
    {
        return ( P[ ( chunkId * chunkSize  + i + ny * j + nx * ny * k ) ] );
    }
    // note that at rotated to y-dir, ny corresponds to i, nx corresponds to j
    else if ( dir == 4 )
    {
        return ( P[( chunkId * chunkSize  + k + nz * j + nz * ny * i ) ] );
    }
    else
    {
        cout<<"incorrect index in RealChunkedArray operator " <<endl;
        exit(1);  
        return ( P[index] );
    }
 
}

#endif

// intel's compiler complaining about inlining
#if ( PITTPACKACC )
#pragma acc routine seq
inline double &RealChunkedArray::operator()( int i )
{
    /*
       if(i>=arraySize)
    {
      printf("index bigger than array size]n");

    //  exit(0);
    }
    */
    return P[i];
}
#else
double &RealChunkedArray::operator()( int i )
{
    /*
       if(i>=arraySize)
    {
      printf("index bigger than array size]n");

    //  exit(0);
    }
    */
    return P[i];
}

#endif
/*!< the order is not important */
RealChunkedArray::~RealChunkedArray()
{
    if ( P != NULL )
    {
#if ( PITTPACKACC )
#pragma acc exit data delete ( P [0:arraySize] )
#endif
        delete[] P;
    }
#if ( PITTPACKACC )
#pragma acc exit data delete ( this )
#endif
}

void RealChunkedArray::getAddress( double *rt )
{
    rt = P;
    cout << "pointer address " << rt << endl;
}
