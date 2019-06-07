#include "statevector.h"

//**************************************************************
//   This class is defined to overload the paranthesis operator 
//
// 
//
// 
// 
//***********************************************************


Pressure::Pressure( int n )
{
    N = n;
    P=new double[N*N*N];
  
    for(int i=0;i<N*N*N;i++)
    {
    P[i]=0.0;
    }


}

void Pressure::allocate( int n )
{
    N = n;
    P=new double[N*N*N];
  
    for(int i=0;i<N*N*N;i++)
    {
    P[i]=0.0;
    }


}


int Pressure::size( )
{
    return(N);
}



double& Pressure::operator()(int i, int j,int k) 
{
// some assertion to catch bugs 
//
    assert(i >= 0 && i < N);
    assert(j >= 0 && j < N);
    assert(k >= 0 && k < N);
 
    return P[i+N*j+N*N*k];
}




Pressure::~Pressure()
{
    delete[] P;
}

void Pressure::print()
{
/*
    for ( int k = 1; k < N-1; k++ )
    {
        for ( int j = 1; j < N-1; j++ )
        {
            for ( int i = 1; i < N-1; i++ )
            {
                std::cout << P[i+N*j+N*N*k] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << "-----------------------" << std::endl;
    }
*/

    for ( int k = 0; k < N; k++ )
    {
        for ( int j = 0; j < N; j++ )
        {
            for ( int i = 0; i < N; i++ )
            {
                std::cout << P[i+N*j+N*N*k] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << "-----------------------" << std::endl;
    }


}

//====================================================

Uvel::Uvel( int n )
{
    N = n;
    u=new double[(N+1)*N*N];
}

double& Uvel::operator()(int i, int j,int k) 
{
// some assertion to catch bugs 
//
    assert(i >= 0 && i < N+1);
    assert(j >= 0 && j < N);
    assert(k >= 0 && k < N);
 
    return u[ k  * ( N + 1 ) *  N  + j * ( N + 1 ) + i];
}


Uvel::~Uvel()
{
    delete[] u;
}

//====================================================

Vvel::Vvel( int n )
{
    N = n;
    v=new double[(N+1)*N*N];
}


double& Vvel::operator()(int i, int j,int k) 
{
// some assertion to catch bugs 
//
    assert(i >= 0 && i < N);
    assert(j >= 0 && j < N+1);
    assert(k >= 0 && k < N);
 
    return v[k * ( N + 1 ) *  N  + j *  N  + i];
}


Vvel::~Vvel()
{
    delete[] v;
}

//===================================================

Wvel::Wvel( int n )
{
    N = n;
    w=new double[(N+1)*N*N];
}


double& Wvel::operator()(int i, int j,int k) 
{
// some assertion to catch bugs 
//
    assert(i >= 0 && i < N);
    assert(j >= 0 && j < N);
    assert(k >= 0 && k < N+1);
 
    return w[k *  N  *  N  + j * N + i];
}


Wvel::~Wvel()
{
    delete[] w;
}













