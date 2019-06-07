#include "pressure.h"

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


