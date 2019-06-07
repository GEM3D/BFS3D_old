#ifndef _STATEVECTOR_H_
#define _STATEVECTOR_H_
#include <algorithm>
#include <bitset>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stack>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <complex>
#include <assert.h>

/*!
 *
 * \class Pressure to Overload Some Operators
 * This class initiates a cell centered pressure field
 *
 * */

class Pressure
{
    private:
    int N;
    double *P;

    public:
    Pressure(){};  
    Pressure( int n ); /*! Constructor  */
    void allocate(int n);
    int size();
    double &operator()( int i, int j, int k );
    void print();

    ~Pressure();
};

class Uvel
{
    private:
    int N;
    double *u;

    public:
    Uvel( int n ); /*! Constructor  */
    double &operator()( int i, int j, int k );

    ~Uvel();
};

class Vvel
{
    private:
    int N;
    double *v;

    public:
    Vvel( int n ); /*! Constructor  */
    double &operator()( int i, int j, int k );

    ~Vvel();
};

class Wvel
{
    private:
    int N;
    double *w;

    public:
    Wvel( int n ); /*! Constructor  */
    double &operator()( int i, int j, int k );

    ~Wvel();
};

#endif
