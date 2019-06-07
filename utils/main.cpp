#include "staggered.h"
#include "tridiag.h"
#include "fft.h"
#include "statevector.h"
#define PLOT 1
#define SOL 0
#define MMS 0


/*  \mainpage
 *
 *   "Single Block implementation of Incompressible Staggered Grid Solver"
 *
 *            <br> Part of NSF project: GEM3D
 *
 *
 *             Required Libraries:
 *               CMAKE <br>
 *
 *
 *             Usage:
 *              progName  <params.txt
 *
 * @authors     Dr. Jaber J. Hasbestan, Dr. Inanc Senocak <br>
 *  PI:                 Inanc Senocak
 * @date        14 Dec 2017
 * \details
 * Copyright (c)  by the Authors 2017
 * Swanson School of Engineering
 * Department of Mechanical Engineering and Material Sciences
 * University of Pittsburgh
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

double rhs( double &x, double &y, double &z );
void Struct_3D( double Xa, double Xb, double Ya, double Yb, int N, int M, double **X, double **Y, double **Z );
void VTK_out( double *X, double *Y, double *Z, int N, int M, int ID );

int main()
{

    int N;
    // N is the number if elements
    cout << "enter N" << endl;
    cin >> N;

    // include the ghost inside initilization

    double a[4] = {0.0, 0.0, 0.0, 0.0};

    uint size;

    double *d = nullptr;
    double *X = nullptr;
    double *Y = nullptr;
    double *Z = nullptr;

    double c1 = -0.5;
    double c2 = -c1;
    /*
       double xshrink,yshrink,zshrink;

        cout<<"This is a Box by default "<<endl;
        cout<< "if you like to shrink at a certain direction enter the coefficinet"<<endl;
        cout<< " if the answer is no type 1.0"<<endl;
        cin>>xshrink;
        cin>>yshrink;
        cin>>zshrink;
     */

    double zshrink = 0.1;
    double yshrink = 1.0;
    double xshrink = 1.0;

    double CC = 0.5;
    double CS = 0.5;

    cout << "shrink coeffs in x,y and z directions"
         << " " << xshrink << " " << yshrink << " " << zshrink << endl;

    double Xa = c1 * xshrink, Xb = c2 * xshrink;
    double Ya = c1 * yshrink, Yb = c2 * yshrink;
    double Za = c1 * zshrink, Zb = c2 * zshrink;

    const int level = 2;

    double val = 0.0;
    double x, y, z;
    double dx = ( Xb - Xa ) / N;
    double dh = dx;
    double C1 = 1. / 6.;
    static const double dh2 = dh * dh;

    Staggered q( N, dx, dx, dx );

    //    q.initialize( N, a );

    // q.wGetFace(1, &size, &d );

    delete[] d;

    double Re = 1000.;
    //   double Re = 1.e15;
    double uold;
    cout << " Re " << Re << endl;

    int count = 0;
    double fv = .5;
    /*
     *
     *  debug tridiagonal
    double a1[4] = {0.0, 1., -0.5,-.9, };
    double b1[4] = { 4,  18,  6,  8 };
    double c4[4] = { -1, 1,  -6,0 };
    double r1[4] = { 5,  5, 10, 23 };
    double u1[4];
    q.tridag( a1, b1, c4, r1, u1, 4);
    for(int i=0;i<4;i++)
    {
    cout<<u1[i]<<endl;
    }
    */

    // q.uSetFace(&fv);
    //
    double res = 1.0;

    double *restrict U = new double[q.getSizeS()];
    double *restrict V = new double[q.getSizeS()];

    double *restrict Rv = new double[3 * q.getSizeS()];
    double *restrict S = new double[3 * q.getSizeS()];

    double *restrict R = new double[3 * q.getSizeS()];
    double *restrict W = new double[q.getSizeS()];
    double *restrict P = new double[q.getSizeP()];
#if ( SOL )

    double *Um = new double[q.getSizeS()];
    double *Vm = new double[q.getSizeS()];
    double *Wm = new double[q.getSizeS()];

    double *R1 = new double[q.getSizeS()];
    double *R2 = new double[q.getSizeS()];
    double *R3 = new double[q.getSizeS()];

    for ( uint i = 0; i < 3 * q.getSizeS(); i++ )
    {
        R[i] = 0.0;
        S[i] = 0.0;
        Rv[i] = 0.0;
    }

    for ( uint i = 0; i < q.getSizeS(); i++ )
    {
        U[i] = 0.0;
        V[i] = 0.0;
        W[i] = 0.0;
        Um[i] = 0.0;
        Vm[i] = 0.0;
        Wm[i] = 0.0;
    }

/**********************************************************************
// keep presssure fixed for MMS analyisis and use it in momentum
//
*/
#if ( MMS )


q.solveManufacturedSolution(Re,Xa,Ya, Za,U,V,W,P );


/********************************************************************/

//    q.setExactMMS( U, V, W, P, Xa, Ya, Za );

// set the boundary for lid driven cavity
#else

#endif

#else
   q.fluxDebug(Re, Xa,Ya,Za);
   
    // this part only to debug the fluxes that given everything exact
    // the l2 norm of the residual should approach zero in the limit
    /**************************************************************/
/*****************************************************************/
//    q.debug( Re, Xa, Ya, Za, U, V, W, P );
#endif
// q.setSolid();

#if ( MMS == 0 )

//q.solveLidDrivenCavity( Re, Xa, Ya,  Za,U,V,W,P );

#endif

// Struct_3D(Xa,Xb,Ya,Yb, N,N,&X, &Y,&Z);
// q.VTK_out(X,Y,Z,N);
#if ( PLOT )

    q.Struct_3D_Ghost( Xa, Xb, Ya, Yb, Za, Zb, &X, &Y, &Z );

    q.VTK_out_with_ghost( X, Y, Z, U, V, W, P );
/*
for(uint i=0;i<q.getSizeS();i++)
{
R1[i]=R[i];
R2[i]=R[i+q.getSizeS()];
R3[i]=R[i+2*q.getSizeS()];

}
q.VTK_out_with_ghost(X,Y,Z,R1,R2,R3,P);
*/

#endif
    delete[] U;
    delete[] V;
    delete[] W;
    delete[] S;
    delete[] R;
    delete[] Rv;
    delete[] P;
    delete[] X;
    delete[] Y;
    delete[] Z;
   
#if ( SOL )
    delete[] Um;
    delete[] Vm;
    delete[] Wm;
    delete[] R1;
    delete[] R2;
    delete[] R3;

#endif

    int mysize = N;
/*
    Tridiag T( mysize );
    T.debug(); 
*/

   Fourier F(mysize);

  // F.debug();

//F.Poisson2DCT1();

//F.Poisson2DCT2();

//F.Poisson2DDir();


//F.Poisson3D();
/*****************************
 *
 *
 * *************************/

#if(1)


Pressure Pn(mysize);

std::vector<std::string>  mybc({"D","D","N","N","D","D"});

F.assignBoundary(mybc);

//F.extractTag(); move to another function

//F.detectTransforms() >> added to assign bc 

/*
std::string str="S";
F.fourierTransform(0,str ,Pn);
*/


//F.detectTridiagDirection(); >> added to assign bc

F.setUpExactFunctionPointers();

F.assignExactFreqs();

F.Poisson3DComplex(Pn);

cout<<Pn.size()<<endl;

double (*f[3])(double);
f[0]=&sin;
cout<<f[0](1.0/2.0)<<endl;

//F.Poisson2DPeriodic();

//double (**f1)()=malloc(3*sizeof *f);

//F.real2Complex();

/*
double *in=new double [mysize];
double *out=new double [2*(mysize+1)];

for(int i=0;i<mysize;i++)
{
in[i]=exp(i);
}

F.extendSignalDST00(in, out);

cout<<"original"<<endl;

for(int i=0;i<mysize;i++)
{
cout<<in[i]<<endl;

}

cout<<"extended"<<endl;
for(int i=0;i<2*(mysize+1);i++)
{
cout<<out[i]<<endl;

}
*/

// DCT00
/*
double *in=new double [mysize];
double *out=new double [2*(mysize-1)];

for(int i=0;i<mysize;i++)
{
in[i]=exp(i);
}

F.extendSignalDCT10(in, out);

cout<<"original"<<endl;

for(int i=0;i<mysize;i++)
{
cout<<in[i]<<endl;

}

cout<<"extended"<<endl;
for(int i=0;i<2*(mysize-1);i++)
{
cout<<out[i]<<endl;

}
*/

// DCT10
/*
double *in=new double [mysize];
double *out=new double [2*(mysize)];

for(int i=0;i<mysize;i++)
{
in[i]=exp(i);
}

F.extendSignalDCT10(in, out);

cout<<"original"<<endl;

for(int i=0;i<mysize;i++)
{
cout<<in[i]<<endl;

}

cout<<"extended"<<endl;
for(int i=0;i<2*(mysize);i++)
{
cout<<out[i]<<endl;

}
*
*/

F.test();




/*
F.perfromAllTransforms( Pn );

F.convertDirichletAtBoundary(Pn);

Pn.print();
*/
//F.setSource(Pn,XX,);


#else
//F.gausSeidel();
#endif

};

void Struct_3D( double Xa, double Xb, double Ya, double Yb, int N, int M, double **X, double **Y, double **Z )
{
    int i;

    double hx = N;
    double hy = N;

    double hz = N;
    double Xh = ( Xb - Xa ) / ( hx );
    double Yh = ( Yb - Ya ) / ( hy );

    cout << "delx " << Xh << endl;

    double Zh = 0.01 * Xh;

    ( *X ) = (double *)new double[N + 1];
    ( *Y ) = (double *)new double[N + 1];
    ( *Z ) = (double *)new double[N + 1];

    for ( i = 0; i < N + 1; i++ )
    {
        ( *X )[i] = Xa + Xh * i;
    }

    for ( i = 0; i < N + 1; i++ )
    {
        ( *Y )[i] = Ya + Yh * i;
    }

    for ( i = 0; i < N + 1; i++ )
    {
        ( *Z )[i] = Xa + Zh * i;
    }
}

double rhs( double &x, double &y, double &z )
{
    static const double pi = atan( 1.0 ) * 4.0;
    static const double pi_2 = pi * pi;
    return ( -( 3.0 * pi_2 * sin( ( pi * x ) / 2.0 ) * sin( ( pi * y ) / 2.0 ) * sin( ( pi * z ) / 2.0 ) ) / 4.0 );
}

