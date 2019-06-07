#ifndef _STAGGERED_H_
#define _STAGGERED_H_
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
using namespace std;


#define restrict __restrict__

/*!    \class Staggered
 *     \brief  This Class encapsulates the methods for incompressible staggered grid calculations
 *   
 *    
 *    
 *    */



// vector<bitset<M>> mesh;



class Staggered
{
    private:
    double dx;
    double dy;
    double dz;
    double dt = 0.001;

    //double boundaryValue[6];
    // every index at its own directon gets an index increment
    double Ubc[6]={0.0,0.0,0.0,1.0,0.0,0.0};
    double Vbc[6]={0.0,0.0,0.0,0.0,0.0,0.0};
    double Wbc[6]={0.0,0.0,0.0,0.0,0.0,0.0};

    bool skew=true;

//    C1=0.5, C2=C1;

   double   C1=0.5;
   double   C2=0.5;


    uint sizeS;
    uint sizeP;

    uint longEnd;
    uint shortEnd;

    public:
    Staggered( int nmax1, double dx, double dy, double dz );

// default projection is P2

    bool P1=true;
// P1 one true means us P2 pressure coorection

    int nmax;

    int index( int i, int j, int k );

    void VTK_out_with_ghost( double *X, double *Y, double *Z, double *U, double *V, double *W, double *P );

    void Struct_3D_Ghost( double Xa, double Xb, double Ya, double Yb, double Za, double Zb, double **X, double **Y, double **Z );

    int uIdx( int i, int j, int k );
    int vIdx( int i, int j, int k );
    int wIdx( int i, int j, int k );
    int pIdx( int i, int j, int k );

    // i=0, face 0; i=1, face 1;
    // j=0, face 2; i=1, face 3;
    // i=0, face 4; i=1, face 5;
    // void uGetFace(uint fId,uint *size , double **face);
    // void vGetFace(uint fId,uint *size , double **face);
    // void wGetFace(uint fId,uint *size , double **face);

    void project( double *U, double *V, double *W, double *P );
    void correct( double *U, double *V, double *W, double *P );

    // void uSetFace( double *val );
    void setNeumanPressure();
    void setSolid();

    void residualInviscidSkew( double *U, double *V, double *W, double *P, double *R ); 
    void residualInviscid( double *U, double *V, double *W, double *P, double *R );
    void sourceMMS( double Re, double *xyz, double *rhs );
    void setExactMMS( double *U, double *V, double *W, double *P, double Xa, double Ya, double Za );
    void exactMMS( double *xyz, double *u1, double *v1, double *w1, double *p1 );
    // void exactMMSPrimitive(double *xyz,double *u1,double *v1,double *w1,double *p1);
    void sourceAdd( double Xa, double Ya, double Za, double Re, double *R );
    void residualViscous( double Re, double *U, double *V, double *W, double *P, double *R );



    void Uxy( double Xa, double Ya, double Za, int i, int j, int k, double *xyz );
    void Vxy( double Xa, double Ya, double Za, int i, int j, int k, double *xyz );
    void Wxy( double Xa, double Ya, double Za, int i, int j, int k, double *xyz );
    void Pxy( double Xa, double Ya, double Za, int i, int j, int k, double *xyz );

    void debug( double Re, double Xa, double Ya, double Za, double *U, double *V, double *W, double *P );
    void getResNorm( double *R, double *del_u );
    void predict( double *U, double *V, double *W, double *Um, double *Vm, double *Wm, double *R );
    void setExactSource( double Re, double *U, double *V, double *W, double Xa, double Ya, double Za );
    void setInteriorToZero( double *U, double *V, double *W, double *P );
    void setBoundaryLidDrivenCavity( double *U, double *V, double *W, double *P );
    void solveManufacturedSolution(double Re,double Xa,double Ya, double Za,double *U,double *V,double *W,double *P ); /*! This function does not project, it is only developed to verify the fluxes */
  
    uint getSizeS();
    uint getSizeP();
    void fluxDebug(double Re, double Xa,double Ya, double Za); /*! set field variables, focus on one point, get the residual, see if converges to zero as grid is refined*/
    
    void solveLidDrivenCavity(double Re,double Xa,double Ya, double Za, double *U,double *V,double *W,double *P);
    void solveLidDrivenCavityCN(double Re,double Xa,double Ya, double Za, double *U,double *V,double *W,double *P);

    void pressureCorrect(double Re,double *Pm, double *phi);

    void residualViscousCrossDerivativeOnly( double Re, double *U, double *V, double *W, double *P, double *R );

    ~Staggered();
};

#endif
