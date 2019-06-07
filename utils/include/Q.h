#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_
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

//vector<bitset<M>> mesh;


class Q
{
private:

double *u=nullptr;
double *v=nullptr;
double *w=nullptr;
double *p=nullptr;

double *up=nullptr;
double *vp=nullptr;
double *wp=nullptr;
double *pp=nullptr;

double *un=nullptr;
double *vn=nullptr;
double *wn=nullptr;
double *pn=nullptr;

double dx;
double dy;
double dz;
double dt=0.01;


uint sizeS;
uint sizeP;

uint longEnd;
uint shortEnd;

public:
Q(int nmax1,double dx,double dy,double dz);

int nmax;

int index(int i,int j,int k);

void initialize(int N,double *a);

void VTK_out( double *X, double *Y, double *Z, int N );

void VTK_out_with_ghost( double *X, double *Y, double *Z,double *U,double *V,double *W, double *P);

void Struct_3D_Ghost(double Xa,double Xb, double Ya, double Yb,double Za,double Zb ,double **X, double **Y, double **Z);

int uIdx( int i,int j,int k);
int vIdx(int i,int j,int k);
int wIdx(int i,int j,int k);
int pIdx(int i,int j,int k);

// i=0, face 0; i=1, face 1;
// j=0, face 2; i=1, face 3;
// i=0, face 4; i=1, face 5;
//void uGetFace(uint fId,uint *size , double **face);
//void vGetFace(uint fId,uint *size , double **face);
//void wGetFace(uint fId,uint *size , double **face);

void project(double *U,double *V,double *W, double *P);
void correct(double *U,double *V,double *W,double *P);
void update(double *del_u);

//void uSetFace( double *val );
void setNeumanPressure();
void setSolid(  );

void residualInviscid(double *U, double *V, double *W,double *P, double *R);
void sourceMMS(double Re,double *xyz,double *rhs);
void setExactMMS(double *U, double *V, double *W,double *P,double Xa,double Ya, double Za);
void exactMMS(double *xyz,double *u1,double *v1,double *w1,double *p1);
//void exactMMSPrimitive(double *xyz,double *u1,double *v1,double *w1,double *p1);
void sourceAdd(double Xa,double Ya,double Za,double Re,double *R);
void residualViscous(double Re, double *U, double *V, double *W,double *P, double *R );
void residualInviscidwithP( double *U, double *V, double *W, double *P, double *R );
void Uxy( double Xa, double Ya,double Za, int i, int j,int k, double *xyz );
void Vxy( double Xa, double Ya,double Za, int i, int j,int k, double *xyz );
void Wxy( double Xa, double Ya,double Za, int i, int j,int k, double *xyz );
void Pxy( double Xa, double Ya,double Za, int i, int j,int k, double *xyz );
void debug(double Re,double Xa,double Ya,double Za,double *U,double *V,double *W,double *P);
void getResNorm(double *R, double *del_u );
void predict(double *U,double *V,double *W, double *Um,double *Vm,double *Wm ,double *R);
void setExactSource(double Re, double *U, double *V, double *W,double Xa, double Ya, double Za );
void setInteriorToZero(double *U,double *V, double *W,double *P);
void setBoundaryLidDrivenCavity(double *U,double *V,double *W,double *P);
void tridag( double * restrict a, double * restrict b, double * restrict c, double * restrict r, double * restrict u, int n );

void residualViscousCrossDerivativeOnly( double Re, double *U, double *V, double *W, double *P, double *R );

uint getSizeS();
uint getSizeP();

~Q();
};





#endif
