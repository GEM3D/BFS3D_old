#include "staggered.h"
#define AVG 1
#define DEBUG 0
#define ZSHRINK 1
#define P_ON 0
#define NITER 10
#define AA 0.1

/*
 *
 *
 *
 *       Method of Manufactured solution can be used to verify the results
 *       We set u,v,w,p to some assumed functions of space and then plug them
 *       into the NS equation, since these functions do not satisfy the NS it will generate
 *       some source terms (S), to verify the fluxes, we initialize the data to exact values,
 *       calculate the fluxes, obtain the l2-norm on the entire domain and monitor the value
 *       every grid refinement the error should get closer to zero depending on the numerical accuracy
 *       of the solution, we monitor residuals for momentum equations,
 *       Residual=R_in-R_visc-S
 *       Notive that not including the prssure at momentum the system becomes very stiff because
 *       of source terms
 *
 *
 *
 *
 */

Staggered::Staggered( int nmax1, double dx1, double dx2, double dx3 )
{
    // nmax1 is the number of elements

    dx = dx1;
    dy = dx2;
    dz = dx3;
    dt = 0.005;

    double volume = dx * dy * dz;

    nmax = nmax1 + 2;

    cout << dx << " " << dy << " " << dz << " " << dt << endl;
    //
    longEnd = nmax;

    shortEnd = nmax - 1;

    // nmax is the number of elements with ghost cell

    cout << "nelem " << nmax1 << endl;

    cout << "nelem with ghost" << nmax << endl;

    sizeS = ( nmax + 1 ) * ( nmax ) * ( nmax );

    // u is u at current time step, un means new values hence at time step n+1 and up is previous value of u means u at n-1
    sizeP = nmax * nmax * nmax;
}

Staggered::~Staggered()
{
}
void Staggered::VTK_out_with_ghost( double *X, double *Y, double *Z, double *U, double *V, double *W, double *P )
{
    unsigned int i, j, k;

    FILE *fp = NULL;
    /*
        i=1;j=2;k=4;

        u[uIdx(i,j,k)]=5.0;
        v[vIdx(i,j,k)]=5.0;
        w[wIdx(i,j,k)]=5.0;
        p[pIdx(i,j,k)]=5.0;

        cout<<" "<<uIdx(i,j,k)<<" "<<vIdx(i,j,k)<<" "<<wIdx(i,j,k)<<endl;
    */

    int M, N;

    // here we get some data into variable data
    char filename[64];
    sprintf( filename, "out%d.vtk", 0 );
    fp = fopen( filename, "w" );

    N = longEnd + 1;
    M = N;

    // fp=fopen("out.vtk","w");
    fprintf( fp, "# vtk DataFile Version 2.0 \n" );
    fprintf( fp, "Grid\n" );
    fprintf( fp, "ASCII\n" );
    fprintf( fp, "DATASET STRUCTURED_GRID\n" );
    fprintf( fp, "DIMENSIONS %d %d %d\n", N, M, N );
    fprintf( fp, "POINTS %d float\n", M * N * N );

    for ( k = 0; k < N; k++ )
    {
        for ( j = 0; j < N; j++ )
        {
            for ( i = 0; i < N; i++ )
            {
                fprintf( fp, "%lf %lf %lf\n", X[i], Y[j], Z[k] );
            }
        }
    }

    fprintf( fp, "CELL_DATA %d\n", ( N - 1 ) * ( N - 1 ) * ( N - 1 ) );

    fprintf( fp, "SCALARS U float 1\n" );
    fprintf( fp, "LOOKUP_TABLE default\n" );

    for ( k = 0; k < longEnd; k++ )
    {

        for ( j = 0; j < longEnd; j++ )
        {
            for ( i = 0; i < longEnd; i++ )
            {
#if ( AVG )
                fprintf( fp, "%lf\n", ( U[uIdx( i, j, k )] + U[uIdx( i + 1, j, k )] ) * 0.5 );
#else
                fprintf( fp, "%lf\n", ( U[uIdx( i, j, k )] ) );
#endif
                //        cout<<q1.Zindex(i,j,0)<<endl;
            }
        }
    }

    fprintf( fp, "SCALARS V float 1\n" );
    fprintf( fp, "LOOKUP_TABLE default\n" );

    for ( k = 0; k < longEnd; k++ )
    {

        for ( j = 0; j < longEnd; j++ )
        {
            for ( i = 0; i < longEnd; i++ )
            {
// fprintf(fp,"%lf\n",(double)q1.Zindex(i,j,k)/((N-1)*(N-1)*(N-1)));
// fprintf( fp, "%lf\n",q[index( i, j, k )] );
//      fprintf(fp,"%lf\n",1.0);

#if ( AVG )
                fprintf( fp, "%lf\n", ( V[vIdx( i, j, k )] + V[vIdx( i, j + 1, k )] ) * 0.5 );
#else
                fprintf( fp, "%lf\n", ( V[vIdx( i, j, k )] ) );
#endif
                //        cout<<q1.Zindex(i,j,0)<<endl;
            }
        }
    }

    fprintf( fp, "SCALARS w float 1\n" );
    fprintf( fp, "LOOKUP_TABLE default\n" );

    for ( k = 0; k < longEnd; k++ )
    {

        for ( j = 0; j < longEnd; j++ )
        {
            for ( i = 0; i < longEnd; i++ )
            {
#if ( AVG )
                fprintf( fp, "%lf\n", ( W[wIdx( i, j, k )] + W[wIdx( i, j, k + 1 )] ) * 0.5 );
#else
                fprintf( fp, "%lf\n", ( W[wIdx( i, j, k )] ) );
#endif
                //        cout<<q1.Zindex(i,j,0)<<endl;
            }
        }
    }

    fprintf( fp, "SCALARS P float 1\n" );
    fprintf( fp, "LOOKUP_TABLE default\n" );

    int cout = 0;

    for ( k = 0; k < longEnd; k++ )
    {

        for ( j = 0; j < longEnd; j++ )
        {
            for ( i = 0; i < longEnd; i++ )
            {
                // fprintf(fp,"%lf\n",(double)q1.Zindex(i,j,k)/((N-1)*(N-1)*(N-1)));
                // fprintf( fp, "%lf\n",q[index( i, j, k )] );
                //      fprintf(fp,"%lf\n",1.0);
                fprintf( fp, " %lf\n", P[pIdx( i, j, k )] );
                //          cout++;
                //        cout<<q1.Zindex(i,j,0)<<endl;
            }
        }
    }

    fclose( fp );
}

// indexes for different vector
//
// i=0:nmax+1, j=0:nmax, k=0:nmax

inline int Staggered::uIdx( int i, int j, int k )
{
#if ( DEBUG )
    {
        if ( i > nmax || j > nmax - 1 || k > nmax - 1 )
        {
            throw std::runtime_error( "index out of bound" );
        }
    }
#endif
    return ( ( k ) * ( nmax + 1 ) * ( nmax ) + j * ( nmax + 1 ) + i );
}

// i=0:nmax, j=0:nmax+1, k=0:nmax
inline int Staggered::vIdx( int i, int j, int k )
{
#if ( DEBUG )
    {
        if ( i > nmax - 1 || j > nmax || k > nmax - 1 )
        {
            throw std::runtime_error( "index out of bound" );
        }
    }

#endif
    return ( k * ( nmax + 1 ) * ( nmax ) + j * ( nmax ) + i );
}

// i=0:nmax, j=0:nmax, k=0:nmax+1
inline int Staggered::wIdx( int i, int j, int k )
{
#if ( DEBUG )
    {
        if ( i > nmax - 1 || j > nmax - 1 || k > nmax )
        {
            throw std::runtime_error( "index out of bound" );
        }
    }

#endif
    return ( k * ( nmax ) * ( nmax ) + j * ( nmax ) + i );
}

// i=0:nmax, j=0:nmax, k=0:nmax
inline int Staggered::pIdx( int i, int j, int k )
{

#if ( DEBUG )
    {
        if ( i > nmax - 1 || j > nmax - 1 || k > nmax - 1 )
        {
            throw std::runtime_error( "index out of bound" );
        }
    }

#endif
    return ( k * ( nmax ) * ( nmax ) + j * ( nmax ) + i );
}

void Staggered::project( double *U, double *V, double *W, double *P )
{
    // Gaus-Seidel Iteration
    double err = 1.0;
    double val;

    double c1 = 2. / dx / dx + 2. / dy / dy + 2. / dz / dz;
    /*
      double *Pm=new double[sizeP];

     for(uint i=0;i<sizeP;i++)
    {
    Pm[i]=P[i];
    }
    */

    //  while (err>1.e-15)
    for ( uint l = 0; l < NITER; l++ )
    {
        err = 0.0;
        for ( uint i = 1; i < shortEnd; i++ )
        {
            for ( uint j = 1; j < shortEnd; j++ )
            {
                for ( uint k = 1; k < shortEnd; k++ )
                {
                    val = -( ( U[uIdx( i + 1, j, k )] - U[uIdx( i, j, k )] ) / dx + ( V[vIdx( i, j + 1, k )] - V[vIdx( i, j, k )] ) / dy
                             + ( W[wIdx( i, j, k + 1 )] - W[wIdx( i, j, k )] ) / dz ) / dt
                          + ( P[pIdx( i + 1, j, k )] + P[pIdx( i - 1, j, k )] ) / dx / dx
                          + ( P[pIdx( i, j + 1, k )] + P[pIdx( i, j - 1, k )] ) / dy / dy
                          + ( P[pIdx( i, j, k + 1 )] + P[pIdx( i, j, k - 1 )] ) / dz / dz;

                    val = val / c1;

                    err += fabs( P[pIdx( i, j, k )] - val );

                    P[pIdx( i, j, k )] = val;
                    //               cout<<"val "<<val<<endl;
                }
            }
        }
        //   cout<<"error "<<err<<endl;
    }

    //  P[pIdx( 1, 1, 1 )] = 0.0;
    // cout<<dx<<"dy "<<dy <<" dz "<<dz<<"dt "<<dt<<endl;
    // cout<<" poisson error "<<err<<endl;
    //

    /*
     for(uint i=0;i<sizeP;i++)
    {
    P[i]=Pm[i]+P[i];
    }


    delete[] Pm;
      */
}

void Staggered::correct( double *U, double *V, double *W, double *P )
{

    for ( uint i = 1; i < longEnd; i++ )
    {

        for ( uint j = 1; j < shortEnd; j++ )
        {
            for ( uint k = 1; k < shortEnd; k++ )
            {

                U[uIdx( i, j, k )] = -dt * ( P[pIdx( i, j, k )] - P[pIdx( i - 1, j, k )] ) / dx + U[uIdx( i, j, k )];
            }
        }
    }

    for ( uint i = 1; i < shortEnd; i++ )
    {

        for ( uint j = 1; j < longEnd; j++ )
        {
            for ( uint k = 1; k < shortEnd; k++ )
            {

                V[vIdx( i, j, k )] = -dt * ( P[pIdx( i, j, k )] - P[pIdx( i, j - 1, k )] ) / dy + V[vIdx( i, j, k )];
            }
        }
    }

    for ( uint i = 1; i < shortEnd; i++ )
    {

        for ( uint j = 1; j < shortEnd; j++ )
        {
            for ( uint k = 1; k < longEnd; k++ )
            {
                W[wIdx( i, j, k )] = -dt * ( P[pIdx( i, j, k )] - P[pIdx( i, j, k - 1 )] ) / dz + W[wIdx( i, j, k )];
            }
        }
    }
}

void Staggered::Struct_3D_Ghost( double Xa, double Xb, double Ya, double Yb, double Za, double Zb, double **X, double **Y, double **Z )

{
    int i;

    double hx = nmax;
    double hy = nmax;

    double hz = nmax;
    double Xh = ( Xb - Xa ) / ( hx );
    double Yh = ( Yb - Ya ) / ( hy );
    double Zh = ( Zb - Za ) / ( hz );

    cout << "delx " << Xh << endl;

    ( *X ) = (double *)new double[nmax + 1];
    ( *Y ) = (double *)new double[nmax + 1];
    ( *Z ) = (double *)new double[nmax + 1];

    for ( i = 0; i < nmax + 1; i++ )
    {
        ( *X )[i] = Xa + Xh * i;
    }

    for ( i = 0; i < nmax + 1; i++ )
    {
        ( *Y )[i] = Ya + Yh * i;
    }

    for ( i = 0; i < nmax + 1; i++ )
    {
        ( *Z )[i] = Xa + Zh * i;
    }
}

void Staggered::sourceMMS( double Re, double *xyz, double *rhs )
{

    double x = xyz[0];
    double y = xyz[1];
    double z = xyz[2];

    const double p = atan( 1.0 ) * 4.0;

    double A = AA;
    double B = 0.1;
    double C = 0.1;
    double D = 0.1;
    /*
        rhs[0] = -( C * ( p * p ) * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0
       )
                    + B * ( p * p ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y * ( 1.0 / 2.0 ) ) * ( 1.0 /
       4.0 )
                    - A * ( p * p ) * sin( p * x ) * sin( p * y ) * sin( p * z ) * 4.0 ) / Re
                 + D * p * cos( p * x * ( 1.0 / 4.0 ) ) * sin( p * y * ( 1.0 / 4.0 ) ) * sin( p * z * ( 1.0 / 4.0 ) ) * ( 1.0 / 4.0 )
                 + ( A * A ) * p * cos( p * x ) * sin( p * x ) * pow( sin( p * y ), 2.0 ) * pow( sin( p * z ), 2.0 ) * 2.0
                 + A * B * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) )
                   * sin( p * x ) * sin( p * z ) - A * B * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x )
                                                   * sin( p * y ) * sin( p * y * ( 1.0 / 2.0 ) ) * sin( p * z ) * ( 1.0 / 2.0 )
                 + A * C * p * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z ) * sin( p * x ) * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y )
                   * sin( p * z * ( 1.0 / 2.0 ) ) + A * C * p * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x )
                                                    * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y ) * sin( p * z ) * ( 1.0 / 2.0 );

        rhs[1]
        = ( B * ( p * p ) * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) )
            - A * ( p * p ) * cos( p * x ) * cos( p * y ) * sin( p * z )
            + C * ( p * p ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
       Re
          + D * p * cos( p * y * ( 1.0 / 4.0 ) ) * sin( p * x * ( 1.0 / 4.0 ) ) * sin( p * z * ( 1.0 / 4.0 ) ) * ( 1.0 / 4.0 )
          - ( B * B ) * p * pow( cos( p * x * ( 1.0 / 2.0 ) ), 2.0 ) * cos( p * y * ( 1.0 / 2.0 ) ) * pow( cos( p * z * ( 1.0 / 2.0 ) ), 2.0
       )
            * sin( p * y * ( 1.0 / 2.0 ) ) + B * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * pow( cos( p * y * ( 1.0 / 2.0 ) ), 2.0 )
                                             * pow( cos( p * z * ( 1.0 / 2.0 ) ), 2.0 ) * sin( p * x * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 )
          - B * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * pow( cos( p * y * ( 1.0 / 2.0 ) ), 2.0 ) * sin( p * x * ( 1.0 / 2.0 ) )
            * pow( sin( p * z * ( 1.0 / 2.0 ) ), 2.0 ) * ( 1.0 / 2.0 )
          + A * B * p * cos( p * x ) * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p *
       y )
            * sin( p * z ) - A * B * p * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x )
                             * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y ) * sin( p * z ) * ( 1.0 / 2.0 );

        rhs[2] = -( A * ( p * p ) * cos( p * x ) * cos( p * z ) * sin( p * y )
                    + B * ( p * p ) * cos( p * x * ( 1.0 / 2.0 ) ) * sin( p * y * ( 1.0 / 2.0 ) ) * sin( p * z * ( 1.0 / 2.0 ) ) * ( 1.0 /
       4.0 )
                    - C * ( p * p ) * cos( p * y * ( 1.0 / 2.0 ) ) * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * z * ( 1.0 / 2.0 ) ) ) / Re
                 + D * p * cos( p * z * ( 1.0 / 4.0 ) ) * sin( p * x * ( 1.0 / 4.0 ) ) * sin( p * y * ( 1.0 / 4.0 ) ) * ( 1.0 / 4.0 )
                 + ( C * C ) * p * pow( cos( p * y * ( 1.0 / 2.0 ) ), 2.0 ) * cos( p * z * ( 1.0 / 2.0 ) )
                   * pow( sin( p * x * ( 1.0 / 2.0 ) ), 2.0 ) * sin( p * z * ( 1.0 / 2.0 ) )
                 - B * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) )
                   * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y * ( 1.0 / 2.0 ) ) * sin( p * z * ( 1.0 / 2.0 ) )
                 + A * C * p * cos( p * x ) * cos( p * y * ( 1.0 / 2.0 ) ) * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y ) * sin( p * z )
                   * sin( p * z * ( 1.0 / 2.0 ) ) + A * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * sin( p * x )
                                                    * sin( p * y ) * sin( p * z ) * sin( p * z * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 );


        // inviscid fluxes for debug
        //

         rhs[0]=D*p*cos(p*x*(1.0/4.0))*sin(p*y*(1.0/4.0))*sin(p*z*(1.0/4.0))*(1.0/4.0)+(A*A)*p*cos(p*x)*sin(p*x)*pow(sin(p*y),2.0)*pow(sin(p*z),2.0)*2.0+A*B*p*cos(p*x*(1.0/2.0))*cos(p*y)*cos(p*y*(1.0/2.0))*cos(p*z*(1.0/2.0))*sin(p*x)*sin(p*z)-A*B*p*cos(p*x*(1.0/2.0))*cos(p*z*(1.0/2.0))*sin(p*x)*sin(p*y)*sin(p*y*(1.0/2.0))*sin(p*z)*(1.0/2.0)+A*C*p*cos(p*y*(1.0/2.0))*cos(p*z)*sin(p*x)*sin(p*x*(1.0/2.0))*sin(p*y)*sin(p*z*(1.0/2.0))+A*C*p*cos(p*y*(1.0/2.0))*cos(p*z*(1.0/2.0))*sin(p*x)*sin(p*x*(1.0/2.0))*sin(p*y)*sin(p*z)*(1.0/2.0);

        rhs[1]=D*p*cos(p*y*(1.0/4.0))*sin(p*x*(1.0/4.0))*sin(p*z*(1.0/4.0))*(1.0/4.0)-(B*B)*p*pow(cos(p*x*(1.0/2.0)),2.0)*cos(p*y*(1.0/2.0))*pow(cos(p*z*(1.0/2.0)),2.0)*sin(p*y*(1.0/2.0))+B*C*p*cos(p*x*(1.0/2.0))*pow(cos(p*y*(1.0/2.0)),2.0)*pow(cos(p*z*(1.0/2.0)),2.0)*sin(p*x*(1.0/2.0))*(1.0/2.0)-B*C*p*cos(p*x*(1.0/2.0))*pow(cos(p*y*(1.0/2.0)),2.0)*sin(p*x*(1.0/2.0))*pow(sin(p*z*(1.0/2.0)),2.0)*(1.0/2.0)+A*B*p*cos(p*x)*cos(p*x*(1.0/2.0))*cos(p*y*(1.0/2.0))*cos(p*z*(1.0/2.0))*sin(p*y)*sin(p*z)-A*B*p*cos(p*y*(1.0/2.0))*cos(p*z*(1.0/2.0))*sin(p*x)*sin(p*x*(1.0/2.0))*sin(p*y)*sin(p*z)*(1.0/2.0);

        rhs[2]=D*p*cos(p*z*(1.0/4.0))*sin(p*x*(1.0/4.0))*sin(p*y*(1.0/4.0))*(1.0/4.0)+(C*C)*p*pow(cos(p*y*(1.0/2.0)),2.0)*cos(p*z*(1.0/2.0))*pow(sin(p*x*(1.0/2.0)),2.0)*sin(p*z*(1.0/2.0))-B*C*p*cos(p*x*(1.0/2.0))*cos(p*y*(1.0/2.0))*cos(p*z*(1.0/2.0))*sin(p*x*(1.0/2.0))*sin(p*y*(1.0/2.0))*sin(p*z*(1.0/2.0))+A*C*p*cos(p*x)*cos(p*y*(1.0/2.0))*sin(p*x*(1.0/2.0))*sin(p*y)*sin(p*z)*sin(p*z*(1.0/2.0))+A*C*p*cos(p*x*(1.0/2.0))*cos(p*y*(1.0/2.0))*sin(p*x)*sin(p*y)*sin(p*z)*sin(p*z*(1.0/2.0))*(1.0/2.0);


        // visous fluxes for debug, minus sign  already incorporated in here

        rhs[0]=-(C*(p*p)*cos(p*x*(1.0/2.0))*cos(p*y*(1.0/2.0))*cos(p*z*(1.0/2.0))*(1.0/4.0)+B*(p*p)*cos(p*z*(1.0/2.0))*sin(p*x*(1.0/2.0))*sin(p*y*(1.0/2.0))*(1.0/4.0)-A*(p*p)*sin(p*x)*sin(p*y)*sin(p*z)*4.0)/Re;

        rhs[1]=(B*(p*p)*cos(p*x*(1.0/2.0))*cos(p*y*(1.0/2.0))*cos(p*z*(1.0/2.0))-A*(p*p)*cos(p*x)*cos(p*y)*sin(p*z)+C*(p*p)*cos(p*z*(1.0/2.0))*sin(p*x*(1.0/2.0))*sin(p*y*(1.0/2.0))*(1.0/4.0))/Re;


        rhs[2]=
        -(A*(p*p)*cos(p*x)*cos(p*z)*sin(p*y)+B*(p*p)*cos(p*x*(1.0/2.0))*sin(p*y*(1.0/2.0))*sin(p*z*(1.0/2.0))*(1.0/4.0)-C*(p*p)*cos(p*y*(1.0/2.0))*sin(p*x*(1.0/2.0))*sin(p*z*(1.0/2.0)))/Re;
      */

    //
    // for skew symmetric
    //
    //   double C1 = 0.5;
    //   double C2 = 0.5;

    rhs[0] =
    -( C * ( p * p ) * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 )
       + B * ( p * p ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 )
       - A * ( p * p ) * sin( p * x ) * sin( p * y ) * sin( p * z ) * 4.0 ) / Re
    + C2 * ( ( A * A ) * p * cos( p * x ) * sin( p * x ) * pow( sin( p * y ), 2.0 ) * pow( sin( p * z ), 2.0 )
             + A * B * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) )
               * sin( p * x ) * sin( p * z ) + A * C * p * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z ) * sin( p * x )
                                               * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y ) * sin( p * z * ( 1.0 / 2.0 ) ) )
    + C1 * ( ( A * A ) * p * cos( p * x ) * sin( p * x ) * pow( sin( p * y ), 2.0 ) * pow( sin( p * z ), 2.0 ) * 2.0
             + A * B * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) )
               * sin( p * x ) * sin( p * z ) - A * B * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x )
                                               * sin( p * y ) * sin( p * y * ( 1.0 / 2.0 ) ) * sin( p * z ) * ( 1.0 / 2.0 )
             + A * C * p * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z ) * sin( p * x ) * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y )
               * sin( p * z * ( 1.0 / 2.0 ) ) + A * C * p * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x )
                                                * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y ) * sin( p * z ) * ( 1.0 / 2.0 ) )
    + D * p * cos( p * x * ( 1.0 / 4.0 ) ) * sin( p * y * ( 1.0 / 4.0 ) ) * sin( p * z * ( 1.0 / 4.0 ) ) * ( 1.0 / 4.0 );

    rhs[1] =
    -C2 * ( ( B * B ) * p * pow( cos( p * x * ( 1.0 / 2.0 ) ), 2.0 ) * cos( p * y * ( 1.0 / 2.0 ) )
            * pow( cos( p * z * ( 1.0 / 2.0 ) ), 2.0 ) * sin( p * y * ( 1.0 / 2.0 ) )
            * ( 1.0 / 2.0 ) + B * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * pow( cos( p * y * ( 1.0 / 2.0 ) ), 2.0 )
                              * sin( p * x * ( 1.0 / 2.0 ) ) * pow( sin( p * z * ( 1.0 / 2.0 ) ), 2.0 )
                              * ( 1.0 / 2.0 ) + A * B * p * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x )
                                                * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y ) * sin( p * z ) * ( 1.0 / 2.0 ) )
    - C1
      * ( ( B * B ) * p * pow( cos( p * x * ( 1.0 / 2.0 ) ), 2.0 ) * cos( p * y * ( 1.0 / 2.0 ) ) * pow( cos( p * z * ( 1.0 / 2.0 ) ), 2.0 )
          * sin( p * y * ( 1.0 / 2.0 ) ) - B * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * pow( cos( p * y * ( 1.0 / 2.0 ) ), 2.0 )
                                           * pow( cos( p * z * ( 1.0 / 2.0 ) ), 2.0 ) * sin( p * x * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 )
          + B * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * pow( cos( p * y * ( 1.0 / 2.0 ) ), 2.0 ) * sin( p * x * ( 1.0 / 2.0 ) )
            * pow( sin( p * z * ( 1.0 / 2.0 ) ), 2.0 ) * ( 1.0 / 2.0 )
          - A * B * p * cos( p * x ) * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) )
            * sin( p * y ) * sin( p * z ) + A * B * p * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x )
                                            * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y ) * sin( p * z ) * ( 1.0 / 2.0 ) )
    + ( B * ( p * p ) * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) )
        - A * ( p * p ) * cos( p * x ) * cos( p * y ) * sin( p * z )
        + C * ( p * p ) * cos( p * z * ( 1.0 / 2.0 ) ) * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) / Re
    + D * p * cos( p * y * ( 1.0 / 4.0 ) ) * sin( p * x * ( 1.0 / 4.0 ) ) * sin( p * z * ( 1.0 / 4.0 ) ) * ( 1.0 / 4.0 );

    rhs[2] =
    -( A * ( p * p ) * cos( p * x ) * cos( p * z ) * sin( p * y )
       + B * ( p * p ) * cos( p * x * ( 1.0 / 2.0 ) ) * sin( p * y * ( 1.0 / 2.0 ) ) * sin( p * z * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 )
       - C * ( p * p ) * cos( p * y * ( 1.0 / 2.0 ) ) * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * z * ( 1.0 / 2.0 ) ) ) / Re
    + C2 * ( ( C * C ) * p * pow( cos( p * y * ( 1.0 / 2.0 ) ), 2.0 ) * cos( p * z * ( 1.0 / 2.0 ) )
             * pow( sin( p * x * ( 1.0 / 2.0 ) ), 2.0 ) * sin( p * z * ( 1.0 / 2.0 ) )
             * ( 1.0 / 2.0 ) - B * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) )
                               * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y * ( 1.0 / 2.0 ) ) * sin( p * z * ( 1.0 / 2.0 ) )
                               * ( 1.0 / 2.0 ) + A * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * sin( p * x )
                                                 * sin( p * y ) * sin( p * z ) * sin( p * z * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) )
    + C1 * ( ( C * C ) * p * pow( cos( p * y * ( 1.0 / 2.0 ) ), 2.0 ) * cos( p * z * ( 1.0 / 2.0 ) )
             * pow( sin( p * x * ( 1.0 / 2.0 ) ), 2.0 ) * sin( p * z * ( 1.0 / 2.0 ) )
             - B * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * cos( p * z * ( 1.0 / 2.0 ) )
               * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y * ( 1.0 / 2.0 ) ) * sin( p * z * ( 1.0 / 2.0 ) )
             + A * C * p * cos( p * x ) * cos( p * y * ( 1.0 / 2.0 ) ) * sin( p * x * ( 1.0 / 2.0 ) ) * sin( p * y ) * sin( p * z )
               * sin( p * z * ( 1.0 / 2.0 ) ) + A * C * p * cos( p * x * ( 1.0 / 2.0 ) ) * cos( p * y * ( 1.0 / 2.0 ) ) * sin( p * x )
                                                * sin( p * y ) * sin( p * z ) * sin( p * z * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) )
    + D * p * cos( p * z * ( 1.0 / 4.0 ) ) * sin( p * x * ( 1.0 / 4.0 ) ) * sin( p * y * ( 1.0 / 4.0 ) ) * ( 1.0 / 4.0 );
}

//=========================================================================
//
//
/*
void Staggered::FluxInviscid(double *qF,double *fc)
{

    fc[0]=qf[0]*qf[0];
    fc[1]=qf[0]*qf[1];
    fc[2]=qf[0]*qf[2];
}
*/
//=========================================================================

void Staggered::residualInviscid( double *U, double *V, double *W, double *P, double *R )
{
    double qf[2];
    double uave[2];
    double flux[3];
    double vave[2];
    double wave[2];
    // cx=0.0, cy=0.0,cz=0.;
    // cout<<cx<<" "<<cy<<" "<<cz<<"dt "<<dt<<endl;

    const double dv = dx * dy * dz;

    for ( uint i = 1; i < longEnd; i++ )
    {
        for ( uint j = 1; j < shortEnd; j++ )
        {
            for ( uint k = 1; k < shortEnd; k++ )
            {
                // East West, x-direction
                // interpolate primitive variables

                uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i - 1, j, k )] );
                uave[1] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i + 1, j, k )] );

                qf[0] = uave[0] * uave[0];
                qf[1] = uave[1] * uave[1];

                flux[0] = ( qf[1] - qf[0] ) * dy * dz;

                // north south, y-direction
                // u ave in y direction
                uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j - 1, k )] );
                uave[1] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j + 1, k )] );
                // v ave in x direction
                vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i - 1, j, k )] );
                vave[1] = 0.5 * ( V[vIdx( i, j + 1, k )] + V[vIdx( i - 1, j + 1, k )] );

                qf[0] = uave[0] * vave[0];
                qf[1] = uave[1] * vave[1];

                flux[1] = ( qf[1] - qf[0] ) * dx * dz;

                // top-bottom

                // u ave in z direction
                uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j, k - 1 )] );
                uave[1] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j, k + 1 )] );
                // w ave in x dirction
                wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i - 1, j, k )] );
                wave[1] = 0.5 * ( W[wIdx( i, j, k + 1 )] + W[wIdx( i - 1, j, k + 1 )] );

                qf[0] = uave[0] * wave[0];
                qf[1] = uave[1] * wave[1];

                flux[2] = ( qf[1] - qf[0] ) * dy * dx;

                if ( P1 == true )
                {
                    R[uIdx( i, j, k )] = C1 * ( flux[0] + flux[1] + flux[2] );
                }
                else
                {
                    R[uIdx( i, j, k )] = C1 * ( flux[0] + flux[1] + flux[2] ) + ( P[pIdx( i, j, k )] - P[pIdx( i - 1, j, k )] ) * dz * dy;
                }
                R[uIdx( i, j, k )] = R[uIdx( i, j, k )] / dv;
            }
        }
    }

    for ( uint i = 1; i < shortEnd; i++ )
    {
        for ( uint j = 1; j < longEnd; j++ )
        {
            for ( uint k = 1; k < shortEnd; k++ )
            {
                // East West, x-direction
                // interpolate primitive variables
                // average in y direction
                //

                // u ave in  y
                uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j - 1, k )] );
                uave[1] = 0.5 * ( U[uIdx( i + 1, j, k )] + U[uIdx( i + 1, j - 1, k )] );
                // v ave in  x
                vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i - 1, j, k )] );
                vave[1] = 0.5 * ( V[vIdx( i + 1, j, k )] + V[vIdx( i, j, k )] );

                qf[0] = uave[0] * vave[0];
                qf[1] = uave[1] * vave[1];

                flux[0] = ( qf[1] - qf[0] ) * dy * dz;

                // north south, y-direction

                // v ave in  y
                vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j - 1, k )] );
                vave[1] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j + 1, k )] );

                qf[0] = vave[0] * vave[0];
                qf[1] = vave[1] * vave[1];

                flux[1] = ( qf[1] - qf[0] ) * dx * dz;

                // top-bottom
                // w ave in y direction
                wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j - 1, k )] );
                wave[1] = 0.5 * ( W[wIdx( i, j, k + 1 )] + W[wIdx( i, j - 1, k + 1 )] );

                // v ave in  z
                vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j, k - 1 )] );
                vave[1] = 0.5 * ( V[vIdx( i, j, k + 1 )] + V[vIdx( i, j, k )] );

                qf[0] = wave[0] * vave[0];
                qf[1] = wave[1] * vave[1];

                flux[2] = ( qf[1] - qf[0] ) * dy * dx;

                if ( P1 == true )
                {
                    R[sizeS + vIdx( i, j, k )] = C1 * ( flux[0] + flux[1] + flux[2] );
                }
                else
                {

                    R[sizeS + vIdx( i, j, k )] = C1 * ( flux[0] + flux[1] + flux[2] ) + ( P[pIdx( i, j, k )] - P[pIdx( i, j - 1, k )] ) * dx
                                                                                        * dz;
                }
                R[sizeS + vIdx( i, j, k )] = R[sizeS + vIdx( i, j, k )] / dv;
            }
        }
    }

    for ( uint i = 1; i < shortEnd; i++ )
    {
        for ( uint j = 1; j < shortEnd; j++ )
        {
            for ( uint k = 1; k < longEnd; k++ )
            {

                // ************************************************************************
                // u ave is z direction
                uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j, k - 1 )] );
                uave[1] = 0.5 * ( U[uIdx( i + 1, j, k )] + U[uIdx( i + 1, j, k - 1 )] );
                // w ave in x dire
                wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i - 1, j, k )] );
                wave[1] = 0.5 * ( W[wIdx( i + 1, j, k )] + W[wIdx( i, j, k )] );

                qf[0] = uave[0] * wave[0];
                qf[1] = uave[1] * wave[1];

                flux[0] = ( qf[1] - qf[0] ) * dy * dz;

                // z-direction

                // top-bottom
                // average in y direction
                //
                wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j - 1, k )] );
                wave[1] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j + 1, k )] );

                // v ave is z direction
                vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j, k - 1 )] );
                vave[1] = 0.5 * ( V[vIdx( i, j + 1, k )] + V[vIdx( i, j + 1, k - 1 )] );

                qf[0] = vave[0] * wave[0];
                qf[1] = vave[1] * wave[1];
                flux[1] = ( qf[1] - qf[0] ) * dz * dx;

                // top-bottom
                //
                wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j, k - 1 )] );
                wave[1] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j, k + 1 )] );

                qf[0] = wave[0] * wave[0];
                qf[1] = wave[1] * wave[1];

                flux[2] = ( qf[1] - qf[0] ) * dx * dy;

                if ( P1 == true )
                {
                    R[2 * sizeS + wIdx( i, j, k )] = C1 * ( flux[0] + flux[1] + flux[2] );
                }
                else
                {
                    R[2 * sizeS + wIdx( i, j, k )] = C1 * ( flux[0] + flux[1] + flux[2] ) + ( P[pIdx( i, j, k )] - P[pIdx( i, j, k - 1 )] )
                                                                                            * dy * dx;
                }

                R[2 * sizeS + wIdx( i, j, k )] = R[2 * sizeS + wIdx( i, j, k )] / dv;
            }
        }
    }
}

/****************************************************************************/

void Staggered::residualInviscidSkew( double *U, double *V, double *W, double *P, double *R )
{
    double s[2];
    double uave[2];
    double flux[3];
    double vave[2];
    double wave[2];
    double dudx[2], dudy[2], dudz[2];
    double dvdx[2], dvdy[2], dvdz[2];
    double dwdx[2], dwdy[2], dwdz[2];

    // cx=0.0, cy=0.0,cz=0.;
    // cout<<cx<<" "<<cy<<" "<<cz<<"dt "<<dt<<endl;

    double r;

    for ( uint i = 1; i < longEnd; i++ )
    {
        for ( uint j = 1; j < shortEnd; j++ )
        {
            for ( uint k = 1; k < shortEnd; k++ )
            {
                // East West, x-direction
                // interpolate primitive variables

                uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i - 1, j, k )] );
                uave[1] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i + 1, j, k )] );

                dudx[0] = ( U[uIdx( i, j, k )] - U[uIdx( i - 1, j, k )] ) / dx;
                dudx[1] = ( U[uIdx( i + 1, j, k )] - U[uIdx( i, j, k )] ) / dx;

                s[0] = dudx[0] * uave[0];
                s[1] = dudx[1] * uave[1];

                flux[0] = 0.5 * ( s[1] + s[0] );

                // north south, y-direction
                // u ave in y direction
                dudy[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j - 1, k )] ) / dy;
                dudy[1] = ( -U[uIdx( i, j, k )] + U[uIdx( i, j + 1, k )] ) / dy;
                // v ave in x direction
                vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i - 1, j, k )] );
                vave[1] = 0.5 * ( V[vIdx( i, j + 1, k )] + V[vIdx( i - 1, j + 1, k )] );

                s[0] = dudy[0] * vave[0];
                s[1] = dudy[1] * vave[1];

                flux[1] = 0.5 * ( s[1] + s[0] );

                // top-bottom

                // u ave in z direction
                dudz[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j, k - 1 )] ) / dz;
                dudz[1] = ( -U[uIdx( i, j, k )] + U[uIdx( i, j, k + 1 )] ) / dz;
                // w ave in x dirction
                wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i - 1, j, k )] );
                wave[1] = 0.5 * ( W[wIdx( i, j, k + 1 )] + W[wIdx( i - 1, j, k + 1 )] );

                s[0] = dudz[0] * wave[0];
                s[1] = dudz[1] * wave[1];

                flux[2] = 0.5 * ( s[1] + s[0] );

                r = flux[0] + flux[1] + flux[2];

                R[uIdx( i, j, k )] = R[uIdx( i, j, k )] + r * C2;
            }
        }
    }

    for ( uint i = 1; i < shortEnd; i++ )
    {
        for ( uint j = 1; j < longEnd; j++ )
        {
            for ( uint k = 1; k < shortEnd; k++ )
            {
                // East West, x-direction
                // interpolate primitive variables
                // average in y direction
                //

                // u ave in  y
                uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j - 1, k )] );
                uave[1] = 0.5 * ( U[uIdx( i + 1, j, k )] + U[uIdx( i + 1, j - 1, k )] );
                // v ave in  x
                dvdx[0] = ( V[vIdx( i, j, k )] - V[vIdx( i - 1, j, k )] ) / dx;
                dvdx[1] = ( V[vIdx( i + 1, j, k )] - V[vIdx( i, j, k )] ) / dx;

                s[0] = uave[0] * dvdx[0];
                s[1] = uave[1] * dvdx[1];

                flux[0] = 0.5 * ( s[1] + s[0] );

                // north south, y-direction

                // v ave in  y
                vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j - 1, k )] );
                vave[1] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j + 1, k )] );

                dvdy[0] = ( V[vIdx( i, j, k )] - V[vIdx( i, j - 1, k )] ) / dy;
                dvdy[1] = ( -V[vIdx( i, j, k )] + V[vIdx( i, j + 1, k )] ) / dy;

                s[0] = vave[0] * dvdy[0];
                s[1] = vave[1] * dvdy[1];

                flux[1] = 0.5 * ( s[1] + s[0] );

                // top-bottom
                // w ave in y direction
                wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j - 1, k )] );
                wave[1] = 0.5 * ( W[wIdx( i, j, k + 1 )] + W[wIdx( i, j - 1, k + 1 )] );

                // v ave in  z
                dvdz[0] = ( V[vIdx( i, j, k )] - V[vIdx( i, j, k - 1 )] ) / dz;
                dvdz[1] = ( V[vIdx( i, j, k + 1 )] - V[vIdx( i, j, k )] ) / dz;

                s[0] = wave[0] * dvdz[0];
                s[1] = wave[1] * dvdz[1];

                flux[2] = ( s[1] + s[0] ) * 0.5;

                r = flux[0] + flux[1] + flux[2];

                R[sizeS + vIdx( i, j, k )] = R[sizeS + vIdx( i, j, k )] + r * C2;
            }
        }
    }

    for ( uint i = 1; i < shortEnd; i++ )
    {
        for ( uint j = 1; j < shortEnd; j++ )
        {
            for ( uint k = 1; k < longEnd; k++ )
            {

                // ************************************************************************
                // u ave is z direction
                uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j, k - 1 )] );
                uave[1] = 0.5 * ( U[uIdx( i + 1, j, k )] + U[uIdx( i + 1, j, k - 1 )] );
                // w ave in x dire
                dwdx[0] = ( W[wIdx( i, j, k )] - W[wIdx( i - 1, j, k )] ) / dx;
                dwdx[1] = ( W[wIdx( i + 1, j, k )] - W[wIdx( i, j, k )] ) / dx;

                s[0] = uave[0] * dwdx[0];
                s[1] = uave[1] * dwdx[1];

                flux[0] = 0.5 * ( s[1] + s[0] );

                // z-direction

                // top-bottom
                // average in y direction
                //
                dwdy[0] = ( W[wIdx( i, j, k )] - W[wIdx( i, j - 1, k )] ) / dy;
                dwdy[1] = ( -W[wIdx( i, j, k )] + W[wIdx( i, j + 1, k )] ) / dy;

                // v ave is z direction
                vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j, k - 1 )] );
                vave[1] = 0.5 * ( V[vIdx( i, j + 1, k )] + V[vIdx( i, j + 1, k - 1 )] );

                s[0] = vave[0] * dwdy[0];
                s[1] = vave[1] * dwdy[1];
                flux[1] = ( s[1] + s[0] ) * 0.5;

                // top-bottom
                //
                wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j, k - 1 )] );
                wave[1] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j, k + 1 )] );

                dwdz[0] = ( W[wIdx( i, j, k )] - W[wIdx( i, j, k - 1 )] ) / dz;
                dwdz[1] = ( -W[wIdx( i, j, k )] + W[wIdx( i, j, k + 1 )] ) / dz;

                s[0] = dwdz[0] * wave[0];
                s[1] = dwdz[1] * wave[1];

                flux[2] = ( s[1] + s[0] ) * 0.5;

                r = flux[0] + flux[1] + flux[2];

                R[2 * sizeS + wIdx( i, j, k )] = R[2 * sizeS + wIdx( i, j, k )] + r * C2;
            }
        }
    }
}

void Staggered::Uxy( double Xa, double Ya, double Za, int i, int j, int k, double *xyz )
{
    double X1 = Xa - dx;
    xyz[0] = X1 + i * dx;
    double Y1 = Ya - dy * 0.5;
    xyz[1] = Y1 + j * dy;
    double Z1 = Za - dz * 0.5;
    xyz[2] = Y1 + k * dz;
}

void Staggered::Vxy( double Xa, double Ya, double Za, int i, int j, int k, double *xyz )
{

    double X1 = Xa - dx * 0.5;
    xyz[0] = X1 + i * dx;
    double Y1 = Ya - dy;
    xyz[1] = Y1 + j * dy;
    double Z1 = Za - dz * 0.5;
    xyz[2] = Z1 + k * dz;
}

void Staggered::Wxy( double Xa, double Ya, double Za, int i, int j, int k, double *xyz )
{

    double X1 = Xa - dx * 0.5;
    xyz[0] = X1 + i * dx;
    double Y1 = Ya - dy * 0.5;
    xyz[1] = Y1 + j * dy;
    double Z1 = Za - dz;
    xyz[2] = Z1 + k * dz;
}

void Staggered::Pxy( double Xa, double Ya, double Za, int i, int j, int k, double *xyz )
{

    double X1 = Xa - dx * 0.5;
    xyz[0] = X1 + i * dx;
    double Y1 = Ya - 0.5 * dy;
    xyz[1] = Y1 + j * dy;
    double Z1 = Za - 0.5 * dz;
    xyz[2] = Z1 + k * dz;
}

void Staggered::exactMMS( double *xyz, double *u1, double *v1, double *w1, double *p1 )
{

    double x = xyz[0];
    double y = xyz[1];
    double z = xyz[2];

    double pi = 4.0 * atan( 1.0 );

    double A = AA;
    double B = 0.1;
    double C = 0.1;
    double D = 0.1;

    *u1 = A * sin( pi * x ) * sin( pi * y ) * sin( pi * z );
    *v1 = B * cos( 0.5 * pi * x ) * cos( 0.5 * pi * y ) * cos( 0.5 * pi * z );
    *w1 = C * sin( 0.5 * pi * x ) * cos( 0.5 * pi * y ) * sin( 0.5 * pi * z );
    *p1 = D * sin( 0.25 * pi * x ) * sin( 0.25 * pi * y ) * sin( 0.25 * pi * z );
}

void Staggered::setExactMMS( double *U, double *V, double *W, double *P, double Xa, double Ya, double Za )
{

    // assign including the entire ghosts
    double x1, y1, z1, u1, v1, w1, p1;
    double xyz[3];

    for ( int i = 0; i < longEnd + 1; i++ )

    {
        for ( int j = 0; j < shortEnd + 1; j++ )
        {

            for ( int k = 0; k < shortEnd + 1; k++ )
            {

                Uxy( Xa, Ya, Za, i, j, k, xyz );
                exactMMS( xyz, &u1, &v1, &w1, &p1 );
                U[uIdx( i, j, k )] = u1;
            }
        }
    }

    for ( int i = 0; i < shortEnd + 1; i++ )

    {
        for ( int j = 0; j < longEnd + 1; j++ )
        {

            for ( int k = 0; k < shortEnd + 1; k++ )
            {

                Vxy( Xa, Ya, Za, i, j, k, xyz );
                exactMMS( xyz, &u1, &v1, &w1, &p1 );
                V[vIdx( i, j, k )] = v1;
            }
        }
    }
    for ( int i = 0; i < shortEnd + 1; i++ )

    {
        for ( int j = 0; j < shortEnd + 1; j++ )
        {

            for ( int k = 0; k < longEnd + 1; k++ )
            {

                Wxy( Xa, Ya, Za, i, j, k, xyz );
                exactMMS( xyz, &u1, &v1, &w1, &p1 );
                W[wIdx( i, j, k )] = w1;
            }
        }
    }

    for ( int i = 0; i < shortEnd + 1; i++ )

    {
        for ( int j = 0; j < shortEnd + 1; j++ )
        {

            for ( int k = 0; k < shortEnd + 1; k++ )
            {

                Pxy( Xa, Ya, Za, i, j, k, xyz );
                exactMMS( xyz, &u1, &v1, &w1, &p1 );
                P[pIdx( i, j, k )] = p1;
            }
        }
    }
}

uint Staggered::getSizeS()
{
    return ( sizeS );
}

uint Staggered::getSizeP()
{
    return ( sizeP );
}
void Staggered::sourceAdd( double Xa, double Ya, double Za, double Re, double *R )
{
    const double dv = dx * dy * dz;

    double x1, y1, z1, u1, v1, w1, p1;
    double xyz[3];
    double rhs[3];

    for ( int i = 1; i < longEnd; i++ )

    {
        for ( int j = 1; j < shortEnd; j++ )
        {

            for ( int k = 1; k < shortEnd; k++ )
            {

                Uxy( Xa, Ya, Za, i, j, k, xyz );
                sourceMMS( Re, xyz, rhs );
                R[uIdx( i, j, k )] = rhs[0];
            }
        }
    }
    for ( int i = 1; i < shortEnd; i++ )

    {
        for ( int j = 1; j < longEnd; j++ )
        {

            for ( int k = 1; k < shortEnd; k++ )
            {

                Vxy( Xa, Ya, Za, i, j, k, xyz );
                sourceMMS( Re, xyz, rhs );
                R[sizeS + vIdx( i, j, k )] = rhs[1];
            }
        }
    }

    for ( int i = 1; i < shortEnd; i++ )

    {
        for ( int j = 1; j < shortEnd; j++ )
        {

            for ( int k = 1; k < longEnd; k++ )
            {

                Wxy( Xa, Ya, Za, i, j, k, xyz );
                sourceMMS( Re, xyz, rhs );
                R[2 * sizeS + wIdx( i, j, k )] = rhs[2];
            }
        }
    }
}

void Staggered::getResNorm( double *R, double *del_u )
{

    double res[3] = {0.0, 0.0, 0.0};

    for ( uint i = 0; i < sizeS; i++ )
    {
        res[0] += R[i] * R[i];
        res[1] += R[sizeS + i] * R[sizeS + i];
        res[2] += R[2 * sizeS + i] * R[2 * sizeS + i];
    }

    /*
            for ( uint i = 0; i < sizeS; i++ )
            {
                res[0] = res[0] + ( un[i] - up[i] ) * ( un[i] - up[i] );
                res[1] = res[1] + ( vn[i] - vp[i] ) * ( vn[i] - vp[i] );
                res[2] = res[2] + ( wn[i] - vp[i] ) * ( vn[i] - vp[i] );
            }
    */
    //      cout << "res " << sqrt( res[0] / sizeS ) << " " << sqrt( res[1] / sizeS ) << " " << sqrt( res[2] / sizeS ) << endl;

    *del_u = sqrt( ( res[0] + res[1] + res[2] ) / 3. / sizeS );

    //  *del_u=sqrt((res[1])/3./sizeS);
}

void Staggered::debug( double Re, double Xa, double Ya, double Za, double *U, double *V, double *W, double *P )
{
    int i, j, k;

    i = 4;

    j = 5;
    k = 6;

    double rhs[3];
    double xyz[3];
    double qf[2];
    double uave[2];
    double flux[3];
    double vave[2];
    double wave[2];

    double R;
    //
    double a = 0.5;
    double b = 0.5;
    double c = 0.5;
    double xyz1[3];
    int idx[3];
    for ( int i = 0; i < shortEnd; i++ )
    {
        Pxy( Xa, Ya, Za, i, j, k, xyz );
        Pxy( Xa, Ya, Za, i + 1, j, k, xyz1 );
        if ( xyz[0] < a && a < xyz1[0] )
        {
            idx[0] = i;
            break;
        }
    }
    for ( int j = 0; j < shortEnd; j++ )
    {
        Pxy( Xa, Ya, Za, i, j, k, xyz );
        Pxy( Xa, Ya, Za, i, j + 1, k, xyz1 );
        if ( xyz[1] < b && b < xyz1[1] )
        {
            idx[1] = j;
            break;
        }
    }

    for ( int k = 0; k < shortEnd; k++ )
    {
        Pxy( Xa, Ya, Za, i, j, k, xyz );
        Pxy( Xa, Ya, Za, i, j, k + 1, xyz1 );
        if ( xyz[2] < c && c < xyz1[2] )
        {
            idx[2] = k;
            break;
        }
    }

    i = idx[0];
    j = idx[1];
    k = idx[2];

    Uxy( Xa, Ya, Za, i, j, k, xyz );

    cout << "U(i,j,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << U[uIdx( i, j, k )] << endl;

    Uxy( Xa, Ya, Za, i + 1, j, k, xyz );
    cout << "U(i+1,j,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << U[uIdx( i + 1, j, k )] << endl;
    Uxy( Xa, Ya, Za, i, j + 1, k, xyz );
    cout << "U(i,j+1,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << U[uIdx( i, j + 1, k )] << endl;
    Uxy( Xa, Ya, Za, i, j, k + 1, xyz );
    cout << "U(i,j,k+1) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << U[uIdx( i, j, k + 1 )] << endl;

    Vxy( Xa, Ya, Za, i, j, k, xyz );
    cout << "V(i,j,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << V[vIdx( i, j, k )] << endl;
    Vxy( Xa, Ya, Za, i + 1, j, k, xyz );
    cout << "V(i+1,j,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << V[vIdx( i + 1, j, k )] << endl;
    Vxy( Xa, Ya, Za, i, j + 1, k, xyz );
    cout << "V(i,j+1,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << V[vIdx( i, j + 1, k )] << endl;
    Vxy( Xa, Ya, Za, i, j, k + 1, xyz );
    cout << "V(i,j,k+1) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << V[vIdx( i, j, k + 1 )] << endl;

    Wxy( Xa, Ya, Za, i, j, k, xyz );
    cout << "W(i,j,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << W[wIdx( i, j, k )] << endl;
    Wxy( Xa, Ya, Za, i + 1, j, k, xyz );
    cout << "W(i+1,j,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << W[wIdx( i + 1, j, k )] << endl;
    Wxy( Xa, Ya, Za, i, j + 1, k, xyz );
    cout << "W(i,j+1,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << W[wIdx( i, j + 1, k )] << endl;
    Wxy( Xa, Ya, Za, i, j, k + 1, xyz );
    cout << "W(i,j,k+1) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << W[wIdx( i, j, k + 1 )] << endl;

    Pxy( Xa, Ya, Za, i, j, k, xyz );
    cout << "P(i,j,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << P[pIdx( i, j, k )] << endl;
    Pxy( Xa, Ya, Za, i + 1, j, k, xyz );
    cout << "P(i+1,j,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << P[pIdx( i + 1, j, k )] << endl;
    Pxy( Xa, Ya, Za, i, j + 1, k, xyz );
    cout << "P(i,j+1,k) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << P[pIdx( i, j + 1, k )] << endl;
    Pxy( Xa, Ya, Za, i, j, k + 1, xyz );
    cout << "P(i,j,k+1) coord " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " value " << P[pIdx( i, j, k + 1 )] << endl;

    double dudx[2];
    double dvdx[2];
    double dwdx[2];
    double dudy[2];
    double dvdy[2];
    double dwdy[2];
    double dudz[2];
    double dvdz[2];
    double dwdz[2];

// debug U momentum fluxes
//

#if ( 1 )
    uave[1] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i + 1, j, k )] );
    uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i - 1, j, k )] );

    qf[0] = 0.5 * uave[0] * uave[0];
    qf[1] = 0.5 * uave[1] * uave[1];

    flux[0] = ( qf[1] - qf[0] ) * dy * dz;

    // north south, y-direction
    // u ave in y direction
    uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j - 1, k )] );

    uave[1] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j + 1, k )] );
    // v ave in x direction
    vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i - 1, j, k )] );
    vave[1] = 0.5 * ( V[vIdx( i, j + 1, k )] + V[vIdx( i - 1, j + 1, k )] );

    qf[0] = uave[0] * vave[0];
    qf[1] = uave[1] * vave[1];

    flux[1] = ( qf[1] - qf[0] ) * dx * dz;

    // top-bottom

    // u ave in z direction
    uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j, k - 1 )] );
    uave[1] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j, k + 1 )] );
    // w ave in x dirction
    wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i - 1, j, k )] );
    wave[1] = 0.5 * ( W[wIdx( i, j, k + 1 )] + W[wIdx( i - 1, j, k + 1 )] );

    qf[0] = uave[0] * wave[0];
    qf[1] = uave[1] * wave[1];

    flux[2] = ( qf[1] - qf[0] ) * dy * dx;

    //  R[uIdx( i, j, k )] = flux[0] + flux[1] + flux[2];

    R = flux[0] + flux[1] + flux[2] + ( P[pIdx( i, j, k )] - P[pIdx( i - 1, j, k )] ) * dz * dy;
    R = R / dx / dy / dz;

    /*********************************************
    *
    *               viscous flux calculation
    *
    * *******************************************
    */
    dudx[0] = ( U[uIdx( i, j, k )] - U[uIdx( i - 1, j, k )] ) / dx;
    dudx[1] = ( U[uIdx( i + 1, j, k )] - U[uIdx( i, j, k )] ) / dx;

    qf[0] = 2. * dudx[0];
    qf[1] = 2. * dudx[1];

    flux[0] = ( qf[1] - qf[0] ) * dy * dz;

    dudy[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j - 1, k )] ) / dy;
    dudy[1] = ( U[uIdx( i, j + 1, k )] - U[uIdx( i, j, k )] ) / dy;
    // v ave in x direction
    dvdx[0] = ( V[vIdx( i, j, k )] - V[vIdx( i - 1, j, k )] ) / dx;
    dvdx[1] = ( V[vIdx( i, j + 1, k )] - V[vIdx( i - 1, j + 1, k )] ) / dx;

    qf[0] = dudy[0] + dvdx[0];
    qf[1] = dudy[1] + dvdx[1];

    flux[1] = ( qf[1] - qf[0] ) * dx * dz;

    // top-bottom

    // u ave in z direction
    dudz[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j, k - 1 )] ) / dz;
    dudz[1] = ( U[uIdx( i, j, k + 1 )] - U[uIdx( i, j, k )] ) / dz;
    // w ave in x dirction
    dwdx[0] = ( W[wIdx( i, j, k )] - W[wIdx( i - 1, j, k )] ) / dx;
    dwdx[1] = ( W[wIdx( i, j, k + 1 )] - W[wIdx( i - 1, j, k + 1 )] ) / dx;

    qf[0] = dudz[0] + dwdx[0];
    qf[1] = dudz[1] + dwdx[1];

    flux[2] = ( qf[1] - qf[0] ) * dy * dx;

    //    R[uIdx( i, j, k )] = flux[0] + flux[1] + flux[2];

    cout << "Conv " << R << endl;
    double Rv = 1. / Re * ( flux[0] + flux[1] + flux[2] ) / dx / dy / dz;

    Uxy( Xa, Ya, Za, i, j, k, xyz );

    cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
    sourceMMS( Re, xyz, rhs );
    cout << "source " << rhs[0] << endl;
    cout << "viscous " << Rv << endl;
    cout << " Res X " << R - Rv - rhs[0] << endl;

    cout << "  " << endl;

#endif

    /*****

    *
    *  y-direction
    *
    * **************************************************************************.
    * **************************************************************************.
    * **************************************************************************.
    *
    */

    // u ave in  y
    uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j - 1, k )] );
    uave[1] = 0.5 * ( U[uIdx( i + 1, j, k )] + U[uIdx( i + 1, j - 1, k )] );
    // v ave in  x
    vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i - 1, j, k )] );
    vave[1] = 0.5 * ( V[vIdx( i + 1, j, k )] + V[vIdx( i, j, k )] );

    qf[0] = uave[0] * vave[0];
    qf[1] = uave[1] * vave[1];

    flux[0] = ( qf[1] - qf[0] ) * dy * dz;

    // north south, y-direction

    // v ave in  y
    vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j - 1, k )] );
    vave[1] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j + 1, k )] );

    qf[0] = 0.5 * vave[0] * vave[0];
    qf[1] = 0.5 * vave[1] * vave[1];

    flux[1] = ( qf[1] - qf[0] ) * dx * dz;

    // top-bottom
    // w ave in y direction
    wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j - 1, k )] );
    wave[1] = 0.5 * ( W[wIdx( i, j, k + 1 )] + W[wIdx( i, j - 1, k + 1 )] );

    // v ave in  z
    vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j, k - 1 )] );
    vave[1] = 0.5 * ( V[vIdx( i, j, k + 1 )] + V[vIdx( i, j, k )] );

    qf[0] = wave[0] * vave[0];
    qf[1] = wave[1] * vave[1];

    flux[2] = ( qf[1] - qf[0] ) * dy * dx;

    //                R[sizeS + vIdx( i, j, k )] = flux[0] + flux[1] + flux[2];
    R = flux[0] + flux[1] + flux[2] + ( P[pIdx( i, j, k )] - P[pIdx( i, j - 1, k )] ) * dx * dz;

    R = R / dx / dy / dz;

    Vxy( Xa, Ya, Za, i, j, k, xyz );

    //                cout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<endl;
    // u ave in  y
    //
    //
    //
    dudy[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j - 1, k )] ) / dy;
    dudy[1] = ( U[uIdx( i + 1, j, k )] - U[uIdx( i + 1, j - 1, k )] ) / dy;
    // v ave in  x
    dvdx[0] = ( V[vIdx( i, j, k )] - V[vIdx( i - 1, j, k )] ) / dx;
    dvdx[1] = ( V[vIdx( i + 1, j, k )] - V[vIdx( i, j, k )] ) / dx;

    qf[0] = dudy[0] + dvdx[0];
    qf[1] = dudy[1] + dvdx[1];

    flux[0] = ( qf[1] - qf[0] ) * dy * dz;

    // north south, y-direction

    // v ave in  y
    dvdy[0] = ( V[vIdx( i, j, k )] - V[vIdx( i, j - 1, k )] ) / dy;
    dvdy[1] = ( V[vIdx( i, j + 1, k )] - V[vIdx( i, j, k )] ) / dy;

    qf[0] = 2.0 * dvdy[0];

    qf[1] = 2.0 * dvdy[1];

    flux[1] = ( qf[1] - qf[0] ) * dx * dz;

    // top-bottom
    // w ave in y direction
    dwdy[0] = ( W[wIdx( i, j, k )] - W[wIdx( i, j - 1, k )] ) / dx;
    dwdy[1] = ( W[wIdx( i, j, k + 1 )] - W[wIdx( i, j - 1, k + 1 )] ) / dx;

    // v ave in  z
    dvdz[0] = ( V[vIdx( i, j, k )] - V[vIdx( i, j, k - 1 )] ) / dz;
    dvdz[1] = ( V[vIdx( i, j, k + 1 )] - V[vIdx( i, j, k )] ) / dz;

    qf[0] = dwdy[0] + dvdz[0];
    qf[1] = dwdy[1] + dvdz[1];

    flux[2] = ( qf[1] - qf[0] ) * dy * dx;

    Rv = 1. / Re * ( flux[0] + flux[1] + flux[2] ) / dx / dy / dz;
    sourceMMS( Re, xyz, rhs );
    cout << "source Y " << rhs[1] << endl;
    cout << "conv Y " << R << endl;
    cout << "viscous " << Rv << endl;
    cout << " Res Y " << R - Rv - rhs[1] << endl;

    cout << "  " << endl;
    // Z
    //
    //
    //
    //

    /*****

    *
    *  z-direction
    *
    * **************************************************************************.
    * **************************************************************************.
    * **************************************************************************.
    *
    */
    // ************************************************************************

    // u ave is z direction
    uave[0] = 0.5 * ( U[uIdx( i, j, k )] + U[uIdx( i, j, k - 1 )] );
    uave[1] = 0.5 * ( U[uIdx( i + 1, j, k )] + U[uIdx( i + 1, j, k - 1 )] );
    // w ave in x dire
    wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i - 1, j, k )] );
    wave[1] = 0.5 * ( W[wIdx( i + 1, j, k )] + W[wIdx( i, j, k )] );

    qf[0] = uave[0] * wave[0];
    qf[1] = uave[1] * wave[1];

    flux[0] = ( qf[1] - qf[0] ) * dy * dz;

    // z-direction

    // top-bottom
    // average in y direction
    //
    wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j - 1, k )] );
    wave[1] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j + 1, k )] );

    // v ave is z direction
    vave[0] = 0.5 * ( V[vIdx( i, j, k )] + V[vIdx( i, j, k - 1 )] );
    vave[1] = 0.5 * ( V[vIdx( i, j + 1, k )] + V[vIdx( i, j + 1, k - 1 )] );

    qf[0] = vave[0] * wave[0];
    qf[1] = vave[1] * wave[1];
    flux[1] = ( qf[1] - qf[0] ) * dz * dx;

    // top-bottom
    //
    wave[0] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j, k - 1 )] );
    wave[1] = 0.5 * ( W[wIdx( i, j, k )] + W[wIdx( i, j, k + 1 )] );

    qf[0] = 0.5 * wave[0] * wave[0];
    qf[1] = 0.5 * wave[1] * wave[1];

    flux[2] = ( qf[1] - qf[0] ) * dx * dy;

    R = flux[0] + flux[1] + flux[2] + ( P[pIdx( i, j, k )] - P[pIdx( i, j, k - 1 )] ) * dx * dy;
    R = R / dx / dy / dz;

    dudz[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j, k - 1 )] ) / dz;
    dudz[1] = ( U[uIdx( i + 1, j, k )] - U[uIdx( i + 1, j, k - 1 )] ) / dz;
    // w ave in x dire
    dwdx[0] = ( W[wIdx( i, j, k )] - W[wIdx( i - 1, j, k )] ) / dx;
    dwdx[1] = ( W[wIdx( i + 1, j, k )] - W[wIdx( i, j, k )] ) / dx;

    qf[0] = dudz[0] + dwdx[0];
    qf[1] = dudz[1] + dwdx[1];

    flux[0] = ( qf[1] - qf[0] ) * dy * dz;

    // z-direction

    // top-bottom
    // average in y direction
    //
    dwdy[0] = ( W[wIdx( i, j, k )] - W[wIdx( i, j - 1, k )] ) / dy;
    dwdy[1] = ( W[wIdx( i, j + 1, k )] - W[wIdx( i, j, k )] ) / dy;

    // v ave is z direction
    dvdz[0] = ( V[vIdx( i, j, k )] - V[vIdx( i, j, k - 1 )] ) / dz;
    dvdz[1] = ( V[vIdx( i, j + 1, k )] - V[vIdx( i, j + 1, k - 1 )] ) / dz;

    qf[0] = dvdz[0] + dwdy[0];
    qf[1] = dvdz[1] + dwdy[1];
    flux[1] = ( qf[1] - qf[0] ) * dz * dx;

    // top-bottom
    //
    dwdz[0] = ( W[wIdx( i, j, k )] - W[wIdx( i, j, k - 1 )] ) / dz;
    dwdz[1] = ( W[wIdx( i, j, k + 1 )] - W[wIdx( i, j, k )] ) / dz;

    qf[0] = 2.0 * dwdz[0];
    qf[1] = 2.0 * dwdz[1];

    flux[2] = ( qf[1] - qf[0] ) * dx * dy;

    //                cout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<endl;
    cout << " Res-Z " << R - rhs[2] << endl;
    Rv = 1. / Re * ( flux[0] + flux[1] + flux[2] ) / dx / dy / dz;
    Wxy( Xa, Ya, Za, i, j, k, xyz );
    sourceMMS( Re, xyz, rhs );
    cout << "source Z " << rhs[2] << endl;
    cout << "conv Z " << R << endl;
    cout << "viscous " << Rv << endl;
    cout << " Res Z " << R - Rv - rhs[2] << endl;
}

void Staggered::residualViscous( double Re, double *U, double *V, double *W, double *P, double *R )
{
    /*
            double cx = 1. / Re / dx / dx;
            double cy = 1. / Re / dy / dy;
            double cz = 1. / Re / dz / dz;
    */

    double qf[2];
    double dudx[2];
    double dvdx[2];
    double dwdx[2];
    double dudy[2];
    double dvdy[2];
    double dwdy[2];
    double dudz[2];
    double dvdz[2];
    double dwdz[2];

    double flux[3];
    double vave[2];
    double wave[2];
    // cx=0.0, cy=0.0,cz=0.;
    // cout<<cx<<" "<<cy<<" "<<cz<<"dt "<<dt<<endl;

    static const double dv = dx * dy * dz;

    for ( uint i = 1; i < longEnd; i++ )
    {
        for ( uint j = 1; j < shortEnd; j++ )
        {
            for ( uint k = 1; k < shortEnd; k++ )
            {
                // East West, x-direction
                // interpolate primitive variables

                dudx[0] = ( U[uIdx( i, j, k )] - U[uIdx( i - 1, j, k )] ) / dx;
                dudx[1] = ( U[uIdx( i + 1, j, k )] - U[uIdx( i, j, k )] ) / dx;

                qf[0] = 2. * dudx[0];
                qf[1] = 2. * dudx[1];

                flux[0] = ( qf[1] - qf[0] ) * dy * dz;

                dudy[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j - 1, k )] ) / dy;
                dudy[1] = ( U[uIdx( i, j + 1, k )] - U[uIdx( i, j, k )] ) / dy;
                // v ave in x direction
                dvdx[0] = ( V[vIdx( i, j, k )] - V[vIdx( i - 1, j, k )] ) / dx;
                dvdx[1] = ( V[vIdx( i, j + 1, k )] - V[vIdx( i - 1, j + 1, k )] ) / dx;

                qf[0] = dudy[0] + dvdx[0];
                qf[1] = dudy[1] + dvdx[1];

                flux[1] = ( qf[1] - qf[0] ) * dx * dz;

                // top-bottom

                // u ave in z direction
                dudz[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j, k - 1 )] ) / dz;
                dudz[1] = ( U[uIdx( i, j, k + 1 )] - U[uIdx( i, j, k )] ) / dz;
                // w ave in x dirction
                dwdx[0] = ( W[wIdx( i, j, k )] - W[wIdx( i - 1, j, k )] ) / dx;
                dwdx[1] = ( W[wIdx( i, j, k + 1 )] - W[wIdx( i - 1, j, k + 1 )] ) / dx;

                qf[0] = dudz[0] + dwdx[0];
                qf[1] = dudz[1] + dwdx[1];

                flux[2] = ( qf[1] - qf[0] ) * dy * dx;

                //    R[uIdx( i, j, k )] = flux[0] + flux[1] + flux[2];

                R[uIdx( i, j, k )] = flux[0] + flux[1] + flux[2];

                R[uIdx( i, j, k )] = 1. / Re * R[uIdx( i, j, k )] / dv;
            }
        }
    }

    for ( uint i = 1; i < shortEnd; i++ )
    {
        for ( uint j = 1; j < longEnd; j++ )
        {
            for ( uint k = 1; k < shortEnd; k++ )
            {
                // East West, x-direction
                // interpolate primitive variables
                // average in y direction
                //

                dudy[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j - 1, k )] ) / dy;
                dudy[1] = ( U[uIdx( i + 1, j, k )] - U[uIdx( i + 1, j - 1, k )] ) / dy;
                // v ave in  x
                dvdx[0] = ( V[vIdx( i, j, k )] - V[vIdx( i - 1, j, k )] ) / dx;
                dvdx[1] = ( V[vIdx( i + 1, j, k )] - V[vIdx( i, j, k )] ) / dx;

                qf[0] = dudy[0] + dvdx[0];
                qf[1] = dudy[1] + dvdx[1];

                flux[0] = ( qf[1] - qf[0] ) * dy * dz;

                // north south, y-direction

                // v ave in  y
                dvdy[0] = ( V[vIdx( i, j, k )] - V[vIdx( i, j - 1, k )] ) / dy;
                dvdy[1] = ( V[vIdx( i, j + 1, k )] - V[vIdx( i, j, k )] ) / dy;

                qf[0] = 2.0 * dvdy[0];

                qf[1] = 2.0 * dvdy[1];

                flux[1] = ( qf[1] - qf[0] ) * dx * dz;

                // top-bottom
                // w ave in y direction
                dwdy[0] = ( W[wIdx( i, j, k )] - W[wIdx( i, j - 1, k )] ) / dx;
                dwdy[1] = ( W[wIdx( i, j, k + 1 )] - W[wIdx( i, j - 1, k + 1 )] ) / dx;

                // v ave in  z
                dvdz[0] = ( V[vIdx( i, j, k )] - V[vIdx( i, j, k - 1 )] ) / dz;
                dvdz[1] = ( V[vIdx( i, j, k + 1 )] - V[vIdx( i, j, k )] ) / dz;

                qf[0] = dwdy[0] + dvdz[0];
                qf[1] = dwdy[1] + dvdz[1];

                flux[2] = ( qf[1] - qf[0] ) * dy * dx;

                R[sizeS + vIdx( i, j, k )] = flux[0] + flux[1] + flux[2];

                R[sizeS + vIdx( i, j, k )] = 1. / Re * R[sizeS + vIdx( i, j, k )] / dv;
            }
        }
    }

    for ( uint i = 1; i < shortEnd; i++ )
    {
        for ( uint j = 1; j < shortEnd; j++ )
        {
            for ( uint k = 1; k < longEnd; k++ )
            {

                dudz[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j, k - 1 )] ) / dz;
                dudz[1] = ( U[uIdx( i + 1, j, k )] - U[uIdx( i + 1, j, k - 1 )] ) / dz;
                // w ave in x dire
                dwdx[0] = ( W[wIdx( i, j, k )] - W[wIdx( i - 1, j, k )] ) / dx;
                dwdx[1] = ( W[wIdx( i + 1, j, k )] - W[wIdx( i, j, k )] ) / dx;

                qf[0] = dudz[0] + dwdx[0];
                qf[1] = dudz[1] + dwdx[1];

                flux[0] = ( qf[1] - qf[0] ) * dy * dz;

                // z-direction

                // top-bottom
                // average in y direction
                //
                dwdy[0] = ( W[wIdx( i, j, k )] - W[wIdx( i, j - 1, k )] ) / dy;
                dwdy[1] = ( W[wIdx( i, j + 1, k )] - W[wIdx( i, j, k )] ) / dy;

                // v ave is z direction
                dvdz[0] = ( V[vIdx( i, j, k )] - V[vIdx( i, j, k - 1 )] ) / dz;
                dvdz[1] = ( V[vIdx( i, j + 1, k )] - V[vIdx( i, j + 1, k - 1 )] ) / dz;

                qf[0] = dvdz[0] + dwdy[0];
                qf[1] = dvdz[1] + dwdy[1];
                flux[1] = ( qf[1] - qf[0] ) * dz * dx;

                // top-bottom
                //
                dwdz[0] = ( W[wIdx( i, j, k )] - W[wIdx( i, j, k - 1 )] ) / dz;
                dwdz[1] = ( W[wIdx( i, j, k + 1 )] - W[wIdx( i, j, k )] ) / dz;

                qf[0] = 2.0 * dwdz[0];
                qf[1] = 2.0 * dwdz[1];

                flux[2] = ( qf[1] - qf[0] ) * dx * dy;

                R[2 * sizeS + wIdx( i, j, k )] = flux[0] + flux[1] + flux[2];

                R[2 * sizeS + wIdx( i, j, k )] = 1. / Re * R[2 * sizeS + wIdx( i, j, k )] / dv;
            }
        }
    }
}

void Staggered::residualViscousCrossDerivativeOnly( double Re, double *U, double *V, double *W, double *P, double *R )
{

    double qf[2];
    double dudx[2];
    double dvdx[2];
    double dwdx[2];
    double dudy[2];
    double dvdy[2];
    double dwdy[2];
    double dudz[2];
    double dvdz[2];
    double dwdz[2];

    double flux[3];
    double vave[2];
    double wave[2];
    // cx=0.0, cy=0.0,cz=0.;
    // cout<<cx<<" "<<cy<<" "<<cz<<"dt "<<dt<<endl;

    static const double dv = dx * dy * dz;

    // in x-direction

    for ( uint i = 1; i < longEnd; i++ )
    {
        for ( uint j = 1; j < shortEnd; j++ )
        {
            for ( uint k = 1; k < shortEnd; k++ )
            {
                // East West, x-direction
                // interpolate primitive variables

                flux[0] = 0.0;

                // y - direction

                dudy[0] = 0.0;
                dudy[1] = 0.0;

                dvdx[0] = ( V[vIdx( i, j, k )] - V[vIdx( i - 1, j, k )] ) / dx;
                dvdx[1] = ( V[vIdx( i, j + 1, k )] - V[vIdx( i - 1, j + 1, k )] ) / dx;

                qf[0] = dudy[0] + dvdx[0];
                qf[1] = dudy[1] + dvdx[1];

                flux[1] = ( qf[1] - qf[0] ) * dx * dz;

                // top-bottom

                // u ave in z direction
                dudz[0] = 0.0;
                dudz[1] = 0.0;
                // w ave in x dirction
                dwdx[0] = ( W[wIdx( i, j, k )] - W[wIdx( i - 1, j, k )] ) / dx;
                dwdx[1] = ( W[wIdx( i, j, k + 1 )] - W[wIdx( i - 1, j, k + 1 )] ) / dx;

                qf[0] = dudz[0] + dwdx[0];
                qf[1] = dudz[1] + dwdx[1];

                flux[2] = ( qf[1] - qf[0] ) * dy * dx;

                //    R[uIdx( i, j, k )] = flux[0] + flux[1] + flux[2];

                R[uIdx( i, j, k )] = flux[0] + flux[1] + flux[2];

                R[uIdx( i, j, k )] = 1. / Re * R[uIdx( i, j, k )] / dv;
            }
        }
    }

    for ( uint i = 1; i < shortEnd; i++ )
    {
        for ( uint j = 1; j < longEnd; j++ )
        {
            for ( uint k = 1; k < shortEnd; k++ )
            {
                // East West, x-direction
                // interpolate primitive variables
                // average in y direction
                //

                dudy[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j - 1, k )] ) / dy;
                dudy[1] = ( U[uIdx( i + 1, j, k )] - U[uIdx( i + 1, j - 1, k )] ) / dy;
                // v ave in  x
                dvdx[0] = 0.0;
                dvdx[1] = 0.0;

                qf[0] = dudy[0] + dvdx[0];
                qf[1] = dudy[1] + dvdx[1];

                flux[0] = ( qf[1] - qf[0] ) * dy * dz;

                // north south, y-direction

                // v ave in  y

                flux[1] = 0.0;

                // top-bottom
                // w ave in y direction

                dwdy[0] = ( W[wIdx( i, j, k )] - W[wIdx( i, j - 1, k )] ) / dx;
                dwdy[1] = ( W[wIdx( i, j, k + 1 )] - W[wIdx( i, j - 1, k + 1 )] ) / dx;

                // v ave in  z
                dvdz[0] = 0.0;
                dvdz[1] = 0.0;

                qf[0] = dwdy[0] + dvdz[0];
                qf[1] = dwdy[1] + dvdz[1];

                flux[2] = ( qf[1] - qf[0] ) * dy * dx;

                R[sizeS + vIdx( i, j, k )] = flux[0] + flux[1] + flux[2];

                R[sizeS + vIdx( i, j, k )] = 1. / Re * R[sizeS + vIdx( i, j, k )] / dv;
            }
        }
    }

    for ( uint i = 1; i < shortEnd; i++ )
    {
        for ( uint j = 1; j < shortEnd; j++ )
        {
            for ( uint k = 1; k < longEnd; k++ )
            {

                dudz[0] = ( U[uIdx( i, j, k )] - U[uIdx( i, j, k - 1 )] ) / dz;
                dudz[1] = ( U[uIdx( i + 1, j, k )] - U[uIdx( i + 1, j, k - 1 )] ) / dz;
                // w ave in x dire
                dwdx[0] = 0.0;
                dwdx[1] = 0.0;

                qf[0] = dudz[0] + dwdx[0];
                qf[1] = dudz[1] + dwdx[1];

                flux[0] = ( qf[1] - qf[0] ) * dy * dz;

                // z-direction

                // top-bottom
                // average in y direction
                //
                dwdy[0] = 0.0;
                dwdy[1] = 0.0;

                // v ave is z direction
                dvdz[0] = ( V[vIdx( i, j, k )] - V[vIdx( i, j, k - 1 )] ) / dz;
                dvdz[1] = ( V[vIdx( i, j + 1, k )] - V[vIdx( i, j + 1, k - 1 )] ) / dz;

                qf[0] = dvdz[0] + dwdy[0];
                qf[1] = dvdz[1] + dwdy[1];
                flux[1] = ( qf[1] - qf[0] ) * dz * dx;

                // top-bottom
                //

                flux[2] = 0.0;

                R[2 * sizeS + wIdx( i, j, k )] = flux[0] + flux[1] + flux[2];

                R[2 * sizeS + wIdx( i, j, k )] = 1. / Re * R[2 * sizeS + wIdx( i, j, k )] / dv;
            }
        }
    }
}

void Staggered::predict( double *U, double *V, double *W, double *Um, double *Vm, double *Wm, double *R )
{

    for ( int i = 0; i < sizeS; i++ )
    {
        U[i] = Um[i] - dt * ( R[i] );
    }

    for ( int i = 0; i < sizeS; i++ )
    {
        V[i] = Vm[i] - dt * ( R[sizeS + i] );
    }

    for ( int i = 0; i < sizeS; i++ )
    {
        W[i] = Wm[i] - dt * ( R[2 * sizeS + i] );
    }
}

void Staggered::setExactSource( double Re, double *U, double *V, double *W, double Xa, double Ya, double Za )
{

    // assign including the entire ghosts
    double x1, y1, z1, u1, v1, w1, p1;
    double xyz[3];
    double rhs[3];
    for ( int i = 0; i < longEnd + 1; i++ )

    {
        for ( int j = 0; j < shortEnd + 1; j++ )
        {

            for ( int k = 0; k < shortEnd + 1; k++ )
            {

                Uxy( Xa, Ya, Za, i, j, k, xyz );
                sourceMMS( Re, xyz, rhs );
                U[uIdx( i, j, k )] = rhs[0];
            }
        }
    }

    for ( int i = 0; i < shortEnd + 1; i++ )

    {
        for ( int j = 0; j < longEnd + 1; j++ )
        {

            for ( int k = 0; k < shortEnd + 1; k++ )
            {

                Vxy( Xa, Ya, Za, i, j, k, xyz );
                sourceMMS( Re, xyz, rhs );
                V[vIdx( i, j, k )] = rhs[1];
            }
        }

        for ( int i = 0; i < shortEnd + 1; i++ )

        {
            for ( int j = 0; j < shortEnd + 1; j++ )
            {

                for ( int k = 0; k < longEnd + 1; k++ )
                {

                    Wxy( Xa, Ya, Za, i, j, k, xyz );
                    sourceMMS( Re, xyz, rhs );
                    W[wIdx( i, j, k )] = rhs[2];
                }
            }
        }
    }
}

void Staggered::setInteriorToZero( double *U, double *V, double *W, double *P )
{

    for ( int i = 1; i < longEnd; i++ )
    {
        for ( int j = 1; j < shortEnd; j++ )
        {
            for ( int k = 1; k < shortEnd; k++ )
            {

                U[uIdx( i, j, k )] = 0.0;
            }
        }
    }

    for ( int i = 1; i < shortEnd; i++ )
    {
        for ( int j = 1; j < longEnd; j++ )
        {
            for ( int k = 1; k < shortEnd; k++ )
            {

                V[vIdx( i, j, k )] = 0.0;
            }
        }
    }

    for ( int i = 1; i < shortEnd; i++ )
    {
        for ( int j = 1; j < shortEnd; j++ )
        {
            for ( int k = 1; k < longEnd; k++ )
            {

                W[wIdx( i, j, k )] = 0.0;
            }
        }
    }

    if ( P1 == true )
    {
        for ( int i = 1; i < shortEnd; i++ )
        {
            for ( int j = 1; j < shortEnd; j++ )
            {
                for ( int k = 1; k < shortEnd; k++ )
                {
                    P[pIdx( i, j, k )] = 0.0;
                }
            }
        }
        cout << "setting interior pressure to zero! " << endl;
    }
}

void Staggered::setBoundaryLidDrivenCavity( double *U, double *V, double *W, double *P )
{

    // size nmax and nmax+1
    // 0: nmax-1 , 0:nmax
    //*****************************************************
    // set surfaces at xmin to xmax faces, Dirichlet
    //****************************************************

    for ( int j = 1; j < shortEnd; j++ )
    {
        for ( int k = 1; k < shortEnd; k++ )
        {
            U[uIdx( 1, j, k )] =Ubc[0];
            U[uIdx( nmax - 1, j, k )] =Ubc[1];
        }
    }

    for ( int j = 1; j < longEnd; j++ )
    {
        for ( int k = 1; k < shortEnd; k++ )
        {
            V[vIdx( 0, j, k )] = 2.*Vbc[0]-V[vIdx( 1, j, k )];
            V[vIdx( nmax - 1, j, k )] =2.*Vbc[1] -V[vIdx( nmax - 2, j, k )];
        }
    }

    for ( int j = 1; j < shortEnd; j++ )
    {
        for ( int k = 1; k < shortEnd; k++ )
        {
            W[wIdx( 0, j, k )] =2.*Wbc[0] -W[wIdx( 1, j, k )];
            W[wIdx( nmax - 1, j, k )] =2.*Wbc[1] -W[wIdx( nmax - 2, j, k )];
        }
    }

    // Neuman

    for ( int j = 1; j < shortEnd; j++ )
    {
        for ( int k = 1; k < shortEnd; k++ )
        {
            P[pIdx( 0, j, k )] = P[pIdx( 1, j, k )];
            P[pIdx( nmax - 1, j, k )] = P[pIdx( nmax - 2, j, k )];
        }
    }

    //*****************************************************
    // set surfaces at zmin to zmax faces, Periodic
    //****************************************************
    // ghost values are updated in a periodic fashion

    for ( int i = 1; i < longEnd; i++ )
    {
        for ( int j = 1; j < shortEnd; j++ )
        {
            U[uIdx( i, j, nmax - 1 )] = U[uIdx( i, j, 1 )];
            U[uIdx( i, j, 0 )] = U[uIdx( i, j, nmax - 2 )];
        }
    }

    for ( int i = 1; i < shortEnd; i++ )
    {
        for ( int j = 1; j < longEnd; j++ )
        {
            V[vIdx( i, j, nmax - 1 )] = V[vIdx( i, j, 1 )];
            V[vIdx( i, j, 0 )] = V[vIdx( i, j, nmax - 2 )];
        }
    }

    for ( int i = 1; i < shortEnd; i++ )
    {
        for ( int j = 1; j < shortEnd; j++ )
        {
            W[wIdx( i, j, nmax )] = W[wIdx( i, j, 1 )];
            W[wIdx( i, j, 0 )] = W[wIdx( i, j, nmax - 1 )];
        }
    }

    for ( int i = 1; i < shortEnd; i++ )
    {
        for ( int j = 1; j < shortEnd; j++ )
        {
            P[pIdx( i, j, nmax - 1 )] = P[pIdx( i, j, 1 )];
            P[pIdx( i, j, 0 )] = P[pIdx( i, j, nmax - 2 )];
        }
    }

    //*****************************************************
    // set surfaces at ymin to ymax faces, Periodic
    //****************************************************
    // pressure to be calculated at the top surface

    for ( int i = 1; i < longEnd; i++ )
    {
        for ( int k = 1; k < shortEnd; k++ )
        {
            U[uIdx( i, 0, k )] = 2.0*Ubc[2]-U[uIdx( i, 1, k )];
            U[uIdx(i, nmax-1,k)]=2.0*Ubc[3]-U[uIdx(i, nmax-2,k)];
//            U[uIdx( i, nmax - 1, k )] = 1.0;
        }
    }

    for ( int i = 1; i < shortEnd; i++ )
    {
        for ( int k = 1; k < shortEnd; k++ )
        {
            V[vIdx( i, 1, k )] = Vbc[2];
            V[vIdx( i, nmax - 1, k )] = Vbc[3];
            // V[vIdx(i, nmax-1,k)]=0.0;
        }
    }

    for ( int j = 1; j < shortEnd; j++ )
    {
        for ( int k = 1; k < longEnd; k++ )
        {
            W[wIdx( 0, j, k )] = 2.*Wbc[2]-W[wIdx( 1, j, k )];
            W[wIdx( nmax - 1, j, k )] =2.*Wbc[3] -W[wIdx( nmax - 2, j, k )];
        }
    }

    // Neuman bottom surface
    for ( int i = 1; i < shortEnd; i++ )
    {
        for ( int k = 1; k < shortEnd; k++ )
        {
            P[pIdx( i, 0, k )] = P[pIdx( i, 1, k )];
        }
    }
}

/***************************

*****************************/
//**************************************************************
//
//    use tagginng for the boundary
//
//***********************************************************
/*
void setBoundaryByTag()
{






}
*/
//
//****************************************************************
void Staggered::solveManufacturedSolution( double Re, double Xa, double Ya, double Za, double *U, double *V, double *W, double *P )
{

    double res = 1.0;

    //    double *restrict U = new double[getSizeS()];
    //    double *restrict V = new double[getSizeS()];

    double *restrict Rv = new double[3 * getSizeS()];
    double *restrict S = new double[3 * getSizeS()];

    double *restrict R = new double[3 * getSizeS()];
    //    double *restrict W = new double[getSizeS()];
    //    double *restrict P = new double[getSizeP()];

    double *Um = new double[getSizeS()];
    double *Vm = new double[getSizeS()];
    double *Wm = new double[getSizeS()];

    for ( uint i = 0; i < 3 * getSizeS(); i++ )
    {
        R[i] = 0.0;
        S[i] = 0.0;
        Rv[i] = 0.0;
    }

    for ( uint i = 0; i < getSizeS(); i++ )
    {
        U[i] = 0.0;
        V[i] = 0.0;
        W[i] = 0.0;
        Um[i] = 0.0;
        Vm[i] = 0.0;
        Wm[i] = 0.0;
    }

    setExactMMS( U, V, W, P, Xa, Ya, Za );
    setExactMMS( Um, Vm, Wm, P, Xa, Ya, Za );

    residualInviscid( U, V, W, P, R );
    residualInviscidSkew( U, V, W, P, R );
    residualViscous( Re, U, V, W, P, Rv );
    /*
    for(uint i=0;i<3*q.getSizeS();i++)
    {
    S[i]=(R[i]-Rv[i]);
    }
    */
    setInteriorToZero( U, V, W, P );
    setInteriorToZero( Um, Vm, Wm, P );

    while ( res > 1e-12 )
    {
        residualInviscid( U, V, W, P, R );

        residualInviscidSkew( U, V, W, P, R );

        residualViscous( Re, U, V, W, P, Rv );

        sourceAdd( Xa, Ya, Za, Re, S );

        for ( uint i = 0; i < 3 * getSizeS(); i++ )
        {
            R[i] = R[i] - Rv[i] - S[i];
        }

        predict( U, V, W, Um, Vm, Wm, R );

        /*
        project(U,V,W,P);

        correct(U,V,W,P);

        res=0.0;

        residualInviscidwithP( U, V, W,P,R );

        for(uint i=0;i<3*q.getSizeS();i++)
        {
        R[i]=R[i]-Rv[i]-S[i];
        }
        */

        for ( uint i = 0; i < getSizeS(); i++ )
        {
            // res+=pow(U[i]-Um[i],2)+pow(W[i]-Wm[i],2)+pow(V[i]-Vm[i],2);
            Um[i] = U[i];
            Vm[i] = V[i];
            Wm[i] = W[i];
        }

        getResNorm( R, &res );

        {
            cout << " res " << res << endl;
        }
    }
    delete[] S;
    delete[] R;
    delete[] Rv;
    delete[] Um;
    delete[] Vm;
    delete[] Wm;
}

void Staggered::fluxDebug( double Re, double Xa, double Ya, double Za )
{
    double res = 0.0;

    double *restrict U = new double[getSizeS()];
    double *restrict V = new double[getSizeS()];

    double *restrict Rv = new double[3 * getSizeS()];
    double *restrict S = new double[3 * getSizeS()];

    double *restrict R = new double[3 * getSizeS()];
    double *restrict W = new double[getSizeS()];
    double *restrict P = new double[getSizeP()];

    setExactMMS( U, V, W, P, Xa, Ya, Za );

    residualInviscid( U, V, W, P, R );

    residualInviscidSkew( U, V, W, P, R );

    residualViscous( Re, U, V, W, P, Rv );

    sourceAdd( Xa, Ya, Za, Re, S );

    for ( uint i = 0; i < 3 * getSizeS(); i++ )
    {
        // cout<<R[i]-S[i]<<endl;
        R[i] = R[i] - S[i] - Rv[i];
    }

    getResNorm( R, &res );
    cout << " Res = " << res << endl;
    cout << " if res is not satifactory make sure the P_ON is one" << endl;

    delete[] U;
    delete[] V;
    delete[] P;
    delete[] W;
    delete[] R;
    delete[] Rv;
    delete[] S;
}

void Staggered::solveLidDrivenCavity( double Re, double Xa, double Ya, double Za, double *U, double *V, double *W, double *P )
{

    double res;
    int co = 0;

    //    double *restrict U = new double[getSizeS()];
    //    double *restrict V = new double[getSizeS()];

    double *restrict Rv = new double[3 * getSizeS()];

    double *restrict R = new double[3 * getSizeS()];
    //    double *restrict W = new double[getSizeS()];
    //    double *restrict P = new double[getSizeP()];

    double *Um = new double[getSizeS()];
    double *Vm = new double[getSizeS()];
    double *Wm = new double[getSizeS()];

    double *Pm = new double[getSizeS()];

    for ( uint i = 0; i < 3 * getSizeS(); i++ )
    {
        R[i] = 0.0;
        Rv[i] = 0.0;
    }

    for ( uint i = 0; i < getSizeS(); i++ )
    {
        U[i] = 0.0;
        V[i] = 0.0;
        W[i] = 0.0;
        Um[i] = 0.0;
        Vm[i] = 0.0;
        Wm[i] = 0.0;
    }

    setBoundaryLidDrivenCavity( U, V, W, P );

    setBoundaryLidDrivenCavity( Um, Vm, Wm, P );

    res = 1.0;

    //  cout << "in here" << endl;

    //  for(uint i=0;i<10;i++)
    while ( res > 1e-13 )
    {
        residualInviscid( U, V, W, P, R );

        residualInviscidSkew( U, V, W, P, R );

        residualViscous( Re, U, V, W, P, Rv );

        // q.sourceAdd(Xa,Ya,Za,Re,S);

        for ( uint i = 0; i < 3 * getSizeS(); i++ )
        {
            // R[i]=R[i]-Rv[i]-S[i];
            R[i] = R[i] - Rv[i];
        }

        predict( U, V, W, Um, Vm, Wm, R );

       // setBoundaryLidDrivenCavity(U,V,W,P);

        project( U, V, W, P );

        //        setBoundaryLidDrivenCavity( U, V, W, P );

        correct( U, V, W, P );

        setBoundaryLidDrivenCavity( U, V, W, P );
        res = 0.0;

        for ( uint i = 0; i < getSizeS(); i++ )
        {
            res += pow( U[i] - Um[i], 2 ) + pow( W[i] - Wm[i], 2 ) + pow( V[i] - Vm[i], 2 );
            Um[i] = U[i];
            Vm[i] = V[i];
            Wm[i] = W[i];
        }
        co++;
        //  cout << "in here" << endl;

        //    if ( co % 100 == 0 )
        {
            cout << " res " << res << endl;
        }
    }

    delete[] R;
    delete[] Rv;
    delete[] Um;
    delete[] Vm;
    delete[] Wm;
}

void Staggered::pressureCorrect( double Re, double *Pm, double *phi )
{

    double *restrict Pn = new double[sizeP];

    for ( uint i = 1; i < shortEnd; i++ )
    {
        for ( uint j = 0; j < shortEnd; j++ )
        {
            for ( uint k = 0; k < shortEnd; k++ )
            {
                Pn[pIdx( i, j, k )] = Pm[pIdx( i, j, k )] + phi[pIdx( i, j, k )]
                                      - 0.0 * dt / Re * ( ( phi[pIdx( i + 1, j, k )] + phi[pIdx( i - 1, j, k )] ) / dx / dx
                                                          + ( phi[pIdx( i, j + 1, k )] + phi[pIdx( i, j - 1, k )] ) / dy / dy
                                                          + ( phi[pIdx( i, j, k + 1 )] + phi[pIdx( i, j, k - 1 )] ) / dz / dz );
            }
        }
    }

    delete[] Pn;
}

