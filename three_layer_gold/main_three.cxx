/************************************************************************************************/
/*  3-layer Green's function for electric field                                                 */
/*                                                                                              */
/*                                                                                              */
/*  Ref : 1. Fast Evaluation of Sommerfeld Integrals for EM scattering and radiation by         */
/*        threee-dimensional buried objects                                                     */
/*                                                                                              */
/*           by Tie Jun Cui and Weng Cho Chew                                                   */
/*                                                                                              */
/*           IEEE Trans. on Geoscience and Remote Sensing, 37, 2, p.887, 1999                   */
/*                                                                                              */
/*                                                                                              */
/*        2. Waves and fields in inhomogeneous media                                            */
/*           by W.C. Chew                                                                       */
/*                                                                                              */
/*                                                                                              */
/*                       G_E(r,r')                                                              */
/*                                                                                              */
/*                          * r = (x,y,z) Target points                                         */
/*                                                                                              */
/*                   epsilon[0], mu_1=1                                                         */
/*                   ____________________________________________________  z = 0                */
/*                                                                                              */
/*                   epsilon[1], mu_2=1                                                         */
/*                                                                                              */
/*                                 * r' = (x',y',z')                                            */
/*                                   Dipole source oriented by (alpha_x, alpha_y, alpha_z)      */
/*                                                                                              */
/*                   ____________________________________________________  z = -thickness = -d  */
/*                                                                                              */
/*                   epsilon[2], mu_3=1                                                         */
/*                                                                                              */
/*                                                                                              */
/*     Note      : There are some typos and errors in graph in Ref. 1                           */
/*     Note      : Source (r) and target(r') points can be anywhere now                         */
/*                                                                                              */
/*     mu_0      : permeability in vacuum  (mu_0 = 1)                                           */
/*     epsilon_0 : permittivity in vacuum  (epsilon_0 = 1)                                      */
/*     lambda    : wavelength                                                                   */
/*     k         : wavenumber = (2*pi/lambda)*sqrt(epsilon_0*mu_0)                              */
/*                                                                                              */
/*     Warning : The first layer interface MUST be located at z = 0                             */
/*                                                                                              */
/*     Compile the code by typing                                                               */
/*     >> g++ main_three.cxx complex.cxx -O3                                                    */
/*                                                                                              */
/*     Run the code by typing                                                                   */
/*     >> ./a.out                                                                               */
/*                                                                                              */
/*                                                                                              */
/*     Plot output of the code in Matlab                                                        */
/*     >> field_plot                                                                            */
/*                                                                                              */
/*                                                                                              */
/*                                                                                              */
/*  10/11/2016                                                                                  */
/*  Min Hyung Cho                                                                               */
/*  UMass Lowell                                                                                */
/*                                                                                              */
/************************************************************************************************/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "defs.h"
#include "complex.h"

//Bessel functions
complex J0(complex z);
complex J1(complex z);
complex J2(complex z);

//for the free-space Dyadic Green's function
complex f1(complex x);
complex f2(complex x);

//Quadrature for Sommerfeld integral
void sommerfeld_quadrature_parabola_contour(complex *ks, complex *weight, int n1, int n2);
void sommerfeld_quadrature_complex_contour(complex *ks, complex *weight, int n1, int n2, int n3, int n4);
void sommerfeld_quadrature_gauss_adaptive_two_pole(complex *ks, complex *weight, complex k1, complex k2, complex k3, int n1, int n2, int n3, int n4, int n5, int n6, int n7);
void generalized_gauss_quadrature(int ngq, double *xi, double *wi);
void gauss_quadrature(int n, double *x, double *w);

//new version
complex g5r(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g5t(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g6r(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g6t(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g7r(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g7t(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g8r(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g8t(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g9r(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g9t(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);

// when the source is in the 2nd layer
complex g5r_u(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g5r_d(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g6r_u(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g6r_d(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g7r_u(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g7r_d(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g8r_u(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g8r_d(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g8r_u_v(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);
complex g8r_d_v(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda);


//Each component of singular part of Green's function
complex Gxx_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gxy_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gxz_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gyx_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gyy_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gyz_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gzx_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gzy_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gzz_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);

//Each component of regular part of Green's function
complex Gxx_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gxy_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gxz_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gyx_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gyy_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gyz_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gzx_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gzy_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);
complex Gzz_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness);


int main(void)
{
    int i, j, res;
    double *x, y, *z;
    double xp, yp, zp;
    double lambda, epsilon_0;
    complex epsilon[3];
    double alpha_x, alpha_y, alpha_z;
    complex Gxx, Gyy, Gzz, Gxy, Gyx, Gxz, Gzx, Gyz, Gzy;
    complex Ex, Ey, Ez;
    double thickness;

    struct timeval tv1, tv2;

    FILE *st1, *st2, *st3;

    st1 = fopen("Ex.dat","w");
    st2 = fopen("Ey.dat","w");
    st3 = fopen("Ez.dat","w");

    res       = 100;
    lambda    = 1.0;
    thickness = 1;

    //Field points (x,y,z)
    x = new double[res];
    for(i = 0 ; i < res ; i++){
        x[i] = -5+(double)i*(10.0)/(res-1);
    }
    z = new double[res];
    for(i = 0 ; i < res ; i++){
        z[i] = (-3)+(double)i*(5.0)/(res-1);
    }
    y = 1.0;

    //Dipole source location (xp, yp, zp) and orientation (alpha_x, alpha_y, alphz_z)
    xp =   0;
    yp =  0.0;
    zp =  -0.5;

    //Dipole orientation
    // alpha_x = sin(PI/4.0)*cos(PI/4.0);
    // alpha_y = sin(PI/4.0)*sin(PI/4.0);
    // alpha_z = cos(PI/4.0);
    alpha_x = 1;
    alpha_y = 0;
    alpha_z = 0;


    //epsilon0 in the free-space
    epsilon_0 = 1.0;

    //epsilon in each layer
    epsilon[0] = complex(1,0);
    epsilon[1] = complex(10.9,0);
    epsilon[2] = complex(1,0);

    gettimeofday(&tv1, NULL);

    //Calling the half-space Green's function for all the points (x,y,z) due to the source at (xp, yp, zp)
    for(i = 0 ; i < res ; i++){
        for(j = 0; j < res ; j++){
            printf("Computing i = %d j = %d point\n",i, j);
            Gxx = Gxx_three_singular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness)+
                   Gxx_three_regular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness);
            Gxy = Gxy_three_singular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness)+
                   Gxy_three_regular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness);
            Gxz = Gxz_three_singular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness)+
                   Gxz_three_regular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness);

            Gyx = Gyx_three_singular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness)+
                   Gyx_three_regular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness);
            Gyy = Gyy_three_singular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness)+
                   Gyy_three_regular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness);
            Gyz = Gyz_three_singular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness)+
                   Gyz_three_regular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness);

            Gzx = Gzx_three_singular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness)+
                   Gzx_three_regular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness);
            Gzy = Gzy_three_singular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness)+
                   Gzy_three_regular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness);
            Gzz = Gzz_three_singular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness)+
                   Gzz_three_regular(lambda, epsilon_0, epsilon, x[i], y, z[j], xp, yp, zp, thickness);

            Ex = Gxx*alpha_x+Gxy*alpha_y+Gxz*alpha_z;
            Ey = Gyx*alpha_x+Gyy*alpha_y+Gyz*alpha_z;
            Ez = Gzx*alpha_x+Gzy*alpha_y+Gzz*alpha_z;

            fprintf(st1, "%f %f %15.12f %15.12f\n ",x[i], z[j], Ex.re, Ex.im);
            fprintf(st2, "%f %f %15.12f %15.12f\n ",x[i], z[j], Ey.re, Ey.im);
            fprintf(st3, "%f %f %15.12f %15.12f\n ",x[i], z[j], Ez.re, Ez.im);
        }
    }

    gettimeofday(&tv2, NULL);
    printf("\n Total time: %f seconds \n\n", ((double)(tv2.tv_usec-tv1.tv_usec)/1000000+(double)(tv2.tv_sec-tv1.tv_sec)));

    fclose(st1);
    fclose(st2);
    fclose(st3);

    return 0;
}


complex Gxx_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gxx_r, Gxx_p, Gxx_t, Gxx_u, Gxx_d;
    complex k, k1, k2, k3, ci;
    double omega, rho;;
    complex epsilon_symmetry[3];

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp-0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z-0.0) < 10e-13 ){
            Gxx_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))
                    *(-0.5*g5r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    +( 0.5*rho*rho-(y-yp)*(y-yp))*g6r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxx_r;
        }
        //When the field point is in the 2nd layer
        else if( z < 0.0 && z > -thickness){
            Gxx_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                    *(-0.5*g5r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    +( 0.5*rho*rho-(y-yp)*(y-yp))*g6r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gxx_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                    *(0.5*g5t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    -(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxx_r+Gxx_t;
        }
        //When the field point is in the 3rd layer
        else{
            Gxx_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))
                    *(0.5*g5t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    -(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxx_t;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0.0 && zp > -thickness){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z-0.0) < 10e-13 ){
            Gxx_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))
                    *(0.5*g5t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    -(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxx_t;
        }
        //When the field point is in the 2nd layer
        else if(z < 0.0 && z > -thickness){
            Gxx_u = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                    *(-0.5*g5r_u(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    +( 0.5*rho*rho-(y-yp)*(y-yp))*g6r_u(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gxx_d = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                    *(-0.5*g5r_d(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    +(0.5*rho*rho-(y-yp)*(y-yp))*g6r_d(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxx_u+Gxx_d;
        }
        //When the field point is in the 3rd layer
        else{
            Gxx_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))
                    *(0.5*g5t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    -(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxx_t;
        }
    }
    //When the source point is in the 3rd layer (use symmetry)
    else{
        zp = fabs(zp)-thickness;
        z = -z-thickness;
        epsilon_symmetry[0] = epsilon[2];
        epsilon_symmetry[1] = epsilon[1];
        epsilon_symmetry[2] = epsilon[0];

        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z-0.0) < 10e-13 ){
            Gxx_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))
                      *(-0.5*g5r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda)
                      +( 0.5*rho*rho-(y-yp)*(y-yp))*g6r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gxx_r;
        }
        //When the field point is in the 2nd layer
        else if(z < 0.0 && z > -thickness){
            Gxx_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                      *(-0.5*g5r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda)
                      +( 0.5*rho*rho-(y-yp)*(y-yp))*g6r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));

            Gxx_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                      *(0.5*g5t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda)
                      -(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gxx_r+Gxx_t;
        }
        //When the field point is in the 3rd layer
        else{
            Gxx_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))
                      *(0.5*g5t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda)
                      -(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gxx_t;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gxy_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gxy_r, Gxy_p, Gxy_t, Gxy_u, Gxy_d;
    complex k, k1, k2, k3, ci;
    double omega, rho;
    complex epsilon_symmetry[3];

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp - 0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z-0.0) < 10e-13){
            Gxy_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(x-xp)*(y-yp)*g6r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda);
            Gtemp = Gxy_r;
        }
        //When the field point is in the 2nd layer
        else if(z < 0.0 && z > -thickness){
            Gxy_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(x-xp)*(y-yp)*(g6r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gxy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(x-xp)*(y-yp)*(-g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxy_r+Gxy_t;
        }
        //When the field point is in the 3rd layer
        else{
            Gxy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(x-xp)*(y-yp)*(-g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxy_t;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0.0 && zp > -thickness){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z-0.0) < 10e-13){
            Gxy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(-(x-xp)*(y-yp)*g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxy_t;
        }
        //When the field point is in the 2nd layer
        else if(z < 0.0 && z > -thickness){
            Gxy_u = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(x-xp)*(y-yp)*(g6r_u(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gxy_d = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(x-xp)*(y-yp)*(g6r_d(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxy_u+Gxy_d;
        }
        //When the field point is in the 3rd layer
        else{
            Gxy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(-(x-xp)*(y-yp)*g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxy_t;
        }
    }
    //When the source point is in the 3rd layer (use symmetry)
    else{
        zp = fabs(zp)-thickness;
        z = -z-thickness;
        epsilon_symmetry[0] = epsilon[2];
        epsilon_symmetry[1] = epsilon[1];
        epsilon_symmetry[2] = epsilon[0];

        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z-0.0) < 10e-13){
            Gxy_r = (1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(x-xp)*(y-yp)*g6r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda);
            Gtemp = Gxy_r;
        }
        //When the field point is in the 2nd layer
        else if(z < 0.0 && z > -thickness){
            Gxy_r = (1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(x-xp)*(y-yp)*(g6r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gxy_t = (1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(x-xp)*(y-yp)*(-g6t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gxy_r+Gxy_t;
        }
        else{
            Gxy_t = (1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(x-xp)*(y-yp)*(-g6t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gxy_t;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gxz_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gxz_r, Gxz_p, Gxz_t, Gxz_r2;
    complex Gxz_d, Gxz_u;
    complex k, k1, k2, k3, ci;
    double omega, rho;
    complex epsilon_symmetry[3];

    ci = complex(0.0,1.0);

    omega = 2.0*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp - 0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gxz_r  = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(-ci*(x-xp)*g8r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxz_r;
        }
        //When the field point is in the 2nd layer
        else if(z < 0 && z > -thickness){
            Gxz_r  = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(-ci*(x-xp)*g9r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gxz_t  = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(x-xp)*g8t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxz_r+Gxz_t;;
        }
        //When the field point is in the 3rd layer
        else{
            Gxz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(ci*(x-xp)*g8t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxz_t;
        }
    }
    //When the source point is in the 2nd layer
    else if (zp < 0.0 && zp > -thickness){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gxz_t  = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(-ci*(x-xp)*g8t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxz_t;
        }
        //When the field point is in the 2nd layer
        else if(z < 0 && z > -thickness){
            Gxz_u  = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(-ci*(x-xp)*g8r_u_v(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gxz_d  = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(x-xp)*g8r_d_v(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxz_u+Gxz_d;;
        }
        //When the field point is in the 3rd layer
        else{
            Gxz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(ci*(x-xp)*g8t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gxz_t;
        }
    }
    //When the source point is in the 3rd layer (use symmetry)
    else{
        zp = fabs(zp)-thickness;
        z = -z-thickness;
        epsilon_symmetry[0] = epsilon[2];
        epsilon_symmetry[1] = epsilon[1];
        epsilon_symmetry[2] = epsilon[0];

        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gxz_r  = (1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(-ci*(x-xp)*g8r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gxz_r;
        }
        //When the field point is in the 2nd layer
        else if(z < 0 && z > -thickness){
            Gxz_r  = (1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(-ci*(x-xp)*g9r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gxz_t  = (1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(x-xp)*g8t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gxz_r+Gxz_t;;
        }
        //When the field point is in the 3rd layer
        else{
            Gxz_t = (1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(ci*(x-xp)*g8t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gxz_t;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}




//By symmetry Gyx = Gxy
complex Gyx_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    return Gxy_three_regular(lambda, epsilon_0, epsilon, x, y, z, xp, yp, zp, thickness);
}





complex Gyy_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gyy_r, Gyy_p, Gyy_t, Gyy_u, Gyy_d;
    complex k, k1, k2, k3, ci;
    double omega, rho;
    complex epsilon_symmetry[3];

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp - 0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gyy_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))
                    *(-0.5*g5r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    -(0.5*rho*rho-(y-yp)*(y-yp))*g6r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyy_r;
        }
        //When the field point is in the 2nd layer
        else if(z < 0.0 && z > -thickness){
            Gyy_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                    *(-0.5*g5r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    -(0.5*rho*rho-(y-yp)*(y-yp))*g6r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gyy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                    *(0.5*g5t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    +(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyy_r+Gyy_t;
        }
        //When the field point is in the 3rd layer
        else {
            Gyy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))
                    *(0.5*g5t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    +(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyy_t;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0.0 && zp > -thickness){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z-0.0) < 10e-13 ){
            Gyy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))
                    *(0.5*g5t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    +(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyy_t;
        }
        //When the field point is in the 2nd layer
        else if(z < 0.0 && z > -thickness){
            Gyy_u = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                    *(-0.5*g5r_u(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                      -( 0.5*rho*rho-(y-yp)*(y-yp))*g6r_u(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));

            Gyy_d = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                    *(-0.5*g5r_d(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                      -(0.5*rho*rho-(y-yp)*(y-yp))*g6r_d(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyy_u+Gyy_d;
        }
        //When the field point is in the 3rd layer
        else{
            Gyy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))
                    *(0.5*g5t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda)
                    +(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyy_t;
        }
    }
    //When the source point is in the 3rd layer (use symmetry)
    else{
        zp = fabs(zp)-thickness;
        z = -z-thickness;
        epsilon_symmetry[0] = epsilon[2];
        epsilon_symmetry[1] = epsilon[1];
        epsilon_symmetry[2] = epsilon[0];

        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gyy_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))
                    *(-0.5*g5r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda)
                    -(0.5*rho*rho-(y-yp)*(y-yp))*g6r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gyy_r;
        }
        //When the field point is in the 2nd layer
        else if(z < 0.0 && z > -thickness){
            Gyy_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                    *(-0.5*g5r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda)
                    -(0.5*rho*rho-(y-yp)*(y-yp))*g6r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gyy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))
                    *(0.5*g5t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda)
                    +(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gyy_r+Gyy_t;
        }
        //When the field point is in the 3rd layer
        else {
            Gyy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))
                    *(0.5*g5t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda)
                    +(0.5*rho*rho-(y-yp)*(y-yp))*g6t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gyy_t;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gyz_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gyz_r, Gyz_p, Gyz_t, Gyz_u, Gyz_d;
    complex k, k1, k2, k3, ci;
    double omega, rho;
    complex epsilon_symmetry[3];

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));

    //When the source point is in the 1st layer  (  z' >= 0  )
    if(zp > 0 || fabs(zp-0.0)< 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z - 0.0)< 10e-13 ){
            Gyz_r  = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(-1*ci*(y-yp)*g8r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyz_r;
        }
        //When the field point is in the 2nd layer  (z < 0 && z >= -thickness)
        else if(z < 0  &&  z > -thickness){
            Gyz_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(-ci*(y-yp)*g9r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gyz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(y-yp)*g8t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyz_r+Gyz_t;
        }
        else{
            Gyz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(ci*(y-yp)*g8t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyz_t;
        }
    }

    //When the source point is in the 2nd layer  (  z' >= 0  )
    else if(zp < 0.0 && zp > -thickness){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z - 0.0)< 10e-13 ){
            Gyz_t  = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(-1*ci*(y-yp)*g8t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyz_t;
        }
        //When the field point is in the 2nd layer  (z < 0 && z >= -thickness)
        else if(z < 0  &&  z > -thickness){
            Gyz_u = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(-ci*(y-yp)*g8r_u_v(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gyz_d = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(y-yp)*g8r_d_v(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyz_u+Gyz_d;
        }
        else{
            Gyz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(ci*(y-yp)*g8t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gyz_t;
        }
    }
    //When the source point is in the 3rd layer (use symmetry)
    else{
        zp = fabs(zp)-thickness;
        z = -z-thickness;
        epsilon_symmetry[0] = epsilon[2];
        epsilon_symmetry[1] = epsilon[1];
        epsilon_symmetry[2] = epsilon[0];
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z - 0.0)< 10e-13 ){
            Gyz_r  = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(-1*ci*(y-yp)*g8r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gyz_r;
        }
        //When the field point is in the 2nd layer  (z < 0 && z >= -thickness)
        else if(z < 0  &&  z > -thickness){
            Gyz_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(-ci*(y-yp)*g9r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gyz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(y-yp)*g8t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gyz_r+Gyz_t;
        }
        else{
            Gyz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(ci*(y-yp)*g8t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gyz_t;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}






complex Gzx_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gzx_r, Gzx_p, Gzx_t, Gzx_u, Gzx_d;
    complex k, k1, k2, k3, ci;
    double omega, rho;
    complex epsilon_symmetry[3];

    ci = complex(0.0,1.0);

    omega = 2.0*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp-0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gzx_r  = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(ci*(x-xp)*g8r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzx_r;
        }
        //When the field point is in the 2nd layer
        else if (z < 0.0 && z > -thickness){
            Gzx_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(x-xp)*g8r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gzx_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(x-xp)*g9t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzx_r+Gzx_t;
        }
        //When the field point is in the 3rd layer
        else {
            Gzx_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(ci*(x-xp)*g9t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzx_t;
        }
    }
    //When the source point is in the 2nd layer
    else if( zp < 0.0 && zp > -thickness){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gzx_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(ci*(x-xp)*g9t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzx_t;
        }
        //When the field point is in the 2nd layer
        else if (z < 0.0 && z > -thickness) {
            Gzx_u = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(x-xp)*g8r_u(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gzx_d =  (1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(x-xp)*g8r_d(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzx_u+Gzx_d;
        }
        //When the field point is in the 3rd layer
        else {
            Gzx_t = (1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(ci*(x-xp)*g9t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzx_t;
        }
    }
    //When the source point is in the 3rd layer (use symmetry)
    else{
        zp = fabs(zp)-thickness;
        z = -z-thickness;
        epsilon_symmetry[0] = epsilon[2];
        epsilon_symmetry[1] = epsilon[1];
        epsilon_symmetry[2] = epsilon[0];

        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gzx_r  = (1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(ci*(x-xp)*g8r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gzx_r;
        }
        //When the field point is in the 2nd layer
        else if (z < 0.0 && z > -thickness){
            Gzx_r = (1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(x-xp)*g8r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gzx_t = (1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(x-xp)*g9t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));

            Gtemp = Gzx_r+Gzx_t;
        }
        //When the field point is in the 3rd layer
        else {
            Gzx_t = (1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(ci*(x-xp)*g9t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gzx_t;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gzy_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gzy_r, Gzy_p, Gzy_t, Gzy_u, Gzy_d;
    complex k, k1, k2, k3, ci;
    double omega, rho;
    complex epsilon_symmetry[3];

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp-0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gzy_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(ci*(y-yp)*g8r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzy_r;
        }
        //When the field point is in the 2nd layer
        else if (z < 0.0 && z > -thickness) {
            Gzy_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(y-yp)*g8r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gzy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(y-yp)*g9t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzy_r+Gzy_t;
        }
        //When the field point is in the 3rd layer
        else {
            Gzy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(ci*(y-yp)*g9t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzy_t;
        }
    }
    //When the source point is in the 2nd layer
    else if ( zp < 0.0 && zp > -thickness){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gzy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(ci*(y-yp)*g9t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzy_t;
        }
        //When the field point is in the 2nd layer
        else if (z < 0.0 && z > -thickness) {
            Gzy_u = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(y-yp)*g8r_u(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gzy_d =  (1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(y-yp)*g8r_d(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzy_u+Gzy_d;
        }
        //When the field point is in the 3rd layer
        else {
            Gzy_t = (1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(ci*(y-yp)*g9t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzy_t;
        }
    }
    //When the source point is in the 3rd layer (use symmetry)
    else{
        zp = fabs(zp)-thickness;
        z = -z-thickness;
        epsilon_symmetry[0] = epsilon[2];
        epsilon_symmetry[1] = epsilon[1];
        epsilon_symmetry[2] = epsilon[0];

        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gzy_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(ci*(y-yp)*g8r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gzy_r;
        }
        //When the field point is in the 2nd layer
        else if (z < 0.0 && z > -thickness) {
            Gzy_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(y-yp)*g8r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gzy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*(ci*(y-yp)*g9t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gzy_r+Gzy_t;
        }
        //When the field point is in the 3rd layer
        else {
            Gzy_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(ci*(y-yp)*g9t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gzy_t;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}


complex Gzz_three_regular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gzz_r, Gzz_p, Gzz_t, Gzz_u, Gzz_d;
    complex k, k1, k2, k3, ci;
    double omega, rho;
    complex epsilon_symmetry[3];

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));
    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp-0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z -0.0) < 10e-13 ){
            Gzz_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(g7r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzz_r;
        }
        //When the field point is in the 2nd layer
        else if(z < 0  &&  z > -thickness){
            Gzz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*g7t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda);
            Gzz_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*g7r(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda);
            Gtemp =  Gzz_r+Gzz_t;
        }
        //When the field point is in the 3rd layer
        else{
            Gzz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*g7t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda);
            Gtemp =  Gzz_t;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0.0 && zp > -thickness){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z - 0.0)< 10e-13 ){
            Gzz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*(g7t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda));
            Gtemp = Gzz_t;
        }
        //When the field point is in the 2nd layer
        else if(z < 0  &&  z > -thickness){
            Gzz_u = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*g7r_d(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda);
            Gzz_d = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*g7r_u(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda);

            Gtemp =  Gzz_u+Gzz_d;
        }
        //When the field point is in the 3rd layer
        else{
            Gzz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*g7t(k, rho, z, zp, thickness, epsilon, k1, k2, k3, lambda);
            Gtemp =  Gzz_t;
        }
    }
    //When the source point is in the 3rd layer (use symmetry)
    else{
        zp = fabs(zp)-thickness;
        z = -z-thickness;
        epsilon_symmetry[0] = epsilon[2];
        epsilon_symmetry[1] = epsilon[1];
        epsilon_symmetry[2] = epsilon[0];

        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z - 0.0) < 10e-13 ){
            Gzz_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[2]))*(g7r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda));
            Gtemp = Gzz_r;
        }
        //When the field point is in the 2nd layer
        else if(z < 0  &&  z > -thickness){
            Gzz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*g7t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda);
            Gzz_r = -(1/(8*PI*PI*omega*epsilon_0*epsilon[1]))*g7r(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda);
            Gtemp =  Gzz_r+Gzz_t;
        }
        //When the field point is in the 3rd layer
        else{
            Gzz_t = -(1/(8*PI*PI*omega*epsilon_0*epsilon[0]))*g7t(k, rho, z, zp, thickness, epsilon_symmetry, k3, k2, k1, lambda);
            Gtemp =  Gzz_t;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gxx_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gxx_r, Gxx_p, Gxx_t;
    complex k, k1, k2, k3, ci;
    double omega, rho, r0;

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));
    r0 = sqrt(rho*rho+(z-zp)*(z-zp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp-0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gxx_p = (1/(ci*4*PI*omega*epsilon_0*epsilon[0]))*(exp(ci*k1*r0)/pow(r0,5))*((x-xp)*(x-xp)*f1(k1*r0) - r0*r0*f2(k1*r0));
            Gtemp = Gxx_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0 &&  zp > -thickness){
        //When the field point is in the 2nd layer
        if(z < 0 &&  z > -thickness){
            Gxx_p = (1/(ci*4*PI*omega*epsilon_0*epsilon[1]))*(exp(ci*k2*r0)/pow(r0,5))*((x-xp)*(x-xp)*f1(k2*r0) - r0*r0*f2(k2*r0));
            Gtemp = Gxx_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 3rd layer
    else{
        //When the field point is in the 3rd layer
        if(z < -thickness ){
            Gxx_p = (1/(ci*4*PI*omega*epsilon_0*epsilon[2]))*(exp(ci*k3*r0)/pow(r0,5))*((x-xp)*(x-xp)*f1(k3*r0) - r0*r0*f2(k3*r0));
            Gtemp = Gxx_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gxy_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gxy_r, Gxy_p, Gxy_t;
    complex k, k1, k2, k3, ci;
    double omega, rho, r0;

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));
    r0 = sqrt(rho*rho+(z-zp)*(z-zp));

    //When the source point is in the 1st layer
    if(zp > 0  || fabs(zp-0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z-0.0) < 10e-13){
            Gxy_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[0]))*(exp(ci*k1*r0)/pow(r0,5))*(x-xp)*(y-yp)*f1(k1*r0);
            Gtemp = Gxy_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0  &&  zp > -thickness){
        //When the field point is in the 2nd layer
        if(z < 0 &&  z > -thickness){
            Gxy_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[1]))*(exp(ci*k2*r0)/pow(r0,5))*(x-xp)*(y-yp)*f1(k2*r0);
            Gtemp = Gxy_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 3rd layer
    else{
        //When the field point is in the 3rd layer
        if(z < -thickness){
            Gxy_p = -(1/(4*ci*PI*omega*epsilon_0*epsilon[2]))*(exp(ci*k3*r0)/pow(r0,5))*(x-xp)*(y-yp)*f1(k3*r0);
            Gtemp = Gxy_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gxz_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gxz_r, Gxz_p, Gxz_t,Gxz_r2;
    complex k, k1, k2, k3, ci;
    double omega, rho, r0;

    ci = complex(0.0,1.0);

    omega = 2.0*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));
    r0  = sqrt(rho*rho+(z-zp)*(z-zp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp-0.0)< 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gxz_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[0]))*(exp(ci*k1*r0)/pow(r0,5))*(x-xp)*(z-zp)*f1(k1*r0);
            Gtemp = Gxz_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0 &&  zp > -thickness){
        //When the field point is in the 2nd layer
        if(z < 0  &&  z > -thickness){
            Gxz_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[1]))*(exp(ci*k2*r0)/pow(r0,5))*(x-xp)*(z-zp)*f1(k2*r0);
            Gtemp = Gxz_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 3rd layer
    else{
        //When the field point is in the 3rd layer
        if (z < -thickness){
            Gxz_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[2]))*(exp(ci*k3*r0)/pow(r0,5))*(x-xp)*(z-zp)*f1(k3*r0);
            Gtemp = Gxz_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gyx_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    return Gxy_three_singular(lambda, epsilon_0, epsilon, x, y, z, xp, yp, zp, thickness);
}





complex Gyy_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gyy_r, Gyy_p, Gyy_t;
    complex k, k1, k2, k3, ci;
    double omega, rho, r0;

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));
    r0 = sqrt(rho*rho+(z-zp)*(z-zp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp-0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z-0.0) < 10e-13){
            Gyy_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[0]))*(exp(ci*k1*r0)/pow(r0,5))*((y-yp)*(y-yp)*f1(k1*r0) - r0*r0*f2(k1*r0));
            Gtemp = Gyy_p;
        }
        else {
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0 &&  zp > -thickness){
        //When the field point is in the 2nd layer
        if( z < 0 && z > -thickness){
            Gyy_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[1]))*(exp(ci*k2*r0)/pow(r0,5))*((y-yp)*(y-yp)*f1(k2*r0) - r0*r0*f2(k2*r0));
            Gtemp = Gyy_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 3rd layer
    else{
        //When the field point is in the 3rd layer
        if(z < -thickness){
            Gyy_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[2]))*(exp(ci*k3*r0)/pow(r0,5))*((y-yp)*(y-yp)*f1(k3*r0) - r0*r0*f2(k3*r0));
            Gtemp = Gyy_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gzx_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gzx_r, Gzx_p, Gzx_t;
    complex k, k1, k2, k3, ci;
    double omega, rho, r0;

    ci = complex(0.0,1.0);

    omega = 2.0*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));
    r0  = sqrt(rho*rho+(z-zp)*(z-zp));


    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp-0.0)< 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gzx_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[0]))*(exp(ci*k1*r0)/pow(r0,5))*(x-xp)*(z-zp)*f1(k1*r0);
            Gtemp = Gzx_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0 &&  zp > -thickness){
        //When the field point is in the 2nd layer
        if(z < 0 &&  z > -thickness){
            Gzx_p =  (1/(4*ci*PI*omega*epsilon_0*epsilon[1]))*(exp(ci*k2*r0)/pow(r0,5))*(x-xp)*(z-zp)*f1(k2*r0);
            Gtemp = Gzx_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 3rd layer
    else{
        //When the field point is in the 3rd layer
        if(z < -thickness){
            Gzx_p =  (1/(4*ci*PI*omega*epsilon_0*epsilon[2]))*(exp(ci*k3*r0)/pow(r0,5))*(x-xp)*(z-zp)*f1(k3*r0);
            Gtemp = Gzx_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gyz_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gyz_r, Gyz_p, Gyz_t;
    complex k, k1, k2, k3, ci;
    double omega, rho, r0;

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));
    r0 = sqrt(rho*rho+(z-zp)*(z-zp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp-0.0)< 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gyz_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[0]))*(exp(ci*k1*r0)/pow(r0,5))*(y-yp)*(z-zp)*f1(k1*r0);
            Gtemp = Gyz_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0 &&  zp > -thickness){
        //When the field point is in the 2nd layer
        if(z < 0 &&  z > -thickness){
            Gyz_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[1]))*(exp(ci*k2*r0)/pow(r0,5))*(y-yp)*(z-zp)*f1(k2*r0);
            Gtemp = Gyz_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 3rd layer
    else{
        //When the field point is in the 3rd layer
        if(z < -thickness){
            Gyz_p = -(1/(4*ci*PI*omega*epsilon_0*epsilon[2]))*(exp(ci*k3*r0)/pow(r0,5))*(y-yp)*(z-zp)*f1(k3*r0);
            Gtemp = Gyz_p;
        }
        else {
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gzy_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gzy_r, Gzy_p, Gzy_t;
    complex k, k1, k2, k3, ci;
    double omega, rho, r0;

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));
    r0  = sqrt(rho*rho+(z-zp)*(z-zp));

    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp-0.0)< 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0  || fabs(z-0.0) < 10e-13){
            Gzy_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[0]))*(exp(ci*k1*r0)/pow(r0,5))*(y-yp)*(z-zp)*f1(k1*r0);
            Gtemp = Gzy_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 2nd layer
    else if(zp < 0 &&  zp > -thickness ){
        //When the field point is in the 2nd layer
        if(z < 0 &&  z > -thickness){
            Gzy_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[1]))*(exp(ci*k2*r0)/pow(r0,5))*(y-yp)*(z-zp)*f1(k2*r0);
            Gtemp = Gzy_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 3rd layer
    else{
        //When the field point is in the 3rd layer
        if(z < -thickness ){
            Gzy_p = -(1/(4*ci*PI*omega*epsilon_0*epsilon[2]))*(exp(ci*k3*r0)/pow(r0,5))*(y-yp)*(z-zp)*f1(k3*r0);
            Gtemp = Gzy_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





complex Gzz_three_singular(double lambda, double epsilon_0, complex *epsilon, double x, double y, double z, double xp, double yp, double zp, double thickness)
{
    complex Gtemp, G;
    complex Gzz_r, Gzz_p, Gzz_t;
    complex k, k1, k2, k3, ci;
    double omega, rho, r0;

    ci = complex(0.0,1.0);

    omega = 2*PI/lambda;
    k = omega*sqrt(epsilon_0);

    k1 = k*sqrt(epsilon[0]);
    k2 = k*sqrt(epsilon[1]);
    k3 = k*sqrt(epsilon[2]);

    rho = sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp));
    r0 = sqrt(rho*rho+(z-zp)*(z-zp));


    //When the source point is in the 1st layer
    if(zp > 0 || fabs(zp - 0.0) < 10e-13){
        //When the field point is in the 1st layer
        if(z > 0.0 || fabs(z-0.0) < 10e-13){
            Gzz_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[0]))*(exp(ci*k1*r0)/pow(r0,5))*((z-zp)*(z-zp)*f1(k1*r0)-r0*r0*f2(k1*r0));
            Gtemp = Gzz_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 2nd layer
    else if( zp < 0  &&  zp > -thickness ){
        //When the field point is in the 2nd layer
        if(z < 0  &&  z > -thickness ){
            Gzz_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[1]))*(exp(ci*k2*r0)/pow(r0,5))*((z-zp)*(z-zp)*f1(k2*r0)-r0*r0*f2(k2*r0));
            Gtemp = Gzz_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    //When the source point is in the 3rd layer
    else{
        //When the field point is in the 3rd layer
        if(z < -thickness){
            Gzz_p = (1/(4*ci*PI*omega*epsilon_0*epsilon[2]))*(exp(ci*k3*r0)/pow(r0,5))*((z-zp)*(z-zp)*f1(k3*r0)-r0*r0*f2(k3*r0));
            Gtemp = Gzz_p;
        }
        else{
            Gtemp.re = 0.0;
            Gtemp.im = 0.0;
        }
    }
    G.re = Gtemp.re;
    G.im = Gtemp.im;

    return G;
}





//For free-space Dyadic Green's function
complex f1(complex x)
{
    complex ci(0.0, 1.0);
    return x*x+ci*3*x-3;
}
//For free-space Dyadic Green's function
complex f2(complex x)
{
    complex ci(0.0, 1.0);
    return x*x+ci*x-1;
}



//Sommerfeld integrals
//Only used for when the field point is in 2nd layer
complex g7r_u(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex spectral_g7_r;
    complex R_12_TM, R_23_TM, R_21_TM, general_R_12_TM, A_2_TM, U_TM;

    sum2 = complex(0.0, 0.0);
    ci  = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    for(i = 0 ; i < n ; i++){
        ks2 = ks[i]*ks[i];
        k1z = sqrt(k1*k1-ks2);
        k2z = sqrt(k2*k2-ks2);
        k3z = sqrt(k3*k3-ks2);

        R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
        R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
        R_21_TM = -R_12_TM;


        U_TM = R_23_TM*exp(2*ci*k2z*thickness)*(exp(ci*k2z*zp)+R_21_TM*exp(-ci*k2z*zp))/(1-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

        spectral_g7_r = ks2*U_TM/k2z;
        sum[i] = ks[i]*weight[i]*spectral_g7_r*J0(ks[i]*rho)*exp(ci*k2z*z);
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}

//Only used for when the field point is in 2nd layer
complex g7r_d(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex spectral_g7_r;
    complex R_12_TM, R_23_TM, general_R_12_TM, R_21_TM, A_2_TM, D_TM;

    sum2 = complex(0.0, 0.0);
    ci  = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    for(i = 0 ; i < n ; i++){
        ks2 = ks[i]*ks[i];
        k1z = sqrt(k1*k1-ks2);
        k2z = sqrt(k2*k2-ks2);
        k3z = sqrt(k3*k3-ks2);

        R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
        R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
        R_21_TM = -R_12_TM;

        D_TM = R_21_TM*(exp(-ci*k2z*zp)+R_23_TM*exp(ci*k2z*(zp+2*thickness)))/(1-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

        spectral_g7_r = ks2*D_TM/k2z;
        sum[i] = ks[i]*weight[i]*spectral_g7_r*J0(ks[i]*rho)*exp(-ci*k2z*z);
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}


complex g7r(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex spectral_g7_r;
    complex R_12_TM, R_23_TM, A_2_TM, general_R_12_TM;

    sum2 = complex(0.0, 0.0);
    ci  = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    //when the field point z is in the 1st layer
    if( z > 0 || fabs(z - 0.0)<10e-13 ){
        for(i = 0 ; i < n ; i++){
            ks2 = ks[i]*ks[i];
            k1z = sqrt(k1*k1-ks2);
            k2z = sqrt(k2*k2-ks2);
            k3z = sqrt(k3*k3-ks2);

            R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
            R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

            general_R_12_TM = (R_12_TM+R_23_TM*exp(2*ci*k2z*thickness))/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));

            spectral_g7_r = ks2*general_R_12_TM/k1z;
            sum[i] = ks[i]*weight[i]*spectral_g7_r*J0(ks[i]*rho)*exp(ci*k1z*(z+zp));
        }
    }
    //When z is in the 2nd layer
    else{
        for(i = 0 ; i < n ; i++){
            ks2 = ks[i]*ks[i];
            k1z = sqrt(k1*k1-ks2);
            k2z = sqrt(k2*k2-ks2);
            k3z = sqrt(k3*k3-ks2);

            R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
            R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

            //R_21 = -R_12 and T_12 = 1+R_12 are used
            A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));

            spectral_g7_r = ks2*A_2_TM*R_23_TM/k1z;
            sum[i] = ks[i]*weight[i]*spectral_g7_r*J0(ks[i]*rho)*exp(ci*k2z*z+ci*k1z*zp+2*ci*k2z*thickness);
        }
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}


complex g7t(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_23_TM, R_21_TM, A_1_TM, A_2_TM, A_3_TM;
    complex spectral_g7_t;

    sum2 = complex(0.0, 0.0);
    ci   = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    //When zp is in the 1st layer
    if(zp > 0.0 || fabs(zp-0.0)<10e-13){
        //When z is in the 2nd layer
        if( z < 0.0 && z > -thickness){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

                A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));

                spectral_g7_t =ks2*A_2_TM/k1z;
                sum[i] = ks[i]*weight[i]*spectral_g7_t*J0(ks[i]*rho)*exp(ci*(-k2z*z+k1z*zp));
            }
        }
        //When z is in the 3rd layer
        else if (z < -thickness || fabs(z+thickness)<10e-13){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;

                A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));

                A_3_TM = A_2_TM*(1.0+R_23_TM)*exp(ci*(k2z-k3z)*thickness);
                spectral_g7_t =ks2*A_3_TM/k1z;
                sum[i] = ks[i]*weight[i]*spectral_g7_t*J0(ks[i]*rho)*exp(ci*(-k3z*z+k1z*zp));
            }
        }
        // (This case never happens)
        else{
        }
    }
    //When zp is in the 2nd layer
    else if( zp < 0.0 && zp > -thickness){
        //when z is in the 1st layer
        if(z > 0 || fabs(z-0.0) < 10e-13 ){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;


                A_1_TM = (1+R_21_TM)*(exp(-ci*k2z*zp)+R_23_TM*exp(ci*k2z*(zp+2*thickness)))/(1.0-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

                spectral_g7_t =ks2*A_1_TM/k2z;
                sum[i] = ks[i]*weight[i]*spectral_g7_t*J0(ks[i]*rho)*exp(ci*k1z*z);
            }
        }
        else if (z < 0.0 && z > -thickness){
            //This case never happen
        }
        //when z is in the 3rd layer
        else{
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;

                A_3_TM = (1+R_23_TM)*exp(ci*(k2z-k3z)*thickness)*(exp(ci*k2z*zp)+R_21_TM*exp(-ci*k2z*zp))/(1.0-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

                spectral_g7_t =ks2*A_3_TM/k2z;
                sum[i] = ks[i]*weight[i]*spectral_g7_t*J0(ks[i]*rho)*exp(-ci*k3z*z);
            }
        }
    }

    //When zp is in the 3rd layer (Use symmetry)
    else{

    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}





complex g9t(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex pole1, pole2, pole3, pole4;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_23_TM, R_21_TM, T_12_TM, T_21_TM, T_23_TM, A_1_TM, A_2_TM, A_3_TM;
    complex spectral_g9_t;

    sum2 = complex(0.0,0.0);
    ci = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    //When zp is in the first layer
    if(zp > 0 || fabs(zp - 0.0)<10e-13){
        if( z < 0.0 && z > -thickness){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                T_12_TM = 2.0*epsilon[1]*k1z/(epsilon[1]*k1z+epsilon[0]*k2z);

                A_2_TM = T_12_TM/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));

                spectral_g9_t = A_2_TM;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*weight[i]*spectral_g9_t*(ks[i]*0.5)*exp(ci*(-k2z*z+k1z*zp));
                }
                else{
                    sum[i] = ks2*weight[i]*spectral_g9_t*(J1(ks[i]*rho)/rho)*exp(ci*(-k2z*z+k1z*zp));
                }
            }
        }
        else {
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

                A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
                A_3_TM = A_2_TM*(1+R_23_TM)*exp(ci*(k2z-k3z)*thickness);
                spectral_g9_t = A_3_TM;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*weight[i]*spectral_g9_t*(ks[i]*0.5)*exp(ci*(-k3z*z+k1z*zp));
                }
                else{
                    sum[i] = ks2*weight[i]*spectral_g9_t*(J1(ks[i]*rho)/rho)*exp(ci*(-k3z*z+k1z*zp));
                }
            }
        }
    }
    //When zp is in the second layer
    else if (zp < 0.0 && zp > -thickness){
        if( z > 0.0 || fabs(z -0.0)<10e-13){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;
                T_21_TM = (2*epsilon[0]*k2z)/(epsilon[0]*k2z+epsilon[1]*k1z);
                //                T_21_TM = (k2z/k1z)*((1+R_12_TM)/(1-R_21_TM*R_23_TM*exp(2*ci*k2z*thickness)));

                //                A_1_TM = (T_21_TM)*(exp(-ci*k2z*zp)+R_23_TM*exp(ci*k2z*zp+2*ci*k2z*thickness))/(1.0+R_23_TM*R_12_TM*exp(2*ci*k2z*thickness));
                A_1_TM = (T_21_TM)*(exp(-ci*k2z*zp)-R_23_TM*exp(ci*k2z*zp+2*ci*k2z*thickness))/(1.0+R_23_TM*R_12_TM*exp(2*ci*k2z*thickness));

                spectral_g9_t = -A_1_TM;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*weight[i]*spectral_g9_t*(ks[i]*0.5)*exp(ci*k1z*z);
                }
                else{
                    sum[i] = ks2*weight[i]*spectral_g9_t*(J1(ks[i]*rho)/rho)*exp(ci*k1z*z);
                }
            }
        }
        else if (z < 0.0 && z > -thickness){

        }
        else {
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;
                T_23_TM = (2*epsilon[2]*k2z)/(epsilon[2]*k2z+epsilon[1]*k3z);

//                A_3_TM = T_23_TM*exp(ci*(k2z-k3z)*thickness)*(exp(ci*k2z*zp)+R_21_TM*exp(-ci*k2z*zp))/(1.0-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

                A_3_TM = T_23_TM*exp(ci*(k2z-k3z)*thickness)*(exp(ci*k2z*zp)-R_21_TM*exp(-ci*k2z*zp))/(1.0-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

                spectral_g9_t = -A_3_TM;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*weight[i]*spectral_g9_t*(ks[i]*0.5)*exp(-ci*k3z*z);
                }
                else{
                    sum[i] = ks2*weight[i]*spectral_g9_t*(J1(ks[i]*rho)/rho)*exp(-ci*k3z*z);
                }
            }
        }

    }
    else{

    }


    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}


complex g6r_u(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_21_TM, R_23_TM;
    complex R_12_TE, R_21_TE, R_23_TE;
    complex general_R_12_TM, general_R_12_TE;
    complex U_TM, U_TE;
    complex spectral_g6_r;

    sum2 = complex(0.0,0.0);
    ci = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);


    //when the field point z is in the 1st layer
    for(i = 0 ; i < n ; i++){
        ks2 = ks[i]*ks[i];
        k1z = sqrt(k1*k1-ks2);
        k2z = sqrt(k2*k2-ks2);
        k3z = sqrt(k3*k3-ks2);

        R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
        R_21_TM = -R_12_TM;
        R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
        R_12_TE = (k1z-k2z)/(k1z+k2z);
        R_21_TE = -R_12_TE;
        R_23_TE = (k2z-k3z)/(k2z+k3z);

        U_TM =   R_23_TM*exp(2*ci*k2z*thickness)*(exp(ci*k2z*zp)-R_21_TM*exp(-ci*k2z*zp))/(1-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));
        U_TE =   R_23_TE*exp(2*ci*k2z*thickness)*(exp(ci*k2z*zp)+R_21_TE*exp(-ci*k2z*zp))/(1-R_23_TE*R_21_TE*exp(2*ci*k2z*thickness));

        spectral_g6_r = (k2z/ks2)*U_TM+(k2*k2/(ks2*k2z))*U_TE;

        if (fabs(rho-0.0)<10e-13){
            sum[i] = ks2*ks[i]*weight[i]*spectral_g6_r*(ks2/8.0)*exp(ci*k2z*z);
        }
        else{
            sum[i] = ks2*ks[i]*weight[i]*spectral_g6_r*J2(ks[i]*rho)/(rho*rho)*exp(ci*k2z*z);
        }
    }


    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}


complex g6r_d(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_21_TM, R_23_TM;
    complex R_12_TE, R_21_TE, R_23_TE;
    complex general_R_12_TM, general_R_12_TE;
    complex D_TM, D_TE;
    complex spectral_g6_r;

    sum2 = complex(0.0,0.0);
    ci = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);


    for(i = 0 ; i < n ; i++){
        ks2 = ks[i]*ks[i];
        k1z = sqrt(k1*k1-ks2);
        k2z = sqrt(k2*k2-ks2);
        k3z = sqrt(k3*k3-ks2);

        R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
        R_21_TM = -R_12_TM;
        R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
        R_12_TE = (k1z-k2z)/(k1z+k2z);
        R_21_TE = -R_12_TE;
        R_23_TE = (k2z-k3z)/(k2z+k3z);

        D_TM   =  R_21_TM*(exp(-ci*k2z*zp)-R_23_TM*exp(ci*k2z*(zp+2*thickness)))/(1-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));
        D_TE   =  R_21_TE*(exp(-ci*k2z*zp)+R_23_TE*exp(ci*k2z*(zp+2*thickness)))/(1-R_23_TE*R_21_TE*exp(2*ci*k2z*thickness));


        spectral_g6_r = (k2z/ks2)*D_TM+(k2*k2/(ks2*k2z))*D_TE;

        if (fabs(rho-0.0)<10e-13){
            sum[i] = ks2*ks[i]*weight[i]*spectral_g6_r*(ks2/8.0)*exp(-ci*k2z*z);
        }
        else{
            sum[i] = ks2*ks[i]*weight[i]*spectral_g6_r*J2(ks[i]*rho)/(rho*rho)*exp(-ci*k2z*z);
        }
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}





complex g6r(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_23_TM;
    complex R_12_TE, R_23_TE;
    complex general_R_12_TM, general_R_12_TE;
    complex A_2_TM, A_2_TE;
    complex spectral_g6_r;

    sum2 = complex(0.0,0.0);
    ci = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);


//    n1 = 200;
//    n2 = 20;
//    n3 = 20;
//    n4 = 200;
//    n5 = 20;
//    n6 = 20;
//    n7 = 200;
//    n8 = 20;
//    n9 = 20;
//    n10 = 200;
//
//    n = n1+n2+n3+n4+n5+n6+n7+n8+n9+n10;
//
//    ks     = new complex[n];
//    weight = new complex[n];
//    sum    = new complex[n];
//
//    sommerfeld_quadrature_gauss_adaptive_three_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10);


    //when the field point z is in the 1st layer
    if( z > 0 || fabs(z-0.0)<10e-13 ){
        for(i = 0 ; i < n ; i++){
            ks2 = ks[i]*ks[i];
            k1z = sqrt(k1*k1-ks2);
            k2z = sqrt(k2*k2-ks2);
            k3z = sqrt(k3*k3-ks2);

            R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
            R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
            R_12_TE = (k1z-k2z)/(k1z+k2z);
            R_23_TE = (k2z-k3z)/(k2z+k3z);


            general_R_12_TM = (R_12_TM+R_23_TM*exp(2*ci*k2z*thickness))/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
            general_R_12_TE = (R_12_TE+R_23_TE*exp(2*ci*k2z*thickness))/(1.0+R_12_TE*R_23_TE*exp(2*ci*k2z*thickness));

            spectral_g6_r = (k1z/ks2)*general_R_12_TM+(k1*k1/(ks2*k1z))*general_R_12_TE;

            if (fabs(rho-0.0)<10e-13){
                sum[i] = ks2*ks[i]*weight[i]*spectral_g6_r*(ks2/8.0)*exp(ci*k1z*(z+zp));
            }
            else{
                sum[i] = ks2*ks[i]*weight[i]*spectral_g6_r*J2(ks[i]*rho)/(rho*rho)*exp(ci*k1z*(z+zp));
            }
        }
    }
    //when the field point z is in the 2nd layer
    else{
        for(i = 0 ; i < n ; i++){
            ks2 = ks[i]*ks[i];
            k1z = sqrt(k1*k1-ks2);
            k2z = sqrt(k2*k2-ks2);
            k3z = sqrt(k3*k3-ks2);

            R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
            R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
            R_12_TE = (k1z-k2z)/(k1z+k2z);
            R_23_TE = (k2z-k3z)/(k2z+k3z);

            A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
            A_2_TE = (1.0+R_12_TE)/(1.0+R_12_TE*R_23_TE*exp(2*ci*k2z*thickness));

            spectral_g6_r = (k2z/ks2)*A_2_TM*R_23_TM+(k2*k2/(ks2*k1z))*A_2_TE*R_23_TE;

            if (fabs(rho-0.0)<10e-13){
                sum[i] = ks2*ks[i]*weight[i]*spectral_g6_r*(ks2/8.0)*exp(ci*k2z*z+ci*k1z*zp+2*ci*k2z*thickness);
            }
            else{
                sum[i] = ks2*ks[i]*weight[i]*spectral_g6_r*J2(ks[i]*rho)/(rho*rho)*exp(ci*k2z*z+ci*k1z*zp+2*ci*k2z*thickness);
            }
        }
    }


    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}






complex g6t(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_21_TM, R_12_TE, T_21_TM, T_23_TM;
    complex R_23_TM, R_23_TE, T_21_TE, T_23_TE, R_21_TE;
    complex A_1_TM, A_1_TE;
    complex A_2_TM, A_2_TE;
    complex A_3_TM, A_3_TE;
    complex spectral_g6_t;

    sum2 = complex(0.0, 0.0);
    ci   = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    //When zp is in the 1st layer
    if(zp > 0 || fabs(zp - 0.0)<10e-13){
        //When z is in the 2nd layer
        if( z < 0.0 && z > -thickness){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_12_TE = (k1z-k2z)/(k1z+k2z);
                R_23_TE = (k2z-k3z)/(k2z+k3z);

                A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
                A_2_TE = (1.0+R_12_TE)/(1.0+R_12_TE*R_23_TE*exp(2*ci*k2z*thickness));

                spectral_g6_t = (k2z/ks2)*A_2_TM-(k2*k2/(k1z*ks2))*A_2_TE;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*ks[i]*weight[i]*spectral_g6_t*(ks2/8.0)*exp(ci*(-k2z*z+k1z*zp));
                }
                else{
                    sum[i] = ks2*ks[i]*weight[i]*spectral_g6_t*(J2(ks[i]*rho)/(rho*rho))*exp(ci*(-k2z*z+k1z*zp));
                }
            }
        }
        //When z is in the 3rd layer
        else{
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_12_TE = (k1z-k2z)/(k1z+k2z);
                R_23_TE = (k2z-k3z)/(k2z+k3z);

                A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
                A_2_TE = (1.0+R_12_TE)/(1.0+R_12_TE*R_23_TE*exp(2*ci*k2z*thickness));

                A_3_TM = A_2_TM*(1.0+R_23_TM)*exp(ci*(k2z-k3z)*thickness);
                A_3_TE = A_2_TE*(1.0+R_23_TE)*exp(ci*(k2z-k3z)*thickness);

                spectral_g6_t = (k3z/ks2)*A_3_TM-(k3*k3/(k1z*ks2))*A_3_TE;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*ks[i]*weight[i]*spectral_g6_t*(ks2/8.0)*exp(ci*(-k3z*z+k1z*zp));
                }
                else{
                    sum[i] = ks2*ks[i]*weight[i]*spectral_g6_t*(J2(ks[i]*rho)/(rho*rho))*exp(ci*(-k3z*z+k1z*zp));
                }
            }
        }
    }
    //When zp is in the 2nd layer
    else if(zp < 0.0 && zp > -thickness){
        if(z > 0.0 || fabs(z-0.0)<10e-13){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;
                T_21_TM = (2*epsilon[0]*k2z)/(epsilon[0]*k2z+epsilon[1]*k1z);

                R_12_TE = (k1z-k2z)/(k1z+k2z);
                R_21_TE = -R_12_TE;
                R_23_TE = (k2z-k3z)/(k2z+k3z);
                T_21_TE = (2*k2z)/(k2z+k1z);


                A_1_TM =  (T_21_TM)*(exp(-ci*k2z*zp)-R_23_TM*exp(ci*k2z*zp+2*ci*k2z*thickness))/(1.0+R_23_TM*R_12_TM*exp(2*ci*k2z*thickness));
                A_1_TE =  (T_21_TE)*(exp(-ci*k2z*zp)+R_23_TE*exp(ci*k2z*zp+2*ci*k2z*thickness))/(1.0+R_23_TE*R_12_TE*exp(2*ci*k2z*thickness));

                spectral_g6_t = (k1z/ks2)*A_1_TM-(k1*k1/(k2z*ks2))*A_1_TE;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*ks[i]*weight[i]*spectral_g6_t*(ks2/8.0)*exp(ci*k1z*z);
                }
                else{
                    sum[i] = ks2*ks[i]*weight[i]*spectral_g6_t*(J2(ks[i]*rho)/(rho*rho))*exp(ci*k1z*z);
                }
            }
        }
        else if (z < 0.0 && z > -thickness){

        }
        else{
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;
                T_23_TM = (2*epsilon[2]*k2z)/(epsilon[2]*k2z+epsilon[1]*k3z);

                R_12_TE = (k1z-k2z)/(k1z+k2z);
                R_21_TE = -R_12_TE;
                R_23_TE = (k2z-k3z)/(k2z+k3z);
                T_23_TE = (2*k2z)/(k2z+k3z);

                A_3_TM =  T_23_TM*exp(ci*(k2z-k3z)*thickness)*(exp(ci*k2z*zp)-R_21_TM*exp(-ci*k2z*zp))/(1.0+R_23_TM*R_12_TM*exp(2*ci*k2z*thickness));
                A_3_TE =  T_23_TE*exp(ci*(k2z-k3z)*thickness)*(exp(ci*k2z*zp)+R_21_TE*exp(-ci*k2z*zp))/(1.0+R_23_TE*R_12_TE*exp(2*ci*k2z*thickness));

                spectral_g6_t = (k3z/ks2)*A_3_TM-(k3*k3/(k2z*ks2))*A_3_TE;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*ks[i]*weight[i]*spectral_g6_t*(ks2/8.0)*exp(-ci*k3z*z);
                }
                else{
                    sum[i] = ks2*ks[i]*weight[i]*spectral_g6_t*(J2(ks[i]*rho)/(rho*rho))*exp(-ci*k3z*z);
                }
            }
        }


    }
    //When zp is in the 3rd layer
    else{
    }


    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}

complex g5r_u(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex spectral_g5_r;
    complex general_R_12_TM, general_R_12_TE;
    complex R_12_TM, R_12_TE, R_23_TM, R_23_TE;
    complex R_21_TM, R_21_TE;
    complex A_2_TM, A_2_TE, U_TM, U_TE;

    sum2 = complex(0.0,0.0);
    ci = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    for(i = 0 ; i < n ; i++){
        ks2 = ks[i]*ks[i];
        k1z = sqrt(k1*k1-ks2);
        k2z = sqrt(k2*k2-ks2);
        k3z = sqrt(k3*k3-ks2);

        R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
        R_21_TM = -R_12_TM;
        R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
        R_12_TE = (k1z-k2z)/(k1z+k2z);
        R_21_TE = -R_12_TE;
        R_23_TE = (k2z-k3z)/(k2z+k3z);

        U_TM =   R_23_TM*exp(2*ci*k2z*thickness)*(exp(ci*k2z*zp)-R_21_TM*exp(-ci*k2z*zp))/(1-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));
        U_TE =   R_23_TE*exp(2*ci*k2z*thickness)*(exp(ci*k2z*zp)+R_21_TE*exp(-ci*k2z*zp))/(1-R_23_TE*R_21_TE*exp(2*ci*k2z*thickness));

        spectral_g5_r = k2z*U_TM-(k2*k2/k2z)*U_TE;

        sum[i] = ks[i]*weight[i]*spectral_g5_r*J0(ks[i]*rho)*exp(ci*k2z*z);
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}



complex g5r_d(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex spectral_g5_r;
    complex general_R_12_TM, general_R_12_TE;
    complex R_12_TM, R_12_TE, R_23_TM, R_23_TE;
    complex R_21_TM, R_21_TE;
    complex A_2_TM, A_2_TE, D_TM, D_TE;

    sum2 = complex(0.0,0.0);
    ci = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    for(i = 0 ; i < n ; i++){
        ks2 = ks[i]*ks[i];
        k1z = sqrt(k1*k1-ks2);
        k2z = sqrt(k2*k2-ks2);
        k3z = sqrt(k3*k3-ks2);

        R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
        R_21_TM = -R_12_TM;
        R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
        R_12_TE = (k1z-k2z)/(k1z+k2z);
        R_21_TE = -R_12_TE;
        R_23_TE = (k2z-k3z)/(k2z+k3z);

        D_TM   =  R_21_TM*(exp(-ci*k2z*zp)-R_23_TM*exp(ci*k2z*(zp+2*thickness)))/(1-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));
        D_TE   =  R_21_TE*(exp(-ci*k2z*zp)+R_23_TE*exp(ci*k2z*(zp+2*thickness)))/(1-R_23_TE*R_21_TE*exp(2*ci*k2z*thickness));

        spectral_g5_r = k2z*D_TM-(k2*k2/k2z)*D_TE;

        sum[i] = ks[i]*weight[i]*spectral_g5_r*J0(ks[i]*rho)*exp(-ci*k2z*z);
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}



complex g5r(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex spectral_g5_r;
    complex general_R_12_TM, general_R_12_TE;
    complex R_12_TM, R_12_TE, R_23_TM, R_23_TE;
    complex A_2_TM, A_2_TE;

    sum2 = complex(0.0,0.0);
    ci = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    //when the field point z is in the 1st layer
    if( z > 0 || fabs(z-0.0)<10e-13 ){
        for(i = 0 ; i < n ; i++){
            ks2 = ks[i]*ks[i];
            k1z = sqrt(k1*k1-ks2);
            k2z = sqrt(k2*k2-ks2);
            k3z = sqrt(k3*k3-ks2);

            R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
            R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
            R_12_TE = (k1z-k2z)/(k1z+k2z);
            R_23_TE = (k2z-k3z)/(k2z+k3z);

            general_R_12_TM = (R_12_TM+R_23_TM*exp(2*ci*k2z*thickness))/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
            general_R_12_TE = (R_12_TE+R_23_TE*exp(2*ci*k2z*thickness))/(1.0+R_12_TE*R_23_TE*exp(2*ci*k2z*thickness));

            spectral_g5_r = k1z*general_R_12_TM-(k1*k1/k1z)*general_R_12_TE;

            sum[i] = ks[i]*weight[i]*spectral_g5_r*J0(ks[i]*rho)*exp(ci*k1z*(z+zp));
        }
    }
    //When z is in the 2nd layer
    else {
        for(i = 0 ; i < n ; i++){
            ks2 = ks[i]*ks[i];
            k1z = sqrt(k1*k1-ks2);
            k2z = sqrt(k2*k2-ks2);
            k3z = sqrt(k3*k3-ks2);

            R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
            R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
            R_12_TE = (k1z-k2z)/(k1z+k2z);
            R_23_TE = (k2z-k3z)/(k2z+k3z);

            A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
            A_2_TE = (1.0+R_12_TE)/(1.0+R_12_TE*R_23_TE*exp(2*ci*k2z*thickness));

            spectral_g5_r = k2z*A_2_TM*R_23_TM-(k2*k2/k1z)*A_2_TE*R_23_TE;

            sum[i] = ks[i]*weight[i]*spectral_g5_r*J0(ks[i]*rho)*exp(ci*k2z*z+ci*k1z*zp+2*ci*k2z*thickness);
        }
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}



complex g5t(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_12_TE, R_23_TM, R_23_TE;
    complex A_1_TM, A_1_TE, A_2_TM, A_2_TE, A_3_TM, A_3_TE;
    complex R_21_TM, T_21_TM, R_21_TE, T_21_TE, T_23_TM, T_23_TE;
    complex spectral_g5_t;

    sum2 = complex(0.0, 0.0);
    ci   = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    //When zp is in the 1st layer
    if(zp > 0 || fabs(zp-0.0)<10e-13){
        //When z is in the 2nd layer
        if( z < 0.0 && z > -thickness){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_12_TE = (k1z-k2z)/(k1z+k2z);
                R_23_TE = (k2z-k3z)/(k2z+k3z);

                A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
                A_2_TE = (1.0+R_12_TE)/(1.0+R_12_TE*R_23_TE*exp(2*ci*k2z*thickness));

                spectral_g5_t = k2z*A_2_TM+(k2*k2/k1z)*A_2_TE;
                sum[i] = ks[i]*weight[i]*spectral_g5_t*J0(ks[i]*rho)*exp(ci*(-k2z*z+k1z*zp));
            }
        }
        //When z is in the 3rd layer
        else{
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_12_TE = (k1z-k2z)/(k1z+k2z);
                R_23_TE = (k2z-k3z)/(k2z+k3z);

                A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
                A_2_TE = (1.0+R_12_TE)/(1.0+R_12_TE*R_23_TE*exp(2*ci*k2z*thickness));


                A_3_TM = A_2_TM*(1.0+R_23_TM)*exp(ci*(k2z-k3z)*thickness);
                A_3_TE = A_2_TE*(1.0+R_23_TE)*exp(ci*(k2z-k3z)*thickness);

                spectral_g5_t = k3z*A_3_TM+(k3*k3/k1z)*A_3_TE;
                sum[i] = ks[i]*weight[i]*spectral_g5_t*J0(ks[i]*rho)*exp(ci*(-k3z*z+k1z*zp));
            }
        }
    }
    //When zp is in the 2nd layer
    else if(zp < 0.0 && zp > -thickness){
        if(z > 0.0 || fabs(z-0.0)<10e-13){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;
                T_21_TM = (2*epsilon[0]*k2z)/(epsilon[0]*k2z+epsilon[1]*k1z);

                R_12_TE = (k1z-k2z)/(k1z+k2z);
                R_21_TE = -R_12_TE;
                R_23_TE = (k2z-k3z)/(k2z+k3z);
                T_21_TE = (2*k2z)/(k2z+k1z);


                A_1_TM =  (T_21_TM)*(exp(-ci*k2z*zp)-R_23_TM*exp(ci*k2z*zp+2*ci*k2z*thickness))/(1.0+R_23_TM*R_12_TM*exp(2*ci*k2z*thickness));
                A_1_TE =  (T_21_TE)*(exp(-ci*k2z*zp)+R_23_TE*exp(ci*k2z*zp+2*ci*k2z*thickness))/(1.0+R_23_TE*R_12_TE*exp(2*ci*k2z*thickness));

                spectral_g5_t = k1z*A_1_TM+(k1*k1/k2z)*A_1_TE;

                sum[i] = ks[i]*weight[i]*spectral_g5_t*J0(ks[i]*rho)*exp(ci*k1z*z);
            }
        }
        else if (z < 0.0 && z > -thickness){

        }
        else{
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;
                T_23_TM = (2*epsilon[2]*k2z)/(epsilon[2]*k2z+epsilon[1]*k3z);

                R_12_TE = (k1z-k2z)/(k1z+k2z);
                R_21_TE = -R_12_TE;
                R_23_TE = (k2z-k3z)/(k2z+k3z);
                T_23_TE = (2*k2z)/(k2z+k3z);

                A_3_TM =  T_23_TM*exp(ci*(k2z-k3z)*thickness)*(exp(ci*k2z*zp)-R_21_TM*exp(-ci*k2z*zp))/(1.0+R_23_TM*R_12_TM*exp(2*ci*k2z*thickness));
                A_3_TE =  T_23_TE*exp(ci*(k2z-k3z)*thickness)*(exp(ci*k2z*zp)+R_21_TE*exp(-ci*k2z*zp))/(1.0+R_23_TE*R_12_TE*exp(2*ci*k2z*thickness));

                spectral_g5_t = k3z*A_3_TM+(k3*k3/k2z)*A_3_TE;
                sum[i] = ks[i]*weight[i]*spectral_g5_t*J0(ks[i]*rho)*exp(-ci*k3z*z);
            }
        }
    }
    //When zp is in the 3rd layer
    else{
    }


    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}



complex g8r_u_v(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_23_TM, R_21_TM, U_TM;
    complex general_R_12_TM;
    complex spectral_g8_r;

    sum2 = complex(0.0,0.0);
    ci   = complex(0.0,1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    for(i = 0 ; i < n ; i++){
        ks2 = ks[i]*ks[i];
        k1z = sqrt(k1*k1-ks2);
        k2z = sqrt(k2*k2-ks2);
        k3z = sqrt(k3*k3-ks2);

        R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
        R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

        general_R_12_TM = R_12_TM+((1+R_12_TM)*R_23_TM*(1-R_12_TM))/(1+R_12_TM*R_23_TM);
        R_21_TM = -R_12_TM;

        U_TM =   R_23_TM*exp(2*ci*k2z*thickness)*(exp(ci*k2z*zp)+R_21_TM*exp(-ci*k2z*zp))/(1-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

        spectral_g8_r = U_TM;
        if (fabs(rho-0.0)<10e-13){
            sum[i] = ks2*weight[i]*spectral_g8_r*(ks[i]*0.5)*exp(ci*k2z*z);
        }
        else{
            sum[i] = ks2*weight[i]*spectral_g8_r*(J1(ks[i]*rho)/rho)*exp(ci*k2z*z);
        }
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }


    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}


complex g8r_d_v(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_23_TM, R_21_TM, D_TM;
    complex general_R_12_TM;
    complex spectral_g8_r;

    sum2 = complex(0.0,0.0);
    ci   = complex(0.0,1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    for(i = 0 ; i < n ; i++){
        ks2 = ks[i]*ks[i];
        k1z = sqrt(k1*k1-ks2);
        k2z = sqrt(k2*k2-ks2);
        k3z = sqrt(k3*k3-ks2);

        R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
        R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
        R_21_TM = -R_12_TM;

        D_TM   = R_21_TM*(exp(-ci*k2z*zp)+R_23_TM*exp(ci*k2z*(zp+2*thickness)))/(1-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

        spectral_g8_r = D_TM;

        if (fabs(rho-0.0)<10e-13){
            sum[i] = ks2*weight[i]*spectral_g8_r*(ks[i]*0.5)*exp(-ci*k2z*z);
        }
        else{
            sum[i] = ks2*weight[i]*spectral_g8_r*(J1(ks[i]*rho)/rho)*exp(-ci*k2z*z);
        }
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}


complex g8r_u(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_23_TM, R_21_TM, U_TM;
    complex general_R_12_TM;
    complex spectral_g8_r;

    sum2 = complex(0.0,0.0);
    ci   = complex(0.0,1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    for(i = 0 ; i < n ; i++){
        ks2 = ks[i]*ks[i];
        k1z = sqrt(k1*k1-ks2);
        k2z = sqrt(k2*k2-ks2);
        k3z = sqrt(k3*k3-ks2);

        R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
        R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

        general_R_12_TM = R_12_TM+((1+R_12_TM)*R_23_TM*(1-R_12_TM))/(1+R_12_TM*R_23_TM);
        R_21_TM = -R_12_TM;

        U_TM =   R_23_TM*exp(2*ci*k2z*thickness)*(exp(ci*k2z*zp)-R_21_TM*exp(-ci*k2z*zp))/(1-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

        spectral_g8_r = U_TM;
        if (fabs(rho-0.0)<10e-13){
            sum[i] = ks2*weight[i]*spectral_g8_r*(ks[i]*0.5)*exp(ci*k2z*z);
        }
        else{
            sum[i] = ks2*weight[i]*spectral_g8_r*(J1(ks[i]*rho)/rho)*exp(ci*k2z*z);
        }
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }


    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}


complex g8r_d(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_23_TM, R_21_TM, D_TM;
    complex general_R_12_TM;
    complex spectral_g8_r;

    sum2 = complex(0.0,0.0);
    ci   = complex(0.0,1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    for(i = 0 ; i < n ; i++){
        ks2 = ks[i]*ks[i];
        k1z = sqrt(k1*k1-ks2);
        k2z = sqrt(k2*k2-ks2);
        k3z = sqrt(k3*k3-ks2);

        R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
        R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
        R_21_TM = -R_12_TM;

        D_TM   = R_21_TM*(exp(-ci*k2z*zp)-R_23_TM*exp(ci*k2z*(zp+2*thickness)))/(1-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

        spectral_g8_r = D_TM;

        if (fabs(rho-0.0)<10e-13){
            sum[i] = ks2*weight[i]*spectral_g8_r*(ks[i]*0.5)*exp(-ci*k2z*z);
        }
        else{
            sum[i] = ks2*weight[i]*spectral_g8_r*(J1(ks[i]*rho)/rho)*exp(-ci*k2z*z);
        }
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}



complex g8r(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_23_TM, A_2_TM;
    complex general_R_12_TM;
    complex spectral_g8_r;

    sum2 = complex(0.0,0.0);
    ci   = complex(0.0,1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    //when the field point z is in the 1st layer
    if( z > 0 || fabs(z-0.0)<10e-13 ){
        for(i = 0 ; i < n ; i++){
            ks2 = ks[i]*ks[i];
            k1z = sqrt(k1*k1-ks2);
            k2z = sqrt(k2*k2-ks2);
            k3z = sqrt(k3*k3-ks2);

            R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
            R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

//            general_R_12_TM = (R_12_TM-R_23_TM*exp(2*ci*k2z*thickness))/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
            general_R_12_TM = (R_12_TM+R_23_TM*exp(2*ci*k2z*thickness))/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));

            spectral_g8_r = general_R_12_TM;
            if (fabs(rho-0.0)<10e-13){
                sum[i] = ks2*weight[i]*spectral_g8_r*(ks[i]*0.5)*exp(ci*k1z*(z+zp));
            }
            else{
                sum[i] = ks2*weight[i]*spectral_g8_r*(J1(ks[i]*rho)/rho)*exp(ci*k1z*(z+zp));
            }
        }
    }
    //when the field point z is in the 2nd layer
    else {
        for(i = 0 ; i < n ; i++){
            ks2 = ks[i]*ks[i];
            k1z = sqrt(k1*k1-ks2);
            k2z = sqrt(k2*k2-ks2);
            k3z = sqrt(k3*k3-ks2);

            R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
            R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

            A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));

            spectral_g8_r = A_2_TM*R_23_TM;

            if (fabs(rho-0.0)<10e-13){
                sum[i] = ks2*weight[i]*spectral_g8_r*(ks[i]*0.5)*exp(ci*k2z*z+ci*k1z*zp+2*ci*k2z*thickness);
            }
            else{
                sum[i] = ks2*weight[i]*spectral_g8_r*(J1(ks[i]*rho)/rho)*exp(ci*k2z*z+ci*k1z*zp+2*ci*k2z*thickness);
            }
        }
    }



    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }


    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}


complex g9r(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_23_TM, A_2_TM;
    complex spectral_g9_r;

    sum2 = complex(0.0,0.0);
    ci   = complex(0.0,1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    //when the field point z is in the 1st layer
    if( z > 0 || fabs(z-0.0)<10e-13 ){
        for(i = 0 ; i < n ; i++){
            sum[i] = 0.0;
        }
    }
    //when the field point z is in the 2nd layer
    else{
        for(i = 0 ; i < n ; i++){
            ks2 = ks[i]*ks[i];
            k1z = sqrt(k1*k1-ks2);
            k2z = sqrt(k2*k2-ks2);
            k3z = sqrt(k3*k3-ks2);

            R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
            R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

            A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
            spectral_g9_r = (k2z/k1z)*A_2_TM*R_23_TM;
            if (fabs(rho-0.0)<10e-13){
                sum[i] = ks2*weight[i]*spectral_g9_r*(ks[i]*0.5)*exp(ci*k2z*z+ci*k1z*zp+2*ci*k2z*thickness);
            }
            else{
                sum[i] = ks2*weight[i]*spectral_g9_r*(J1(ks[i]*rho)/rho)*exp(ci*k2z*z+ci*k1z*zp+2*ci*k2z*thickness);
            }
        }
    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }


    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}





complex g8t(complex k, double rho, double z, double zp, double thickness, complex *epsilon, complex k1, complex k2, complex k3, double lambda)
{
    int i;
    int n1, n2, n3, n4, n5, n6, n7, n;
    complex *ks, *weight;
    complex *sum, sum2, ci;
    complex k1z, k2z, k3z, ks2;
    complex R_12_TM, R_21_TM, T_21_TM, R_23_TM, T_23_TM;
    complex spectral_g8_t;
    complex A_1_TM, A_2_TM, A_3_TM;

    sum2 = complex(0.0,0.0);
    ci = complex(0.0, 1.0);

    n1 = 200;
    n2 = 20;
    n3 = 20;
    n4 = 200;
    n5 = 20;
    n6 = 20;
    n7 = 200;
    n = n1+n2+n3+n4+n5+n6+n7;

    ks     = new complex[n];
    weight = new complex[n];
    sum    = new complex[n];

    sommerfeld_quadrature_gauss_adaptive_two_pole(ks, weight, k1, k2, k3, n1, n2, n3, n4, n5, n6, n7);

    //When the source point is in the 1st layer  (  z' >= 0  )
    if(zp > 0 || fabs(z-zp)< 10e-13){
        //When z is in the 2nd layer
        if( z < 0.0 && z > -thickness){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

                A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));

                spectral_g8_t = k2z*A_2_TM/k1z;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*weight[i]*spectral_g8_t*(ks[i]*0.5)*exp(ci*(-k2z*z+k1z*zp));
                }
                else{
                    sum[i] = ks2*weight[i]*spectral_g8_t*(J1(ks[i]*rho)/rho)*exp(ci*(-k2z*z+k1z*zp));
                }
            }
        }
        //When z is in the 3rd layer
        else if(  z < -thickness || fabs(z+thickness)<10e-13  ){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);

                A_2_TM = (1.0+R_12_TM)/(1.0+R_12_TM*R_23_TM*exp(2*ci*k2z*thickness));
                A_3_TM = A_2_TM*(1.0+R_23_TM)*exp(ci*(k2z-k3z)*thickness);

                spectral_g8_t = k3z*A_3_TM/k1z;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*weight[i]*spectral_g8_t*(ks[i]*0.5)*exp(ci*(-k3z*z+k1z*zp));
                }
                else{
                    sum[i] = ks2*weight[i]*spectral_g8_t*(J1(ks[i]*rho)/rho)*exp(ci*(-k3z*z+k1z*zp));
                }
            }
        }
        // (This case never happens)
        else{

        }
    }
    //When the source point is in the 2nd layer  (  z' >= 0  )
    else if( zp < 0.0 && zp > -thickness){
        if( z > 0.0 || fabs(z-0.0)<10e-13){
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;
                T_21_TM = (2*epsilon[0]*k2z)/(epsilon[0]*k2z+epsilon[1]*k1z);

                A_1_TM = (T_21_TM)*(exp(-ci*k2z*zp)+R_23_TM*exp(ci*k2z*zp+2*ci*k2z*thickness))/(1.0+R_23_TM*R_12_TM*exp(2*ci*k2z*thickness));

                spectral_g8_t = k1z*A_1_TM/k2z;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*weight[i]*spectral_g8_t*(ks[i]*0.5)*exp(ci*k1z*z);
                }
                else{
                    sum[i] = ks2*weight[i]*spectral_g8_t*(J1(ks[i]*rho)/rho)*exp(ci*k1z*z);
                }
            }
        }
        else if( z < 0.0 && z > -thickness){

        }
        else{
            for(i = 0 ; i < n ; i++){
                ks2 = ks[i]*ks[i];
                k1z = sqrt(k1*k1-ks2);
                k2z = sqrt(k2*k2-ks2);
                k3z = sqrt(k3*k3-ks2);

                R_12_TM = (epsilon[1]*k1z-epsilon[0]*k2z)/(epsilon[1]*k1z+epsilon[0]*k2z);
                R_23_TM = (epsilon[2]*k2z-epsilon[1]*k3z)/(epsilon[2]*k2z+epsilon[1]*k3z);
                R_21_TM = -R_12_TM;
                T_23_TM = (2*epsilon[2]*k2z)/(epsilon[2]*k2z+epsilon[1]*k3z);

                A_3_TM = T_23_TM*exp(ci*(k2z-k3z)*thickness)*(exp(ci*k2z*zp)+R_21_TM*exp(-ci*k2z*zp))/(1.0-R_23_TM*R_21_TM*exp(2*ci*k2z*thickness));

                spectral_g8_t = k3z*A_3_TM/k2z;

                if(fabs(rho-0.0)<10e-13){
                    sum[i] = ks2*weight[i]*spectral_g8_t*(ks[i]*0.5)*exp(-ci*k3z*z);
                }
                else{
                    sum[i] = ks2*weight[i]*spectral_g8_t*(J1(ks[i]*rho)/rho)*exp(-ci*k3z*z);
                }
            }
        }
    }
    //When the source point is in the 3rd layer
    else{

    }

    for(i = 0 ; i < n ; i++){
        sum2 = sum2+sum[i];
    }

    delete[] ks;
    delete[] weight;
    delete[] sum;

    return 2.0*PI*sum2;
}



//Adaptive contour with 2 poles (correct works for both ka < kb and kb < ka) located at k1.re and k1.re
void sommerfeld_quadrature_gauss_adaptive_two_pole(complex *ks, complex *weight, complex k1, complex k2, complex k3, int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{
    int i;
    double L; //contour truncation
    double L1;
    double L2;
    double L3;
    double L4;

    double *x1, *w1;
    double *x2, *w2;
    double *x3, *w3;
    double *x4, *w4;
    double *x5, *w5;
    double *x6, *w6;
    double *x7, *w7;

    complex *ks1, *ks2, *ks3, *ks4, *ks5, *ks6, *ks7;
    complex *weight1, *weight2, *weight3, *weight4, *weight5, *weight6, *weight7;

    if(fabs(k1.re - k2.re)<10e-13 && fabs(k1.im - k2.im)<10e-13){
        k2 =  k3;
    }

    if(k1.re < k2.re){
        L1 = k1.re-0.01;
        L2 = k1.re+0.01;
        L3 = k2.re-0.01;
        L4 = k2.re+0.01;
        L  = k2.re+40.0;
    }
    else{
        L1 = k2.re-0.01;
        L2 = k2.re+0.01;
        L3 = k1.re-0.01;
        L4 = k1.re+0.01;
        L  = k1.re+40.0;
    }

    ks1 = new complex[n1];
    ks2 = new complex[n2];
    ks3 = new complex[n3];
    ks4 = new complex[n4];
    ks5 = new complex[n5];
    ks6 = new complex[n6];
    ks7 = new complex[n7];

    weight1 = new complex[n1];
    weight2 = new complex[n2];
    weight3 = new complex[n3];
    weight4 = new complex[n4];
    weight5 = new complex[n5];
    weight6 = new complex[n6];
    weight7 = new complex[n7];


    x1 = new double[n1];
    w1 = new double[n1];
    gauss_quadrature(n1, x1, w1);
    for(i = 0 ; i < n1 ; i++){
        ks1[i].re = L1*(x1[i]+1.0)*0.5;
        ks1[i].im = 0.0;
        weight1[i].re = (w1[i]*0.5)*L1;
        weight1[i].im = 0.0;
    }

    if(k1.re < k2.re){
        x2 = new double[n2];
        w2 = new double[n2];
        generalized_gauss_quadrature(n2, x2, w2);
        for(i = 0 ; i < n2 ; i++){
            ks2[i] = -(k1-L1)*x2[i]+k1;
            weight2[i] = w2[i]*(k1-L1);
        }
    }
    else{
        x2 = new double[n2];
        w2 = new double[n2];
        generalized_gauss_quadrature(n2, x2, w2);
        for(i = 0 ; i < n2 ; i++){
            ks2[i] = -(k2-L1)*x2[i]+k2;
            weight2[i] = w2[i]*(k2-L1);
        }
    }

    if( k1.re < k2.re){
        x3 = new double[n3];
        w3 = new double[n3];
        generalized_gauss_quadrature(n3, x3, w3);
        for(i = 0 ; i < n3 ; i++){
            ks3[i] = (L2-k1)*x3[i]+k1;
            weight3[i] = w2[i]*(L2-k1);
        }
    }
    else{
        x3 = new double[n3];
        w3 = new double[n3];
        generalized_gauss_quadrature(n3, x3, w3);
        for(i = 0 ; i < n3 ; i++){
            ks3[i] = (L2-k2)*x3[i]+k2;
            weight3[i] = w2[i]*(L2-k2);
        }
    }

    x4 = new double[n4];
    w4 = new double[n4];
    gauss_quadrature(n4, x4, w4);
    for(i = 0 ; i < n4 ; i++){
        ks4[i].re = (L3-L2)*(x4[i]+1.0)*0.5+L2;
        ks4[i].im = 0.0;
        weight4[i].re = (w4[i]*0.5)*(L3-L2);
        weight4[i].im = 0.0;
    }

    if(k1.re < k2.re){
        x5 = new double[n5];
        w5 = new double[n5];
        generalized_gauss_quadrature(n5, x5, w5);
        for(i = 0 ; i < n5 ; i++){
            ks5[i] = -(k2-L3)*x5[i]+k2;
            weight5[i] = w5[i]*(k2-L3);
        }
    }
    else{
        x5 = new double[n5];
        w5 = new double[n5];
        generalized_gauss_quadrature(n5, x5, w5);
        for(i = 0 ; i < n5 ; i++){
            ks5[i] = -(k1-L3)*x5[i]+k1;
            weight5[i] = w5[i]*(k1-L3);
        }
    }

    if(k1.re < k2.re){
        x6 = new double[n6];
        w6 = new double[n6];
        generalized_gauss_quadrature(n6, x6, w6);
        for(i = 0 ; i < n6 ; i++){
            ks6[i] = (L4-k2)*x6[i]+k2;
            weight6[i] = w6[i]*(L4-k2);
        }
    }
    else{
        x6 = new double[n6];
        w6 = new double[n6];
        generalized_gauss_quadrature(n6, x6, w6);
        for(i = 0 ; i < n6 ; i++){
            ks6[i] = (L4-k1)*x6[i]+k1;
            weight6[i] = w6[i]*(L4-k1);
        }
    }

    x7 = new double[n7];
    w7 = new double[n7];
    gauss_quadrature(n7, x7, w7);
    for(i = 0 ; i < n7 ; i++){
        ks7[i].re = (L-L4)*(x7[i]+1.0)*0.5+L4;
        ks7[i].im = 0.0;
        weight7[i].re = (w7[i]*0.5)*(L-L4);
        weight7[i].im = 0.0;
    }


    for(i = 0 ; i < n1 ; i++){
        ks[i] = ks1[i];
        weight[i] = weight1[i];
    }
    for(i = n1 ; i < n1+n2 ; i++){
        ks[i] = ks2[i-n1];
        weight[i] = weight2[i-n1];
    }
    for(i = n1+n2 ; i < n1+n2+n3 ; i++){
        ks[i] = ks3[i-n1-n2];
        weight[i] = weight3[i-n1-n2];
    }
    for(i = n1+n2+n3 ; i < n1+n2+n3+n4 ; i++){
        ks[i] = ks4[i-n1-n2-n3];
        weight[i] = weight4[i-n1-n2-n3];
    }
    for(i = n1+n2+n3+n4 ; i < n1+n2+n3+n4+n5 ; i++){
        ks[i] = ks5[i-n1-n2-n3-n4];
        weight[i] = weight5[i-n1-n2-n3-n4];
    }
    for(i = n1+n2+n3+n4+n5 ; i < n1+n2+n3+n4+n5+n6 ; i++){
        ks[i] = ks6[i-n1-n2-n3-n4-n5];
        weight[i] = weight6[i-n1-n2-n3-n4-n5];
    }
    for(i = n1+n2+n3+n4+n5+n6 ; i < n1+n2+n3+n4+n5+n6+n7 ; i++){
        ks[i] = ks7[i-n1-n2-n3-n4-n5-n6];
        weight[i] = weight7[i-n1-n2-n3-n4-n5-n6];
    }

    delete[] ks1;
    delete[] ks2;
    delete[] ks3;
    delete[] ks4;
    delete[] ks5;
    delete[] ks6;
    delete[] ks7;

    delete[] x1;
    delete[] x2;
    delete[] x3;
    delete[] x4;
    delete[] x5;
    delete[] x6;
    delete[] x7;

    delete[] weight1;
    delete[] weight2;
    delete[] weight3;
    delete[] weight4;
    delete[] weight5;
    delete[] weight6;
    delete[] weight7;

    delete[] w1;
    delete[] w2;
    delete[] w3;
    delete[] w4;
    delete[] w5;
    delete[] w6;
    delete[] w7;
}


//Gaussian quadrature for [-1,1]
void gauss_quadrature(int n, double *x, double *w)
{
    int i, j, m;
    double xm, xl, z, zl;
    double p1, p2, p3, pp;
    double *xtemp, *wtemp, x1, x2;

    x1 = -1.0;
    x2 =  1.0;

    xtemp  = (double *)calloc(sizeof(double), n+1);
    wtemp  = (double *)calloc(sizeof(double), n+1);

    for(i = 0 ; i <= n ; i++){
        xtemp[i] = 0.0;
        wtemp[i] = 0.0;
    }

    m = (int)( n + 1 )/ (int)2.0;

    xm = ( x2 + x1 ) * 0.5;
    xl = ( x2 - x1 ) * 0.5;

    for(i = 1; i <= m ; i++){
        z = cos(PI*(i-0.25)/((double)n+0.5));
        do{
            p1 = 1.0;
            p2 = 0.0;
            for(j = 1 ; j <= n ; j++){
                p3 = p2;
                p2 = p1;
                p1 = ((2.0*(double)j-1.0)*z*p2-((double)j-1.0)*p3)/(double)j;
            }
            pp = (double)n*(z*p1-p2)/(z*z-1.0);
            zl = z;
            z = zl- p1/pp;

        } while( fabs( z - zl ) > 10e-14);
        xtemp[i] = xm-xl*z;
        xtemp[n+1-i]=xm+xl*z;
        wtemp[i] = 2.0*xl/((1.0-z*z)*pp*pp);
        wtemp[n+1-i] = wtemp[i];
    }

    for(i = 0 ; i < n ; i++){
        x[i] = xtemp[i+1];
        w[i] = wtemp[i+1];
    }

    free(xtemp);
    free(wtemp);
}


//Quadrature for x^-0.5 in [0,1] interval
void generalized_gauss_quadrature(int ngq, double *xi, double *wi)
{
    if(ngq ==5){
        xi[0]=0.220055532702321e-2;
        xi[1]=0.532526444285811e-1;
        xi[2]=0.25;
        xi[3]=0.591721954534264;
        xi[4]=0.908380401265687;

        wi[0]=0.111142584286221e-1;
        wi[1]=0.110450910249386;
        wi[2]=0.284444444444444;
        wi[3]=0.368177760249980;
        wi[4]=0.225812626627567;
    }
    else if(ngq ==10){
        xi[0]=0.170217313506295e-3;
        xi[1]=0.455197375232787e-2;
        xi[2]=0.256945562245545e-1;
        xi[3]=0.802601948484878e-1;
        xi[4]=0.181103722710989;
        xi[5]=0.329978061692620;
        xi[6]=0.513655588977735;
        xi[7]=0.705104124523579;
        xi[8]=0.869615340441312;
        xi[9]=0.974076745830678;

        wi[0]=0.869843410720295e-3;
        wi[1]=0.100832309490842e-1;
        wi[2]=0.351184957693976e-1;
        wi[3]=0.762838816843756e-1;
        wi[4]=0.125764125553643;
        wi[5]=0.169760099161110;
        wi[6]=0.192982837625621;
        wi[7]=0.183967866746584;
        wi[8]=0.139368118201496;
        wi[9]=0.658015008979678e-1;
    }
    else if(ngq == 15){
        xi[0] = 0.360449058720925e-4;
        xi[1] = 0.983656825228959e-3;
        xi[2] = 0.576031032998391e-2;
        xi[3] = 0.189863966971689e-1;
        xi[4] = 0.460162191690594e-1;
        xi[5] = 0.917631475619828e-1;
        xi[6] = 0.159522718866145;
        xi[7] = 0.250000000000000;
        xi[8] = 0.360716812863579;
        xi[9] = 0.485914494639546;
        xi[10] = 0.616988391777598;
        xi[11] = 0.743404128057339;
        xi[12] = 0.853966893740411;
        xi[13] = 0.938257049225935;
        xi[14] = 0.988028562926358;

        wi[0] = 0.184634499540016e-3;
        wi[1] = 0.220691172454993e-2;
        wi[2] = 0.813303209689366e-2;
        wi[3] = 0.192316020292443e-1;
        wi[4] = 0.356670580668843e-1;
        wi[5] = 0.563926955430610e-1;
        wi[6] = 0.792541212080791e-1;
        wi[7] = 0.101289120962781;
        wi[8] = 0.119177364119033;
        wi[9] = 0.129768304472501;
        wi[10] = 0.130602147750110;
        wi[11] = 0.120339075896910;
        wi[12] = 0.990261883702783e-1;
        wi[13] = 0.681591357635582e-1;
        wi[14] = 0.305686074965773e-1;
    }
    else if(ngq == 20){
        xi[0] = 0.118040372897702e-4;
        xi[1] = 0.324505506016988e-3;
        xi[2] = 0.192569889609293e-2;
        xi[3] = 0.647083718891321e-2;
        xi[4] = 0.160868754200355e-1;
        xi[5] = 0.331142308281794e-1;
        xi[6] = 0.598127724451431e-1;
        xi[7] = 0.980610158280346e-1;
        xi[8] = 0.149078672924258;
        xi[9] = 0.213200816542450;
        xi[10] = 0.289727337675948;
        xi[11] = 0.376864524065904;
        xi[12] = 0.471767104543454;
        xi[13] = 0.570679774395970;
        xi[14] = 0.669167911554694;
        xi[15] = 0.762418781880186;
        xi[16] = 0.845587809011132;
        xi[17] = 0.914160127147419;
        xi[18] = 0.964296432783931;
        xi[19] = 0.993140403222385;


        wi[0] = 0.605164515048587e-4;
        wi[1] = 0.731395632734516e-3;
        wi[2] = 0.275022407735182e-2;
        wi[3] = 0.669890718081941e-2;
        wi[4] = 0.129282095841639e-1;
        wi[5] = 0.215082324328235e-1;
        wi[6] = 0.322066292668297e-1;
        wi[7] = 0.444969640416654e-1;
        wi[8] = 0.575967453908003e-1;
        wi[9] = 0.705318509111265e-1;
        wi[10] = 0.822215362195993e-1;
        wi[11] = 0.915762410818034e-1;
        wi[12] = 0.975991452767166e-1;
        wi[13] = 0.994820091823469e-1;
        wi[14] = 0.966862995286949e-1;
        wi[15] = 0.890019102330765e-1;
        wi[16] = 0.765778343958853e-1;
        wi[17] = 0.599218242567572e-1;
        wi[18] = 0.398700341676524e-1;
        wi[19] = 0.175534906876473e-1;
    }

    else{
        printf("Error : use 5, 10, 15, 20 quadrature points\n\n");
    }
}


//Complex argument 0th order the first kind Bessel function
complex J0(complex z)
{
    int k, k0;
    double a0;
    complex j0;
    complex z1, z2, rp2, cbj0, cr, ct1, cp0, cq0, cu;
    complex nquarter, arg, temp;
    double a[12], b[12];

    a[0] = -7.03125e-2;
    a[1] =  0.112152099609375;
    a[2] = -0.5725014209747314;
    a[3] =  6.074042001273483;
    a[4] = -1.100171402692467e+2;
    a[5] =  3.038090510922384e+3;
    a[6] = -1.188384262567832e+5;
    a[7] =  6.252951493434797e+6;
    a[8] = -4.259392165047669e+8;
    a[9] =  3.646840080706556e+10;
    a[10]= -3.833534661393944e+12;
    a[11]=  4.854014686852901e+14;


    b[0] =	7.32421875e-2;
    b[1] = -0.2271080017089844;
    b[2] =	1.727727502584457;
    b[3] = -2.438052969955606e+1;
    b[4] =	5.513358961220206e+2;
    b[5] = -1.825775547429318e+4;
    b[6] =	8.328593040162893e+5;
    b[7] = -5.006958953198893e+7;
    b[8] =	3.836255180230433e+9;
    b[9] = -3.649010818849833e+11;
    b[10]=	4.218971570284096e+13;
    b[11]= -5.827244631566907e+15;

    nquarter.re = -0.25;
    nquarter.im = -0.0;

    z1 = z;
    z2 = z*z;
    rp2.re = 2.0/PI;
    rp2.im = 0.0;

    a0 = abs(z);

    if(fabs(a0-0.0)<10e-13){
        j0.re = 1.0;
        j0.im = 0.0;
    }

    if(z.re < 0.0){
        z1.re = -z.re;
        z1.im = -z.im;
    }

    //For small argument (|z|<=12)
    if(a0 < 12.0 || fabs(a0-12.0)<10e-13){
        cbj0.re = 1.0;
        cbj0.im = 0.0;
        cr.re = 1.0;
        cr.im = 0.0;
        for( k = 1 ; k <= 40 ; k++){
            arg.re = (double)k*(double)k;
            arg.im = 0.0;
            cr = nquarter*cr*z2/arg;
            cbj0 = cbj0+cr;
            if(abs(cr) < abs(cbj0)*10e-13){
                break;
            }
        }
        j0 = cbj0;
    }
    //For large argument (|z|>12)
    else{
        k0 = 12;
        if(a0 > 35.0 || fabs(a0 - 35.0)<10e-13 ){
            k0 = 10;
        }
        if(a0 > 50.0 || fabs(a0 - 50.0)<10e-13){
            k0 = 8;
        }

        temp.re = 0.25*PI;
        temp.im = 0.0;
        ct1 = z1-temp;
        cp0.re = 1.0;
        cp0.im = 0.0;
        for(k = 0 ; k < k0 ; k++){
            arg.re = a[k];
            arg.im = 0.0;
            cp0 = cp0+arg*pow(z1, -2.0*(double)k-2.0);
        }
        temp.re = -0.125;
        temp.im =  0.0;

        cq0 = temp/z1;
        for(k = 0 ; k < k0 ; k++){
            arg.re = b[k];
            arg.im = 0.0;
            cq0 = cq0 + arg*pow(z1, -2.0*(double)k-3.0);
        }
        cu = sqrt(rp2/z1);
        j0 = cu*(cp0*cos(ct1) - cq0*sin(ct1));
    }
    return j0;
}



//Complex argument 1st order the first kind Bessel function
complex J1(complex z)
{
    int k, k0;
    double a0;
    complex z1, z2, cbj1, cr, cp1, ct2, cq1, cu, rp2, arg, j1;
    complex nquarter, half, temp;
    double a1[12], b1[12];

    a1[0] =  0.1171875;
    a1[1] = -0.1441955566406250;
    a1[2] =  0.6765925884246826;
    a1[3] = -6.883914268109947;
    a1[4] =	 1.215978918765359e2;
    a1[5] = -3.302272294480852e3;
    a1[6] =	 1.276412726461746e5;
    a1[7] = -6.656367718817688e6;
    a1[8] =	 4.502786003050393e8;
    a1[9] = -3.833857520742790e10;
    a1[10]=	 4.011838599133198e12;
    a1[11]= -5.060568503314727e14;

    b1[0] =	-0.1025390625;
    b1[1] =  0.2775764465332031;
    b1[2] = -1.993531733751297;
    b1[3] =  2.724882731126854e1;
    b1[4] = -6.038440767050702e2;
    b1[5] =  1.971837591223663e4;
    b1[6] = -8.902978767070678e5;
    b1[7] =  5.310411010968522e7;
    b1[8] = -4.043620325107754e9;
    b1[9] =  3.827011346598605e11;
    b1[10]= -4.406481417852278e13;
    b1[11]=  6.065091351222699e15;

    nquarter.re = -0.25;
    nquarter.im = -0.0;
    half.re = 0.5;
    half.im = 0.0;
    a0 = abs(z);
    z1 = z;
    z2 = z*z;
    rp2.re = 2.0/PI;
    rp2.im = 0.0;


    if( fabs(a0-0.0)<10e-13 ){
        j1.re = 0.0;
        j1.im = 0.0;
    }

    if(z.re < 0.0){
        z1.re = -z.re;
        z1.im = -z.im;
    }

    //For small argument (|z|<=12)
    if( a0 < 12.0 || fabs(a0-0.0) < 10e-13 ){
        cbj1.re = 1.0;
        cbj1.im = 0.0;
        cr.re = 1.0;
        cr.im = 0.0;
        for(k = 1 ; k <= 40 ; k++){
            arg.re = (double)k*((double)k+1.0);
            arg.im = 0.0;
            cr = (nquarter*cr*z2)/arg;
            cbj1 = cbj1+cr;
            if(abs(cr) < abs(cbj1)*10e-13){
                break;
            }
        }
        cbj1 = 0.5*z1*cbj1;
    }
    //For large argument (|z|>12)
    else{
        k0 = 12;
        if(a0 > 35.0 || fabs(a0 - 35.0)<10e-13 ){
            k0 = 10;
        }
        if(a0 > 50.0 || fabs(a0 - 50.0)<10e-13){
            k0 = 8;
        }

        temp.re = 0.75*PI;
        temp.im = 0.0;
        ct2 = z1-temp;

        cp1.re = 1.0;
        cp1.im = 0.0;
        for(k = 0 ; k < k0 ; k++){
            arg.re = a1[k];
            arg.im = 0.0;
            cp1 = cp1+ arg*pow(z1, -2.0*(double)k-2.0);
        }
        temp.re = 0.375;
        temp.im = 0.0;
        cq1 = temp/z1;
        for(k = 0 ; k < k0 ; k++){
            arg.re = b1[k];
            arg.im = 0.0;
            cq1 = cq1+arg*pow(z1, -2.0*(double)k-3.0);
        }
        cu = sqrt(rp2/z1);
        cbj1 = cu*(cp1*cos(ct2)-cq1*sin(ct2));
    }

    if(z.re < 0.0){
        j1.re= -cbj1.re;
        j1.im= -cbj1.im;
    }
    else{
        j1 = cbj1;
    }

    return j1;
}

//Complex argument 2nd order the first kind Bessel function
complex J2(complex z)
{
    return (2.0/z)*J1(z)-J0(z);
}
