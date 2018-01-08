#ifndef _MBPT_H_
#define _MBPT_H_

#include<stdio.h>
#include<math.h>

#include<gsl/gsl_odeiv.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_errno.h>

#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>

#include<gsl/gsl_randist.h>
#include<gsl/gsl_rng.h>

#include<gsl/gsl_statistics_double.h>

#include<gsl/gsl_integration.h>
/*================================================================================================*/
/*  Cosmological Parameters     */
#define omegab  0.044
#define omegac  0.236
#define omegam  (omegab + omegac)
#define omegav  (1. - omegam)
#define omegag  1E-4
#define h       0.71
#define H0      100                     /*  km/s/Mpc    */
#define Tcmb    2.726                   /*  K           */

/*================================================================================================*/
/*  Moving background velocity  */
#define VBC     30                      /*  km/s        */

/*================================================================================================*/
/*  Some Physical Constatns     */
#define gamma   (5./3.)
#define MU      1.22
#define mH      510998.902              /*  ev/c^2      */
#define kB      8.617342E-5             /*  ev/K        */
#define c_light 2.99792458E5            /*  km/s        */

/*================================================================================================*/
/*  Some Mathematical Constatns     */
#define Pi      3.1415926535

/*================================================================================================*/
/*  Primordial PowerSpectrum    */
#define P_zeta2 2.42E-9

/*================================================================================================*/
/*  Data Structure  */



/*================================================================================================*/
/*  intial condition array size */
#define     ic_size 582

/*================================================================================================*/
/*  background equations    */
double H(double a);                     /*  km/s/Mpc    */
double Tb(double a);                    /*  K           */
double cs(double a);                    /*  km/s        */

/*================================================================================================*/
/*  Initial conditions      */
void setup_ode_initial( char fname[], 
                        double lnk[], 
                        double deltac[], double deltab[],
                        double thetac[], double thetab[] );
void set_ode_ini_k_mode(double deltac[], double deltab[], 
                        double thetac[], double thetab[], 
                        double lnk[], double ktemp,
                        double delta_theta[] );

/*================================================================================================*/
/*  ODEs                    */
int MBPT_ode_system(double t, const double y[], double f[], void *params);


/*================================================================================================*/
/*  Interface to transfer functions     */
void transfer_solver(double kmin, double kmax, double mu, double vbc);
void transfer_solver_k_mode(double k, double mu, double vbc, double tfk[]);


/*================================================================================================*/
/*  Power Spectrum  */
double power_of_k_mu(double mu, void *params);
double power_k_averaged_over_mu_vbc(double k, double N_ensemble);

/*================================================================================================*/
/*  GSL tools       */
/*  interpolation   */
double gsl_interpolation(double x[], double y[], int size, double xtemp);
/*  Integration : QUAG  */
double gsl_integration(double (*f)(), double *params, double x1, double x2);
/*  Use GSL Gaussian Random Number Generator    */
double vbc_generator(double sigma_vbc);

#endif
