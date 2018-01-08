#include"../include/mbpt.h"
double H(double a)      /*  km/s/Mpc    */
{   
    double a3   = a*a*a;
    double a4   = a3*a;
    return H0*h*sqrt(omegam/a3 + omegav + omegag/a4);
}

double Tb(double a)     /*  K       */
{
    double a1   = 1./119.;
    double a2   = 1./115.;
    
    return Tcmb/a/(1. + (a/a1)/(1. + pow(a2/a, 1.5)));
}

double cs(double a)     /*  km/s    */
{
    return sqrt(gamma*kB*Tb(a)/MU/mH)*c_light;
}

double vbc_generator(double sigma_vbc)
{
    srand(time(NULL));
    const gsl_rng_type *T;
    T   = gsl_rng_taus;
    gsl_rng *r  = gsl_rng_alloc(T);
    gsl_rng_set(r,rand());

    double vbc = gsl_ran_gaussian(r, sigma_vbc);

    gsl_rng_free(r);

    return vbc;
}
