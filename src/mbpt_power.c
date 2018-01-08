/***************************************************************************************************
    The power spectrum we need is averaged over angles and vbcs(gaussian distribution).
***************************************************************************************************/
#include"../include/mbpt.h"

#define N_MU    25

double power_k_averaged_over_mu_vbc(double k, double N_ensemble)
{
    int i,j;
    double *Pk      = malloc(sizeof(double)*N_ensemble);
    double *Pk_mu   = malloc(sizeof(double)*N_MU);
//    double *Pk_mu   = malloc(sizeof(double)*100);
    double Pk_averaged = 0.;

    double deltac_real, deltac_image;
    double deltab_real, deltab_image;
    double omegacm, omegabm;

    double mu, vbc;
    double transfer_function[8];

    omegacm = omegac/omegam;
    omegabm = omegab/omegam;

/*  make ensembles of Pk    */

    for(i=0; i<N_ensemble; i++)
    {
//printf("## %5d th ensemble ....\n",i);
        /*  generate a vbc  */
        vbc = vbc_generator(10.0);
        transfer_solver_k_mode(k, 1., vbc, transfer_function);
        deltac_real     = transfer_function[0];
        deltac_image    = transfer_function[1];
        deltab_real     = transfer_function[4];
        deltab_image    = transfer_function[5];

        Pk[i]           = omegacm*omegacm*(deltac_real*deltac_real + deltac_image*deltac_image)
                        + omegabm*omegabm*(deltab_real*deltab_real + deltab_image*deltab_image);

//        Pk_mu[j]        = omegacm*omegacm*(deltac_real*deltac_real + deltac_image*deltac_image)
//                        + omegabm*omegabm*(deltab_real*deltab_real + deltab_image*deltab_image);

        /*  Initialize the array stores transfer functions as a function of mu    */
/*
        for(j=0; j<N_MU; j++)
        {
//printf("J = %5d\n",j);
            mu              = -1.0 + j*2.0/(N_MU - 1);

            transfer_solver_k_mode(k, mu, vbc, transfer_function);

            deltac_real     = transfer_function[0];
            deltac_image    = transfer_function[1];
            deltab_real     = transfer_function[4];
            deltab_image    = transfer_function[5];

            Pk_mu[j]        = omegacm*omegacm*(deltac_real*deltac_real + deltac_image*deltac_image)
                            + omegabm*omegabm*(deltab_real*deltab_real + deltab_image*deltab_image);
        }
        Pk[i]   = gsl_integration(power_of_k_mu, Pk_mu, -1.0, 1.0);
*/
    }

/*  Take ensemble average of the "ENSEMBLES"    */
    Pk_averaged = gsl_stats_mean(Pk, 1, N_ensemble);

    free(Pk);
    free(Pk_mu);
}

double power_of_k_mu(double mu, void *params)
{
    int i;
    double *power_mu    = (double *)params;
    double mu_temp[N_MU];
    double power;

    for(i=0; i<N_MU; i++)
        mu_temp[i]   = -1.0 + i*2.0/(N_MU - 1);

    return gsl_interpolation(mu_temp, power_mu, N_MU, mu);
}
