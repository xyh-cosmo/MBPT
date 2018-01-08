#include"../include/mbpt.h"

extern double   lnk[ic_size];
extern double   deltac[ic_size];
extern double   deltab[ic_size];
extern double   thetac[ic_size];
extern double   thetab[ic_size];

/*  the ODE system  */
int MBPT_ode_system(double t, const double y[], double f[], void *params)
{
    double a        = t;
    double *p       = (double *)params;
    double k        = p[0];
    double mu       = p[1];
    double vbc      = p[2];
    double aH       = a*H(a);
    double cs2      = cs(a)*cs(a);
//printf("k = %10.5E\t mu = %10.5E\t vbc = %10.5E\n",k, mu, vbc);exit(0);

    double deltac_r = y[0];
    double deltac_i = y[1];
    double thetac_r = y[2];
    double thetac_i = y[3];

    double deltab_r = y[4];
    double deltab_i = y[5];
    double thetab_r = y[6];
    double thetab_i = y[7];

    /*  vbc scales as a^-1  */
    double a1020    = 1./1021.;
    vbc = vbc*a1020/a;
    /*  dy/dt   */
    double deltac_r_dot =-vbc*mu*k*deltac_i/a/aH - thetac_r/aH;
    double deltac_i_dot = vbc*mu*k*deltac_r/a/aH - thetac_i/aH;
    double thetac_r_dot =-vbc*mu*k*thetac_i/a/aH 
                        - 1.5*(omegac*deltac_r + omegab*deltab_r)*aH/a/a
                        - 2.0*thetac_r/a;
    double thetac_i_dot = vbc*mu*k*thetac_r/a/aH
                        - 1.5*(omegac*deltac_i + omegab*deltab_i)*aH/a/a
                        - 2.0*thetac_i/a;

    double deltab_r_dot =-thetab_r/aH;
    double deltab_i_dot =-thetab_i/aH;
    double thetab_r_dot =-1.5*(omegac*deltac_r + omegab*deltab_r)*aH/a/a
                        - 2.0*thetab_r/a
                        + cs2*k*k*deltab_r/aH/a;
    double thetab_i_dot =-1.5*(omegac*deltac_i + omegab*deltab_i)*aH/a/a
                        - 2.0*thetab_i/a
                        + cs2*k*k*deltab_i/aH/a;

    /*  put dy/dt into f[]  */
    f[0]    = deltac_r_dot;
    f[1]    = deltac_i_dot;
    f[2]    = thetac_r_dot;
    f[3]    = thetac_i_dot;
    f[4]    = deltab_r_dot;
    f[5]    = deltab_i_dot;
    f[6]    = thetab_r_dot;
    f[7]    = thetab_i_dot;

    return GSL_SUCCESS;
}


void transfer_solver(double kmin, double kmax, double mu, double vbc)
{
    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(T, 8);
    gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(1E-10, 0.0);
    gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(8);

//    double k, kmin = 1E-3, kmax = 10;
    double k;
    double dlnk = (log(kmax) - log(kmin))/100.;

    double ai = 1./(1020. + 1.), ae = 1./(40. + 1.);
    double da = 1E-8;

    double delta_theta[8];

    /*  angle & vbc */
    double params[3] = {0, mu, vbc};
    int i;
    for(i=0; i<100; i++)
    {
        k   = exp(log(kmin) + i*dlnk);
        params[0]    = k;
        gsl_odeiv_system sys = {MBPT_ode_system, NULL, 8, params};
        
        ai  = 1./(1020. + 1.);
        da  = 1E-6;

        set_ode_ini_k_mode( deltac, deltab, 
                            thetac, thetab, 
                            lnk, k,
                            delta_theta );

        while(ai < ae){
            int status = gsl_odeiv_evolve_apply(e, c, s, 
                                                &sys, 
                                                &ai, ae, 
                                                &da, delta_theta);
            if(status != GSL_SUCCESS)
                break;
        }

        printf("k = %10.5E\t|deltac| = %10.5E\t|deltab| = %10.5E\n",
                k,
                sqrt(delta_theta[0]*delta_theta[0] + delta_theta[1]*delta_theta[1]),
                sqrt(delta_theta[4]*delta_theta[4] + delta_theta[5]*delta_theta[5]) );

    }

    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
}


/*================================================================================================*/
/***************************************************************************************************
    The following functions are used to check (output) the transfer function(s)
***************************************************************************************************/
/*================================================================================================*/
void transfer_solver_k_mode(double k, double mu, double vbc, double tfk[])
{
    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
//    const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
//    const gsl_odeiv_step_type *T = gsl_odeiv_step_rkck;
//    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4imp;
    gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(T, 8);
    gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(1E-5, 0);
    gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(8);

    double ai = 1./(1020. + 1.), ae = 1./(40. + 1.);
    double da = 1E-6;

    double delta_theta[8];

    //  angle & vbc 
    double params[3] = {k, mu, vbc};

    gsl_odeiv_system sys = {MBPT_ode_system, NULL, 8, params};

    set_ode_ini_k_mode( deltac, deltab, 
                        thetac, thetab, 
                        lnk, k,
                        delta_theta );

    while(ai < ae){
        int status = gsl_odeiv_evolve_apply(e, c, s, 
                                            &sys, 
                                            &ai, ae, 
                                            &da, delta_theta);
        if(status != GSL_SUCCESS)
            break;
    }

    tfk[0]  = delta_theta[0];
    tfk[1]  = delta_theta[1];
    tfk[2]  = delta_theta[2];
    tfk[3]  = delta_theta[3];
    tfk[4]  = delta_theta[4];
    tfk[5]  = delta_theta[5];
    tfk[6]  = delta_theta[6];
    tfk[7]  = delta_theta[7];

    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
}

/*================================================================================================*/
void write_transfer_function(   double kmin, double kmax, int nk, 
                                double mu, double vbc, char fname[])
{
    int i,j;
    double k, dlnk = (log(kmax) - log(kmin))/nk;
    double tfk[8];

    FILE *fp = fopen(fname,"w");
    for(i=0; i<nk; i++)
    {
        k   = exp(log(kmin) + i*dlnk);
        transfer_solver_k_mode(k, 0, 30, tfk);

        fprintf(fp,"%10.5E\t",k);        
        for(j=0; j<8; j++)
        {
            fprintf(fp,"%10.5E\t",tfk[j]);
        }
        fprintf(fp,"\n");
    }
}


void write_transfer_function_mu_vbc(double kmin, double kmax, int nk, char fname[])
{
    int i,j,l;
    double k, dlnk = (log(kmax) - log(kmin))/(nk - 1.);
    double tfk[8];
    double mu;

    FILE *fp = fopen(fname,"w");
    int nmu = 25;
    for(i=0; i<nk; i++)
    {
        for(j=0; j< nmu; j++)
        {
            k   = exp(log(kmin) + i*dlnk);
            mu = -1.0 + j*2.0/nmu;

            transfer_solver_k_mode(k, mu, 30, tfk);
            /*========================================================*/
            fprintf(stdout,"%6.3E %6.3E  ", k, mu);        
            
            fprintf(stdout,"%10.5E\t%10.5E\t%10.5E\t%10.5E\t%10.5E\t%10.5E\t%10.5E\t%10.5E\n ",
                            tfk[0]/k/k, tfk[1]/k/k,
                            tfk[2]/k, tfk[3]/k,
                            tfk[4]/k/k, tfk[5]/k/k,
                            tfk[6]/k, tfk[7]/k);

            /*========================================================*/
            fprintf(fp,"%10.3E\t%10.3E\t", k, mu);        
            fprintf(fp,"%10.5E\t%10.5E\t%10.5E\t%10.5E\t%10.5E\t%10.5E\t%10.5E\t%10.5E\n ",
                        tfk[0]/k/k, tfk[1]/k/k,
                        tfk[2]/k, tfk[3]/k,
                        tfk[4]/k/k, tfk[5]/k/k,
                        tfk[6]/k, tfk[7]/k);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

