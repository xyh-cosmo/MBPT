/***************************************************************************************************
    For the moment, the intial conditions for each k-mode is set from the "ini_z_1020.dat", later
    I will modify the code to store the IC in a array, and thus lot's of time may be saved.
***************************************************************************************************/

#include"../include/mbpt.h"

//extern int      ic_size;
//extern double   lnk[ic_size];
//extern double   deltac[ic_size];
//extern double   deltab[ic_size];
//extern double   thetac[ic_size];
//extern double   thetab[ic_size];


void setup_ode_initial( char fname[], 
                        double lnk[], 
                        double deltac[], double deltab[],
                        double thetac[], double thetab[] )
/*
    We set up initial condition by interpolating the output of cmbfast at z = 1020.
    The cmbfast-output data file is stored with a name like ini_z_1020.dat
*/
{
    int i,j,count;
    double temp;
    FILE *get;

    if(!(get = fopen(fname,"r"))){
        printf("Can't open \" %s \"\n Exiting ... \n",fname);
        exit(0);
    }

    count = 0;
    while(1 == fscanf(get,"%lf ",&temp)){
        i = count/5;     //  get row number
        j = count%5;     //  get colume number
        if(j==0)        lnk[i]    = log(temp);
        else if(j==1)   deltac[i] = temp;
        else if(j==2)   deltab[i] = temp;
        else if(j==3)   thetac[i] = temp;
        else if(j==4)   thetab[i] = temp;
        count ++;
    }

    fclose(get);
}

void set_ode_ini_k_mode(double deltac[], double deltab[], 
                        double thetac[], double thetab[], 
                        double lnk[], double ktemp,
                        double delta_theta[] )
/*
    delta_theta = { deltac_r, deltac_i,
                    thetac_r, thetac_i,
                    deltab_r, deltab_i,
                    thetab_r, thetab_i  }
*/
{
    double lnk_temp = log(ktemp);
    double phase = 0.0;

    double Deltac = gsl_interpolation(lnk, deltac, ic_size, lnk_temp);
    double Deltab = gsl_interpolation(lnk, deltab, ic_size, lnk_temp);
    double Thetac = gsl_interpolation(lnk, thetac, ic_size, lnk_temp);
    double Thetab = gsl_interpolation(lnk, thetab, ic_size, lnk_temp);

    delta_theta[0]  = Deltac*cos(phase);
    delta_theta[1]  = Deltac*sin(phase);
    delta_theta[2]  = Thetac*cos(phase);
    delta_theta[3]  = Thetac*sin(phase);
    delta_theta[4]  = Deltab*cos(phase);
    delta_theta[5]  = Deltab*sin(phase);
    delta_theta[6]  = Thetab*cos(phase);
    delta_theta[7]  = Thetab*sin(phase);
}

