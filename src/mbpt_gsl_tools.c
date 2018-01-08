#include"../include/mbpt.h"

/*  interpolation   */
double gsl_interpolation(double x[], double y[], int size, double xtemp)
{
    gsl_interp_accel    *acc    = gsl_interp_accel_alloc();
    gsl_spline          *spline = gsl_spline_alloc(gsl_interp_cspline, size);

    gsl_spline_init(spline, x, y, size);

    double ytemp    = gsl_spline_eval(spline, xtemp, acc);

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return ytemp;
}


/*  Integration : QUAG  */
double gsl_integration(double (*f)(), double *params, double x1, double x2)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function  = f;
    F.params    = params;

    gsl_integration_qags(&F, x1, x2, 0, 1E-7, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}
