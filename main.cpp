#include "header_classes.hpp"
#include "header_functions.hpp"

double v_init(double x)
{
  return sin(4 * M_PI * x);
}

double h_init(double x)
{
  return cos(3 * M_PI * x) + 1.5;
}

double v(double x, double t)
{
  return cos(2*M_PI*t)*sin(4*M_PI*x);
}

double h(double x, double t)
{
  return (cos(3*M_PI*x)+1.5)*exp(t);
}

double dro_dt(double x, double t)
{
  return exp(t)*(cos(3*M_PI*x) + 1.5);
}

double dro_u_dx(double x, double t)
{
  return exp(t)*cos(2*M_PI*t)*(4*M_PI*cos(4*M_PI*x)*(cos(3*M_PI*x) + 1.5) -
         3*M_PI*sin(3*M_PI*x)*sin(4*M_PI*x));
}

double du_dt(double x, double t)
{
  return (-2*M_PI)*sin(2*M_PI*t)*sin(4*M_PI*x);
}

double du_dx(double x, double t)
{
  return (4*M_PI)*cos(2*M_PI*t)*cos(4*M_PI*x);
}

double d2u_dx2(double x, double t)
{
  return (-16)*M_PI*M_PI*cos(2*M_PI*t)*sin(4*M_PI*x);
}

double dP_dx(double x, double t, double constant, double gamma)
{
  return (-3)*M_PI*constant*gamma*std::pow(((cos(3*M_PI*x)+1.5)*exp(t)),gamma - 1)
         *exp(t)*sin(3*M_PI*x);
}

double F(double x, double t, double constant, double gamma, double mu)
{
  return h(x,t)*du_dt(x,t) + h(x,t)*v(x,t)*du_dx(x,t) +
         dP_dx(x,t,constant,gamma) - mu*d2u_dx2(x,t);
}

double F0(double x, double t)
{
  return dro_dt(x,t) + dro_u_dx(x,t);
}

  // double v_init(double x)
  //   {
  //     return sin(M_PI*x);
  //   }
  // double h_init(double x)
  //   {
  //     return 1;
  //   }
  // double h(double x, double t)
  //   {
  //     return 1;
  //   }
  // double v(double x, double t)
  //   {
  //     return sin(M_PI*x);
  //   }
  // double F0(double x, double t)
  //   {
  //     return (M_PI*cos(M_PI*x));
  //   }
  // double F(double x, double t, double constant, double gamma, double mu)
  //   {
  //     return M_PI*sin(M_PI*x)*cos(M_PI*x)+mu*M_PI*M_PI*sin(M_PI*x);
  //   }


int main (int argc, char * argv [])
{
    FILE * res_test1 = fopen ("./res_test1.txt", "w");
    P_gas gas;
    P_she shem;
    param_diff(gas);
    param_she(gas, shem);

    std::vector<double> V_old;
    std::vector<double> V;
    std::vector<double> H_old;
    std::vector<double> H;

    fill_h0_v0 (shem, V_old, H_old, V, H, v_init, h_init);
    if (solve (gas, shem, V_old, H_old, V, H, F, F0) < 0)
      {
        printf ("\n\nERROR!!!!\n\n");
        return -1;
      }
    get_residual_V(V, v, shem, gas, shem.n, res_test1);

    fclose(res_test1);
    return 0;
}
