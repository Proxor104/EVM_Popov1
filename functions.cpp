#include "header_classes.hpp"
#include "header_functions.hpp"
#define EPS 1.0e-16

int solve (P_gas & gas, P_she & shem, std::vector<double> & V_old,
          std::vector<double> & H_old, std::vector<double> & V,
          std::vector<double> & H,
          double (*F) (double, double, double, double, double),
          double (*F0) (double, double))
{
  int i = 0; // номер вычисляемого временного слоя
  int N = shem.n;
  int M = shem.m_x;
  std::vector<double> up_diag;
  std::vector<double> mid_diag;
  std::vector<double> down_diag;
  std::vector<double> right_vec;
  std::vector<double> alpha;
  std::vector<double> beta;
  alpha.resize(M + 1);
  beta.resize(M + 1);
  up_diag.resize(M + 1);
  down_diag.resize(M + 1);
  mid_diag.resize(M + 1);
  right_vec.resize(M + 1);
  for(i = 1; i < N + 1; i++)
    {
      if (get_eta(gas, shem, H_old) < 0)
        {
          return -1;
        }
      solve_step(gas, shem, V_old, H_old, V, H, up_diag, mid_diag, down_diag,
                right_vec, alpha, beta, F, F0, i - 1);
    }
  return 0;
}

int solve_step (const P_gas & gas, const P_she & shem,
               std::vector<double> & V_old, std::vector<double> & H_old,
               std::vector<double> & V, std::vector<double> & H,
               std::vector<double> & up_diag, std::vector<double> & mid_diag,
               std::vector<double> & down_diag, std::vector<double> & right_vec,
               std::vector<double> & alpha, std::vector<double> & beta,
               double (*F) (double, double, double, double, double),
               double (*F0) (double, double), int i)
{
  int j = 0;
  int M = shem.m_x;
  int N = shem.n;
  double h = shem.h_x;
  double tau = shem.tau;
  double eta = shem.eta;
  double mu = gas.mu;
  double constant = gas.p_ro;
  double gamma = gas.p_gamma;
  std::fill(up_diag.begin(), up_diag.end(), 0);
  std::fill(mid_diag.begin(), mid_diag.end(), 0);
  std::fill(down_diag.begin(), down_diag.end(), 0);
  std::fill(right_vec.begin(), right_vec.end(), 0);
  std::fill(alpha.begin(), alpha.end(), 0);
  std::fill(beta.begin(), beta.end(), 0);
  // ---------------------------------------------------------------------------
  // вычисляем вектор V
  j = 0;
  // printf("\n === \n %10.3e\n === \n", std::pow(H_old[j - 1], gamma));
  up_diag[j] = 0;
  mid_diag[j] = 1;
  down_diag[j] = 0;
  right_vec[j] = 0;
  for (j = 1; j < M; j++)
    {
      up_diag[j] = (V_old[j] + V_old[j + 1]) / (6 * h) - eta / (h * h);
      down_diag[j] = ((-1) * (V_old[j] + V_old[j - 1])) / (6 * h) - (eta) / (h * h);
      mid_diag[j] = 1 / tau + (2 * eta) / (h * h);
      right_vec[j] = V_old[j] / tau - (constant * (std::pow(H_old[j + 1], gamma) -
                     std::pow(H_old[j - 1], gamma)))/(2 * h * H_old[j]) -
                     (eta - mu / H_old[j]) * (V_old[j - 1] - 2 * V_old[j] +
                     V_old[j + 1]) / (h * h) + F(j * h, i * tau, constant, gamma, mu);
      // printf("\n==H_old[j]==\n%10.3e\n=====\n",std::pow(H_old[j + 1], gamma));
    }
  down_diag[j] = 0;
  mid_diag[j] = 1;
  up_diag[j] = 0;
  right_vec[j] = 0;
  // метод прогонки для V
  run_through_method(up_diag, mid_diag, down_diag, right_vec, alpha, beta, V, M);
  // print_v(V);
  // ---------------------------------------------------------------------------
  // вычисляем вектор H
  j = 0;
  up_diag[j] = (V[1]) / (2 * h);
  mid_diag[j] = (1 / tau) - (V[0]) / (2 * h);
  down_diag[j] = 0;
  right_vec[j] = F0(0, tau * i) + H_old[j]/tau - H_old[j] * (V[j + 1] - V[j])/(2 * h) +
                 (1/(4*h)) * ((-5/2)*H_old[j+1]*V_old[j+1] + 2*H_old[j+2]*V_old[j+2]-
                (1/2)*H_old[j+3]*V_old[j+3] + H_old[j]*(2*V_old[j] + (-5/2)*V_old[j+1] + 
                2*V_old[j+2] - (1/2)*V_old[j+3]));
  // printf("\n==right_vec==\n%10.3e\n=====\n",right_vec[j]);
  for (j = 1; j < M; j++)
   {
    up_diag[j] = (V[j] + V[j + 1]) / (4 * h);
    down_diag[j] = ((-1) * (V[j] + V[j - 1])) / (4 * h);
    mid_diag[j] = 1 / tau;
    right_vec[j] = F0(j * h, i * tau) + H_old[j]/tau - H_old[j] * (V[j+1]-V[j-1])/(2*h);
    // printf("\n==right_vec==\n%10.3e\n=====\n",right_vec[j]);
   }
  up_diag[j] = 0;
  mid_diag[j] = 1 / tau + V[j] / (2*h);
  down_diag[j] = ((-1)*V[j-1])/(2*h);
  right_vec[j] = F0(h * j, tau * i) + (H_old[j])/tau - H_old[j] * (V[j] - V[j - 1])/(2 * h) -
                 (1/(4*h)) * ((-5/2)*H_old[j-1]*V_old[j-1] + 2*H_old[j-2]*V_old[j-2]-
                (1/2)*H_old[j-3]*V_old[j-3] + H_old[j]*(2*V_old[j] + (-5/2)*V_old[j-1] + 
                2*V_old[j-2] - (1/2)*V_old[j-3]));
  // printf("\n==right_vec==\n%10.3e\n=====\n",right_vec[j]);
  //метод прогонки для H
  run_through_method(up_diag, mid_diag, down_diag, right_vec, alpha, beta, H, M);
  // print_v(H);
  // ---------------------------------------------------------------------------
  for(j = 0; j < M + 1; j++)
    {
      V_old[j] = V[j];
    }
  for(j = 0; j < M + 1; j++)
    {
      H_old[j] = H[j];
    }
  // std::vector<double> vec;
  // vec.resize(M + 1);
  // for (j = 0; j < M + 1; j++)
  //  {
  //    up_diag[j] = -1;
  //    down_diag[j] = 0;
  //    mid_diag[j] = 1;
  //    right_vec[j] = 1;
  //  }
  // run_through_method(up_diag, mid_diag, down_diag, right_vec, alpha, beta, vec, M);
  // print_v(vec);
  return 0;
}

void fill_h0_v0 (const P_she & shem, std::vector<double> & V_old,
                std::vector<double> & H_old, std::vector<double> & V,
                std::vector<double> & H, double (* velocity0)(double),
                double (* dencity0)(double))
{
  double step = shem.h_x;
  int M = shem.m_x;
  double v_m = 0, h_m = 0;
  for (int i = 0; i < M + 1; i++)
    {
      v_m = velocity0(i * step);
      V_old.push_back(v_m);
      V.push_back(v_m);
      h_m = dencity0(i * step);
      H_old.push_back(h_m);
      H.push_back(h_m);
    }
}

int get_eta (const P_gas & gas, P_she & shem, std::vector<double> & H)
{
  int i = 0;
  int M = shem.m_x;
  double mu = gas.mu;
  if (fabs(H[i]) < EPS)
    {
      return -1;
    }
  double max = mu / H[0];
  double tmp = 0;
  for (i = 1; i < M + 1; i++)
    {
      tmp = mu / H[i];
      if (fabs(H[i]) < EPS)
        {
          return -1;
        }
      if(tmp > max)
        {
          max = tmp;
        }
    }
  shem.eta = max;
  return 0;
}

int run_through_method(const std::vector<double> &  up, const std::vector<double> &  mid,
                      const std::vector<double> &  down, const std::vector<double> &  f_v, 
                      std::vector<double> &  alpha, std::vector<double> &  beta, 
                      std::vector<double> &  x, int N)
{
  int j = 0;
  alpha[j] = (-1) * (up[j]) / (mid[j]);
  beta[j] = (f_v[j]) / (mid[j]);
  for (j = 1; j < N; j++)
    {
      alpha[j] = (-1) * (up[j])/(mid[j] + down[j] * alpha[j - 1]);
      beta[j] = (f_v[j] - down[j] * beta[j - 1])/(mid[j] + down[j] * alpha[j - 1]);
    }
  x[N] = (f_v[N] - down[N] * beta[N - 1]) / (mid[N] + down[N] * alpha[N - 1]);
  for(j = N - 1; j >= 0; j--)
    {
      x[j] = x[j + 1] * alpha[j] + beta[j];
    }
  return 0;
}

int print_v (std::vector<double> v)
{
  printf ("===//===\n");
  for(const double & elem : v)
    {
      printf("%10.3e ", elem);
    }
  printf ("\n======\n");
  return 0;
}

int get_residual_V(std::vector<double> V, double (*v)(double, double), const P_she & shem, 
                  const P_gas & gas, int i, FILE * test)
{
  int j = 0; int M = shem.m_x;
  double norm_l2 = 0;
  double norm_c = 0;
  double norm_w = 0;
  double h = shem.h_x;
  double tau = shem.tau;
  double tmp = 0;
  double constant = gas.p_ro;
  double gamma = gas.p_gamma;
  double mu = gas.mu;
  for(j = 0; j < M + 1; j++)
    {
      tmp = (V[j] - v(j * h, i * tau));
      norm_l2 += tmp * tmp;
      if(fabs(tmp) > norm_c)
        {
          norm_c = fabs(tmp);
        }
    }
  norm_l2 = sqrt(norm_l2);
  fprintf(test, " ____________ ____________ ____________ ____________ ____________ ____________\n");
  fprintf(test, "|     c      |    gamma   |     mu     |     l2     |     c      |     w      |\n");
  fprintf(test, "|%12.3e|%12.3e|%12.3e|%12.3e|%12.3e|%12.3e|\n", 
          constant, gamma, mu, norm_l2, norm_c, norm_w);
  return 0;
}

