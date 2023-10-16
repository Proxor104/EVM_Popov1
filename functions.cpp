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
  up_diag.resize(M);
  down_diag.resize(M);
  mid_diag.resize(M + 1);
  right_vec.resize(M + 1);
  for(i = 1; i < N + 1; i++)
    {
      if (get_eta(gas, shem, H) < 0)
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
  up_diag[j] = (V_old[0] + V_old[1]) / (6 * h);
  mid_diag[j] = 1 / tau + (2 * eta) / (h * h);
  right_vec[j] = V_old[j] / tau - (constant * (std::pow(H_old[j + 1], gamma) -
                 std::pow(H_old[j - 1], gamma)))/(2 * h * H_old[j]) -
                 (eta - mu / H_old[j]) * (V_old[j - 1] - 2 * V_old[j] +
                 V_old[j + 1]) / (h * h) +
                 F(j * M, i * N, constant, gamma, mu);
  for (j = 1; j < M - 2; j++)
    {
      up_diag[j] = (V_old[j] + V_old[j + 1]) / (6 * h);
      down_diag[j] = ((-1) * (V_old[j] + V_old[j - 1])) / (6 * h) - (eta) / (h * h);
      mid_diag[j] = 1 / tau + (2 * eta) / (h * h);
      right_vec[j] = V_old[j] / tau - (constant * (std::pow(H_old[j + 1], gamma) -
                     std::pow(H_old[j - 1], gamma)))/(2 * h * H_old[j]) -
                     (eta - mu / H_old[j]) * (V_old[j - 1] - 2 * V_old[j] +
                     V_old[j + 1]) / (h * h) +
                     F(j * M, i * N, constant, gamma, mu);
    }
  down_diag[j] = ((-1) * (V_old[j] + V_old[j - 1])) / (6 * h) - (eta) / (h * h);
  mid_diag[j] = 1 / tau + (2 * eta) / (h * h);
  right_vec[j] = V_old[j] / tau - (constant * (std::pow(H_old[j + 1], gamma) -
                 std::pow(H_old[j - 1], gamma)))/(2 * h * H_old[j]) -
                 (eta - mu / H_old[j]) * (V_old[j - 1] - 2 * V_old[j] +
                 V_old[j + 1]) / (h * h) +
                 F(j * h, i * tau, constant, gamma, mu);
  // метод прогонки для V
  alpha[0] = (-1) * (up_diag[0]) / (mid_diag[0]);
  beta[0] = (right_vec[0]) / (mid_diag[0]);
  for (j = 1; j < M - 2; j++)
    {
      alpha[j] = (-1) * (up_diag[j])/(mid_diag[j]
                 + down_diag[j] * alpha[j - 1]);
      beta[j] = (right_vec[j] - down_diag[j] * beta[j - 1])/(mid_diag[j]
                + down_diag[j] * alpha[j - 1]);
    }
  V[M - 1] = (right_vec[M - 2] - down_diag[M - 2] * beta[M - 3]) / (mid_diag[M - 2]
             + down_diag[M - 2] * alpha[M - 3]);
  for(j = M - 2; j >= 1; j--)
    {
      V[j] = V[j + 1] * alpha[j - 1] + beta[j - 1];
    }
  V[M] = 0; V[0] = 0;
  // ---------------------------------------------------------------------------
  // вычисляем вектор H
  j = 0;
  up_diag[j] = (V[1]) / (2 * h);
  mid_diag[j] = 1 / tau - (V[0]) / (2 * h);
  down_diag[j] = 0;
  right_vec[j] = F0(0, tau * i) + H_old[j]/tau - H_old[j] * (V[j + 1] - V[j])/(2 * h) +
                 (h/4) * ((H_old[j]*V_old[j] - 2*H_old[j+1]*V_old[j+1] +
                 H_old[j+2]*V_old[j+2])/(h * h) - (1/2)*(H_old[j+1]*V_old[j+1] -
                 2*H_old[j+2]*V_old[j+2] + H_old[j+3]*V_old[j+3])/(h * h) +
                 H_old[j]*((V_old[j] - 2*V_old[j+1] + V_old[j+2])/(h * h) -
                 (1/2)*(V_old[j+1] - 2*V_old[j+2] + V_old[j+3])/(h * h)));
  for (j = 1; j < M; j++)
   {
     up_diag[j] = (V[j] + V[j + 1]) / (4 * h);
     down_diag[j] = ((-1) * (V[j] + V[j - 1])) / (4 * h);
     mid_diag[j] = 1 / tau;
     right_vec[j] = F0(j * h, i * tau) + H_old[j]/tau - H_old[j] * (V[j+1]-V[j-1])/(2*h);
   }
  up_diag[j] = 0;
  mid_diag[j] = 1 / tau;
  down_diag[j] = ((-1)*V[j-1])/(2*h);
  right_vec[j] = F0(h * j, tau * i) + (H_old[j])/tau - H_old[j] * (V[j] - V[j - 1])/(2 * h) -
                 (h/2) * ((H_old[j-2]*V_old[j-2] - 2*H_old[j-1]*V_old[j-1] +
                 H_old[j]*V_old[j])/(h * h) - (1/2)*(H_old[j-3]*V_old[j-3] -
                 2*H_old[j-2]*V_old[j-2] + H_old[j-1]*V_old[j-1])/(h * h) +
                 H_old[j]*((V_old[j-2] - 2*V_old[j-1] + V_old[j])/(h * h) -
                 (1/2)*(V_old[j-3] - 2*V_old[j-1] + V_old[j])/(h * h)));
  //метод прогонки для H
  alpha[0] = (-1) * (up_diag[0]) / (mid_diag[0]);
  beta[0] = (right_vec[0]) / (mid_diag[0]);
  for (j = 1; j < M; j++)
    {
      alpha[j] = (-1) * (up_diag[j])/(mid_diag[j]
                 + down_diag[j] * alpha[j - 1]);
      beta[j] = (right_vec[j] - down_diag[j] * beta[j - 1])/(mid_diag[j]
                + down_diag[j] * alpha[j - 1]);
    }
  H[M] = (right_vec[M] - down_diag[M] * beta[M - 1]) / (mid_diag[M]
         + down_diag[M] * alpha[M - 1]);
  for(j = M - 1; j >= 0; j--)
    {
      H[j] = H[j + 1] * alpha[j] + beta[j];
    }
  // ---------------------------------------------------------------------------
  j = 0;
  for(const double & elem : V)
    {
      V_old[j] = elem; j++;
    }
  j = 0;
  for(const double & elem : H)
    {
      H_old[j] = elem; j++;
    }
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
  double max = H[0];
  for (i = 1; i < M + 1; i++)
    {
      if(H[i] > max)
        {
          max = H[i];
        }
    }
  if (fabs(max) < EPS)
    {
      return -1;
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
      tmp = (V[j] - v(j * h, j * tau));
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

