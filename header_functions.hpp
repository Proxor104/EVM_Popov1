#ifndef HEADER_FUNCTIONS
#define HEADER_FUNCTIONS
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>

int solve (P_gas & gas, P_she & shem, std::vector<double> & V_old,
          std::vector<double> & H_old, std::vector<double> & V,
          std::vector<double> & H,
          double (*F) (double, double, double, double, double),
          double (*F0) (double, double));

int solve_step (const P_gas & gas, const P_she & shem,
               std::vector<double> & V_old, std::vector<double> & H_old,
               std::vector<double> & V, std::vector<double> & H,
               std::vector<double> & up_diag, std::vector<double> & mid_diag,
               std::vector<double> & down_diag, std::vector<double> & right_vec,
               std::vector<double> & alpha, std::vector<double> & beta,
               double (*F) (double, double, double, double, double),
               double (*F0) (double, double), int i);

void fill_h0_v0 (const P_she & shem, std::vector<double> & V_old,
                std::vector<double> & H_old, std::vector<double> & V,
                std::vector<double> & H, double (* velocity0)(double),
                double (* dencity0)(double));

int get_eta (const P_gas & gas, P_she & shem, std::vector<double> & H);

int print_v (std::vector<double> v);

#endif //HEADER_FUNCTIONS
