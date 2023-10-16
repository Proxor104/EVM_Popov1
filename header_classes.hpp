#ifndef HEADER_CLASSES
#define HEADER_CLASSES
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>

// параметры дифференциальной задачи
typedef struct
{
    double segm_T; // задаёт временной отрезок [0, segm_T]
    double segm_X; // задаёт пространственный отрезок [0, segm_X]
    double mu; // вязкость газа в исходной задаче
    double p_ro; // задаёт константу C для давления газа в случае p = p_ro * H
    double p_gamma; // задает степень для давления газа в случае p = H ^ p_gamma
}P_gas;

void param_diff (P_gas & gas);

//параметры сетки
typedef struct
{
    int m_x; // количество пространственных точек
    int n; // количество временных точек
    int dim; // m_x + 1
    double h_x; // шаг h
    double tau; // шаг tau
    double eta; // искусственная вязкость
}P_she;

void param_she (P_gas & gas, P_she & shem);

#endif// HEADER_CLASSES
