#ifndef CLASSES_H
#define CLASSES_H
#endif
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include<sys/time.h>
#include<omp.h>
#include "vector"
using namespace std;

#define PI 3.14159265359
#define d0 10.0 //background density

class matrix
{ 
  int dim;
  double **value;
  double h;
public:
    matrix(int n,double hi);
    ~matrix();
    void display();
    void Error_Message(const matrix &b);
    double Error(const matrix& b);

    void SOR_smoothing(const matrix& rho,double omega,int steps);
    void SOR_smoothing(const matrix& rho, const matrix& ans, double omega, double  err_threshold, int num_thread);
    double averaging(int i,int j);
    matrix Restriction();
    double insertion(int i,int j,int dim_in);
    matrix Interpolation(int);
    matrix Residual(const matrix& rho);

    void init_constant_dens();
    void init_sin_dens();

    double get_h();
    double get_dim();
    void input_answer(int i,int j,double ans);

    matrix operator+(const matrix&);
};
