#ifndef HEADER_H
#define HEADER_H

double F(double x, double y);
void InitRepository(int _UNumber, int _FNumber);
double U(double x, double y);
void CreateMatrix();
void InitMatrix();
double rightDerivativeByY(double x, double y);
double rightDerivativeByX(double x, double y);
double leftDerivativeByY(double x, double y);
double leftDerivativeByX(double x, double y);

void InitSolving(int _iter, double _eps);
void Solve();
double Norm(double* vector);

extern double* di, * au1, * au2, * al1, * al2, * f, * x;
extern int N;

extern int nX, nY, gX, gY, UNumber, FNumber, al_shift, au_shift, di_shift;
extern double lambda, gamma, beta;
extern double* X, * Y;
extern int* boundary;

const double h = 1e-9;

#endif
