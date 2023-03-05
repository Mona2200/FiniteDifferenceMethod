#include <iostream>
#include "Repository.cpp"
using namespace std;

double* di, * au1, * au2, * al1, * al2, * f;
int N;

int nX, nY, gX, gY, UNumber, FNumber, al_shift, au_shift, di_shift;
double lambda, gamma, beta;
double* X, * Y;
int* boundary;

const double h = 1e-9;

void CreateMatrix();
void InitSolving(int _nX, int _nY, int _gX, int _gY, double* _X, double* _Y, int* _boundary, int _UNumber, int _FNumber, double _lambda, double _gamma);
double rightDerivativeByY(double x, double y);
double rightDerivativeByX(double x, double y);
double leftDerivativeByY(double x, double y);
double leftDerivativeByX(double x, double y);

void InitSolving()
{
   al_shift = gX, au_shift = nX;

   N = gX * (nY - gY) + nX * gY;

   di = new double[N];
   au1 = new double[N];
   au2 = new double[N];
   al1 = new double[N];
   al2 = new double[N];
   f = new double[N];

   for (int i = 0; i < N; i++)
      di[i] = 1;

   InitRepository(UNumber, FNumber);
   CreateMatrix();
}

void CreateMatrix() {
   int index = 0;
   double hx1, hx2, hY_1, hY_2;

   f[0] = U(X[0], Y[0]); //левая нижняя точка

   // Нижняя линия области "Г"
   hY_1 = abs(Y[0] - Y[1]);
   for (int i = 1; i < gX - 1; i++, index++) {
      f[index] = U(X[i], Y[i]);
      if (boundary[0] == 3) {
         di[index] = lambda / hY_1 + beta;
         au1[index] = -lambda / hY_1;
         f[index] = -lambda * rightDerivativeByY(X[i], Y[i]) + beta * U(X[i], Y[i]);
      }
   }

   f[index++] = U(X[gX - 1], Y[0]); //правая нижняя точка

   // До линии между шапкой и ножкой
   for (int i = 1; i < gY - 1; i++, index++) {
      hY_1 = abs(Y[i] - Y[i - 1]);
      hY_2 = abs(Y[i + 1] - Y[i]);
      f[index] = U(X[0], Y[i]);
      if (boundary[1] == 3) {
         hx1 = abs(X[0] - X[1]);
         di[index] = -lambda / hx1;
         au1[index] = lambda / hx1;
      }
      index++;

      for (int j = 1; j < gX - 1; j++, index++) {
         hx1 = abs(X[j] - X[j - 1]);
         hx2 = abs(X[j + 1] - X[j]);
         f[index] = F(X[j], Y[i]);
         al1[index - 1] = -2 * lambda / (hx1 * (hx2 + hx1));
         al2[index - al_shift] = -2 * lambda / (hY_1 * (hY_2 + hY_1));
         au1[index] = -2 * lambda / (hx2 * (hx2 + hx1));
         au2[index] = -2 * lambda / (hY_2 * (hY_2 + hY_1));
         di[index] = lambda * (2 / (hx1 * hx2) + 2 / (hY_1 * hY_2)) + gamma;
      }

      f[index] = U(X[gX - 1], Y[i]);
      if (boundary[2] == 3) {
         hx1 = abs(X[gX - 1] - X[gX - 2]);
         di[index] = lambda / hx1;
         al1[index - 1] = -lambda / hx1;
      }
   }

   di_shift = index;

   // Между шляпкой и ножкой 
   hY_1 = abs(Y[gY - 1] - Y[gY - 2]);
   hY_2 = abs(Y[gY] - Y[gY - 1]);
   f[index] = U(X[0], Y[gY - 1]);
   if (boundary[3] == 3) {
      hx1 = abs(X[0] - X[1]);
      di[index] = -lambda / hx1;
      au1[index] = lambda / hx1;
   }
   index++;

   for (int i = 1; i < gX; i++, index++) {
      hx1 = abs(X[i] - X[i - 1]);
      hx2 = abs(X[i + 1] - X[i]);
      f[index] = F(X[i], Y[gY - 1]);
      al1[index - 1] = -2 * lambda / (hx1 * (hx2 + hx1));
      al2[index - al_shift] = -2 * lambda / (hY_1 * (hY_2 + hY_1));
      au1[index] = -2 * lambda / (hx2 * (hx2 + hx1));
      au2[index] = -2 * lambda / (hY_2 * (hY_2 + hY_1));
      di[index] = lambda * (2 / (hx1 * hx2) + 2 / (hY_1 * hY_2)) + gamma;
   }

   for (int i = gX; i < nX - 1; i++, index++) {
      f[index] = U(X[i], Y[gY - 1]);
      if (boundary[4] == 3) {
         di[index] = -lambda / hY_2;
         au2[index] = lambda / hY_2;
      }
   }
   f[index] = U(X[nX - 1], Y[gY - 1]);
   index++;

   // Шляпка 
   for (int i = gY; i < nY - 1; i++) {
      hY_1 = abs(Y[i] - Y[i - 1]);
      hY_2 = abs(Y[i + 1] - Y[i]);
      f[index] = U(X[0], Y[i]);
      if (boundary[5] == 3) {
         hx1 = abs(X[0] - X[1]);
         di[index] = -lambda / hx1;
         au1[index] = lambda / hx1;
      }
      index++;

      for (int j = 1; j < nX - 1; j++, index++) {
         hx1 = abs(X[j] - X[j - 1]);
         hx2 = abs(X[j + 1] - X[j]);
         f[index] = F(X[j], Y[i]);
         al1[index - 1] = -2 * lambda / (hx1 * (hx2 + hx1));
         al2[index - al_shift] = -2 * lambda / (hY_1 * (hY_2 + hY_1));
         au1[index] = -2 * lambda / (hx2 * (hx2 + hx1));
         au2[index] = -2 * lambda / (hY_2 * (hY_2 + hY_1));
         di[index] = lambda * (2 / (hx1 * hx2) + 2 / (hY_1 * hY_2)) + gamma;
      }

      f[index] = U(X[nX - 1], Y[i]);
      if (boundary[6] == 3) {
         hx1 = abs(X[nX - 1] - X[nX - 2]);
         di[index] = lambda / hx1;
         al1[index - 1] = -lambda / hx1;
      }
      index++;
   }

   // Верхушка шляпки
   f[index] = U(X[0], Y[nY - 1]);
   index++;
   hY_1 = abs(Y[nY - 1] - Y[nY - 2]);
   for (int i = 1; i < nX - 1; i++, index++) {
      f[index] = U(X[i], Y[nY - 1]);
      if (boundary[7] == 3) {
         di[index] = lambda / hY_1;
         al2[index - al_shift] = -lambda / hY_1;
      }
   }
   f[index] = U(X[nX - 1], Y[nY - 1]);
}

double rightDerivativeByY(double x, double y)
{
   return (U(x, y + h) - U(x, y)) / h;
}

double rightDerivativeByX(double x, double y)
{
   return (U(x + h, y) - U(x, y)) / h;
}

double leftDerivativeByY(double x, double y)
{
   return (U(x, y) - U(x, y - h)) / h;
}

double leftDerivativeByX(double x, double y)
{
   return (U(x, y) - U(x - h, y)) / h;
}