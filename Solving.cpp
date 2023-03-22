#include <iostream>
#include "Header.h"
using namespace std;

double omega, iter, eps;  /// Параметр релаксации

void InitSolving(int _iter, double _eps, double _omega) {
   omega = _omega;
   iter = _iter;
   eps = _eps;

   Solve();
}

//* Решение СЛАУ
void Solve() {
   double sum, Nev = 0, norm_f;
   int Iter = 0;
   norm_f = Norm(f);

   do {
      Nev = 0;
      for (int i = 0; i < N; i++) {
         sum = di[i] * x[i];

         if (i < N - 1)
            sum += au1[i] * x[i + 1];
         if (i >= 2)
            sum += al1[i - 1] * x[i - 1];

         if (i < dshift)
            sum += au2[i] * x[shift1 + i];
         else if (i < N - shift2)
            sum += au2[i] * x[shift2 + i];

         if (i >= shift1 + dshift)
            sum += al2[i - shift1] * x[i - shift2];
         else if (i >= shift1)
            sum += al2[i - shift1] * x[i - shift1];

         Nev += (f[i] - sum) * (f[i] - sum);
         x[i] += omega / di[i] * (f[i] - sum);
      }

      Nev = sqrt(Nev) / norm_f; /// Относительная невязка
      Iter++;
      
   } while (Nev > eps &&
      Iter <= iter);

   cout << "Iter: " << Iter << " Nev: " << Nev;
}

double OtnNev(double* _f, double* Ax)
{
   double* f_ = new double[N];
   for (int i = 0; i < N; i++)
   {
      f_[i] = _f[i] - Ax[i];
   }
   return Norm(f_) / Norm(_f);
}

double Norm(double* vector)
{
   double sum = 0;
   for (int i = 0; i < N; i++)
   {
      sum += pow(vector[i], 2);
   }
   return sqrt(sum);
}