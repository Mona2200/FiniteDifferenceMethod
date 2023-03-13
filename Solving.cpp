#include <iostream>
#include "Header.h"
using namespace std;

double omega, iter, eps;  /// Параметр релаксации

void InitSolving(int _iter, double _eps) {
   omega = 1;
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

         if (i < di_shift)
            sum += au2[i] * x[al_shift + i];
         else if (i < N - au_shift)
            sum += au2[i] * x[au_shift + i];

         if (i >= al_shift + di_shift)
            sum += al2[i - al_shift] * x[i - au_shift];
         else if (i >= al_shift)
            sum += al2[i - al_shift] * x[i - al_shift];

         Nev += (f[i] - sum) * (f[i] - sum);
         x[i] += omega / di[i] * (f[i] - sum);
      }

      Nev = sqrt(Nev) / norm_f; /// Относительная невязка
      Iter++;
      cout << "Iter: " << Iter << " Nev: " << Nev;
   } while (Nev > eps &&
      Iter <= iter);
}

double Norm(double* vector) {
   double norm = 0;
   for (int i = 0; i < N; i++)
      norm += vector[i] * vector[i];
   return sqrt(norm);
}