#include <fstream>
#include <iostream>
#include "Header.h"
using namespace std;

void main()
{
   ifstream readerMesh("Mesh.txt");
   readerMesh >> nX >> nY;

   X = new double[nX];
   Y = new double[nY];

   for (int i = 0; i < nX; i++)
      readerMesh >> X[i];
   for (int i = 0; i < nY; i++)
      readerMesh >> Y[i];

   readerMesh >> gX >> gY;

   readerMesh.close();


   boundary = new int[8];
   ifstream readerCond("BoundaryConditions.txt");

   for (int i = 0; i < 8; i++)
      readerCond >> boundary[i];

   readerCond.close();

   ifstream readerFunc("Functions.txt");
   readerFunc >> UNumber >> FNumber >> lambda >> gamma >> beta;
   readerFunc.close();

   InitMatrix();
   InitSolving(1000, 1e-14, 1);

   double* _u = new double[N];

   ofstream writer("Result.txt");
   writer.setf(ios::scientific);

   int ind = 0;

   for (int i = 0; i < gY - 1; i++)
   {
      for (int j = 0; j < gX; j++)
      {
         _u[i] = U(X[j], Y[i]);
         writer << X[j] << "\t";
         writer << Y[i] << "\t";
         writer << x[ind++] << "\t";
         writer << U(X[j], Y[i]) << "\n";
         /*printf("%.14f\t", X[j]);
         printf("%.14f\t", Y[i]);
         printf("%.14f\t", x[ind++]);
         printf("%.14f\n", U(X[j], Y[i]));*/
      }
         
   }

   for (int i = gY - 1; i < nY; i++)
   {
      for (int j = 0; j < nX; j++)
      {
         _u[i] = U(X[j], Y[i]);
         writer << X[j] << "\t";
         writer << Y[i] << "\t";
         writer << x[ind++] << "\t";
         writer << U(X[j], Y[i]) << "\n";
         /*printf("%.14f\t", X[j]);
         printf("%.14f\t", Y[i]);
         printf("%.14f\t", x[ind++]);
         printf("%.14f\n", U(X[j], Y[i]));*/
      }

   }
   writer.close();

   cout.setf(ios::scientific);
   //cout << OtnNev(x, _u);
}