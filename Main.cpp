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
   InitSolving(1000, 1e-14);

   for (int i = 0; i < N; i++)
   {
      printf("%.14f\t", x[i]);
   }
      cout << endl;
}