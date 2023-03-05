#include <fstream>
#include <iostream>
#include "Solving.cpp"
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


   boundary = new int[5];
   ifstream readerCond("BoundaryConditions.txt");

   for (int i = 0; i < 5; i++)
      readerCond >> boundary[i];

   readerCond.close();

   ifstream readerFunc("Functions.txt");
   readerFunc >> UNumber >> FNumber >> lambda >> gamma;
   readerFunc.close();

   InitSolving();
}