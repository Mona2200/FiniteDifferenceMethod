#include <iostream>
#include "Header.h"
using namespace std;

double* di, * au1, * au2, * al1, * al2, * f, * x;
int N;

int nX, nY, gX, gY, UNumber, FNumber, shift1, shift2, dshift;
double lambda, gamma, beta;
double* X, * Y;
int* boundary;

void InitMatrix()
{
   shift1 = gX, shift2 = nX;

   N = gX * (nY - gY) + nX * gY;

   di = new double[N];
   au1 = new double[N];
   au2 = new double[N];
   al1 = new double[N];
   al2 = new double[N];
   f = new double[N];
   x = new double[N];

   for (int i = 0; i < N; i++)
   {
      di[i] = 1;
      x[i] = 0;
      au1[i] = 0;
      au2[i] = 0;
      al1[i] = 0;
      al2[i] = 0;
   }

   InitRepository(UNumber, FNumber);
   CreateMatrix();
}

void CreateMatrix()
{
   int index = 0; //: ������ ����
   double hx1, hx2, hy1, hy2, ubeta; //: h-�� ���������� ����������
   // ������ ����� ������� "�"
   f[0] = U(X[0], Y[0]); // ����� ������ ����
   index++;
   hy1 = abs(Y[0] - Y[1]);
   for (int i = 1; i < gX - 1; i++, index++) {
      f[i] = U(X[i], Y[0]);
      if (boundary[0] == 3) {
         ubeta = lambda * rightDerivativeByY(X[i], Y[0]) / beta + U(X[i], Y[0]);
         di[index] = lambda / hy1 + beta;
         au1[index] = -lambda / hy1;
         f[index] = -lambda * rightDerivativeByY(X[i], Y[0]) + beta * (U(X[i], Y[0]) - ubeta) + beta * ubeta;
      }
   }
   f[index] = U(X[gX - 1], Y[0]); // ������ ������ ����
   index++;
   // �� ����� ����� ������ � ������
   for (int i = 1; i < gY - 1; i++, index++) {
      hy1 = abs(Y[i] - Y[i - 1]);
      hy2 = abs(Y[i + 1] - Y[i]);
      f[index] = U(X[0], Y[i]);
      if (boundary[1] == 3) { // ����� �������
         ubeta = lambda * rightDerivativeByX(X[0], Y[i]) / beta + U(X[0], Y[i]);
         hx1 = abs(X[0] - X[1]);
         di[index] = lambda / hx1 + beta;
         au2[index] = -lambda / hx1;
         f[index] = -lambda * rightDerivativeByX(X[0], Y[i]) + beta * (U(X[0], Y[i]) - ubeta) + beta * ubeta;
      }
      index++;
      for (int j = 1; j < gX - 1; j++, index++) {
         hx1 = abs(X[j] - X[j - 1]);
         hx2 = abs(X[j + 1] - X[j]);
         f[index] = F(X[j], Y[i]);
         al1[index - 1] = -2 * lambda / (hx1 * (hx2 + hx1));
         al2[index - shift1] = -2 * lambda / (hy1 * (hy2 + hy1));
         au1[index] = -2 * lambda / (hx2 * (hx2 + hx1));
         au2[index] = -2 * lambda / (hy2 * (hy2 + hy1));
         di[index] = lambda * (2 / (hx1 * hx2) + 2 / (hy1 * hy2)) + gamma;
      }
      f[index] = U(X[gX - 1], Y[i]);
      if (boundary[2] == 3) { // ������ �������
         ubeta = lambda * leftDerivativeByX(X[gX - 1], Y[i]) / beta + U(X[gX - 1], Y[i]);
         hx1 = abs(X[gX] - X[gX - 1]);
         di[index] = lambda / hx1 + beta;
         al2[index - shift1] = -lambda / hx1;
         f[index] = lambda * leftDerivativeByX(X[gX - 1], Y[i]) + beta * (U(X[gX - 1], Y[i]) - ubeta) + beta * ubeta;
      }
   }
   dshift = index;
   // ����� ������� � ������
   hy1 = abs(Y[gY - 1] - Y[gY - 2]);
   hy2 = abs(Y[gY] - Y[gY - 1]);
   f[index] = U(X[0], Y[gY - 1]);
   if (boundary[1] == 3) { // ����� �������
      ubeta = lambda * rightDerivativeByX(X[0], Y[gY - 1]) / beta + U(X[0], Y[gY - 1]);
      hx1 = abs(X[0] - X[1]);
      di[index] = lambda / hx1 + beta;
      au2[index] = -lambda / hx1;
      f[index] = -lambda * rightDerivativeByX(X[0], Y[gY - 1]) + beta * (U(X[0], Y[gY - 1]) - ubeta) + beta * ubeta;
   }
   index++;
   for (int i = 1; i < gX; i++, index++) {
      hx1 = abs(X[i] - X[i - 1]);
      hx2 = abs(X[i + 1] - X[i]);
      f[index] = F(X[i], Y[gY - 1]);
      al1[index - 1] = -2 * lambda / (hx1 * (hx2 + hx1));
      al2[index - shift1] = -2 * lambda / (hy1 * (hy2 + hy1));
      au1[index] = -2 * lambda / (hx2 * (hx2 + hx1));
      au2[index] = -2 * lambda / (hy2 * (hy2 + hy1));
      di[index] = lambda * (2 / (hx1 * hx2) + 2 / (hy1 * hy2)) + gamma;
   }
   for (int i = gX; i < nX - 1; i++, index++) {
      f[index] = U(X[i], Y[gY - 1]);
      if (boundary[3] == 3) {
         ubeta = lambda * rightDerivativeByY(X[i], Y[gY - 1]) / beta + U(X[i], Y[gY - 1]);
         di[index] = lambda / hy2 + beta;
         au1[index] = -lambda / hy2;
         f[index] = -lambda * rightDerivativeByY(X[i], Y[gY - 1]) + beta * (U(X[i], Y[gY - 1]) - ubeta) + beta * ubeta;
      }
   }
   f[index] = U(X[nX - 1], Y[gY - 1]);
   index++;
   // ������
   for (int i = gY; i < nY - 1; i++) {
      hy1 = abs(Y[i] - Y[i - 1]);
      hy2 = abs(Y[i + 1] - Y[i]);
      f[index] = U(X[0], Y[i]);
      if (boundary[1] == 3) { // ����� �������
         ubeta = lambda * rightDerivativeByX(X[0], Y[i]) / beta + U(X[0], Y[i]);
         hx1 = abs(X[0] - X[1]);
         di[index] = lambda / hx1 + beta;
         au2[index] = -lambda / hx1;
         f[index] = -lambda * rightDerivativeByX(X[0], Y[i]) + beta * (U(X[0], Y[i]) - ubeta) + beta * ubeta;
      }
      index++;
      for (int j = 1; j < nX - 1; j++, index++) {
         hx1 = abs(X[j] - X[j - 1]);
         hx2 = abs(X[j + 1] - X[j]);
         f[index] = F(X[j], Y[i]);
         al1[index - 1] = -2 * lambda / (hx1 * (hx2 + hx1));
         al2[index - shift1] = -2 * lambda / (hy1 * (hy2 + hy1));
         au1[index] = -2 * lambda / (hx2 * (hx2 + hx1));
         au2[index] = -2 * lambda / (hy2 * (hy2 + hy1));
         di[index] = lambda * (2 / (hx1 * hx2) + 2 / (hy1 * hy2)) + gamma;
      }
      f[index] = U(X[nX - 1], Y[i]);
      if (boundary[4] == 3) { // ������ �������
         ubeta = lambda * leftDerivativeByX(X[nX - 1], Y[i]) / beta + U(X[nX - 1], Y[i]);
         hx1 = abs(X[nX - 1] - X[nX - 2]);
         di[index] = lambda / hx1 + beta;
         al2[index - shift1] = -lambda / hx1;
         f[index] = lambda * leftDerivativeByX(X[nX - 1], Y[i]) + beta * (U(X[nX - 1], Y[i]) - ubeta) + beta * ubeta;
      }
      index++;
   }
   // �������� ������
   f[index] = U(X[0], Y[nY - 1]);
   index++;
   hy1 = abs(Y[nY - 1] - Y[nY - 2]);
   for (int i = 1; i < nX - 1; i++, index++) {
      f[index] = U(X[i], Y[nY - 1]);
      if (boundary[5] == 3) {
         ubeta = lambda * leftDerivativeByY(X[i], Y[nY - 1]) / beta + U(X[i], Y[nY - 1]);
         di[index] = lambda / hy1 + beta;
         al1[index - 1] = -lambda / hy1;
         f[index] = lambda * leftDerivativeByY(X[i], Y[nY - 1]) + beta * (U(X[i], Y[nY - 1]) - ubeta) + beta * ubeta;
      }
   }
   f[index] = U(X[nX - 1], Y[nY - 1]);


   //int index = 0;
   //double hx1, hx2, hY_1, hY_2, ubeta;

   //f[0] = U(X[0], Y[0]); //����� ������ �����
   //index++;

   //// ������ ����� ������� "�"
   //hY_1 = abs(Y[0] - Y[1]);
   //for (int i = 1; i < gX - 1; i++, index++) {
   //   f[index] = U(X[i], Y[0]);
   //   if (boundary[0] == 3) {
   //      ubeta = lambda * rightDerivativeByY(X[i], Y[0]) / beta + U(X[i], Y[0]);
   //      di[index] = lambda / hY_1 + beta;
   //      au1[index] = -lambda / hY_1;
   //      f[index] = -lambda * rightDerivativeByY(X[i], Y[0]) + beta * (U(X[i], Y[0]) - ubeta) + beta * ubeta;
   //   }
   //}

   //f[index++] = U(X[gX - 1], Y[0]); //������ ������ �����

   //// �� ����� ����� ������ � ������
   //for (int i = 1; i < gY - 1; i++, index++) {
   //   hY_1 = abs(Y[i] - Y[i - 1]);
   //   hY_2 = abs(Y[i + 1] - Y[i]);
   //   f[index] = U(X[0], Y[i]);
   //   if (boundary[1] == 3) { // ����� �������
   //      ubeta = lambda * rightDerivativeByX(X[0], Y[i]) / beta + U(X[0], Y[i]);
   //      hx1 = abs(X[0] - X[1]);
   //      di[index] = lambda / hx1 + beta;
   //      au2[index] = -lambda / hx1;
   //      f[index] = -lambda * rightDerivativeByX(X[0], Y[i]) + beta * (U(X[0], Y[i]) - ubeta) + beta * ubeta;
   //   }
   //   index++;

   //   for (int j = 1; j < gX - 1; j++, index++) {
   //      hx1 = abs(X[j] - X[j - 1]);
   //      hx2 = abs(X[j + 1] - X[j]);
   //      f[index] = F(X[j], Y[i]);
   //      al1[index - 1] = -2 * lambda / (hx1 * (hx2 + hx1));
   //      al2[index - shift1] = -2 * lambda / (hY_1 * (hY_2 + hY_1));
   //      au1[index] = -2 * lambda / (hx2 * (hx2 + hx1));
   //      au2[index] = -2 * lambda / (hY_2 * (hY_2 + hY_1));
   //      di[index] = lambda * (2 / (hx1 * hx2) + 2 / (hY_1 * hY_2)) + gamma;
   //   }

   //   f[index] = U(X[gX - 1], Y[i]);
   //   if (boundary[2] == 3) { // ������ �������
   //      ubeta = lambda * leftDerivativeByX(X[gX - 1], Y[i]) / beta + U(X[gX - 1], Y[i]);
   //      hx1 = abs(X[gX] - X[gX - 1]);
   //      di[index] = lambda / hx1 + beta;
   //      al2[index - shift1] = -lambda / hx1;
   //      f[index] = lambda * leftDerivativeByX(X[gX - 1], Y[i]) + beta * (U(X[gX - 1], Y[i]) - ubeta) + beta * ubeta;
   //   }
   //}

   //dshift = index;

   //// ����� ������� � ������ 
   //hY_1 = abs(Y[gY - 1] - Y[gY - 2]);
   //hY_2 = abs(Y[gY] - Y[gY - 1]);
   //f[index] = U(X[0], Y[gY - 1]);
   //if (boundary[1] == 3) { // ����� �������
   //   ubeta = lambda * rightDerivativeByX(X[0], Y[gY - 1]) / beta + U(X[0], Y[gY - 1]);
   //   hx1 = abs(X[0] - X[1]);
   //   di[index] = lambda / hx1 + beta;
   //   au2[index] = -lambda / hx1;
   //   f[index] = -lambda * rightDerivativeByX(X[0], Y[gY - 1]) + beta * (U(X[0], Y[gY - 1]) - ubeta) + beta * ubeta;
   //}
   //index++;

   //for (int i = 1; i < gX; i++, index++) {
   //   hx1 = abs(X[i] - X[i - 1]);
   //   hx2 = abs(X[i + 1] - X[i]);
   //   f[index] = F(X[i], Y[gY - 1]);
   //   al1[index - 1] = -2 * lambda / (hx1 * (hx2 + hx1));
   //   al2[index - shift1] = -2 * lambda / (hY_1 * (hY_2 + hY_1));
   //   au1[index] = -2 * lambda / (hx2 * (hx2 + hx1));
   //   au2[index] = -2 * lambda / (hY_2 * (hY_2 + hY_1));
   //   di[index] = lambda * (2 / (hx1 * hx2) + 2 / (hY_1 * hY_2)) + gamma;
   //}

   //for (int i = gX; i < nX - 1; i++, index++) {
   //   f[index] = U(X[i], Y[gY - 1]);
   //   if (boundary[3] == 3) {
   //      ubeta = lambda * rightDerivativeByY(X[i], Y[gY - 1]) / beta + U(X[i], Y[gY - 1]);
   //      di[index] = lambda / hY_2 + beta;
   //      au1[index] = -lambda / hY_2;
   //      f[index] = -lambda * rightDerivativeByY(X[i], Y[gY - 1]) + beta * (U(X[i], Y[gY - 1]) - ubeta) + beta * ubeta;
   //   }
   //}
   //f[index] = U(X[nX - 1], Y[gY - 1]);
   //index++;

   //// ������ 
   //for (int i = gY; i < nY - 1; i++, index++) {
   //   hY_1 = abs(Y[i] - Y[i - 1]);
   //   hY_2 = abs(Y[i + 1] - Y[i]);
   //   f[index] = U(X[0], Y[i]);
   //   if (boundary[1] == 3) { // ����� �������
   //      ubeta = lambda * rightDerivativeByX(X[0], Y[i]) / beta + U(X[0], Y[i]);
   //      hx1 = abs(X[0] - X[1]);
   //      di[index] = lambda / hx1 + beta;
   //      au2[index] = -lambda / hx1;
   //      f[index] = -lambda * rightDerivativeByX(X[0], Y[i]) + beta * (U(X[0], Y[i]) - ubeta) + beta * ubeta;
   //   }
   //   index++;

   //   for (int j = 1; j < nX - 1; j++, index++) {
   //      hx1 = abs(X[j] - X[j - 1]);
   //      hx2 = abs(X[j + 1] - X[j]);
   //      f[index] = F(X[j], Y[i]);
   //      al1[index - 1] = -2 * lambda / (hx1 * (hx2 + hx1));
   //      al2[index - shift1] = -2 * lambda / (hY_1 * (hY_2 + hY_1));
   //      au1[index] = -2 * lambda / (hx2 * (hx2 + hx1));
   //      au2[index] = -2 * lambda / (hY_2 * (hY_2 + hY_1));
   //      di[index] = lambda * (2 / (hx1 * hx2) + 2 / (hY_1 * hY_2)) + gamma;
   //   }

   //   f[index] = U(X[nX - 1], Y[i]);
   //   if (boundary[4] == 3) { // ������ �������
   //      ubeta = lambda * leftDerivativeByX(X[nX - 1], Y[i]) / beta + U(X[nX - 1], Y[i]);
   //      hx1 = abs(X[nX - 1] - X[nX - 2]);
   //      di[index] = lambda / hx1 + beta;
   //      al2[index - shift1] = -lambda / hx1;
   //      f[index] = lambda * leftDerivativeByX(X[nX - 1], Y[i]) + beta * (U(X[nX - 1], Y[i]) - ubeta) + beta * ubeta;
   //   }
   //}

   //// �������� ������
   //f[index] = U(X[0], Y[nY - 1]);
   //index++;
   //hY_1 = abs(Y[nY - 1] - Y[nY - 2]);
   //for (int i = 1; i < nX - 1; i++, index++) {
   //   f[index] = U(X[i], Y[nY - 1]);
   //   if (boundary[5] == 3) {
   //      ubeta = lambda * leftDerivativeByY(X[i], Y[nY - 1]) / beta + U(X[i], Y[nY - 1]);
   //      di[index] = lambda / hY_1 + beta;
   //      al1[index - 1] = -lambda / hY_1;
   //      f[index] = lambda * leftDerivativeByY(X[i], Y[nY - 1]) + beta * (U(X[i], Y[nY - 1]) - ubeta) + beta * ubeta;
   //   }
   //}
   //f[index] = U(X[nX - 1], Y[nY - 1]);
   //index++;
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

