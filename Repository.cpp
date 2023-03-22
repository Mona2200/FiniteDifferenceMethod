#include <iostream>
using namespace std;

int uN, fN;

void InitRepository(int _UNumber, int _FNumber)
{
   uN = _UNumber; fN = _FNumber;
}

double U(double x, double y)
{
   switch (uN)
   {
   case 1:
      return x + y;

   case 2:
      return x * x + y * x;

   case 3:
      return x * x * x + y * y * y;

   case 4:
      return cos(2 * x + 2 * y);

   case 5:
      return sin(x + y);
   }
   return 0;
}

double F(double x, double y)
{
   switch (fN)
   {
   case 1:
      return 2 * x + 2 * y;

   case 2:
      return 2 * x * x + 2 * y * y - 4;

   case 3:
      return 2 * x * x * x + 2 * y * y * y - 6 * x - 6 * y;

   case 4: /// polynom_3
      return 9 * cos(2 * x + 2 * y);

   case 5: /// polynom_4
      return -48 * x * x - 96 * y * y + 6 * x * x * x * x + 12 * y * y * y * y;

   case 6: /// not_polynom
      return 3 * sin(x + y);
   }
   return 0;
}