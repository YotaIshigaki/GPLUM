/* Standard headers */
#include <iostream>

extern "C" void f_main_();

int c_extern;
long myVariable;
struct {float r, s;} com;
float single;

double d;
double D;
double darray[16];
double array2d[4][4];

struct point {
   double x,y,z;
};
struct point p;
struct point parray[4];

int main(int argc, char *argv[]) {

   //* Set global vars
   c_extern = 16;
   myVariable = 18;
   com.r = 4.0f;
   com.s = 6.0f;
   single = 1.0f;
   d = 1.0e0;
   D = 2.0e0;

   for (int i=0; i<16;i++) {
      darray[i] = (double)i;
   }
   for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
         array2d[i][j] = (double)(i*i+j);
         std::cout << i << " " << j << " " << array2d[i][j] << std::endl;
      }
   }

   p.x = 1.0e0;
   p.y = 2.0e0;
   p.z = 3.0e0;
   
   for (int i=0; i<4; i++) {
      parray[i].x = 1.0e0;
      parray[i].y = 2.0e0;
      parray[i].z = 3.0e0;
   }

   //* Call Fortran routine
   f_main_();

   return 0;
}

