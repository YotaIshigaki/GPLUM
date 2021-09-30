#include <iostream>
#include <Eigen/Dense>

int main(int argc, char* argv[]){
   using namespace Eigen;
   //* Local variables
   Matrix<double,3,3> mat;

   mat << 1, 2, 3,
          4, 5, 6,
          7, 8, 9;
   std::cout << mat << std::endl;

   mat = mat.Zero();
   std::cout << mat << std::endl;

   return 0;

}
