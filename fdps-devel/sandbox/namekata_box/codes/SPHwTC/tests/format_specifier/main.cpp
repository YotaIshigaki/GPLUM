#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

int main(int argc, char* argv[]) {
   //* Local variablesa
   const int ndump=5;
   std::ostringstream string_stream;
   std::string filename;

   string_stream << "hydro" << std::setfill('0') << std::setw(5) << ndump;
   filename = string_stream.str();
   std::cout << filename << std::endl; 

   return 0;

}
