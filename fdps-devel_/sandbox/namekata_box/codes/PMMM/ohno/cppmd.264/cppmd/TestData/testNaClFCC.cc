#include "NaClFCC.h"

int main(int argc, char **argv) {
  int natom;
  SpaceVector<int> latticeNum;

  if ( argc == 2 ) {
    latticeNum = SpaceVector<int>(atoi(argv[1]));
  } else if (argc == 4) {
    latticeNum = SpaceVector<int>(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
  } else {
    latticeNum = SpaceVector<int>(1);
  }

  NaClFCC naclfcc(latticeNum);

  //  natom = arfcc.setFccLattice();
  
  std::cerr << "#atoms " << naclfcc.natom << std::endl;
  std::cerr << "Box " << naclfcc.side << std::endl;
  /*
  vector<Particle>::iterator it = pa.begin();
  while(it != pa.end()) {
    cerr << it->position << endl;
    ++it;
  }
  */
  naclfcc.writePDB();
  return 0;
}
