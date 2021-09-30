#include "ArFCC.h"

int main(int argc, char **argv) {
  char aname[5]="AR  ";
  char rname[4]="AR ";
  int natom;
  SpaceVector<int> latticeNum;

  if ( argc == 2 ) {
    latticeNum = SpaceVector<int>(atoi(argv[1]));
  } else if (argc == 4) {
    latticeNum = SpaceVector<int>(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
  } else {
    latticeNum = SpaceVector<int>(1);
  }

  ArFCC arfcc(latticeNum);

  //  natom = arfcc.setFccLattice();
  
  std::cerr << "#atoms " << arfcc.natom << std::endl;
  std::cerr << "Box " << arfcc.side << std::endl;
  /*
  vector<Particle>::iterator it = pa.begin();
  while(it != pa.end()) {
    cerr << it->position << endl;
    ++it;
  }
  */
  arfcc.writePDB();
  return 0;
}
