#include <cstdlib>

#include "Common.h"

template<typename PA>
void qpos_out(const PA &ptcl, const int num, const double len=1.0)
{
  double scale = 1.0/len;
  FILE *fp = fopen("qpos.dat", "w");
  assert(fp);
  fprintf(fp, "%d\n", num);
  for(int i=0; i<num; i++){
    fprintf(fp, "%A %A %A %A\n",
	    getcharge(ptcl,i),
	    getpos(ptcl,i).x*scale,
	    getpos(ptcl,i).y*scale,
	    getpos(ptcl,i).z*scale);
  }
  fclose(fp);
}

