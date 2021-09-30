#include "CubicCell.h"

#include "TestParticle.h"

void
TestParticle::makeparticle(ParticleArray& particlearray,
                           std::vector<int>& particle_setid,
                           std::vector<TypeRange>& typerangearray,
                           int total_num_particle, int total_num_set, 
                           int total_num_unit, int unit_id,
                           PostProcess& postprocess)
{
  int num_particle = total_num_particle/total_num_unit;
  //  int num_set = total_num_set/total_num_unit;
  int npc = (total_num_particle+total_num_set-1)/total_num_set;
  int pnpc;
  if(!power_of_two(npc, pnpc)){
    pnpc++;
  }
  int sp[3];
  split_power_of_two(pnpc,sp);
  int div[3] = {(1<<sp[2]),(1<<sp[1]),(1<<sp[0])};
  //  int n=0;
  int aid0=num_particle*unit_id;
  size_t sindex = 0;
  for(int i=0;i<num_particle;i++){
    int sid = particle_setid[sindex];
    SpaceVector<int> cellpos = postprocess.cell_position[sid];
    SpaceVector<double> pos0 = postprocess.cellsize;
    pos0.x *= (double(cellpos.x));
    pos0.y *= (double(cellpos.y));
    pos0.z *= (double(cellpos.z));
    int ns = typerangearray[sindex].end-typerangearray[sindex].begin;
    int pindex = typerangearray[sindex].lj.end;
    int nsx = ns%(div[0]);
    int nsy = (ns/div[0])%(div[1]);
    int nsz = ns/(div[0]*div[1]);
    pos0.x += postprocess.cellsize.x*((double(nsx)+0.5)/double(div[0]));
    pos0.y += postprocess.cellsize.y*((double(nsy)+0.5)/double(div[1]));
    pos0.z += postprocess.cellsize.z*((double(nsz)+0.5)/double(div[2]));
    particlearray[pindex].position = pos0;
    particlearray[pindex].velocity(0.0,0.0,0.0);
    // override for test
    if(num_particle==1){
      particlearray[0].position.z = postprocess.cellsize.z*double(cellpos.z) - 0.98;
      particlearray[0].velocity.z = -0.1;
    }

    particlearray[pindex].force(0.0,0.0,0.0);
    particlearray[pindex].mass=3.32603;
    particlearray[pindex].inv_mass=1.0/particlearray[pindex].mass;
    particlearray[pindex].atomtype = 0;
    particlearray[pindex].atomid = aid0+i;
    typerangearray[sindex].end++;
    typerangearray[sindex].lj.end++;
    typerangearray[sindex].ljcoulomb.shift(1);
    typerangearray[sindex].coulomb.shift(1);
    sindex++;
    if(sindex>=particle_setid.size())sindex=0;
  }
  if(DebugLog::verbose>1){
    if(unit_id==0){
      std::cout << "0 cell 0 ";
      for(int i=typerangearray[0].begin;i<typerangearray[0].end;i++){
        std::cout << particlearray[i].position;
      }
      std::cout << std::endl;
    }
    std::cout << "generate particle " << num_particle << std::endl;
  }
}
