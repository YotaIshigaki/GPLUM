#include "Log.h"

int DebugLog::verbose;
int DebugLog::particle_dump;

void dump_typerange(std::vector<TypeRange>& typerangearray)
{
  std::cout << "typerange" ;
  for(size_t s=0;s<typerangearray.size();s++){
    std::cout << " " << typerangearray[s].begin << "-" << typerangearray[s].end;
  }
  std::cout << std::endl;
  std::cout << "lj" ;
  for(size_t s=0;s<typerangearray.size();s++){
    std::cout << " " << typerangearray[s].lj.begin << "-" << typerangearray[s].lj.end;
  }
  std::cout << std::endl;
  std::cout << "ljcoulomb" ;
  for(size_t s=0;s<typerangearray.size();s++){
    std::cout << " " << typerangearray[s].ljcoulomb.begin << "-" << typerangearray[s].ljcoulomb.end;
  }
  std::cout << std::endl;
  std::cout << "coulomb" ;
  for(size_t s=0;s<typerangearray.size();s++){
    std::cout << " " << typerangearray[s].coulomb.begin << "-" << typerangearray[s].coulomb.end;
  }
  std::cout << std::endl;
}

template<typename PA>
void dump_particle(PA& particlearray, 
          std::vector<TypeRange>& typerangearray, 
          std::vector<int>& setid,
          double energy)
{
  std::cout << " position ";
  for(size_t s=0;s<typerangearray.size();s++){
    if(typerangearray[s].end>typerangearray[s].begin){
      std::cout << " cell " << setid[s] << " ";
      for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
        std::cout << "index " << i << " ";
        std::cout << "atomID " << getatomid(particlearray,i) << " ";
        std::cout << getpos(particlearray,i);
      }
    }
  }
  std::cout << std::endl;
  std::cout << "velocity ";
  for(size_t s=0;s<typerangearray.size();s++){
    if(typerangearray[s].end>typerangearray[s].begin){
      std::cout << " cell " << setid[s] << " ";
      for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
        std::cout << getvelocity(particlearray,i);
      }
    }
  }
  std::cout << std::endl;
  std::cout << "force ";
  for(size_t s=0;s<typerangearray.size();s++){
    if(typerangearray[s].end>typerangearray[s].begin){
      std::cout << " cell " << setid[s] << " ";
      for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
        std::cout << getforce(particlearray,i);
      }
    }
  }
  std::cout << std::endl;
  std::cout << "energy " << energy  << std::endl;
}
template
void dump_particle(ParticleArray& particlearray, 
                   std::vector<TypeRange>& typerangearray, 
                   std::vector<int>& setid,
                   double energy);
template
void dump_particle(CombinedParticleArray& particlearray, 
                   std::vector<TypeRange>& typerangearray, 
                   std::vector<int>& setid,
                   double energy);

template<typename PA>
void dump_atomid(const PA& particlearray, 
                 const std::vector<TypeRange>& typerangearray, 
                 const std::vector<int>& setid)
{
  std::cout << " atomID ";
  for(size_t s=0;s<typerangearray.size();s++){
    if(typerangearray[s].end>typerangearray[s].begin){
      std::cout << " cell " << setid[s] << " ";
      for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
        std::cout << " " << getatomid(particlearray,i);
      }
    }
  }
  std::cout << std::endl;
}
template
void dump_atomid(const ParticleArray& particlearray, 
                 const std::vector<TypeRange>& typerangearray, 
                 const std::vector<int>& setid);
template
void dump_atomid(const CombinedParticleArray& particlearray, 
                 const std::vector<TypeRange>& typerangearray, 
                 const std::vector<int>& setid);
template
void dump_atomid(const GhostParticleArray& particlearray, 
                 const std::vector<TypeRange>& typerangearray, 
                 const std::vector<int>& setid);

void print_cb(CovalentBondInfo::Bond bond)
{
  std::cout << " (" << bond.id_of_atom[0] << "," << bond.id_of_atom[1] << ")" <<  bond.typeofbond << ">>=" <<  bond.shake;
}
void print_cb(CovalentBondInfo::Angle angle)
{
  std::cout << " (" << angle.id_of_atom[0] << "," << angle.id_of_atom[1] << "," << angle.id_of_atom[2] << ")" <<  angle.typeofangle;
}
void print_cb(CovalentBondInfo::Torsion torsion)
{
  std::cout << " (" << torsion.id_of_atom[0] << "," << torsion.id_of_atom[1] << "," << torsion.id_of_atom[2] << "," << torsion.id_of_atom[3] << ")" <<  torsion.typeoftorsion;
}
void print_cb(CovalentBondInfo::Improper improper)
{
  std::cout << " (" << improper.id_of_atom[0] << "," << improper.id_of_atom[1] << "," << improper.id_of_atom[2] << "," << improper.id_of_atom[3] << ")" <<  improper.typeofimproper;
}

template<typename CB>
void print_bondarray(std::vector<CB> cbarray)
{
  std::cout << cbarray.size();
  for(size_t b=0;b<cbarray.size();b++){
    print_cb(cbarray[b]);
  }
}

void dump_bondlistarray(std::vector<CovalentBondInfo::BondList> bondlistarray)
{
  for(size_t s=0;s<bondlistarray.size();s++){
    if(bondlistarray[s].BondArray.size()>0){
      std::cout << "Bond in set " << s << " ";
      print_bondarray(bondlistarray[s].BondArray);
      std::cout << std::endl;
    }
    if(bondlistarray[s].AngleArray.size()>0){
      std::cout << "Angle in set " << s << " ";
      print_bondarray(bondlistarray[s].AngleArray);
      std::cout << std::endl;
    }
    if(bondlistarray[s].TorsionArray.size()>0){
      std::cout << "Torsion in set " << s << " ";
      print_bondarray(bondlistarray[s].TorsionArray);
      std::cout << std::endl;
    }
    if(bondlistarray[s].ImproperArray.size()>0){
      std::cout << "Improper in set " << s << " ";
      print_bondarray(bondlistarray[s].ImproperArray);
      std::cout << std::endl;
    }
  }
}

template<class PA>
void dump_shakelist(ShakeList& shakelist, PA& particlearray)
{
  std::cout << "Number of Shake " << shakelist.size() << std::endl;
  for(ShakeList::iterator it=shakelist.begin(); it != shakelist.end(); it++){
    int nh1 = it->second.nh1;
    std::cout << "Shake nh1= " << nh1 << "  it->first= " << it->first << " atomid=" << getatomid(particlearray,it->first) << std::endl;
    for(int n1=0; n1<nh1; n1++){
      std::cout << "     " << "n1= " << "it->second.h1[n1]= " << it->second.h1[n1] << " atomid=" << getatomid(particlearray,it->second.h1[n1]) << std::endl;
    }
  }
}
template
void dump_shakelist(ShakeList& shakelist, ParticleArray& particlearray);

template<class PA>
void dump_shakelist_distance(ShakeList& shakelist, PA& particlearray)
{
  std::cout << "Number of Shake " << shakelist.size() << std::endl;
  for(ShakeList::iterator it=shakelist.begin(); it != shakelist.end(); it++){
    int nh1 = it->second.nh1;
    std::cout << "Shake nh1= " << nh1 << "  it->first= " << it->first << " atomid=" << getatomid(particlearray,it->first) << " type=" << getatomtype(particlearray,it->first) << std::endl;
    for(int n1=0; n1<nh1; n1++){
      Position hah1(getpos(particlearray,it->first) - getpos(particlearray,it->second.h1[n1]));
      std::cout << "     " << "n1= " << "it->second.h1[n1]= " << it->second.h1[n1] << " atomid=" << getatomid(particlearray,it->second.h1[n1]) << " type=" << getatomtype(particlearray,it->second.h1[n1]) << " distance=" << sqrt(hah1.norm2()) << std::endl;
    }
  }
}

#ifdef OLDPARTICLE
template
void dump_shakelist_distance(ShakeList& shakelist, ParticleArray& particlearray);
#else
template
void dump_shakelist_distance(ShakeList& shakelist, CombinedParticleArray& particlearray);
#endif
