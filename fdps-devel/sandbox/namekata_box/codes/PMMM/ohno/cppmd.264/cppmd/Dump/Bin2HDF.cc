#include <mpi.h>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include "Bin2HDF.h"
#include "Converter.h"
#include "HDFSnapshotData.h"


void setParameter(HDFDump& hdfdump, double dt=1.0, int step=0, int rc=0, int pc=0,
		  int crdc=0, int rstc=0, AtomID sp=1)
{
  hdfdump.setParameter(dt,step,rc,pc,crdc,rstc,sp);
}

void setInfo(HDFDump& hdfdump, int rank, int num_rank, string* version, AtomID sp=1, HDFSnapshot::HDFSnapshotFile::SnapshotType sst=HDFSnapshot::HDFSnapshotFile::Single)
{
  hdfdump.setInfo(rank,num_rank,version,sp,sst);
}

void setIntegration(HDFDump& hdfdump, 
		    double LJEnergyCorrection,
			     double KineticEnergy,
			     double PotentialEnergy,
			     double TotalEnergy,
			     double Virial,
			     double EtaPosition,
			     double EtaVelocity,
			     double EtaForce,
			     double LogVPosition,
			     double LogVVelocity,
			     double LogVForce,
			     double Volume)
{
  hdfdump.setIntegration(LJEnergyCorrection, KineticEnergy, PotentialEnergy, TotalEnergy, Virial,
			 EtaPosition, EtaVelocity, EtaForce,
			 LogVPosition, LogVVelocity, LogVForce, Volume );
}

void setBoxDatatype(HDFDump& hdfdump, SpaceVector<double> boxsize, int type, SpaceVector<double> angle)
{
  hdfdump.setBoxDatatype(boxsize, type, angle);
}

void Bin2HDF::outputSingleHDF(Converter<CombinedParticleArray>& converter)
{
  converter.gather_size();
  converter.merge_bondlistarray();
  if(with_waterlist>0){
    converter.merge_waterlists_atomid();
  }
  CombinedParticleArray pa(0);
  converter.gather_particle(pa,converter.particlearray);
  std::vector<CovalentBondInfo::BondList> bl(0);
  converter.gather_bondlist(bl,converter.bondlistarray);
  WaterList wl;
  if(with_waterlist>0){
    converter.gather_waterlist(wl,converter.waterlist);
  }
  converter.merge_shakelists_atomid();
  ShakeList sl;
  converter.gather_shakelist(sl,converter.shakelist);

  if(comm_rank==0){
    HDFDump hdfdump(hdfname,HDFDUMP);
    double dt = 1.0;
    int step = converter.timesteps[0];
    int rc = converter.ri;
    int pc = converter.pi;
    int crdc = converter.dci;
    int rstc = converter.dri;
    AtomID size=converter.total_number_of_particle;
    std::string ver("alpha");
    double LJEnergyCorrection = converter.tljcec;
    double KineticEnergy = converter.kenergy;
    double PotentialEnergy = converter.penergy;
    double TotalEnergy = converter.total_penergy;
    double Virial = converter.virial;
    double EtaPosition = converter.eta_pos;
    double EtaVelocity = converter.eta_vel;
    double EtaForce = converter.eta_force;
    double LogVPosition = converter.logv_pos;
    double LogVVelocity = converter.logv_vel;
    double LogVForce = converter.logv_force;
    double Volume = converter.volume;
    //    SpaceVector<double> boxsize(64.0,64.0,64.0);
    int boxtype = 0;
    SpaceVector<double> angle(90.0,90.0,90.0);
    setParameter(hdfdump,dt,step,rc,pc,crdc,rstc,size);
    setInfo(hdfdump,0,1,&ver,size,snapshottype);
    setIntegration(hdfdump,
		   LJEnergyCorrection, KineticEnergy, PotentialEnergy, TotalEnergy, Virial,
		   EtaPosition, EtaVelocity, EtaForce,
		   LogVPosition, LogVVelocity, LogVForce, Volume );
    setBoxDatatype(hdfdump, converter.boxsize, boxtype, angle);
    std::vector<TypeRange> tr(1);
    tr[0].begin = 0;
    tr[0].end = pa.size();
    tr[0].lj.begin = 0;
    tr[0].lj.end = 0;
    tr[0].ljcoulomb.begin = 0;
    tr[0].ljcoulomb.end = tr[0].end;
    tr[0].coulomb.begin = tr[0].ljcoulomb.end;
    tr[0].coulomb.end = tr[0].coulomb.begin;
    hdfdump.setParticle(pa,tr);
    hdfdump.setBondlists(bl);
    if(with_waterlist>0){
      hdfdump.setWaterlist(wl);
    }
    hdfdump.setShakelist(sl);
    hdfdump.dump();
  }  
}

void Bin2HDF::dumpParallelHDF(Converter<CombinedParticleArray>& converter)
{
  std::string namebase;
  std::string sufix(".hdf");
  std::string::size_type posdot = hdfname.find_last_of('.');
  if(posdot!=std::string::npos){
    if(posdot<hdfname.length()){
      char c = hdfname.at(posdot+1);
      if((c=='H')||(c=='h')){
	namebase = hdfname.substr(0,posdot);
	sufix = hdfname.substr(posdot);
      }else{
	namebase = hdfname;
      }
    }else{
      namebase = hdfname.substr(0,posdot);
    }
  }else{
    namebase = hdfname;
  }
  for(int i=0;i<converter.number_of_localrestorefiles;i++){
    int index = converter.restorefileindexes[i];
    char name[256];
    snprintf(name,256,"%s.%05d",namebase.c_str(),index);
    std::string filename(name);
    filename.append(sufix);
    HDFDump hdfdump(filename,HDFDUMP);
    double dt = 1.0;
    int step = converter.timesteps[i];
    int rc = converter.ris[i];
    int pc = converter.pis[i];
    int crdc = converter.dcis[i];
    int rstc = converter.dris[i];
    AtomID size=converter.total_number_of_particle;
    std::string ver("alpha");
    double LJEnergyCorrection = converter.tljcecs[i];
    double KineticEnergy = converter.kenergys[i];
    double PotentialEnergy = converter.penergys[i];
    double TotalEnergy = converter.total_penergys[i];
    double Virial = converter.virials[i];
    double EtaPosition = converter.eta_poss[i];
    double EtaVelocity = converter.eta_vels[i];
    double EtaForce = converter.eta_forces[i];
    double LogVPosition = converter.logv_poss[i];
    double LogVVelocity = converter.logv_vels[i];
    double LogVForce = converter.logv_forces[i];
    double Volume = converter.volumes[i];
    int boxtype = 0;
    SpaceVector<double> angle(90.0,90.0,90.0);
    setParameter(hdfdump,dt,step,rc,pc,crdc,rstc,size);
    setInfo(hdfdump,index,number_of_restorefiles,&ver,size,snapshottype);
    setIntegration(hdfdump,
		   LJEnergyCorrection, KineticEnergy, PotentialEnergy, TotalEnergy, Virial,
		   EtaPosition, EtaVelocity, EtaForce,
		   LogVPosition, LogVVelocity, LogVForce, Volume );
    setBoxDatatype(hdfdump, converter.boxsizes[i], boxtype, angle);
    hdfdump.setParticle(converter.particlearrays[i],converter.typerangearrays[i]);
    hdfdump.setBondlists(converter.bondlistarrays[i]);
    if(with_waterlist>0){
      hdfdump.setWaterlist(converter.waterlists[i]);
    }
    hdfdump.setShakelist(converter.shakelists[i]);
    hdfdump.dump();
  }
}

Bin2HDF::Bin2HDF(int argc, char* argv[], MPI_Comm comm):
  argv(argv), argc(argc), comm(comm), comm_rank(0), comm_size(1), args(""),
  filenamebase("binaryrestorefile"), number_of_restorefiles(1),
  hdfname(""), with_waterlist(0), shrink(0), snapshottype(HDFSnapshot::HDFSnapshotFile::Single)
{
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);
  ParseArguments(argc, argv);
  Converter<CombinedParticleArray> converter(filenamebase,number_of_restorefiles,comm);
  if(!hdfname.empty()) {
    if(snapshottype==HDFSnapshot::HDFSnapshotFile::Single){
      outputSingleHDF(converter);
    }else if(snapshottype==HDFSnapshot::HDFSnapshotFile::ParallelDump){
      dumpParallelHDF(converter);
    }
  }
}

void Bin2HDF::NormalizeString(std::string &s) {
  // convert to lower case
  std::transform(s.begin(), s.end(), s.begin(),
                 static_cast<int(*)(int)>(std::tolower));
  // convert '_' to '-'
  std::replace(s.begin(), s.end(), '_', '-');
}

template<typename T>
void Bin2HDF::Read(T &t) {
  std::istringstream iss(::optarg);
  iss.exceptions(std::ios::badbit | std::ios::failbit);
  iss >> t;
}

void Bin2HDF::ReadSnapshotType(HDFSnapshot::HDFSnapshotFile::SnapshotType &t) {
  std::string s(::optarg);
  NormalizeString(s);
  if (s == "0" || s == "single")  { t = HDFSnapshot::HDFSnapshotFile::Single; return; }
  if (s == "1" || s == "parallel") { t = HDFSnapshot::HDFSnapshotFile::Parallel; return; }
  if (s == "2" || s == "singledump") { t = HDFSnapshot::HDFSnapshotFile::SingleDump; return; }
  if (s == "3" || s == "paralleldump") { t = HDFSnapshot::HDFSnapshotFile::ParallelDump; return; }
  throw s;
}

int Bin2HDF::ParseArguments(int argc, char* argv[])
{
  // record arguments
  args += argv[0];
  for (char **p = argv + 1; *p; ++p) { args += " "; args += *p; }

  enum {
    LONG_OPTION_START = 300,

    HELP = LONG_OPTION_START,

    FILENAME,
    NUMBER,
    HDFSNAPSHOT,
    WITHWATER,
    SHRINK,
    TYPE
  };

  for (;;) {
    static struct option long_opts[] = {
      {"help",     no_argument,       0, 'h'},

      {"filename", required_argument, 0, 'f'},
      {"number",   required_argument, 0, 'n'},
      {"hdf",    required_argument, 0, HDFSNAPSHOT},
      {"with-water", no_argument,   0, 'w'},
      {"shrink", no_argument,   0, 's'},
      {"type",  required_argument, 0, 't'},

      {0, 0, 0, 0}
    };
    int opt_index = 0;
    int opt = ::getopt_long(argc, argv,
			       ":hf:n:t:ws",
			       long_opts, &opt_index);
    if (opt == -1) break;
    try{
      switch (opt) {
      case 'f': Read(filenamebase); break;
      case 'n': Read(number_of_restorefiles); break;
      case HDFSNAPSHOT: Read(hdfname); break;
      case 'w' : with_waterlist = 1; break;
      case 's' : shrink = 1; break;
      case 't' : ReadSnapshotType(snapshottype);	break; 
      case 'h'/*HELP*/: {
	{
	  std::cout <<
	    "options:\n"
	    "  -h,  --help\n"
	    "  -f,  --filename<s>\n"
	    "  -n,  --number<s>\n"
	    "       --hdf=<s>\n"
	    "  -w,  --with-water\n"
	    "  -s,  --shrink\n"
	    "  -t,  --type=<s>\n"
		    << std::flush;
	}
	MPI_Finalize();
	exit(EXIT_SUCCESS);
      }
      case ':':  {  // missing option argument
	std::cout << "option requires an argument";
	if (0 < ::optopt && ::optopt < LONG_OPTION_START)
	  std::cout << " -- '" << static_cast<char>(::optopt) << "'";
	else if (::optopt >= LONG_OPTION_START)
	  std::cout << " -- '"
	       << long_opts[::optopt - LONG_OPTION_START].name << "'";
	std::cout << ".  try '-h' for help\n" << std::flush;
	MPI_Finalize();
	exit(EXIT_FAILURE);
      }
      default:  /* case '?': */ {  // unknown option
	std::cout << "cppmd: unknown option";
	if (0 < ::optopt && ::optopt < LONG_OPTION_START)
	  std::cout << " -- '" << static_cast<char>(::optopt) << "'";
	else
	  std::cout << " -- '" << &argv[::optind - 1][2] << "'";
	std::cout << ".  try '-h' for help\n" << std::flush;
	MPI_Finalize();
	exit(EXIT_FAILURE);
      }
      }
    }
    catch (...) {  // invalid option argument
      {
	std::cout << "invalid option argument";
        if (0 < ::optopt && ::optopt < LONG_OPTION_START)
	  std::cout << " -- '" << static_cast<char>(::optopt) << "'";
        else
	  std::cout << " -- '" << argv[::optind - 1] << "'";
	std::cout << ".  try '-h' for help\n" << std::flush;
      }
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }
}
