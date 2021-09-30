#include <mpi.h>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include "Dump.h"
#include "HDF2Bin.h"
#include "Converter.h"
#include "HDFSnapshotData.h"

void HDF2Bin::dumpBin(const std::string& binfilename,
		      const CombinedParticleArray& particlearray,
		      const std::vector<TypeRange>& typerangearray,
		      const std::vector<CovalentBondInfo::BondList>& bondlistarray,
		      const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
		      const WaterList& waterlist,
		      const ShakeList& shakelist,
		      const double od[15],
		      const int oi[4],
		      const long t)
{
  BinaryDump<CombinedParticleArray> binarydump;
  binarydump.init(binfilename);
  binarydump.dump_binary_basic(particlearray, typerangearray,
			       bondlistarray, bondlistarray_idx,
			       waterlist, shakelist,
			       t);
  binarydump.dump_binary_optional_double(od, 15);
  binarydump.dump_binary_optional_int(oi, 4);
  binarydump.close();
}

void HDF2Bin::restoreHDF(const std::string& filename,
			 size_t step,
			 CombinedParticleArray& pa, std::vector<TypeRange>& tr,
			 std::vector<CovalentBondInfo::BondList>& bondlistarray,
			 std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
			 WaterList& waterlist,
			 ShakeList& shakelist,
			 long& t,
			 double od[15],
			 int oi[4])
{
  HDFDump hdfdump(filename,HDFRESTORE);
  hdfdump.restore(step);

  SpaceVector<double> boxsize;
  int type;
  SpaceVector<double> angle;
  hdfdump.getBoxDatatype(boxsize,type,angle);
  std::cout << "boxsize " << boxsize << std::endl;

  hdfdump.getParticle(pa,tr);
  AtomID max = pa.size();

  hdfdump.getBondlists(bondlistarray);
  {
    std::cout <<" bondlist " << bondlistarray.size() << std::endl;
  }
  
  hdfdump.getWaterlist(waterlist);

  hdfdump.getShakelist(shakelist);

  double dt;
  int ti;
  AtomID numatom;
  hdfdump.getIntegration(od[14],od[0],od[1],od[3],od[2],
			 od[4],od[5],od[6],od[7],od[8],od[9],od[13]);
  od[10] = boxsize.x;
  od[11] = boxsize.y;
  od[12] = boxsize.z;
  hdfdump.getParameter(dt,ti,oi[0],oi[1],oi[2],oi[3],numatom);
  t = ti;
}

void HDF2Bin::restoreSingleHDF()
{
  if(comm_rank==0){
    /*
    HDFDump hdfdump(hdfname,HDFRESTORE);
    hdfdump.restore(step);
    HDFSnapshot::HDFSnapshotFile::SnapshotType sst = hdfdump.getSnapshotType();
    if(sst!=snapshottype)
      {
	std::cout << "SnapshotType miss match : requested " << snapshottype << " : file " << sst << std::endl;
      }

      SpaceVector<double> boxsize;
      int type;
      SpaceVector<double> angle;
      hdfdump.getBoxDatatype(boxsize,type,angle);
      std::cout << "boxsize " << boxsize << std::endl;
      CombinedParticleArray pa;
      std::vector<TypeRange> tr;
      hdfdump.getParticle(pa,tr);
      AtomID max = pa.size();
      std::vector<CovalentBondInfo::BondList> bondlistarray;
      std::vector<CovalentBondInfo::BondList> bondlistarray_idx;
      hdfdump.getBondlists(bondlistarray);
      {
	std::cout <<" bondlist " << bondlistarray.size() << std::endl;
      }
      WaterList waterlist;
      hdfdump.getWaterlist(waterlist);
      ShakeList shakelist;
      hdfdump.getShakelist(shakelist);
      long t;
      double od[15];
      int oi[4];
      double dt;
      int ti;
      AtomID numatom;
      hdfdump.getIntegration(od[14],od[0],od[1],od[3],od[2],
			      od[4],od[5],od[6],od[7],od[8],od[9],od[13]);
      od[10] = boxsize.x;
      od[11] = boxsize.y;
      od[12] = boxsize.z;
      hdfdump.getParameter(dt,ti,oi[0],oi[1],oi[2],oi[3],numatom);
      t = ti;
    */
    CombinedParticleArray pa;
    std::vector<TypeRange> tr;
    std::vector<CovalentBondInfo::BondList> bondlistarray;
    std::vector<CovalentBondInfo::BondList> bondlistarray_idx;
    WaterList waterlist;
    ShakeList shakelist;
    long t;
    double od[15];
    int oi[4];
    restoreHDF(hdfname, step,
	       pa, tr,
	       bondlistarray, bondlistarray_idx,
	       waterlist, shakelist,
	       t, od, oi);
    dumpBin(binname,pa,tr,bondlistarray,bondlistarray_idx,waterlist,shakelist,od,oi,t);
  }
}

void HDF2Bin::restoreParallelDumpHDF()
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

  std::cout << namebase << sufix << std::endl;

  number_of_localrestorefiles=0;
  restorefileindexes.clear();
  if(number_of_restorefiles>0){
    int m = number_of_restorefiles/comm_size;
    int l = number_of_restorefiles%comm_size;
    int first_index, last_index;
    if(comm_rank<l){
      number_of_localrestorefiles = m+1;
      first_index = (m+1)*comm_rank;
      last_index = (m+1)*(comm_rank+1)-1;
    }else{
      number_of_localrestorefiles = m;
      first_index = (m+1)*l + m*(comm_rank-l);
      last_index = (m+1)*l + m*(comm_rank-l+1)-1;
    }
    for(int i=first_index;i<=last_index;i++){
      restorefileindexes.push_back(i);
    }
    std::cout << "Rank " << comm_rank << " : number of file " << number_of_localrestorefiles << std::endl;
  }

  SpaceVector<double> boxsize;
  int n=0;
  int t=0;
  int i;
  int li=0;
  local_maxid=0;
  for(li=0;li<number_of_localrestorefiles;li++){
    i = restorefileindexes[li];
    char name[256];
    snprintf(name,256,"%s.%05d",namebase.c_str(),i);
    std::string filename(name);
    filename.append(sufix);
    snprintf(name,256,"%s.%05d",binname.c_str(),i);
    std::string binfilename(name);
    {
      CombinedParticleArray pa;
      std::vector<TypeRange> tr;
      std::vector<CovalentBondInfo::BondList> bondlistarray;
      std::vector<CovalentBondInfo::BondList> bondlistarray_idx;
      WaterList waterlist;
      ShakeList shakelist;
      long t;
      double od[15];
      int oi[4];
      restoreHDF(filename, step,
		 pa, tr,
		 bondlistarray, bondlistarray_idx,
		 waterlist, shakelist,
		 t, od, oi);
      dumpBin(binfilename,pa,tr,bondlistarray,bondlistarray_idx,waterlist,shakelist,od,oi,t);
    }
  }
}

HDF2Bin::HDF2Bin(int argc, char* argv[], MPI_Comm comm):
    argv(argv), argc(argc), comm(comm), comm_rank(0), comm_size(1), args(""),
    binname("md.bin"),
    hdfname(""), step(0), snapshottype(HDFSnapshot::HDFSnapshotFile::Single)
{
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);
  ParseArguments(argc, argv);

  if(number_of_restorefiles>1)snapshottype=HDFSnapshot::HDFSnapshotFile::ParallelDump;
  if(!hdfname.empty()){
    if(snapshottype==HDFSnapshot::HDFSnapshotFile::Single){
      restoreSingleHDF();
    }else if(snapshottype==HDFSnapshot::HDFSnapshotFile::ParallelDump){
      restoreParallelDumpHDF();
    }
  }
}

void HDF2Bin::NormalizeString(std::string &s) {
  // convert to lower case
  std::transform(s.begin(), s.end(), s.begin(),
                 static_cast<int(*)(int)>(std::tolower));
  // convert '_' to '-'
  std::replace(s.begin(), s.end(), '_', '-');
}

template<typename T>
void HDF2Bin::Read(T &t) {
  std::istringstream iss(::optarg);
  iss.exceptions(std::ios::badbit | std::ios::failbit);
  iss >> t;
}

void HDF2Bin::ReadSnapshotType(HDFSnapshot::HDFSnapshotFile::SnapshotType &t) {
  std::string s(::optarg);
  NormalizeString(s);
  if (s == "0" || s == "single")  { t = HDFSnapshot::HDFSnapshotFile::Single; return; }
  if (s == "1" || s == "parallel") { t = HDFSnapshot::HDFSnapshotFile::Parallel; return; }
  if (s == "2" || s == "singledump") { t = HDFSnapshot::HDFSnapshotFile::SingleDump; return; }
  if (s == "3" || s == "paralleldump") { t = HDFSnapshot::HDFSnapshotFile::ParallelDump; return; }
  throw s;
}

int HDF2Bin::ParseArguments(int argc, char* argv[])
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
    HDFSTEP,
    TYPE
  };

  for (;;) {
    static struct option long_opts[] = {
      {"help",     no_argument,       0, 'h'},

      {"filename", required_argument, 0, 'f'},
      {"number",   required_argument, 0, 'n'},
      {"hdf",    required_argument, 0, HDFSNAPSHOT},
      {"step",   required_argument, 0, 's'},
      {"type",  required_argument, 0, 't'},

      {0, 0, 0, 0}
    };
    int opt_index = 0;
    int opt = ::getopt_long(argc, argv,
			       ":hf:s:",
			       long_opts, &opt_index);
    if (opt == -1) break;
    try{
      switch (opt) {
      case 'f': Read(binname); break;
      case 'n': Read(number_of_restorefiles); break;
      case HDFSNAPSHOT: Read(hdfname); break;
      case 's' : Read(step); break;
      case 't' : ReadSnapshotType(snapshottype);	break; 
     case 'h'/*HELP*/: {
	{
	  std::cout <<
	    "options:\n"
	    "  -h,  --help\n"
	    "  -n,  --number<s>\n"
	    "  -f,  --filename<s>\n"
	    "       --hdf=<s>\n"
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
