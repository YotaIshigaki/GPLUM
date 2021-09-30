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
#include "HDF2RST.h"
#include "Converter.h"
#include "HDFSnapshotData.h"

void HDF2RST::restoreSingleHDF()
{
  if(comm_rank==0){
    HDFDump hdfdump(hdfname,HDFRESTORE);
    hdfdump.restore(step);
    HDFSnapshot::HDFSnapshotFile::SnapshotType sst = hdfdump.getSnapshotType();
    if(sst!=snapshottype)
      {
	std::cout << "SnapshotType miss match : requested " << snapshottype << " : file " << sst << std::endl;
      }
    {
      SpaceVector<double> boxsize;
      int type;
      SpaceVector<double> angle;
      hdfdump.getBoxDatatype(boxsize,type,angle);
      std::cout << "boxsize " << boxsize << std::endl;
      CombinedParticleArray pa;
      std::vector<TypeRange> tr;
      hdfdump.getParticle(pa,tr);
      AtomID max = pa.size();
      Dump dump_rst(boxsize,max,1,0,0,rstname,1);
      dump_rst.DumpAmberCRD(pa,tr);
    }
  }
}

void HDF2RST::restoreParallelDumpHDF()
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
  CombinedParticleArray pa;
  std::vector<TypeRange> tr;
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
    HDFDump hdfdump(filename,HDFRESTORE);
    hdfdump.restore(step);
    HDFSnapshot::HDFSnapshotFile::SnapshotType sst = hdfdump.getSnapshotType();
    if(sst!=snapshottype)
      {
	std::cout << "SnapshotType miss match : requested " << snapshottype << " : file " << sst << std::endl;
      }
    {
      SpaceVector<double> angle;
      int type;
      hdfdump.getBoxDatatype(boxsize,type,angle);

      CombinedParticleArray ppa;
      std::vector<TypeRange> ptr;
      hdfdump.getParticle(ppa,ptr);
      tr.resize(tr.size()+ptr.size());
      pa.resize(pa.size()+ppa.size());
      for(int r=0;r<ptr.size();r++){
	int index_shift = n-ptr[r].begin;
	tr[t] = ptr[r];
	tr[t].shift(index_shift);
	for(int p=ptr[r].begin;p<ptr[r].end;p++){
	  AtomID aid = getatomid(ppa,p);
	  if(aid>local_maxid)local_maxid=aid;
	  setparticle(pa,getparticle(ppa,p),n);
	  n++;
	}
	t++;
      }
    }
  }
  std::cout << "Number of resotred typerange " << tr.size() << std::endl;
  std::cout << "Reserved size of particle " << pa.size() << std::endl;
  std::cout << "Number of restored atom " << n << std::endl;

  {
    int lnum,tnum;
    int lmid,tmid;
    lnum = n;
    MPI_Allreduce(&lnum,&tnum,1,MPI_INT,MPI_SUM,comm);
    total_number_of_particle = tnum;
    lmid = local_maxid;
    MPI_Allreduce(&lmid,&tmid,1,MPI_INT,MPI_MAX,comm);
    maxid = tmid;
    if(maxid>=total_number_of_particle){
      if(comm_rank==0){
	printf("Max AtomID %d is larger than number of atom %d\n",maxid,total_number_of_particle);
      }
      total_number_of_particle = maxid;
    } 
    if(comm_rank==0){
      std::cout << "boxsize " << boxsize << std::endl;
      std::cout << "Total number of atom " << total_number_of_particle << std::endl;
    }
    Dump dump_rst(boxsize,total_number_of_particle,comm_size,comm_rank,0,rstname,1);
    dump_rst.GatherDumpAmberCRD(pa,tr);
  }
}

HDF2RST::HDF2RST(int argc, char* argv[], MPI_Comm comm):
    argv(argv), argc(argc), comm(comm), comm_rank(0), comm_size(1), args(""),
    rstname("md.rst"),
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

void HDF2RST::NormalizeString(std::string &s) {
  // convert to lower case
  std::transform(s.begin(), s.end(), s.begin(),
                 static_cast<int(*)(int)>(std::tolower));
  // convert '_' to '-'
  std::replace(s.begin(), s.end(), '_', '-');
}

template<typename T>
void HDF2RST::Read(T &t) {
  std::istringstream iss(::optarg);
  iss.exceptions(std::ios::badbit | std::ios::failbit);
  iss >> t;
}

void HDF2RST::ReadSnapshotType(HDFSnapshot::HDFSnapshotFile::SnapshotType &t) {
  std::string s(::optarg);
  NormalizeString(s);
  if (s == "0" || s == "single")  { t = HDFSnapshot::HDFSnapshotFile::Single; return; }
  if (s == "1" || s == "parallel") { t = HDFSnapshot::HDFSnapshotFile::Parallel; return; }
  if (s == "2" || s == "singledump") { t = HDFSnapshot::HDFSnapshotFile::SingleDump; return; }
  if (s == "3" || s == "paralleldump") { t = HDFSnapshot::HDFSnapshotFile::ParallelDump; return; }
  throw s;
}

int HDF2RST::ParseArguments(int argc, char* argv[])
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
      case 'f': Read(rstname); break;
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
