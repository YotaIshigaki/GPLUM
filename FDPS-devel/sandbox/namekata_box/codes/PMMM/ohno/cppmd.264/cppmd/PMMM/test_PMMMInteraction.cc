#ifdef _OPENMP
#include <omp.h>
#endif
#include <mpi.h>

#include <cstdlib>

#include "Common.h"

#include "OperationSelector.h"
#include "CellIndex.h"
#include "CubicCell.h"
#include "CalcForce.h"
#include "MPIParallel.h"

#include "PMMMPreparator.h"

#include "PMMMInteraction.h"

#include "PairListInteraction.h"

#include "PP.h"

#include "ewald.h"

static void print_f(
		    const ForceArray &key,
		    const ForceArray &f,
		    const char * name,
		    const int icut,
		    const int p)
{
  static char fname[256];
  sprintf(fname, "%s.c%dp%d.dat", name, icut, p);
  FILE *fp = fopen(fname, "w");
  assert(fp);
  const int len = std::min(key.size(),f.size());
  
  for(int i=0;i<len;i++){
    fprintf(fp, "%e %e %e   %e %e %e\n", key[i].x, key[i].y, key[i].z, f[i].x, f[i].y, f[i].z);
  }
  fclose(fp);
}

static void print_err(
		      const std::vector<double> &key,
		      const std::vector<double> &err,
		      const char * name,
		      const int icut,
		      const int p)
{
  static char fname[256];
  sprintf(fname, "%s.c%dp%d.dat", name, icut, p);
  FILE *fp = fopen(fname, "w");
  assert(fp);
  const int len = std::min(key.size(),err.size());
  
  for(int i=0;i<len;i++){
    fprintf(fp, "%e %e\n", key[i], err[i]);
  }
  fclose(fp);
}

static Position rand_vec(){
  return Position(drand48(), drand48(), drand48());
}

template<typename PA>
static double gen_rand_dist(
			  const int np,
			  PA &ptcl,
			  const double scale,
			  const long seed = 19810614)
{
  //  printf("gen_rand_dist %d atom\n",np);

  srand48(seed);

  ptcl.resize(np);
  for(int i=0; i<np; i++){
    getpos(ptcl,i) = rand_vec()*scale;
    getatomtype(ptcl,i) = 12;
    getatomid(ptcl,i) = i;
  }
  double msum = 0.0;
  for(int i=0; i<np; i++){
    msum += 
      (getcharge(ptcl,i)  = drand48());
  }
  double coffset = msum/np;
#ifdef NON_CHARGE_NEUTRAL
  coffset *= 0.5;
#endif
  for(int i=0; i<np; i++){
    getcharge(ptcl,i) -= coffset;
  }

  msum = 0.0;
  for(int i=0; i<np; i++){
    msum += getcharge(ptcl,i);
  }
  return msum;
}

template<typename PA>
static double gen_testnp(
			  const int np,
			  PA &ptcl,
			  const double scale)
{
  //  printf("gen_rand_dist %d atom\n",np);

  ptcl.resize(np);
  for(int i=0; i<np; i++){
    getpos(ptcl,i).x = (0.5+(double)(i-(np>>1))*0.1)*scale;
    getpos(ptcl,i).y = 0.0;
    getpos(ptcl,i).z = 0.0;
    getatomtype(ptcl,i) = 12;
    getatomid(ptcl,i) = i;
  }
  double msum = 0.0;
  for(int i=0; i<np; i++){
    msum += 
      (getcharge(ptcl,i)  = 1.0*((i&0x1)*2-1));
  }
  double coffset = msum/np;
#ifdef NON_CHARGE_NEUTRAL
  coffset *= 0.5;
  if(coffset<1.0e-6)coffset=-0.25;
#endif
  for(int i=0; i<np; i++){
    getcharge(ptcl,i) -= coffset;
  }

  msum = 0.0;
  for(int i=0; i<np; i++){
    msum += getcharge(ptcl,i);
  }
  return msum;
}


int
main(int argc, char *argv[])
{

  enum{
    NP   = 128,
    NC   = 8,
    NC3  = NC*NC*NC,
    PFMM = 5,
    ICUT = 2,
  };

  MPI_Init(&argc,&argv);

  DebugLog::verbose = 2;
#ifdef _OPENMP
  int num_omp_threads = omp_get_num_threads();
  int max_omp_threads = omp_get_max_threads();
  printf("OpenMP num threads %d  max threads %d\n",num_omp_threads,max_omp_threads);
#endif

  int node_id;
  int num_node;

  MPI_Comm comm = MPI_COMM_WORLD;
  
  MPI_Comm_rank(comm, &node_id);
  MPI_Comm_size(comm, &num_node);

  if(num_node<2){
    printf("MPI_Comm_size smaller than 2\n");
    exit(EXIT_FAILURE);
  }

  int num_pp;
  int num_mm;

  num_pp = 1;
  while(num_pp*2<num_node){
    num_pp *= 2;
  }

  num_mm = num_node-num_pp;

  OperationSelector operation;
  if(node_id<num_pp){
    operation.doShortrangecalculation  = true;
    operation.doLongrangecalculation   = false;
    operation.doShortLongCommunication = true;
  }else{
    operation.doShortrangecalculation  = false;
    operation.doLongrangecalculation   = true;
    operation.doShortLongCommunication = true;
  }
  operation.doEnergycalculation = true;

  MPI_Comm ppmm_comm;
  MPI_Comm pp_comm;
  MPI_Comm mm_comm;

  int ppmm_comm_color=0;
  if(operation.doLongrangecalculation)ppmm_comm_color=1;
  MPI_Comm_split(comm,ppmm_comm_color,node_id,&ppmm_comm);

  int short_id;
  if(operation.doShortrangecalculation){
    MPI_Comm_rank(ppmm_comm, &short_id);
    pp_comm = ppmm_comm;
  }else{
    short_id = NO_SHORT_ID;
    pp_comm = MPI_COMM_NULL;
  }
  if(operation.doLongrangecalculation){
    mm_comm = ppmm_comm;
  }else{
    mm_comm = MPI_COMM_NULL;
  }

  CombinedParticleArray particlearray;
  std::vector<int> particle_setid;
  std::vector<TypeRange> typerangearray;
 
  std::vector<int> sendparticle_target_id;
  std::vector< std::vector<int> > sendparticle_setid;
  std::vector<int> recvparticle_target_id;
  std::vector< std::vector<int> > recvparticle_setid;
  std::vector<int> move_target_id;
  std::vector< std::vector<int> > move_setid;
  CellIndex cellindex;

  double cutoff = 12.0;
  double cellmargin1d = 1.0;
  double box_len = 62.23;
  SpaceVector<double> boxsize(box_len,box_len,box_len);
  SpaceVector<double> cellmargin(cellmargin1d,cellmargin1d,cellmargin1d);
  int cell_div1d = 8;
  SpaceVector<int> node_div(2,2,2);
  SpaceVector<int> cell_div3d(cell_div1d,cell_div1d,cell_div1d);
  SpaceVector<int> cell_div3d_in_node(4,4,4);

  CellSubset cellsubset;

  CellIndexType cellindextype = FullCube;

  std::vector< std::vector<int> > targetset;
  std::vector< std::pair<int,int> > ghostpairs;

  GhostParticleArray ghost;
  std::vector<int> ghostsetid;
  std::vector<TypeRange> ghosttyperange;
  std::vector<ParticleRange> targetparticlerange;
  std::vector<int> recvsetid;
  std::map<int,int> recvsetid_to_index;


  ShortRange::CoulombType cltype = ShortRange::OriginalCoulomb;

  std::vector<CovalentBondInfo::BondList> bondlistarray;
  std::vector<CovalentBondInfo::BondList> bondlistarray_idx;
  std::vector<CovalentBondInfo::BondList> targetbond(0);

  PostProcess postprocess;

  int num_total_set = cell_div3d[0] * cell_div3d[1] * cell_div3d[2];
  int max_particle_in_cell = 200;
  {
    int m[3];
    if(power_of_two_of_cube(num_pp,m)){
      node_div[0] = m[0];
      node_div[1] = m[1];
      node_div[2] = m[2];
    }
    cell_div3d_in_node[0] = cell_div3d[0]/node_div[0];
    cell_div3d_in_node[1] = cell_div3d[1]/node_div[1];
    cell_div3d_in_node[2] = cell_div3d[2]/node_div[2];

    cellindex.set_cell_parameter(boxsize, cell_div3d, cellmargin, cutoff);
    cellindex.calc_volume();

    postprocess.celldiv = cellindex.celldiv;
    postprocess.cell_position = cellindex.cell_position;

    printf("cell size %f %f %f\n",cellindex.cellsize.x,cellindex.cellsize.y,cellindex.cellsize.z);
    printf("cell target range %d %d %d\n",cellindex.target_range.x,cellindex.target_range.y,cellindex.target_range.z);

    int num_cell = cellindex.generate_cell_index();

    if(!cellindex.distribute_cell(node_div)){
      std::cout << "distribute_cell false " << std::endl;
    }

    int num_set_per_node = num_total_set/num_pp;

    bondlistarray.clear();
    bondlistarray_idx.clear();

    particlearray.resize(max_particle_in_cell*num_set_per_node);

    if(operation.doShortrangecalculation){
      particle_setid.resize(num_set_per_node);
      typerangearray.resize(num_set_per_node);  
      for(size_t s=0;s<particle_setid.size();s++){
	particle_setid[s] = cellindex.cell_id_in_node[short_id][s];
	typerangearray[s].begin = s*max_particle_in_cell;
	typerangearray[s].end = typerangearray[s].begin;
	typerangearray[s].lj.begin = typerangearray[s].begin;
	typerangearray[s].lj.end = typerangearray[s].end;
	typerangearray[s].ljcoulomb.begin = typerangearray[s].lj.end;
	typerangearray[s].ljcoulomb.end = typerangearray[s].lj.end;
	typerangearray[s].coulomb.begin = typerangearray[s].lj.end;
	typerangearray[s].coulomb.end = typerangearray[s].lj.end;
      }

      cellsubset = CellSubset(cellindex,short_id);
      cellsubset.makefulltarget(particle_setid,cellindextype,operation);

      postprocess.shift_cell=cellsubset.shift_cell_all;

      cellsubset.makecommtarget(sendparticle_target_id, 
				sendparticle_setid,
				recvparticle_target_id,
				recvparticle_setid,
				move_target_id, move_setid);
      {
	printf("recvparticle_target_id");
	for(int i;i<recvparticle_target_id.size();i++){
	  printf(" %d",recvparticle_target_id[i]);
	}
	printf("\n");
      }
      targetset = cellsubset.short_target_cell;
      ghostpairs = cellsubset.ghost_cell_pairs;
      bondlistarray.resize(num_set_per_node);
      bondlistarray_idx.resize(num_set_per_node);
     }else{
      particle_setid.resize(0);
      typerangearray.resize(0);
      sendparticle_target_id.resize(0);
      sendparticle_setid.resize(0);
      recvparticle_target_id.resize(0);
      recvparticle_setid.resize(0);
      move_target_id.resize(0);
      move_setid.resize(0);
      targetset.resize(0);
      ghostpairs.resize(0);
      bondlistarray.resize(0);
      bondlistarray_idx.resize(0);
    }
  }

  {
    postprocess.max_particle_in_cell = max_particle_in_cell;
    postprocess.margin = cellindex.margin;
    postprocess.fine_margin = cellindex.margin;
    postprocess.boxsize = cellindex.boxsize;
    postprocess.cellsize = cellindex.cellsize;
    postprocess.invcellsize.x = 1.0/postprocess.cellsize.x;
    postprocess.invcellsize.y = 1.0/postprocess.cellsize.y;
    postprocess.invcellsize.z = 1.0/postprocess.cellsize.z;

    postprocess.calc_volume_margin();
  }


  if(node_id==0){
    printf("num_pp %d  num_mm %d\n",num_pp,num_mm);
  }

  if(0)
  {
    int i;
    for(i=0;i<cellsubset.short_target_cell.size();i++){
      printf("%d :",i);
      for(int j=0;j<cellsubset.short_target_cell[i].size();j++){
	printf(" %d",cellsubset.short_target_cell[i][j]);
      }
      printf("\n");
    }
  }

  WaterList waterlist;
  ShakeList shakelist;

  double total_charge = 0.0;
  CombinedParticleArray pa;
  int nump = NP;
#ifdef EWALD_FILE
        {
	  FILE *fp = fopen("qpos.dat", "r");
	  assert(fp);
	  fscanf(fp, "%d", &nump);
	  pa.resize(nump);
	  for(int i=0; i<nump; i++){
	    double m,x,y,z;
	    fscanf(fp, "%lf %lf %lf %lf", &m,&x,&y,&z);
	    getcharge(pa,i) = m;
	    total_charge += m;
	    getpos(pa,i).x = x*box_len;
	    getpos(pa,i).y = y*box_len;
	    getpos(pa,i).z = z*box_len;
	    getatomtype(pa,i) = 12;
	    getatomid(pa,i) = i;
	  }
	  fclose(fp);
        }

#else
# ifdef TEST2P
	nump = 2;
	total_charge = gen_testnp(nump,pa,box_len);
# else
    total_charge = gen_rand_dist(nump,pa,box_len);
# endif
#endif
    printf("%d atom\n",pa.size());
    printf("total_charge %e\n",total_charge);

  if(operation.doShortrangecalculation){
    WaterList wl(pa);
    ShakeList sl(pa);
    std::vector<PotentialModel> pm(nump,LJCoulombPotential);

    {
      int n;
      n = cellsubset.distribute_particle_cell(particlearray,
					      particle_setid,
					      typerangearray,
					      waterlist,
					      shakelist,
					      pa, pm, wl, sl);
      double cs=0.0;
      for(int t=0;t<typerangearray.size();t++){
	for(int i=typerangearray[t].begin;i<typerangearray[t].end;i++)cs+=getcharge(particlearray,i);
      }
      printf("local atom %d / %d , local charge sum %f\n",n,pa.poscharge.size(),cs);
    }
  }

    std::map<int,int> shortidtorank;
    for (int u = 0; u < num_pp; ++u) shortidtorank.insert(std::pair<int,int>(u,u));
    CommPattern short_c_p = Direct;
    CommPattern move_c_p = Direct;
    int max_id = num_total_set;
    MPICommunicator communicator(shortidtorank,
				   pp_comm,
				   sendparticle_target_id,
				   sendparticle_setid,
				   particle_setid,
				   recvparticle_target_id,
				   recvparticle_setid,
				   move_target_id,
				   move_setid,
				   max_particle_in_cell,
				   max_id, short_id,
				   node_div,
				   cellindex.celldiv,
				   cellindex.move_range,
				   cellindex.node_move_surface_ratio,
				   short_c_p,
				   move_c_p
				   );


  double alpha  = 2.4/box_len;
  int order = 5;

  LongRangeParameter longparam;
  {
    longparam.cutoff = cutoff;
    longparam.alpha = alpha;
    longparam.order = order;
    longparam.boxSize = boxsize;
  }
  
  double pairlistmargin = 2.0;
  bool expire_reverse = false;
 
  CalcForce calcforce(particle_setid, targetset, ghostpairs, cltype, longparam, node_id, short_id, cutoff, 
#ifdef USE_PAIRLIST
		      pairlistmargin,
		      int(particlearray.size()),
#endif
		      mm_comm,expire_reverse);
  
  { // CalcPreparator_::construct_pmmm_target_cell
    std::vector< std::vector<int> > pmmm_target_cell;
    std::vector<int> scell;
    std::vector<int> rcell;
    ShiftCellArray shift_cell;
    pmmm_target_cell.resize(particle_setid.size());
    for(int i=0;i<particle_setid.size();i++){
      scell.clear();
      rcell.clear();
      shift_cell.clear();
      pmmm_target_cell[i].clear();
      cubic_target(cellindex.celldiv,cellindex.cell_position[particle_setid[i]],cellindex.target_range.x,scell,rcell,shift_cell);
      pmmm_target_cell[i] = rcell;
    }
    calcforce.clear_pmmmtarget_set();
    calcforce.generate_pmmmtarget_set(pmmm_target_cell,expire_reverse);
    if(0)
    {
      for(int s=0;s<calcforce.self_pmmmtarget_set.size();s++){
	printf("self_pmmmtarget_set %d :",s);
	for(int i=0;i<calcforce.self_pmmmtarget_set[s].size();i++){
	  printf(" %d",calcforce.self_pmmmtarget_set[s][i]);
	}
	printf("\n");
      }
      for(int s=0;s<calcforce.ghost_pmmmtarget_set.size();s++){
	printf("ghost_pmmmtarget_set %d :",s);
	for(int i=0;i<calcforce.ghost_pmmmtarget_set[s].size();i++){
	  printf(" %d",calcforce.ghost_pmmmtarget_set[s][i]);
	}
	printf("\n");
      }

      for(int s=0;s<calcforce.self_pmmmtarget_index.size();s++){
	printf("self_pmmmtarget_index %d :",s);
	for(int i=0;i<calcforce.self_pmmmtarget_index[s].size();i++){
	  printf(" %d",calcforce.self_pmmmtarget_index[s][i]);
	}
	printf("\n");
      }
      for(int s=0;s<calcforce.ghost_pmmmtarget_index.size();s++){
	printf("ghost_pmmmtarget_index %d :",s);
	for(int i=0;i<calcforce.ghost_pmmmtarget_index[s].size();i++){
	  printf(" %d",calcforce.ghost_pmmmtarget_index[s][i]);
	}
	printf("\n");
      }
    }
  }

  std::vector<CellMethodModule::CellRange> cell_range;
  std::vector<Position> cell_center;
#if 0
  CellMethodModule::CellRange my_range;
  {
    my_range.min = cellindex.celldiv;
    my_range.max.x = my_range.max.y = my_range.max.z = 0;
    Position cs = cellindex.cellsize;
    cell_center.resize(cellsubset.num_mycell);
    for(int i=0;i<cellsubset.num_mycell;i++){
      SpaceVector<int> pos = cellindex.cell_position[cellsubset.mycell[i]];
      if(pos.x>my_range.max.x)my_range.max.x=pos.x;
      if(pos.y>my_range.max.y)my_range.max.y=pos.y;
      if(pos.z>my_range.max.z)my_range.max.z=pos.z;
      if(pos.x<my_range.min.x)my_range.min.x=pos.x;
      if(pos.y<my_range.min.y)my_range.min.y=pos.y;
      if(pos.z<my_range.min.z)my_range.min.z=pos.z;
      cell_center[i].x = cs.x*(0.5+pos.x);
      cell_center[i].y = cs.y*(0.5+pos.y);
      cell_center[i].z = cs.z*(0.5+pos.z);
      //      printf(" (%f %f %f)",cell_center[i].x,cell_center[i].y,cell_center[i].z);
    }
    //    printf("\n");
    my_range.max.x++;
    my_range.max.y++;
    my_range.max.z++;
    //    printf("cell_range %d:%d %d:%d %d:%d\n",my_range.min.x,my_range.max.x,my_range.min.y,my_range.max.y,my_range.min.z,my_range.max.z);
  }
  if(operation.doShortrangecalculation){
    printf("Gather cell_range %d\n",node_id);
    MPI_Gather(&(my_range.min[0]), 6, MPI_INTEGER,
	       &(cell_range[0].min[0]), 6, MPI_INTEGER, 0, pp_comm);
  }
  {
    if(node_id==0){
      printf("send cell_range to %d : by %d\n",num_pp,node_id);
      MPI_Send(&(cell_range[0].min[0]), 6*num_pp, MPI_INTEGER, num_pp, 0, comm);
    }
    if(node_id==num_pp){
      MPI_Status stat;
      printf("recv cell_range to %d : by %d\n",0, node_id);
      MPI_Recv(&(cell_range[0].min[0]), 6*num_pp,  MPI_INTEGER, 0, 0, comm, &stat);
    }
  }
  if(operation.doLongrangecalculation){
    printf("Bcast cell_range %d\n",node_id);
    MPI_Bcast(&(cell_range[0].min[0]), 6*num_pp,  MPI_INTEGER, 0, mm_comm);
  }
#else
  PMMMPreparator::set_cell_range(comm, pp_comm, num_pp, mm_comm,
				 cellindex, cellsubset, operation,
				 cell_center, cell_range);
#endif

  std::vector<int> m_boundary;
  std::vector<int> mm_targetid;
  std::vector<int> pp_targetid;
#if 0
  
  mm_decompose(m_boundary,num_mm,order);
  if(0)
  {
    printf("m_boundary[%d]",m_boundary.size());
    for(int m=0;m<m_boundary.size();m++){
      printf(" %d",m_boundary[m]);
    }
    printf("\n");
  }
  mm_targetid.resize(num_mm);
  for(int m=0;m<num_mm;m++)mm_targetid[m]=num_pp+m;
  pp_targetid.resize(num_pp);
  for(int m=0;m<num_pp;m++)pp_targetid[m]=m;
#else
  PMMMPreparator::make_mm_info(order, num_pp, num_mm,
			       m_boundary, pp_targetid, mm_targetid);
			       
#endif

  if(0)
  {
    printf("particle_setid");
    for(int i=0;i<particle_setid.size();i++){
      printf(" %d", particle_setid[i]);
    }
    printf("\n");
  }
  PMMMLocalInteraction pmmm_local;
  PMMMLongRangeInteraction pmmm_long;

#if 0
  if(operation.doShortrangecalculation){
    pmmm_local.settings(node_id, longparam, ppmm_comm);
    pmmm_local.initialize(cell_center, &(cellindex.celldiv_node[0]), 
			  comm, mm_targetid, m_boundary);
  }
  if(operation.doLongrangecalculation){
    pmmm_long.settings(node_id, longparam, ppmm_comm);
    pmmm_long.initialize(&(cellindex.celldiv[0]), m_boundary, cellindex.cellsize.x,
			 3, 5, 2,
			 comm, pp_targetid, cell_range);
    printf("local mm %d %d\n",pmmm_long.pmmm_mm.m_begin,pmmm_long.pmmm_mm.m_end);
  }
#else
  PMMMPreparator::initialize_PMMM_Interaction(comm,pp_comm,mm_comm,
					      cellindex, longparam, operation,
					      m_boundary, pp_targetid, mm_targetid,
					      cell_center, cell_range, total_charge,
					      pmmm_local, pmmm_long);
#endif

  double energy = 0.0;
  double total_energy = 0.0;
  double shortenergy = 0.0;
  double total_shortenergy = 0.0;
  ForceArray shortforce(particlearray.size());
  MPI_Barrier(comm);
  if(operation.doShortrangecalculation){


    printf("calc_local_top_half\n");
    pmmm_local.calc_local_top_half(particlearray, typerangearray);

    communicator.resize_receive_data(targetparticlerange);
    printf("exchangeParticleArraysubset_top_half\n");
    communicator.exchangeParticleArraysubset_top_half(particlearray, 
						      typerangearray,
						      bondlistarray, 
						      ghost,
						      targetparticlerange, 
						      recvsetid, 
						      recvsetid_to_index,
						      ghosttyperange, 
						      targetbond);

    printf("cellpairsinteraction\n");
    cellpairsinteraction(particlearray,typerangearray,
			 particlearray,typerangearray, calcforce.self_pmmmtarget_index,
			 shortforce, shortenergy,true, operation);
    printf("exchangeParticleArraysubset_bottom_half\n");
    communicator.exchangeParticleArraysubset_bottom_half(particlearray, 
							 typerangearray, 
							 bondlistarray, 
							 ghost,
							 targetparticlerange, 
							 recvsetid, 
							 recvsetid_to_index,
							 ghosttyperange, 
							 targetbond);
    postprocess.postprocess_receive(ghost,recvsetid,recvsetid_to_index,ghosttyperange);


    for(size_t s=0;s<calcforce.ghost_pmmmtarget_set.size();s++){
      calcforce.convert_setid_to_index(recvsetid_to_index,calcforce.ghost_pmmmtarget_set[s],
                             calcforce.ghost_pmmmtarget_index[s]);
    }
    cellpairsinteraction(particlearray,typerangearray,
			 ghost, ghosttyperange, calcforce.ghost_pmmmtarget_index,
			 shortforce, shortenergy, false, operation);

    printf("calc_long_progress\n");
    pmmm_local.calc_long_progress();
    printf("calc_local_bottom_half\n");
    pmmm_local.calc_local_bottom_half(particlearray, typerangearray, particlearray.force, energy);
    printf("Done\n");
    shortenergy*=0.5;
    energy*=0.5;
    MPI_Allreduce(&shortenergy, &total_shortenergy, 1, MPI_DOUBLE, MPI_SUM, pp_comm);
    if(node_id==0){
      printf("PMMM total_shortenergy %24.16e\n",total_shortenergy);
    }
    MPI_Allreduce(&energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, pp_comm);
    if(node_id==0){
      printf("PMMM PM total_energy %24.16e\n",total_energy);
    }
    if(node_id==0){
      printf("PMMM total_energy %24.16e\n",total_energy+total_shortenergy);
    }
  }
  if(operation.doLongrangecalculation){
    printf("prep\n");
    pmmm_long.prep();
    printf("calcForce\n");
    pmmm_long.calcForce();
    printf("post\n");
    pmmm_long.post();
    printf("Done\n");
  }

  double inv_box_len = 1.0/box_len;
  double inv_box_len2 = inv_box_len*inv_box_len;

  if(operation.doShortrangecalculation){
    if(node_id==0){
      int npa = particlearray.size()*num_pp;
      int ntr = typerangearray.size()*num_pp;
      CombinedParticleArray fullpa(npa);
      std::vector<TypeRange> ta(ntr);
      CellIndex fullcell;
      std::vector<int> fullps(ntr);
      fullcell.set_cell_parameter(boxsize,cell_div3d, cellmargin, cutoff);
      fullcell.generate_cell_index();
      SpaceVector<int> one_div(1,1,1);
      fullcell.distribute_cell(one_div);
      for(size_t s=0;s<ntr;s++){
	fullps[s] = fullcell.cell_id_in_node[0][s];
	ta[s].begin = s*max_particle_in_cell;
	ta[s].end = ta[s].begin;
	ta[s].lj.begin = ta[s].begin;
	ta[s].lj.end = ta[s].end;
	ta[s].ljcoulomb.begin = ta[s].lj.end;
	ta[s].ljcoulomb.end = ta[s].lj.end;
	ta[s].coulomb.begin = ta[s].lj.end;
	ta[s].coulomb.end = ta[s].lj.end;
      }
      if(0)printf("full cell num %d\n",fullcell.cell_id_in_node[0].size());
      CellSubset fullcellsubset;
      fullcellsubset = CellSubset(fullcell,0);
      WaterList wla;
      ShakeList sla;
      WaterList wl(pa);
      ShakeList sl(pa);
      std::vector<PotentialModel> pm(nump,LJCoulombPotential);
      if(0){
	for(int i=0;i<pa.size();i++){
	  printf(" %d (%f %f %f) %f",i, getpos(pa,i).x,getpos(pa,i).y,getpos(pa,i).z,
		 getcharge(pa,i));
	}
	printf("\n");
      }
      int n;
      n = fullcellsubset.distribute_particle_cell(fullpa,
						  fullps,
						  ta,
						  wla,
						  sla,
						  pa, pm, wl, sl);
      if(0)printf("fullcellsubset.distribute_particle_cell return %d\n",n);
      double inv_box_len = 1.0/box_len;
      ForceArray acc_pp(npa);
      std::vector<double> phi_pp(npa);
      if(0){
	for(int t=0;t<ntr;t++){
	  if(ta[t].begin<ta[t].end){
	    printf(" range %d %d:%d",t,ta[t].begin,ta[t].end);
	  }
	}
	printf("\n");
      }
      n=0;
      for(int t=0;t<ntr;t++){
	for(int i=ta[t].begin;i<ta[t].end;i++){
	  n++;
	  getpos(fullpa,i) *= inv_box_len;
	  acc_pp[i].x = acc_pp[i].y = acc_pp[i].z = 0.0;
	  phi_pp[i] = 0.0;
	}
      }
      if(0){
	for(int t=0;t<ntr;t++){
	  for(int i=ta[t].begin;i<ta[t].end;i++){
	    printf(" %d:%d (%f %f %f) %f",n,i, getpos(fullpa,i).x,getpos(fullpa,i).y,getpos(fullpa,i).z, getcharge(fullpa,i));
	  }
	}
	printf("\n");
      }
      PP_interact_PBC(fullpa.poscharge, ta, acc_pp, phi_pp,
		      cellindex.target_range.x, cellindex.celldiv.x, cellindex.celldiv.y, cellindex.celldiv.z);
      double e=0.0;
      int t;
      //#pragma omp parallel for reduction(+:e)
      for(t=0;t<ta.size();t++){
	for(int i=ta[t].begin;i<ta[t].end;i++){
	  e += 0.5*phi_pp[i]*getcharge(fullpa,i);
	  acc_pp[i] *= getcharge(fullpa,i)*inv_box_len2;
	  //	  printf("PP force[%d] %e %e %e\n",getatomid(fullpa,i),acc_pp[i].x,acc_pp[i].y,acc_pp[i].z);
	}
      }
      e *= inv_box_len;
      printf("energy PP : %24.16e\n", e);
    }
  }

#ifndef PMMM_ONLY
  ForceArray acc_direct(nump);
  PosChargeArray  ptcl(pa.poscharge);
  if(node_id==0){
    std::vector<double> phi_k(nump,0.0);
    std::vector<double> phi_r(nump,0.0);
    std::vector<double> phi_direct(nump,0.0);

    double csum = 0.0;

    for(int i=0;i<nump;i++){
      //      ptcl[i].position *= inv_box_len;
      ptcl[i].position.x *= inv_box_len;
      ptcl[i].position.y *= inv_box_len;
      ptcl[i].position.z *= inv_box_len;
      csum += ptcl[i].charge;
      //      printf(" (%f %f %f)",ptcl[i].position.x,ptcl[i].position.y,ptcl[i].position.z);
    }
    //    printf("\n");
    double ail = alpha*box_len;

#ifdef _OPENMP
  omp_set_num_threads(max_omp_threads);
  num_omp_threads = omp_get_num_threads();
  printf("OpenMP num threads %d  max threads %d\n",num_omp_threads,max_omp_threads);
#endif
    eval_k_space(nump,ail,ptcl,acc_direct,phi_k,5);
    eval_r_space(nump,ail,csum,ptcl,acc_direct,phi_r,3);

    for(int i=0; i<nump; i++){
      phi_direct[i] = phi_k[i] + phi_r[i];
    }
#ifndef EWALD_FILE
    {
      FILE *fp = fopen("qpos.dat", "w");
      assert(fp);
      int n = nump;
      fprintf(fp, "%d\n", n);
      for(int i=0; i<nump; i++){
	fprintf(fp, "%A %A %A %A\n",
		ptcl[i].charge,
		ptcl[i].position.x,
		ptcl[i].position.y,
		ptcl[i].position.z);
      }
      fclose(fp);
    }
    {
      FILE *fp = fopen("ewald.dat", "w");
      assert(fp);
      int n = nump;
      fprintf(fp, "%d\n", n);
      for(int i=0; i<nump; i++){
	fprintf(fp, "%A %A %A %A\n",
	       phi_direct[i],
	       acc_direct[i].x,
	       acc_direct[i].y,
	       acc_direct[i].z);
      }
      fclose(fp);
    }
#endif

    double en_dir=0.0;
    double en_k=0.0, en_r = 0.0;
    for(int i=0; i<nump; i++){
      en_dir += 0.5*phi_direct[i]*ptcl[i].charge;
      en_k += 0.5*phi_k[i]*ptcl[i].charge;
      en_r += 0.5*phi_r[i]*ptcl[i].charge;
      acc_direct[i] *= ptcl[i].charge*inv_box_len2;
    }
    en_dir *= inv_box_len;
    printf("energy k %24.16e + r %24.16e =  %24.16e\n", en_k,en_r,en_dir);

  }

  if(node_id==0){
    printf("scaleing %e\n",inv_box_len);
    ForceArray acc_bn(nump);
    std::vector<double> phi_bn(nump,0.0);
    double e0 = 0.0;
    int nmir  = 3;
    brute_nbody(nump,ptcl,acc_bn,phi_bn,e0,nmir);
    printf("energy no mirror of %d  %24.16e\n",nmir,e0*0.5*inv_box_len);

    double en_bn=0.0;
    for(int i=0; i<nump; i++){
      en_bn += 0.5 * phi_bn[i]*ptcl[i].charge;
      acc_bn[i] *= ptcl[i].charge*inv_box_len2;
    }
    en_bn *= inv_box_len;
    printf("energy : %24.16e\n", en_bn);

    for(int i=0;i<nump;i++){
      printf("%e %e %e\n",acc_bn[i].x,acc_bn[i].y,acc_bn[i].z);
    }

  }
  
  if(operation.doShortrangecalculation){
    int np =  particlearray.size();
    int nt =  typerangearray.size();
    ForceArray fall(nump);
    std::vector<int> atomid(nump+1);
    ForceArray f(nump);
    if(node_id==0){
      for(int t=0;t<nt;t++){
	for(int i=typerangearray[t].begin;i<typerangearray[t].end;i++){
	  fall[getatomid(particlearray,i)] = shortforce[i]+particlearray.force[i];
	  //	  printf("force[%d] short %e %e %e  PM %e %e %e\n",getatomid(particlearray,i),shortforce[i].x,shortforce[i].y,shortforce[i].z, particlearray.force[i].x,particlearray.force[i].y,particlearray.force[i].z);
	}
      }
      for(int t=1;t<num_pp;t++){
	MPI_Recv(&(atomid[0]), nump+1, MPI_INTEGER, t, t, pp_comm, (MPI_Status *)NULL);
	int n = atomid[0];
	MPI_Recv(&(f[0].x), n*3, MPI_DOUBLE, t, t, pp_comm, (MPI_Status *)NULL);
	for(int i=0;i<n;i++){
	  fall[atomid[i+1]] = f[i];
	}
      }
      print_f(acc_direct, fall, "FeF",cellindex.target_range.x,order);
      std::vector<double> fne(nump);
      std::vector<double> fn(nump);
      ForceArray fe(nump);
      std::vector<double> fer(nump);
      for(int i=0;i<nump;i++){
	fne[i] = acc_direct[i].norm();
	fn[i] = fall[i].norm();
	fe[i] = fall[i]-acc_direct[i];
	fer[i] = fe[i].norm()/fne[i];
      }
      print_err(fne,fn,"NFeF",cellindex.target_range.x,order);
      print_err(fne,fer,"NFeEr",cellindex.target_range.x,order);
      
    }else{
      int n = 0;
      for(int t=0;t<nt;t++){
	for(int i=typerangearray[t].begin;i<typerangearray[t].end;i++,n++){
	  f[n] = shortforce[i]+particlearray.force[i];
	  atomid[n+1] = getatomid(particlearray,i);
	  //	  printf("force[%d] short %e %e %e  PM %e %e %e\n",atomid[n+1],shortforce[i].x,shortforce[i].y,shortforce[i].z, particlearray.force[i].x,particlearray.force[i].y,particlearray.force[i].z);
	}
      }
      atomid[0] = n;
      MPI_Send(&(atomid[0]), n+1, MPI_INTEGER, 0, node_id, pp_comm);
      MPI_Send(&(f[0].x), n*3, MPI_DOUBLE, 0, node_id, pp_comm);
    }
  }
#endif

  MPI_Finalize();
}
