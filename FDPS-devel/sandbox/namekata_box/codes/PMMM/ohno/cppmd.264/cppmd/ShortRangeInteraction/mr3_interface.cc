#include "mr3_interface.h"
#include "mr3.h"

template<class PA, class GPA>
void MR3::calcforce(PA& particlearray,
		    std::vector<TypeRange>& typerangearray,
		    std::vector< std::vector<int> >& self_shorttarget_index,
		    GPA& ghost,
		    std::vector<TypeRange>& ghosttyperange, 
		    std::vector< std::vector<int> >& ghost_shorttarget_index,
		    std::vector< std::pair<int,int> >& ghost_pair_index,
		    ForceArray& shortforce, 
		    ForceArray& ghostshortforce,
		    double& shortenergyself,
		    double& shortenergy,
		    const double cutoff2)
{
  int ti;
  int t;
  int i;

  ni = 0;
  for(ti=0;ti<self_shorttarget_index.size();ti++){
    t = self_shorttarget_index[ti];
    for(i=typerangearray[t].begin;i<typerangearray[t].end;i++){
      xi[ni*3  ] = getpos(particlearray,i).x;
      xi[ni*3+1] = getpos(particlearray,i).y;
      xi[ni*3+2] = getpos(particlearray,i).z;
      qi[ni] = getcharge(particlearray,i);
      atypei[ni] = getatomtype(particlearray,i);
      ni++;
    }
  }
  nj = ni;
  for(ti=0;ti<ghost_shorttarget_index.size();ti++){
    t = ghost_shorttarget_index[ti];
    for(i=ghosttyperange[t].begin;i<ghosttyperange[t].end;i++){
      xi[nj*3  ] = getpos(ghost,i).x;
      xi[nj*3+1] = getpos(ghost,i).y;
      xi[nj*3+2] = getpos(ghost,i).z;
      qi[nj] = getcharge(ghost,i);
      atypei[nj] = getatomtype(ghost,i);
      nj++;
    }
  }

  memset((void *)ei,(int)0,sizeof(double)*ni*3);
  memset((void *)force,(int)0,sizeof(double)*ni*3);
  MR3calccoulomb_ij(ni, xi, qi, ei,
		    nj, xi, qi,
		    1.0, 1, xmax, 3);
  MR3calccoulomb_ij(ni, xi, qi, force,
		    nj, xi, qi,
		    1.0, 0, xmax, 3);
  MR3calcvdw_ij(ni, xi, atypei, ei,
		nj, xi, atypei,
		nat, gscale, rscale, 3, xmax, 1);
  MR3calcvdw_ij(ni, xi, atypei, force,
		nj, xi, atypei,
		nat, gscale, rscale, 2, xmax, 1);

  int li=0;
  for(ti=0;ti<self_shorttarget_index.size();ti++){
    t = self_shorttarget_index[ti];
    for(i=typerangearray[t].begin;i<typerangearray[t].end;i++){
      shortforce[i].x += force[li*3  ] ;
      shortforce[i].y += force[li*3+1] ;
      shortforce[i].z += force[li*3+2] ;
      li++;
    }
  }
  for(i=0;i<ni;i++){
    shortenergy += ei[i*3];
  }
}

template<class PA, class GPA>
void MR3::calccoulomb(PA& particlearray,
		      std::vector<TypeRange>& typerangearray,
		      std::vector< std::vector<int> >& self_shorttarget_index,
		      GPA& ghost,
		      std::vector<TypeRange>& ghosttyperange, 
		      std::vector< std::vector<int> >& ghost_shorttarget_index,
		      std::vector< std::pair<int,int> >& ghost_pair_index,
		      ForceArray& shortforce, 
		      ForceArray& ghostshortforce,
		      double& shortenergyself,
		      double& shortenergy,
		      const double cutoff2)
{
  int ti;
  int t;
  int i;

  ni = 0;
  for(ti=0;ti<self_shorttarget_index.size();ti++){
    t = self_shorttarget_index[ti];
    for(i=typerangearray[t].begin;i<typerangearray[t].end;i++){
      xi[ni*3  ] = getpos(particlearray,i).x;
      xi[ni*3+1] = getpos(particlearray,i).y;
      xi[ni*3+2] = getpos(particlearray,i).z;
      qi[ni] = getcharge(particlearray,i);
      ni++;
    }
  }
  nj = ni;
  for(ti=0;ti<ghost_shorttarget_index.size();ti++){
    t = ghost_shorttarget_index[ti];
    for(i=ghosttyperange[t].begin;i<ghosttyperange[t].end;i++){
      xi[nj*3  ] = getpos(ghost,i).x;
      xi[nj*3+1] = getpos(ghost,i).y;
      xi[nj*3+2] = getpos(ghost,i).z;
      qi[nj] = getcharge(ghost,i);
      nj++;
    }
  }

  memset((void *)ei,(int)0,sizeof(double)*ni*3);
  memset((void *)force,(int)0,sizeof(double)*ni*3);
  MR3calccoulomb_ij(ni, xi, qi, ei,
		    nj, xi, qi,
		    1.0, 1, xmax, 3);
  MR3calccoulomb_ij(ni, xi, qi, force,
		    nj, xi, qi,
		    1.0, 0, xmax, 3);

  int li=0;
  for(ti=0;ti<self_shorttarget_index.size();ti++){
    t = self_shorttarget_index[ti];
    for(i=typerangearray[t].begin;i<typerangearray[t].end;i++){
      shortforce[i].x += force[li*3  ] ;
      shortforce[i].y += force[li*3+1] ;
      shortforce[i].z += force[li*3+2] ;
      li++;
    }
  }
  for(i=0;i<ni;i++){
    shortenergy += ei[i*3];
  }
}

template<class PA, class GPA>
void MR3::calcvdw(PA& particlearray,
		  std::vector<TypeRange>& typerangearray,
		  std::vector< std::vector<int> >& self_shorttarget_index,
		  GPA& ghost,
		  std::vector<TypeRange>& ghosttyperange, 
		  std::vector< std::vector<int> >& ghost_shorttarget_index,
		  std::vector< std::pair<int,int> >& ghost_pair_index,
		  ForceArray& shortforce, 
		  ForceArray& ghostshortforce,
		  double& shortenergyself,
		  double& shortenergy,
		  const double cutoff2)
{
  int ti;
  int t;
  int i;

  ni = 0;
  for(ti=0;ti<self_shorttarget_index.size();ti++){
    t = self_shorttarget_index[ti];
    for(i=typerangearray[t].begin;i<typerangearray[t].end;i++){
      xi[ni*3  ] = getpos(particlearray,i).x;
      xi[ni*3+1] = getpos(particlearray,i).y;
      xi[ni*3+2] = getpos(particlearray,i).z;
      atypei[ni] = getatomtype(particlearray,i);
      ni++;
    }
  }
  nj = ni;
  for(ti=0;ti<ghost_shorttarget_index.size();ti++){
    t = ghost_shorttarget_index[ti];
    for(i=ghosttyperange[t].begin;i<ghosttyperange[t].end;i++){
      xi[nj*3  ] = getpos(ghost,i).x;
      xi[nj*3+1] = getpos(ghost,i).y;
      xi[nj*3+2] = getpos(ghost,i).z;
      atypei[nj] = getatomtype(ghost,i);
      nj++;
    }
  }

  memset((void *)ei,(int)0,sizeof(double)*ni*3);
  memset((void *)force,(int)0,sizeof(double)*ni*3);
  MR3calcvdw_ij(ni, xi, atypei, ei,
		nj, xi, atypei,
		nat, gscale, rscale, 3, xmax, 1);
  MR3calcvdw_ij(ni, xi, atypei, force,
		nj, xi, atypei,
		nat, gscale, rscale, 2, xmax, 1);

  int li=0;
  for(ti=0;ti<self_shorttarget_index.size();ti++){
    t = self_shorttarget_index[ti];
    for(i=typerangearray[t].begin;i<typerangearray[t].end;i++){
      shortforce[i].x += force[li*3  ] ;
      shortforce[i].y += force[li*3+1] ;
      shortforce[i].z += force[li*3+2] ;
      li++;
    }
  }
  for(i=0;i<ni;i++){
    shortenergy += ei[i*3];
  }
}

