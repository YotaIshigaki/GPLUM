#include "samplecell.h"

bool
power_of_two(int n, int &p)
{
  int tmp=n;
  p = 0;
  while(tmp>1){
    tmp = tmp>>1;
    p++;
  }
  int p2 = 1<<p;

  if(p2!=n){
    return false;
  }
  return true;
}

void
split_power_of_two(int p, int sp[3])
{
  sp[0] = p/3;
  sp[1] = p/3;
  sp[2] = p/3;
  if((p%3)==2){
    sp[1]++;
    sp[2]++;
  }else if((p%3)==1){
    sp[2]++;
  }
}

size_t
generate_cell_index(int celldiv, 
                    std::vector< SpaceVector<int> >& cell_position)
{
  cell_position.resize(celldiv*celldiv*celldiv);
  size_t i=0;
  for(int cz=0;cz<celldiv;cz++){
    for(int cy=0;cy<celldiv;cy++){
      for(int cx=0;cx<celldiv;cx++){
        cell_position[i] = SpaceVector<int>(cx,cy,cz);
        i++;
      }
    }
  }
  return i;
}

void
distribute_cell(int div[3], int celldiv, 
                std::vector< SpaceVector<int> >& cell_position,
                std::vector< std::vector<int> >& cell_index_node,
                std::map<int,int>& cell_to_node)
{
  size_t num_node = div[0]*div[1]*div[2];
  cell_index_node.resize(num_node);
  std::cout << "resize " << num_node << std::endl;
  cell_index_node.clear();
  for(size_t n=0;n<num_node;n++){
    std::vector<int> tmp;
    cell_index_node.push_back(tmp);
  }
  std::cout << "clear" << std::endl;
  int cellperx = celldiv/div[0];
  int cellpery = celldiv/div[1];
  int cellperz = celldiv/div[2];
  for(size_t i=0;i<cell_position.size();i++){
    int distx = cell_position[i].x/cellperx;
    int disty = cell_position[i].y/cellpery;
    int distz = cell_position[i].z/cellperz;
    int node = distx + disty*div[0] + distz*div[0]*div[1];
    cell_index_node[node].push_back(i);
    cell_to_node.insert(std::pair<int,int>(i,node));
  }
  /*
  for(size_t n=0;n<num_node;n++){
    std::cout << "node " << n << " have ";
    for(size_t i=0;i<cell_index_node[n].size();i++){
      std::cout << cell_index_node[n][i] << cell_position[cell_index_node[n][i]] << " ";
    }
    std::cout  <<  cell_index_node[n].size() << " cell";
    if(cell_index_node[n].size()!=cell_position.size()/num_node){
      std::cout << " != " << cell_position.size()/num_node; 
    }
    std::cout << std::endl;
  }
  */
}

double
min_cell_distance(SpaceVector<int> ci, SpaceVector<int> cj, 
                  SpaceVector<double> scale)
{
  SpaceVector<int> dc = cj - ci;
  int dx_min = std::min(abs(dc.x),0);
  int dy_min = std::min(abs(dc.y),0);
  int dz_min = std::min(abs(dc.z),0);
  double dx = scale.x*dx_min;
  double dy = scale.y*dy_min;
  double dz = scale.z*dz_min;
  return sqrt(dx*dx+dy*dy+dz*dz);
}


size_t
cubic_target(int celldiv, SpaceVector<int> ci, int reach,
             std::vector<int>& target_cell, 
             ShiftCellArray& shift_cell)
{
  int shiftx, shifty, shiftz;
  target_cell.clear();
  for(int rz=-reach;rz<=reach;rz++){
    int z = ci.z + rz;
    shiftz = 0;
    if(z<0){
      z+=celldiv;
      shiftz = -1;
    }
    if(z>=celldiv){
      z-=celldiv;
      shiftz = 1;
    }
    for(int ry=-reach;ry<=reach;ry++){
      int y = ci.y + ry;
      shifty = 0;
      if(y<0){
        y+=celldiv;
        shifty = -1;
      }
      if(y>=celldiv){
        y-=celldiv;
        shifty = 1;
      }
      for(int rx=-reach;rx<=reach;rx++){
        int x = ci.x + rx;
        shiftx = 0;
        if(x<0){
          x+=celldiv;
          shiftx = -1;
        }
        if(x>=celldiv){
          x-=celldiv;
          shiftx = 1;
        }
        if(!((rz==0)&&(ry==0)&&(rx==0))){
          int cj = x+y*celldiv+z*celldiv*celldiv;
          target_cell.push_back(cj);
          {
            if(shiftx==1) shift_cell.plusx.push_back(cj);
            if(shiftx==-1) shift_cell.minusx.push_back(cj);
            if(shifty==1) shift_cell.plusy.push_back(cj);
            if(shifty==-1) shift_cell.minusy.push_back(cj);
            if(shiftz==1) shift_cell.plusz.push_back(cj);
            if(shiftz==-1) shift_cell.minusz.push_back(cj);
          }
        }
      }
    }
  }
  return target_cell.size();
}

inline
void addshift(std::vector<int>& set, SpaceVector<int> shift,
              std::map< int,SpaceVector<int> >& set_to_shift)
{
  for(std::vector<int>::iterator s=set.begin(); s!=set.end(); ++s){
    std::map< int,SpaceVector<int> >::iterator se = set_to_shift.find((*s));
    if(se!=set_to_shift.end()){
      se->second += shift;
    }else{
      set_to_shift.insert(std::pair<int,SpaceVector<int> >((*s),shift));
    }
  }
}

// When number of node is small, different ghosts of a cell will appear in one node.
// This version only warn such different ghosts and ignore.
void
merge_shift_cell(std::vector<ShiftCellArray>& shift_cell,
                 ShiftCellArray& shift_cell_all)
{
  std::vector< std::map< int,SpaceVector<int> > > set_to_shift;

  set_to_shift.resize(shift_cell.size());
  for(int ci=0;ci<shift_cell.size();ci++){
    addshift(shift_cell[ci].plusx, SpaceVector<int>(1,0,0),set_to_shift[ci]);
    addshift(shift_cell[ci].minusx, SpaceVector<int>(-1,0,0),set_to_shift[ci]);
    addshift(shift_cell[ci].plusy, SpaceVector<int>(0,1,0),set_to_shift[ci]);
    addshift(shift_cell[ci].minusy, SpaceVector<int>(0,-1,0),set_to_shift[ci]);
    addshift(shift_cell[ci].plusz, SpaceVector<int>(0,0,1),set_to_shift[ci]);
    addshift(shift_cell[ci].minusz, SpaceVector<int>(0,0,-1),set_to_shift[ci]);
  }
  /*
  for(int ci=0;ci<shift_cell.size();ci++){
    for(std::map< int,SpaceVector<int> >::iterator s=set_to_shift[ci].begin();
        s!=set_to_shift[ci].end();++s){
      std::cout << " " << s->first << s->second;
    }
    std::cout << std::endl;
  }
  */
  for(int ci=0;ci<shift_cell.size()-1;ci++){
    for(int ci2=ci+1;ci2<shift_cell.size();ci2++){
      for(std::map< int,SpaceVector<int> >::iterator s=set_to_shift[ci].begin();
          s!=set_to_shift[ci].end();++s){
        std::map< int,SpaceVector<int> >::iterator se = set_to_shift[ci2].find((s->first));
        if(se!=set_to_shift[ci2].end()){
          if(se->second!=s->second){
            std::cout << s->first << " " << ci << s->second << " " << ci2 << se->second << std::endl;
          }
        }
      }
    }
  }
  
  std::set<int> tmpset;
  for(int ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].plusx.begin(),shift_cell[ci].plusx.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.plusx.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.plusx.begin());
  
  tmpset.clear();
  for(int ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].minusx.begin(),shift_cell[ci].minusx.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.minusx.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.minusx.begin());

  tmpset.clear();
  for(int ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].plusy.begin(),shift_cell[ci].plusy.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.plusy.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.plusy.begin());

  tmpset.clear();
  for(int ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].minusy.begin(),shift_cell[ci].minusy.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.minusy.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.minusy.begin());

  tmpset.clear();
  for(int ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].plusz.begin(),shift_cell[ci].plusz.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.plusz.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.plusz.begin());

  tmpset.clear();
  for(int ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].minusz.begin(),shift_cell[ci].minusz.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.minusz.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.minusz.begin());
}

int
calc_target_gather(std::map<int,int>& cell_to_node,
                   int mynode, 
                   std::vector< std::vector<int> >& target_cell,
                   std::vector<int>& target_node,
                   std::vector< std::vector<int> >& target_dist
                   )
{
  std::map<int,int> target_to_index;
  std::vector< std::set<int> > td;
  std::set<int> tmptarget;
  int num_target=0;
  for(int i=0;i<target_cell.size();i++){
    for(int j=0;j<target_cell[i].size();j++){
      int tcell = target_cell[i][j];
      int n = cell_to_node[tcell];
      if(n!=mynode){
        std::map<int,int>::iterator ti = target_to_index.find(n);
        if(ti==target_to_index.end()){
          tmptarget.insert(n);
          target_to_index.insert(std::pair<int,int>(n,num_target));
          num_target++;
          std::set<int> newdist;
          newdist.insert(tcell);
          td.push_back(newdist);
        }else{
          td[(*ti).second].insert(tcell);
        }
      }
    }
  }
  for(std::set<int>::iterator t=tmptarget.begin();t!=tmptarget.end();++t){
    int target=(*t);
    int ti = target_to_index[target];
    target_node.push_back(target);
    std::vector<int> dist;
    for(std::set<int>::iterator cell=td[ti].begin();cell!=td[ti].end();++cell){
      dist.push_back((*cell));
    }
    target_dist.push_back(dist);
  }
  return target_node.size();
}

int
calc_target_scatter(std::map<int,int>& cell_to_node,
                    int mynode, std::vector<int> mycell,
                    std::vector< std::vector<int> >& target_cell,
                    std::vector<int>& target_node,
                    std::vector< std::vector<int> >& target_dist
                    )
{
  std::map<int,int> target_to_index;
  std::vector< std::set<int> > td;
  std::set<int> tmptarget;
  int num_target=0;
  for(int i=0;i<target_cell.size();i++){
    int icell = mycell[i];
    for(int j=0;j<target_cell[i].size();j++){
      int tcell = target_cell[i][j];
      int n = cell_to_node[tcell];
      if(n!=mynode){
        std::map<int,int>::iterator ti = target_to_index.find(n);
        if(ti==target_to_index.end()){
          tmptarget.insert(n);
          target_to_index.insert(std::pair<int,int>(n,num_target));
          num_target++;
          std::set<int> newdist;
          newdist.insert(icell);
          td.push_back(newdist);
        }else{
          td[(*ti).second].insert(icell);
        }
      }
    }
  }
  for(std::set<int>::iterator t=tmptarget.begin();t!=tmptarget.end();++t){
    int target=(*t);
    int ti = target_to_index[target];
    target_node.push_back(target);
    std::vector<int> dist;
    for(std::set<int>::iterator cell=td[ti].begin();cell!=td[ti].end();++cell){
      dist.push_back((*cell));
    }
    target_dist.push_back(dist);
  }
  return target_node.size();
}

