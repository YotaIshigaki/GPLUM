/*
 * postprocess_vtk.h
 *
 *  Created on: Nov 28, 2017
 *      Author: jchen@AICS
 */

#ifndef POSTPROCESS_VTK_H_
#define POSTPROCESS_VTK_H_

#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream, std::stringbuf

#include "MPI_File.h"
#include "ByteSwap.h"


////////////////////////////////////////////////////////////////////////////////
namespace NS_VIS
{
template<class T> class VTK;
}
////////////////////////////////////////////////////////////////////////////////
// class NS_VIS::VTK                                                          //
////////////////////////////////////////////////////////////////////////////////

template<class T> class NS_VIS::VTK
{
// ---
 private:
  bool   finalized;
  bool   point_data_tag;
 public:
  int    count;
  double time;
// ---
  const MPI_Comm &my_comm;
  int  my_rank;
  int  num_cpu;
  std::vector<MPI_Offset>  num_node_cum_cpu;
  std::vector<MPI_Offset>  num_node_cum_all;
// ---
  std::vector<T>  node;                           // node/point position
  std::map<std::string, std::vector<T> > NSset;   // scalar data associated to node
  std::map<std::string, std::vector<T> > NVset;   // vector data associated to node; perhaps "node" and "displacement" should be put here
// ---
 public:
  VTK( MPI_Comm comm = MPI_COMM_WORLD )
   : finalized(false), point_data_tag(false), count(0), time(0.0),
     my_comm(comm)
  {
    my_rank  = NS_MPI::GetCommRank(comm);
    num_cpu  = NS_MPI::GetCommSize(comm);
    num_node_cum_cpu   .resize( 1, 0 );
  }
 private:
  void InitializeNumNode();
  void FinalizeMesh();
  void ConvertDataEndian();    // .vtk file handled by Paraview expect Big endien...

 public:
  void AppendNode       ( const PS::ParticleSystem<FP>& sph_system );
  void AppendNodeScalar( const std::string &name, const PS::ParticleSystem<FP>& sph_system );
  void AppendNodeVector( const std::string &name, const PS::ParticleSystem<FP>& sph_system );

// ---
 private:
  void WriteHead        ( NS_MPI::File &file );
  void WriteTime        ( NS_MPI::File &file );
  void WriteNode        ( NS_MPI::File &file );
  void WriteNodeScalar  ( NS_MPI::File &file );
  void WriteNodeVector  ( NS_MPI::File &file );
 public:
  void Write() const{}
  void Write( const std::string &directory, const std::string &prefix );
  void WriteMesh           ( NS_MPI::File &file );
  void ClearData();
};

////////////////////////////////////////////////////////////////////////////////
template<class T>
inline void NS_VIS::VTK<T>::FinalizeMesh()
{
// --- make num_cum_all
  InitializeNumNode();
}

template<class T>
inline void NS_VIS::VTK<T>::InitializeNumNode()
{
  std::vector<MPI_Offset> buf( num_cpu );
  {
    MPI_Datatype type_offset;
    MPI_Type_contiguous( sizeof(MPI_Offset), MPI_BYTE, &type_offset );
    MPI_Type_commit( &type_offset );

    // Gather the # of nodes on each process (cpu) and store them in buf;
    MPI_Allgather( &(num_node_cum_cpu.back()), 1, type_offset, &buf.front(), 1, type_offset, my_comm );

    MPI_Type_free( &type_offset );
  }
  num_node_cum_all.resize( num_cpu + 1 );
  num_node_cum_all[0] = 0;


  for( int i = 0; i < num_cpu; ++i ) {
    num_node_cum_all[i+1] = num_node_cum_all[i] + buf[i];   // let every cpu aware of # of nodes on all cpus
    //  Debug output
//    std::cout << "i, num_node_cum_cpu_all[i], buf[i], num_node_cum_cpu: "
//    		  << i << ", "
//			  << num_node_cum_all[i] << ", "
//			  << buf[i] << ", "
//			  << num_node_cum_cpu.back() << std::endl;
  }
}

template<class T>
inline void NS_VIS::VTK<T>::AppendNode( const PS::ParticleSystem<FP>& sph_system )
{
  MPI_Offset offset( node.size() ), size( sph_system.getNumberOfParticleLocal() );
  node.resize( offset + 3*size );
  typename std::vector<T>::iterator idata( node.begin() + offset );
  for( int iptcl = 0; iptcl < size; iptcl++, idata+=3 )
  {
    *idata     = sph_system[iptcl].pos.x;
    *(idata+1) = sph_system[iptcl].pos.y;
    *(idata+2) = sph_system[iptcl].pos.z;
  }
  num_node_cum_cpu.push_back( num_node_cum_cpu[count] + size );

  //  Debug output
  std::cout << "AppendNode: myrank, count, num_node_cum_cpu.back(), num_node_cum_cpu[count], size " << my_rank << ", " << count << ", "
		    << num_node_cum_cpu.back() << "," << num_node_cum_cpu[count] << ", " <<  size << ", "
		    << std::endl;
  std::cout << "num_node_cum_cpu.size(): " << num_node_cum_cpu.size() << std::endl;
  for (int i = 0; i < num_node_cum_cpu.size(); i++ ) std::cout << num_node_cum_cpu[i] << ", ";
  std::cout << std::endl;
//  std::cin.get();
}

template<class T>
inline void NS_VIS::VTK<T>::AppendNodeScalar( const std::string &name, const PS::ParticleSystem<FP>& sph_system )
{
  std::vector<T> &set( NSset[name] );
  set.reserve( sph_system.getNumberOfParticleLocal() );
  MPI_Offset nmb_node( sph_system.getNumberOfParticleLocal() );

  if ( name == "density" ) {
    for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
      set.push_back( sph_system[nodeID].dens );
  } else if ( name == "pressure" ) {
	  for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
	    set.push_back( sph_system[nodeID].pres );
  } else if ( name == "itype" ) {
	  for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
	    set.push_back( sph_system[nodeID].itype );
  } else if  ( name == "energy" ) {
	  for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
	    set.push_back( sph_system[nodeID].eng );
  } else if  ( name == "shearing_strain" ) {
    for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
      set.push_back( sph_system[nodeID].eps.xy );
  } else if  ( name == "state" ) {
    for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
      set.push_back( sph_system[nodeID].istate );
  } else if  ( name == "stress_xx" ) {
    for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
      set.push_back( sph_system[nodeID].sig.xx );
  } else if  ( name == "stress_xy" ) {
    for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
      set.push_back( sph_system[nodeID].sig.xy );
  } else if  ( name == "epsp_xx" ){
    for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
      set.push_back( sph_system[nodeID].epsp.xx );
  } else if  ( name == "epsp_xy"){
    for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
      set.push_back( sph_system[nodeID].epsp.xy );
  } else if  ( name == "rank"){
    for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
      set.push_back( sph_system[nodeID].prank );
  }
  	else {
	  std::cout << "No such data for particle: " << name << std::endl;
  }

}

template<class T>
inline void NS_VIS::VTK<T>::AppendNodeVector( const std::string &name, const PS::ParticleSystem<FP>& sph_system )
{
  std::vector<T> &set( NVset[name] );
  set.reserve( 3 * sph_system.getNumberOfParticleLocal() );
  MPI_Offset nmb_node( sph_system.getNumberOfParticleLocal() );

  if ( name == "velocity" ) {
	  for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
	  {
	    set.push_back( sph_system[nodeID].vel.x );
	    set.push_back( sph_system[nodeID].vel.y );
	    set.push_back( sph_system[nodeID].vel.z );
	  }
  } else if ( name == "displacement" ) {
	  for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
	  {
	    set.push_back( sph_system[nodeID].disp.x );
	    set.push_back( sph_system[nodeID].disp.y );
	    set.push_back( sph_system[nodeID].disp.z );
	  }
  } else if  ( name == "acceleration" ) {
	  for( int nodeID = 0; nodeID < nmb_node; ++nodeID )
	  {
	    set.push_back( sph_system[nodeID].acc.x );
	    set.push_back( sph_system[nodeID].acc.y );
	    set.push_back( sph_system[nodeID].acc.z );
	  }
  }	else {
	  std::cout << "No such data for particle: " << name << std::endl;
  }

}

template<class T>
inline void displayEndian( int ndata=0, T* data=NULL )
{
    unsigned char *p;
    for (int j = 0; j<ndata; j++ ) {
       p = (unsigned char*) (data+j);
	   for (int i = 0; i != sizeof(T); ++i)  printf("%02X ", p[i]);
	   std::cout << std::endl;
    }
}


template<class T>
inline void NS_VIS::VTK<T>::ConvertDataEndian()
{

  SwapBytes<T>( node.size(), &node[0] );

  for( typename std::map<std::string, std::vector<T> >::iterator it = NSset.begin(); it != NSset.end(); ++it )
  {
	  SwapBytes<T>( it->second.size(), &it->second[0] );
  }

  for( typename std::map<std::string, std::vector<T> >::iterator it = NVset.begin(); it != NVset.end(); ++it )
    SwapBytes<T>( it->second.size(), &it->second[0] );
}

template<class T>
inline void NS_VIS::VTK<T>::WriteHead( NS_MPI::File &file )
{
  MPI_Barrier(my_comm); //Wait for the others, if the computation load is not well balanced
  std::stringstream out;
  out << "# vtk DataFile Version 3.0\n";
  out << "visualization file\n";
  out << "BINARY\n";
  out << "DATASET UNSTRUCTURED_GRID\n";
  file.Write( out.str() );
}


template<class T>
inline void NS_VIS::VTK<T>::WriteTime( NS_MPI::File &file )
{
// --- head
  std::stringstream out;
  out << "FIELD FieldData 1" << '\n';
  out << "TIME " << 1 << " " << 1 << " double" << '\n';
  file.Write( out.str() );
  file.Write( 1, 1, &time );
}

template<class T>
inline void NS_VIS::VTK<T>::WriteNode( NS_MPI::File &file )
{
// --- head
  std::stringstream out;
  if ( typeid(T) == typeid(double) ) out << "\nPOINTS " << num_node_cum_all[num_cpu] << " double" << '\n';
  if ( typeid(T) == typeid(float) )  out << "\nPOINTS " << num_node_cum_all[num_cpu] << " float" << '\n';
  file.Write( out.str() );

// --- body
  file.CollectiveWrite( 3, num_node_cum_all, &node[0] );

//  if (my_rank == 0) {
//	  for (int i =0; i<num_node_cum_all.size(); i++ )
//	    std::cout << "i, num_node_cum_all[i]: " << i << ", " << num_node_cum_all[i] << std::endl;
//	    std::cin.get();
//  }
}

template<class T>
inline void NS_VIS::VTK<T>::WriteMesh( NS_MPI::File &file )
{
  if( ! finalized ){ FinalizeMesh(); finalized = true; }
  WriteHead       ( file );
  //std::cout << "WriteHead finished!" << std::endl;
  WriteTime       ( file );
  //std::cout << "WriteTime finished!" << std::endl;
  WriteNode       ( file );
  //std::cout << "WriteNode finished!" << std::endl;
}


template<class T>
inline void NS_VIS::VTK<T>::WriteNodeScalar( NS_MPI::File &file )
{
  for( typename std::map<std::string, std::vector<T> >::const_iterator ci = NSset.begin(); ci != NSset.end(); ++ci )
  {
// --- head
    std::stringstream out;
    if (!point_data_tag) {
      out << "\nPOINT_DATA " << num_node_cum_all.back() << '\n';
      point_data_tag = true;
    }
    else out << '\n';
    if ( typeid(T) == typeid(double) )  out << "SCALARS " << ci->first << " double" << '\n';
    if ( typeid(T) == typeid(float) )   out << "SCALARS " << ci->first << " float" << '\n';
    out << "LOOKUP_TABLE default" << '\n';
    file.Write( out.str() );

    std::cout << out.str() << ci->second.size() << std::endl;
    //std::cin.get();

// --- body
    file.CollectiveWrite( 1, num_node_cum_all, (T *) &(ci->second[0]) );

  }
}

template<class T>
inline void NS_VIS::VTK<T>::WriteNodeVector( NS_MPI::File &file )
{
  for( typename std::map<std::string, std::vector<T> >::const_iterator ci = NVset.begin(); ci != NVset.end(); ++ci )
  {
// --- head
    std::stringstream out;
    if (!point_data_tag) {
      out << "\nPOINT_DATA " << num_node_cum_all.back() << '\n';
      point_data_tag = true;
    }
    else out << '\n';
    if ( typeid(T) == typeid(double) )  out << "VECTORS " << ci->first << " double\n";
    if ( typeid(T) == typeid(float) ) out << "VECTORS " << ci->first << " float\n";
    file.Write( out.str() );

// --- body
    file.CollectiveWrite( 3, num_node_cum_all, (T *)&(ci->second[0]) );
  }
}

inline std::string AppendFileSeparator( const std::string &directory_ )
{
  std::string directory( directory_ );
  if( directory.empty() ) directory.assign("./");
  char last( directory[directory.size() - 1] );
  if( last != '/' && last != '\\' )
  {
    if( directory.find("\\") == std::string::npos )
      directory.append("/");
    else
      directory.append("\\");
  }
  return directory;
}

template<class T>
inline void NS_VIS::VTK<T>::Write( const std::string &dirname, const std::string &prefix )
{

#if BYTE_ORDER == LITTLE_ENDIAN
  ConvertDataEndian();
  if(my_rank == 0) std::cout << "LITTLE_ENDIAN: ConvertDataEndian" << std::endl;
#endif
  NS_MPI::File file;
  file.Open( AppendFileSeparator(dirname) + prefix + ".vtk", false, my_comm );

  //MPI_Barrier(my_comm);
  WriteMesh( file );
  //std::cout << "WriteMesh finished!" << std::endl;
  //std::cin.get();

  WriteNodeScalar( file );
  //std::cout << "WriteNodeScalar finished!" << std::endl;
  //std::cin.get();

  WriteNodeVector( file );
  //std::cout << "WriteNodeVector finished!" << std::endl;
  //std::cin.get();

  file.Close();
}

template<class T>
inline void NS_VIS::VTK<T>::ClearData()
{
  NSset.clear();
  NVset.clear();
}



#endif /* POSTPROCESS_VTK_H_ */
