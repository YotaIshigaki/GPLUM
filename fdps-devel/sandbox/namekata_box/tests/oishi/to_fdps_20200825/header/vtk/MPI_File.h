/*
 * File:   MPI_File.h
 * Ver1:   by Lalith @ U. of Tokyo
 * Ver2:   Copied from IES_SRA and simplified for FDPS-SPH by J. Chen @AICS Nov 29, 2017
 */
#ifndef NS_MPI_FILE_H
#define NS_MPI_FILE_H
////////////////////////////////////////////////////////////////////////////////
#include <mpi.h>
#include "MPI_Datatypes.h"

namespace NS_MPI
{
  class File;
  inline int GetCommRank(const MPI_Comm &comm = MPI_COMM_WORLD) {
     int my_rank;  MPI_Comm_rank(comm, &my_rank); return my_rank;
  }
  inline int GetCommSize(const MPI_Comm &comm = MPI_COMM_WORLD) {
     int no_cpus;   MPI_Comm_size(comm, &no_cpus); return no_cpus;
  }
}

////////////////////////////////////////////////////////////////////////////////
//   NS_MPI::File                                                             //
////////////////////////////////////////////////////////////////////////////////

class NS_MPI::File
{
 private:
  File( const File & );            // copy is prohibited
  File operator=( const File & );  // substitution is prohibited
 protected:
  MPI_File fh;
  MPI_Offset disp;                // offset for mpi-file output
  MPI_Comm comm;
  std::string directory;
  std::string name;
  int num_cpu;
  int my_rank;
 public:
  File() : disp(0) {}
  File( std::string name_ ) : disp(0), name(name_) {}
  File( std::string directory_, std::string name_ ) : disp(0), directory(directory_), name(name_) {}
  void Open( bool append = false, const MPI_Comm &comm = MPI_COMM_WORLD );
  void Open( const std::string &filename, bool append = false, const MPI_Comm &comm = MPI_COMM_WORLD );
  void Open( const char *       filename, bool append = false, const MPI_Comm &comm = MPI_COMM_WORLD );
  void Close();
  template< class T >
  void Write( MPI_Offset size_chunk, MPI_Offset num_chunk, T *data )
  {
    MPI_Datatype type;
    MPI_Type_contiguous( size_chunk, mpi_type<T>::get(), &type );
    MPI_Type_commit( &type );
    
    if( my_rank == 0 )
    {
      MPI_File_seek( fh, disp, MPI_SEEK_SET );
      MPI_File_write( fh, data, num_chunk, type, MPI_STATUS_IGNORE );
    }
    disp += size_chunk * num_chunk * sizeof(T);
    
    MPI_Type_free( &type );
  }
  void Write( std::string data )
  {
    Write( data.size(), 1, &data[0] );
  }

  template< class T >
  void CollectiveWrite( int size_chunk, std::vector<MPI_Offset> num_chunk_cum, MPI_Offset num_segment, T *data )
  {
    char data_rep[] = "native"; //external32";
    MPI_Offset fsize;
    MPI_Datatype ftype;
    MPI_Aint extent_in_file;
    
    int globalarrsize = num_chunk_cum[num_cpu] * sizeof(T);
    int start = num_chunk_cum[my_rank] * sizeof(T);
    int localarrsize = (num_chunk_cum[my_rank+1] - num_chunk_cum[my_rank]) * sizeof(T);
    MPI_Type_create_subarray( 1, &globalarrsize, &localarrsize, &start, MPI_ORDER_C, MPI_BYTE, &ftype);
    MPI_Type_commit( &ftype );
    
    // only for preallocation
    // view of memory and file are different, especially when external32 is used
    // use extent_in_file to calculate disp and other byte lengths for portability
    MPI_File_set_view( fh, disp, MPI_BYTE, MPI_BYTE, data_rep, MPI_INFO_NULL );
    MPI_File_get_type_extent( fh, MPI_BYTE, &extent_in_file );
    fsize = disp + extent_in_file * globalarrsize;
    
    MPI_File_preallocate( fh, fsize );
    
    // set file view and collective write
    MPI_File_set_view( fh, disp, MPI_BYTE, ftype, data_rep, MPI_INFO_NULL );
    MPI_File_write_all( fh, data, num_segment * localarrsize, MPI_BYTE, MPI_STATUS_IGNORE );
    disp = fsize * num_segment;
    
    MPI_Type_free( &ftype );
  }
  template< class T >
  void CollectiveWrite( int size_chunk, std::vector<MPI_Offset> num_chunk_cum, MPI_Offset num_segment, const std::vector< std::vector<T> > &dataset )
  {
    char data_rep[] = "native"; //external32";
    MPI_Offset fsize;
    MPI_Datatype ftype;
    MPI_Aint extent_in_file;
    
    int globalarrsize = num_chunk_cum[num_cpu] * sizeof(T);
    int start = num_chunk_cum[my_rank] * sizeof(T);
    int localarrsize = (num_chunk_cum[my_rank+1] - num_chunk_cum[my_rank]) * sizeof(T);
    MPI_Type_create_subarray( 1, &globalarrsize, &localarrsize, &start, MPI_ORDER_C, MPI_BYTE, &ftype);
    MPI_Type_commit( &ftype );
    
    // only for preallocation
    // view of memory and file are different, especially when external32 is used
    // use extent_in_file to calculate disp and other byte lengths for portability
    MPI_File_set_view( fh, disp, MPI_BYTE, MPI_BYTE, data_rep, MPI_INFO_NULL );
    MPI_File_get_type_extent( fh, MPI_BYTE, &extent_in_file );
    fsize = disp + extent_in_file * globalarrsize;
    
    MPI_File_preallocate( fh, fsize );
    
    // set file view and collective write
    MPI_File_set_view( fh, disp, MPI_BYTE, ftype, data_rep, MPI_INFO_NULL );
    for( MPI_Offset i = 0; i < num_segment; ++i )
      MPI_File_write_all( fh, (T*)&dataset[i][0], localarrsize, MPI_BYTE, MPI_STATUS_IGNORE );
    disp = fsize * num_segment;
    
    MPI_Type_free( &ftype );
  }
  template< class T >
  void CollectiveWrite( MPI_Offset size_chunk, MPI_Offset num_chunk, T *data )
  {
    std::vector<MPI_Offset>  num_chunk_cum( num_cpu + 1, 0 );
    std::vector<MPI_Offset>  buf( num_cpu );
    {
      MPI_Datatype type_offset;
      MPI_Type_contiguous( sizeof(MPI_Offset), MPI_BYTE, &type_offset );
      MPI_Type_commit( &type_offset );
      MPI_Allgather( &num_chunk, 1, type_offset, &buf.front(), 1, type_offset, comm );
      MPI_Type_free( &type_offset );
    }
    for( int i = 0 ; i < num_cpu; ++i )
      num_chunk_cum[i + 1] = num_chunk_cum[i] + buf[i];
    
    MPI_Datatype type;
    MPI_Type_contiguous( size_chunk, mpi_type<T>::get(), &type );
    MPI_Type_commit( &type );
    
    MPI_File_seek( fh, disp + sizeof(T) * size_chunk * num_chunk_cum[my_rank], MPI_SEEK_SET );
    MPI_File_write_all( fh, data, num_chunk, type, MPI_STATUS_IGNORE );
    disp += sizeof(T) * size_chunk * num_chunk_cum[num_cpu];
    
    MPI_Type_free( &type );
  }
  template< class T >
  void CollectiveWrite( MPI_Offset size_chunk, const std::vector<MPI_Offset> &num_chunk_cum, T *data )  //mostly used
  {
    MPI_Offset num_chunk = num_chunk_cum[my_rank+1] - num_chunk_cum[my_rank];
    
    MPI_Datatype type;
    MPI_Type_contiguous( size_chunk, mpi_type<T>::get(), &type );
    MPI_Type_commit( &type );
    
    MPI_File_seek( fh, disp + sizeof(T) * size_chunk * num_chunk_cum[my_rank], MPI_SEEK_SET );
    MPI_File_write_all( fh, data, num_chunk, type, MPI_STATUS_IGNORE );
    disp += sizeof(T) * size_chunk * num_chunk_cum[num_cpu];

    MPI_Type_free( &type );
  }
};

inline void NS_MPI::File::Open( const std::string &filename, bool append, const MPI_Comm &comm_ )
{
  Open( filename.c_str(), append, comm_ );
}

inline void NS_MPI::File::Open( const char *filename, bool append, const MPI_Comm &comm_ )
{
  name = filename;
  comm = comm_;
  num_cpu = NS_MPI::GetCommSize( comm );
  my_rank = NS_MPI::GetCommRank( comm );
  
  MPI_File_open( comm, (char *)name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh );
  MPI_File_set_errhandler( fh, MPI_ERRORS_ARE_FATAL );
  if( fh == MPI_FILE_NULL )
  {
    if( NS_MPI::GetCommRank(comm) == 0 )
      std::cerr << "Faile to open " << name << std::endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  if( append )
    MPI_File_get_size( fh, &disp );
  else
    MPI_File_set_size( fh, 0 );
}

inline void NS_MPI::File::Open ( bool append, const MPI_Comm &comm_ )
{
  Open( name, append, comm_ );
}

inline void NS_MPI::File::Close()
{
  MPI_File_sync  ( fh );
  MPI_File_close ( &fh );
}

////////////////////////////////////////////////////////////////////////////////
#endif
