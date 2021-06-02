/*
 * File:   MPI_Datatypes.h
 * Ver1:   by Lalith @ U. of Tokyo
 * Ver2:   Copied from IES_SRA and simplified for FDPS-SPH by J. Chen @AICS Nov 29, 2017
 */

#ifndef MPI_DATATYPES_H
#define	MPI_DATATYPES_H

#include <mpi.h>
#include <iostream>
#include <vector>
#include <utility>


//include only the frequently used data types to NS_MPI. Keep the other in-frequent types local to the related data structure

   /******************************************************************************
    * Show an error message and abort MPI                                        *
    ******************************************************************************/
   template <typename T> void MPI_Err_Abort(T err) {
      std::cerr << "Error : " << err << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   // Templates for finding equivelent MPI_Datatype of POD(Plain Old Data) and
   // some C++ STL data types
   //specialize mpi_type for POD types

   template<class T> struct mpi_type {
      static MPI_Datatype get() {
         return MPI_DATATYPE_NULL;
      };
   };

   template <>MPI_Datatype  mpi_type <char>::get()           { return MPI_CHAR; }
   template <>MPI_Datatype  mpi_type <short>::get()          { return MPI_SHORT; }
   template <>MPI_Datatype  mpi_type <int>::get()            { return MPI_INT;}
   //template <>MPI_Datatype  mpi_type <int8_t>::get()         { return MPI_INT8_T;}
   template <>MPI_Datatype  mpi_type <long>::get()           { return MPI_LONG;}
   template <>MPI_Datatype  mpi_type <signed char>::get()    { return MPI_CHAR;}
   template <>MPI_Datatype  mpi_type <unsigned char>::get()  { return MPI_UNSIGNED_CHAR;}
   template <>MPI_Datatype  mpi_type <unsigned short>::get() { return MPI_UNSIGNED_SHORT;}
   template <>MPI_Datatype  mpi_type <unsigned int>::get()   { return MPI_UNSIGNED;}
   template <>MPI_Datatype  mpi_type <unsigned long>::get()  { return MPI_UNSIGNED_LONG;}

   template <>MPI_Datatype  mpi_type <float>::get()          { return MPI_FLOAT;}
   template <>MPI_Datatype  mpi_type <double>::get()         { return MPI_DOUBLE;}
   template <>MPI_Datatype  mpi_type <long double>::get()    { return MPI_LONG_DOUBLE;}

   //some useful functions
   inline bool Is_Predefined_MPIType(const MPI_Datatype &mpitype);
   inline bool Is_UserDefined(const MPI_Datatype &mpitype);

   //use the following functions to Free MPI_Datatypes. These can prevent user defined
   //static MPI_Datatype being freed
   inline void Free_MPIType(MPI_Datatype &mpitype);
   template<template<class, class> class STLC>void Free_MPIType(STLC<MPI_Datatype, std::allocator<MPI_Datatype> > &mpi_types);
   inline void Free_MPIType(const unsigned int N, MPI_Datatype mpi_types[]);

   //not yet being used
   inline int Decode_MPI_Type(MPI_Datatype &datatype);


/*******************************************************************************
     * Uncomment part with FORTRAN data types, if this is to be used with FORTRAN
     * MPI_Datatype
     *******************************************************************************/
    inline bool Is_Predefined_MPIType(const MPI_Datatype &mpitype) {


        if (mpitype == MPI_DATATYPE_NULL) return true;

        //Data types for C language bindings
        if (mpitype == MPI_INT) return true; //32-bit integer
        if (mpitype == MPI_DOUBLE) return true; //64-bit floating point
        if (mpitype == MPI_FLOAT) return true; //32-bit floating point
        if (mpitype == MPI_CHAR) return true; //8-bit character
        if (mpitype == MPI_LONG) return true; //32-bit integer
        if (mpitype == MPI_LONG_DOUBLE) return true; //64-bit floating point
        if (mpitype == MPI_LONG_LONG) return true; //64-bit integer
        if (mpitype == MPI_LONG_LONG_INT) return true; //64-bit integer
        if (mpitype == MPI_SHORT) return true; //16-bit integer
        if (mpitype == MPI_SIGNED_CHAR) return true; //8-bit signed character
        if (mpitype == MPI_UNSIGNED) return true; //32-bit unsigned integer
        if (mpitype == MPI_UNSIGNED_CHAR) return true; //8-bit unsigned character
        if (mpitype == MPI_UNSIGNED_LONG) return true; //32-bit unsigned integer
        if (mpitype == MPI_UNSIGNED_LONG_LONG) return true; //64-bit unsigned integer
        if (mpitype == MPI_UNSIGNED_SHORT) return true; //16-bit unsigned integer
        if (mpitype == MPI_WCHAR) return true; //Wide (16-bit) unsigned character

        //Data types for reduction functions (Fortran reduction types)
        /*
        if (mpitype == MPI_2COMPLEX) return true; //MPI_COMPLEX, MPI_COMPLEX}
        if (mpitype == MPI_2DOUBLE_PRECISION) return true; //{MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION}
        if (mpitype == MPI_2INTEGER) return true; //{MPI_INTEGER, MPI_INTEGER}
        if (mpitype == MPI_2REAL) return true; //{MPI_REAL, MPI_REAL}
        */

        //Special purpose data types
        if (mpitype == MPI_BYTE) return true; //Untyped byte data
//        if (mpitype == MPI_LB) return true; //Explicit lower bound marker
        if (mpitype == MPI_PACKED) return true; //Packed data (byte)
//        if (mpitype == MPI_UB) return true; //Explicit upper bound marker


        //Data types for Fortran language bindings
        /*
        if(mpitype == MPI_CHARACTER) return true; //8-bit character
        if(mpitype == MPI_COMPLEX) return true; //32-bit floating point real, 32-bit floating point imaginary
        if(mpitype == MPI_COMPLEX8) return true; //32-bit floating point real, 32-bit floating point imaginary
        if(mpitype == MPI_COMPLEX16) return true; //64-bit floating point real, 64-bit floating point imaginary
        if(mpitype == MPI_COMPLEX32) return true; //128-bit floating point real, 128-bit floating point imaginary
        if(mpitype == MPI_DOUBLE_COMPLEX) return true; //64-bit floating point real, 64-bit floating point imaginary
        if(mpitype == MPI_DOUBLE_PRECISION) return true; //64-bit floating point
        if(mpitype == MPI_INTEGER) return true; //32-bit integer
        if(mpitype == MPI_INTEGER1) return true; //8-bit integer
        if(mpitype == MPI_INTEGER2) return true; //16-bit integer
        if(mpitype == MPI_INTEGER4) return true; //32-bit integer
        if(mpitype == MPI_INTEGER8) return true; //64-bit integer
        if(mpitype == MPI_LOGICAL) return true; //32-bit logical
        if(mpitype == MPI_LOGICAL1) return true; //8-bit logical
        if(mpitype == MPI_LOGICAL2) return true; //16-bit logical
        if(mpitype == MPI_LOGICAL4) return true; //32-bit logical
        if(mpitype == MPI_LOGICAL8) return true; //64-bit logical
        if(mpitype == MPI_REAL) return true; //32-bit floating point
        if(mpitype == MPI_REAL4) return true; //32-bit floating point
        if(mpitype == MPI_REAL8) return true; //64-bit floating point
        if(mpitype == MPI_REAL16) return true; //128-bit floating point

        //Data types for reduction functions (C reduction types)
        if(mpitype == MPI_DOUBLE_INT) return true; //{MPI_DOUBLE, MPI_INT}
        if(mpitype == MPI_FLOAT_INT) return true; //{MPI_FLOAT, MPI_INT}
        if(mpitype == MPI_LONG_DOUBLE_INT) return true; //{MPI_LONG_DOUBLE, MPI_INT}
        if(mpitype == MPI_LONG_INT) return true; //{MPI_LONG, MPI_INT}
        if(mpitype == MPI_SHORT_INT) return true; //{MPI_SHORT, MPI_INT}
        if(mpitype == MPI_2INT) return true; //{MPI_INT, MPI_INT}
         */

        return false;
    }

    /******************************************************************************
     * test whether a MPI_Type is predefined or user defined
     * Same as  !Is_Predefined_MPIType()
    ******************************************************************************/
   inline bool Is_UserDefined(const MPI_Datatype &mpitype) {
      int num_ints, num_adds, num_dtypes, combiner;

      if (mpitype == MPI_DATATYPE_NULL) return false;

      MPI_Type_get_envelope(mpitype, &num_ints, &num_adds, &num_dtypes, &combiner);

      if (combiner == MPI_COMBINER_NAMED) return false;

      return false;
   }

   /******************************************************************************
    * Free a MPI type, if it is user defined
    ******************************************************************************/
   inline void Free_MPIType(MPI_Datatype &mpitype) {
      if (!Is_Predefined_MPIType(mpitype)) {
	  MPI_Type_free(&mpitype);
	  mpitype=MPI_DATATYPE_NULL;
      }
   }


   /******************************************************************************
    * Free a STL list or vector of user defined MPI_Data_type
    ******************************************************************************/
   template<template<class, class> class STLC>void Free_MPIType(STLC<MPI_Datatype, std::allocator<MPI_Datatype> > &mpi_types) {

      for (unsigned int i = 0; i < mpi_types.size(); ++i)
         if (!Is_Predefined_MPIType(mpi_types[i])) {
	     MPI_Type_free(&mpi_types[i]);
	     mpi_types[i]=MPI_DATATYPE_NULL;
	 }
      //if (mpi_types[i] != MPI_DATATYPE_NULL) MPI_Type_free(&mpi_types[i]);

      mpi_types.resize(0);
   }

   /******************************************************************************
    * Free an array of user defined MPI_Data_type
    * N = array size
    ******************************************************************************/
   inline void Free_MPIType(const unsigned int N, MPI_Datatype mpi_types[]) {

      for (unsigned int i = 0; i < N; ++i)
         if (!Is_Predefined_MPIType(mpi_types[i])){
	     MPI_Type_free(&mpi_types[i]);
	     mpi_types[i]=MPI_DATATYPE_NULL;
	 }
   }

   /******************************************************************************
    * Free a STL list or vector of pair<int,MPI_Datatype> where int is the
    * CPU id and MPI_Data_type is the user defind MPI data type
    ******************************************************************************/
   template<template<class, class> class STLC>void Free_MPIType(STLC< std::pair<int, MPI_Datatype>, std::allocator<std::pair<int, MPI_Datatype> > > &mpi_types) {

      for (unsigned int i = 0; i < mpi_types.size(); ++i)
         Free_MPIType(mpi_types[i].second);

      mpi_types.resize(0);
   }

  /******************************************************************************
    * decode a MPI_datatype and free the input MPI_datatype including the nested
    * This cannot be used to free user defined data types, sine the function
    * MPI_Type_get_contents() return a copy to a user defined data types, which
    * itself should be freed
    * This function not complete
    ******************************************************************************/
   inline int Decode_MPI_Type(MPI_Datatype &datatype) {

      std::vector<int> ints;
      std::vector<MPI_Aint> adds;
      std::vector<MPI_Datatype> dtypes;
      int num_ints, num_adds, num_dtypes, combiner;

      MPI_Type_get_envelope(datatype, &num_ints, &num_adds, &num_dtypes, &combiner);

      ints.resize(num_ints);
      adds.resize(num_adds);
      dtypes.resize(num_dtypes);
      MPI_Type_get_contents(datatype, num_ints, num_adds, num_dtypes, &ints[0], &adds[0], &dtypes[0]);

      switch (combiner) {
         case MPI_COMBINER_NAMED: //a named predened datatype
            if (datatype == MPI_INT) std::cout << "MPI_INT\n";
            else if (datatype == MPI_DOUBLE) std::cout << "MPI_DOUBLE\n";
            return 0;
//            break; //does not make sense

         case MPI_COMBINER_STRUCT: //MPI_TYPE_CREATE_STRUCT
         case MPI_COMBINER_STRUCT_INTEGER:
            for (int i = 0; i < ints[0]; ++i)
               if (Decode_MPI_Type(dtypes[i])) MPI_Type_free(&dtypes[i]);
            break;
         case MPI_COMBINER_DUP: //MPI_TYPE_DUP

            break;
//         case MPI_COMBINER_CONTIGUOUS: //MPI_TYPE_CONTIGUOUS
//         case MPI_COMBINER_VECTOR: //MPI_TYPE_VECTOR
//         case MPI_COMBINER_HVECTOR: //MPI_TYPE_CREATE_HVECTOR
//         case MPI_COMBINER_INDEXED: //MPI_TYPE_INDEXED
//         case MPI_COMBINER_HINDEXED: //MPI_TYPE_CREATE_HINDEXED
//         case MPI_COMBINER_INDEXED_BLOCK: //MPI_TYPE_CREATE_INDEXED_BLOCK
//         case MPI_COMBINER_HINDEXED_BLOCK: //MPI_TYPE_CREATE_HINDEXED_BLOCK
//         case MPI_COMBINER_SUBARRAY: //MPI_TYPE_CREATE_SUBARRAY
//         case MPI_COMBINER_DARRAY: //MPI_TYPE_CREATE_DARRAY
//         case MPI_COMBINER_F90_REAL: //MPI_TYPE_CREATE_F90_REAL
//         case MPI_COMBINER_F90_COMPLEX: //MPI_TYPE_CREATE_F90_COMPLEX
//         case MPI_COMBINER_F90_INTEGER: //MPI_TYPE_CREATE_F90_INTEGER
//         case MPI_COMBINER_RESIZED: //MPI_TYPE_CREATE_RESIZED

         default:
            MPI_Err_Abort("Unrecognized MPI_Data combiner type\n");
      }
      return 1;

   }


   /******************************************************************************
    * A class which makes the creation of MPI_Type_create_struct simple
    ******************************************************************************/
   template <typename T, const unsigned int N = 1 > struct mpi_struct_basic {
      T dumy[2];
      MPI_Datatype type[N], mpitype;
      int blocklen[N];
      unsigned int p;
      MPI_Aint disp[N], origin, extent;

      //virtual void set()=0;

      mpi_struct_basic() : mpitype(MPI_DATATYPE_NULL), p(0) { }

      MPI_Datatype& get_mpi_type() {
         MPI_Get_address(&dumy[0], &origin);
         for (int i = 0; i < N; ++i) disp[i] -= origin;

         MPI_Datatype tmp_mpi_str;
         MPI_Type_create_struct(N, blocklen, disp, type, &tmp_mpi_str); // build datatype describing structure
         MPI_Get_address(&dumy[1], &extent);
         extent -= origin;
         MPI_Type_create_resized(tmp_mpi_str, (MPI_Aint) 0, extent, &mpitype); // Take care of padding
         MPI_Type_commit(&mpitype); // Commit the data type

	 //free any user defined MPI_datatypes
         Free_MPIType(N, type);

	 return mpitype;
      }

      //

      template <class U> void include(U* begin, const int &size) {
         assert(p < N);
         type[p] = mpi_type<U>::get();
         blocklen[p] = size;
         MPI_Get_address(begin, &disp[p++]);
      }

   };



   ////////////////////////////////////////////////////////////////////////////////
   //   User defined frequently used data types
   ////////////////////////////////////////////////////////////////////////////////

   //modify these later such that the MPI_Datatype is static. Then care has to be
   //taken not to free these MPI_Datatype's when calling MPI_Type_free. Adding
   //these to Is_Predefined_MPIType is a good way to prevent these user defined data types free,

//
//   /******************************************************************************
//    * Functors for generating MPI_Datatypes of IntVector2/3Ds
//    ******************************************************************************/
//   template < > struct mpi_type < IntVector2D > {
//
//      static MPI_Datatype get() {
//         mpi_struct_basic< IntVector2D, 2> sd; //structure data
//
//         sd.include(&sd.dumy[0].i, 1);
//         sd.include(&sd.dumy[0].j, 1);
//
//         return sd.get_mpi_type();
//      }
//   };
//
//   /******************************************************************************
//    *
//    ******************************************************************************/
//   template < > struct mpi_type < IntVector3D > {
//
//      static MPI_Datatype get() {
//         mpi_struct_basic< IntVector3D, 3> sd; //structure data
//
//         sd.include(&sd.dumy[0].i, 1);
//         sd.include(&sd.dumy[0].j, 1);
//         sd.include(&sd.dumy[0].k, 1);
//
//         return sd.get_mpi_type();
//      }
//   };
//
//   /******************************************************************************
//    *
//    ******************************************************************************/
//   template < > struct mpi_type < Vector2D > {
//
//      static MPI_Datatype get() {
//         mpi_struct_basic< Vector2D, 2> sd; //structure data
//
//         sd.include(&sd.dumy[0].x, 1);
//         sd.include(&sd.dumy[0].y, 1);
//
//         return sd.get_mpi_type();
//      }
//   };
//
//   /******************************************************************************
//    *
//    ******************************************************************************/
//   template < > struct mpi_type < Vector3D > {
//
//      static MPI_Datatype get() {
//         mpi_struct_basic< Vector3D, 3> sd; //structure data
//
//         sd.include(&sd.dumy[0].x, 1);
//         sd.include(&sd.dumy[0].y, 1);
//         sd.include(&sd.dumy[0].z, 1);
//
//         return sd.get_mpi_type();
//      }
//   };
//
//
//   ////////////////////////////////////////////////////////////////////////////////
//   //   std::pair<U,V>
//   ////////////////////////////////////////////////////////////////////////////////
//
//   template <typename U, typename V> struct mpi_type <std::pair<U, V> > {
//
//      static MPI_Datatype get() {
//
//         mpi_struct_basic< std::pair<U, V>, 2> sd; //structure data
//
//         sd.type[0] = mpi_type<U>::get();
//         sd.blocklen[0] = 1;
//         MPI_Get_address(&sd.dumy[0].first, &sd.disp[0]);
//         sd.type[1] = mpi_type<V>::get();
//         sd.blocklen[1] = 1;
//         MPI_Get_address(&sd.dumy[0].second, &sd.disp[1]);
//
//         return sd.get_mpi_type();
//      }
//   };
//
//}


#endif	/* MPI_DATATYPES_H */

