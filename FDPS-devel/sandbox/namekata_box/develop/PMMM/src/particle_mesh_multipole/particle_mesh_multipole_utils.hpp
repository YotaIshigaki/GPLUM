#pragma once
#include "../ps_macro_defs.h"
#include "../ps_defs.hpp"

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

        class ParticleMeshMultipoleParameters {
        public:
            S32 p; // the maximum order of multipole expansion
            S32 icut; // the minimum cell separation 
            S32 bc; // the boundary condition 
            F64ort pos_unit_cell; // the size of unit cell (real image of the computational box).
            S32vec n_cell; // the number of particle mesh cells
            S32 nmax; // this parameter determines how far virtual image boxes
                      // are to be considered in the calculation of the Ewald
                      // summation in the real space.
            S32 mmax; // this parameter determines how large wavenumber vectors
                      // are to be considered in the calculation of the Ewald
                      // summation in the wavenumber space.
            F64 alpha; // A parameter of Ewald summation. 

            ParticleMeshMultipoleParameters() {
                p = FDPS_DFLT_VAL_PMMM_P;
                icut = FDPS_DFLT_VAL_PMMM_ICUT;
                bc = BOUNDARY_CONDITION_PERIODIC_XYZ;
                pos_unit_cell.init();
                n_cell = S32vec(FDPS_DFLT_VAL_PMMM_N_CELL_1D,
                                FDPS_DFLT_VAL_PMMM_N_CELL_1D,
                                FDPS_DFLT_VAL_PMMM_N_CELL_1D);
                nmax = FDPS_DFLT_VAL_PMMM_NMAX;
                mmax = FDPS_DFLT_VAL_PMMM_MMAX;
                alpha = FDPS_DFLT_VAL_PMMM_ALPHA;
            }

            const ParticleMeshMultipoleParameters & operator = (const ParticleMeshMultipoleParameters & rhs) {
                p             = rhs.p;
                icut          = rhs.icut;
                bc            = rhs.bc;
                pos_unit_cell = rhs.pos_unit_cell; 
                n_cell        = rhs.n_cell;
                nmax          = rhs.nmax;
                mmax          = rhs.mmax;
                alpha         = rhs.alpha;
                return (*this);
            }

            bool operator == (const ParticleMeshMultipoleParameters & rhs) const {
                return ((p == rhs.p) &&
                        (icut == rhs.icut) && 
                        (bc == rhs.bc) &&
                        (pos_unit_cell.low_ == rhs.pos_unit_cell.low_) &&
                        (pos_unit_cell.high_ == rhs.pos_unit_cell.high_) &&
                        (n_cell == rhs.n_cell) &&
                        (nmax == rhs.nmax) &&
                        (mmax == rhs.mmax) &&
                        (alpha == rhs.alpha));
            }

            S32 sanity_check(const S32 check_level, std::ostream & fout=std::cout) const {
                // check_level = 0 : required by initialize() of TreeForForce.
                // check_level = 1 : required by getPMMParam() of TreeForForce.
                // check_level = 2 : required by PMM module.
                S32 ret = 0;
                if (!((0 <= check_level) && (check_level <= 2))) {
                    std::string errmsg;
                    errmsg =  "the value of argument check_level of ";
                    errmsg += "class ParticleMeshMultipoleParameters is invalid.";
                    PARTICLE_SIMULATOR_PRINT_ERROR(errmsg);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    Abort(-1);
#endif
                    std::exit(-1);
                    
                }
                if ((check_level >= 0) && (p <= 0)) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("The order of multipole expansion must be > 0");
                    fout << "p = " << p << std::endl;
                }
                if ((check_level >= 0) && (icut < 1)) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("The cell separation must be >= 1");
                    fout << "icut = " << icut << std::endl;
                }
                if ((check_level >= 1) &&
                    (bc != BOUNDARY_CONDITION_OPEN) && (bc != BOUNDARY_CONDITION_PERIODIC_XYZ)) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("This boundary condition is not supported");
                    fout << "bc = " << bc << std::endl;
                }
                if ((check_level >= 1) && (!pos_unit_cell.isValid())) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("The size of unit cell is invalid");
                    fout << "pos_unit_cell = " << pos_unit_cell << std::endl;
                }
                if ((check_level >= 0) && 
                    ((n_cell.x <= 0) || (n_cell.y <= 0) || (n_cell.z <= 0))) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("The number of the cells in each dimension must be > 0");
                    fout << "n_cell = " << n_cell << std::endl;
                }
                if ((check_level >= 2) && ((nmax <= 0) || (mmax <= 0) || (alpha <= 0.0))) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("The parameters for Ewald summation method are invalid");
                    fout << "nmax  = " << nmax << std::endl;
                    fout << "mmax  = " << mmax << std::endl;
                    fout << "alpha = " << alpha << std::endl;
                }
                return ret;
            }

            void dump(std::ostream & fout=std::cout) const {
                fout << "------------------------------------------------" << std::endl;
                fout << "Parameters of Particle Mesh Multipole method"     << std::endl;
                fout << "    p      = " << p << std::endl;
                fout << "    icut   = " << icut << std::endl;
                fout << "    bc     = " << bc << std::endl;
                fout << "    n_cell = " << n_cell << std::endl;
                fout << "    nmax   = " << nmax << std::endl;
                fout << "    mmax   = " << mmax << std::endl;
                fout << "    alpha  = " << alpha << std::endl;
                fout << "------------------------------------------------" << std::endl;
            }

        };


        //---- GetWidthOfParticleMeshCell
        inline F64vec GetWidthOfParticleMeshCell(const F64ort & pos_unit_cell,
                                                 const S32vec & n_cell) {
            F64vec width;
            width.x = pos_unit_cell.getFullLength().x / n_cell.x;
            width.y = pos_unit_cell.getFullLength().y / n_cell.y;
            width.z = pos_unit_cell.getFullLength().z / n_cell.z;
            return width; 
        }
    
        //---- GetCenterOfParticleMeshCell
        inline F64vec GetCenterOfParticleMeshCell(const F64ort & pos_unit_cell,
                                                  const F64vec & width_cell,
                                                  const S32vec & idx) {
            F64vec pos;
            pos.x = pos_unit_cell.low_.x + (idx.x + 0.5) * width_cell.x;
            pos.y = pos_unit_cell.low_.y + (idx.y + 0.5) * width_cell.y;
            pos.z = pos_unit_cell.low_.z + (idx.z + 0.5) * width_cell.z;
            return pos;
        }
    
        inline F64vec GetCenterOfParticleMeshCell(const F64ort & pos_unit_cell,
                                                  const S32vec & n_cell,
                                                  const S32vec & idx) {
            F64vec width = GetWidthOfParticleMeshCell(pos_unit_cell, n_cell);
            return GetCenterOfParticleMeshCell(pos_unit_cell, width, idx);
        }
    
        //---- GetCellIDMeasuredInUnitCell
        inline S32vec GetCellIDMeasuredInUnitCell(const F64vec & pos_unit_cell,
                                                  const F64vec & width_cell,
                                                  const F64vec & pos) {
            // This function returns a cell index of a PM cell to which
            // a particle whose position is `pos` belongs.
            // The cell index is calculated so that the index (0,0,0)
            // corresponds to the PM cell located at the lower limits
            // of the unit cell (computational box).
            F64vec idx_f;
            idx_f.x = (pos.x - pos_unit_cell.x) / width_cell.x;
            idx_f.y = (pos.y - pos_unit_cell.y) / width_cell.y;
            idx_f.z = (pos.z - pos_unit_cell.z) / width_cell.z;
            assert(fabs(idx_f.x) < (F64)std::numeric_limits<S32>::max());
            assert(fabs(idx_f.y) < (F64)std::numeric_limits<S32>::max());
            assert(fabs(idx_f.z) < (F64)std::numeric_limits<S32>::max());
            S32vec idx;
            if (idx_f.x >= 0.0) idx.x = (S32) idx_f.x;
            else idx.x = ((S32) idx_f.x) - 1;
            if (idx_f.y >= 0.0) idx.y = (S32) idx_f.y;
            else idx.y = ((S32) idx_f.y) - 1;
            if (idx_f.z >= 0.0) idx.z = (S32) idx_f.z;
            else idx.z = ((S32) idx_f.z) - 1;
            return idx;
        }
    
        inline S32vec GetCellIDMeasuredInUnitCell(const F64ort & pos_unit_cell,
                                                  const F64vec & width_cell,
                                                  const F64vec & pos) {
            return GetCellIDMeasuredInUnitCell(pos_unit_cell.low_, width_cell, pos);
        }
    
        //---- GetCellIDMeasuredInTreeRootCell
        inline S32vec GetCellIDMeasuredInTreeRootCell(const F64vec pos_unit_cell,
                                                      const F64vec width_cell,
                                                      const S32vec idx_offset,
                                                      const F64vec pos) {
            // This function returns a cell index of a PM cell to which
            // a particle whose position is `pos` belongs.
            // The cell index is calculated so that the index (0,0,0)
            // corresponds to the PM cell located at the lower limits
            // of the root cell of a tree.
            S32vec idx;
            idx = GetCellIDMeasuredInUnitCell(pos_unit_cell, width_cell, pos);
            idx.x += idx_offset.x;
            idx.y += idx_offset.y;
            idx.z += idx_offset.z; 
            return idx;
        }
        inline S32vec GetCellIDMeasuredInTreeRootCell(const F64ort pos_unit_cell,
                                                      const F64vec width_cell,
                                                      const S32vec idx_offset,
                                                      const F64vec pos) {
            return GetCellIDMeasuredInTreeRootCell(pos_unit_cell.low_,
                                                   width_cell, idx_offset, pos);
        }
    
    
        //---- GetMinBoxRoundedUpToParticleMeshCell
        inline F64ort GetMinBoxRoundedUpToParticleMeshCell(const F64vec & pos,
                                                           const F64vec & pos_unit_cell,
                                                           const F64vec & width_cell,
                                                           const S32 icut = 0) {
            S32vec idx = GetCellIDMeasuredInUnitCell(pos_unit_cell, width_cell, pos);
            F64ort box;
            box.low_.x  = pos_unit_cell.x + width_cell.x * (idx.x - icut);
            box.low_.y  = pos_unit_cell.y + width_cell.y * (idx.y - icut);
            box.low_.z  = pos_unit_cell.z + width_cell.z * (idx.z - icut);
            box.high_.x = pos_unit_cell.x + width_cell.x * (idx.x + icut + 1);
            box.high_.y = pos_unit_cell.y + width_cell.y * (idx.y + icut + 1);
            box.high_.z = pos_unit_cell.z + width_cell.z * (idx.z + icut + 1);
            return box;
        }
    
        inline F64ort GetMinBoxRoundedUpToParticleMeshCell(const F64ort & box_in,
                                                           const F64vec & pos_unit_cell,
                                                           const F64vec & width_cell,
                                                           const S32 icut = 0) {
            assert(box_in.isValid());
    
            bool is_infinity[DIMENSION_LIMIT];
            F64vec pos_tmp;
            F64ort box_out;
            S32vec idx;
    
            // Calculate box_out.low_
            for (S32 cid=0; cid < DIMENSION_LIMIT; cid++) {
                if (box_in.low_[cid] == -LARGE_FLOAT) is_infinity[cid] = true;
                else is_infinity[cid] = false;
            }
            pos_tmp = 0.0;
            for (S32 cid=0; cid < DIMENSION_LIMIT; cid++) {
                if (is_infinity[cid]==false) pos_tmp[cid] = box_in.low_[cid];
                else pos_tmp[cid] = pos_unit_cell[cid]; // dummy value to prevent assert.
            }
            idx = GetCellIDMeasuredInUnitCell(pos_unit_cell, width_cell, pos_tmp);
            for (S32 cid=0; cid < DIMENSION_LIMIT; cid++) {
                if (is_infinity[cid]==false) {
                    box_out.low_[cid] = pos_unit_cell[cid] + width_cell[cid] * (idx[cid] - icut);
                } else {
                    box_out.low_[cid] = box_in.low_[cid];
                }
            }
    
            // Calculate box_out.high_
            for (S32 cid=0; cid < DIMENSION_LIMIT; cid++) {
                if (box_in.high_[cid] == LARGE_FLOAT) is_infinity[cid] = true;
                else is_infinity[cid] = false;
            }
            pos_tmp = 0.0;
            for (S32 cid=0; cid < DIMENSION_LIMIT; cid++) {
                if (is_infinity[cid]==false) pos_tmp[cid] = box_in.high_[cid];
                else pos_tmp[cid] = pos_unit_cell[cid]; // dummy value to prevent assert.
            }
            idx = GetCellIDMeasuredInUnitCell(pos_unit_cell, width_cell, pos_tmp);
            for (S32 cid=0; cid < DIMENSION_LIMIT; cid++) {
                if (is_infinity[cid]==false) {
                    box_out.high_[cid] = pos_unit_cell[cid] + width_cell[cid] * (idx[cid] + icut + 1);
                } else {
                    box_out.high_[cid] = box_in.high_[cid];
                }
            }
    
            return box_out;
        }

    } // END of namespace of ParticleMeshMultipole
    namespace PMM = ParticleMeshMultipole;
} // END of namespace of ParticleSimulator
