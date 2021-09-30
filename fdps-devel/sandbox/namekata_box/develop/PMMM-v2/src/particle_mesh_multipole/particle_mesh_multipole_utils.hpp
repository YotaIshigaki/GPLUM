#pragma once
#include "../ps_macro_defs.h"
#include "../ps_defs.hpp"

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

        template <class Tforce, class Tepi, class Tepj, class Tspj>
        class CellInfo {
        public:
            bool epj_is_shared;

            S32 n_epi;
            S32 adr_epi_first;
            Tepi * epi_first;
            Tforce * force_first;

            std::vector<S32> n_epj;
            std::vector<Tepj * > epj_first;
            std::vector<S32> n_spj;
            std::vector<Tspj * > spj_first;

            CellInfo() {
                epj_is_shared = false;

                n_epi = 0;
                adr_epi_first = -1;
                epi_first = nullptr;
                force_first = nullptr;

                n_epj.clear();
                epj_first.clear();
                n_spj.clear();
                spj_first.clear();
            }

            friend std::ostream & operator <<(std::ostream & c, const CellInfo & info){
                c << info.n_epi << " " 
                  << info.adr_epi_first << " "
                  << info.epi_first << " "
                  << info.force_first << " "
                  << info.epj_is_shared << " "
                  << info.n_epj.size() << " "
                  << info.epj_first.size() << " "
                  << info.n_spj.size() << " "
                  << info.spj_first.size();
                return c;
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
   
        //---- GetBoxOfParticleMeshCell
        inline F64ort GetBoxOfParticleMeshCell(const S32vec & idx,
                                               const F64vec & pos_unit_cell,
                                               const F64vec & width_cell)  {
            F64ort box;
            box.low_.x  = pos_unit_cell.x + width_cell.x * idx.x;
            box.low_.y  = pos_unit_cell.y + width_cell.y * idx.y;
            box.low_.z  = pos_unit_cell.z + width_cell.z * idx.z;
            box.high_.x = pos_unit_cell.x + width_cell.x * (idx.x + 1);
            box.high_.y = pos_unit_cell.y + width_cell.y * (idx.y + 1);
            box.high_.z = pos_unit_cell.z + width_cell.z * (idx.z + 1);
            return box;
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
