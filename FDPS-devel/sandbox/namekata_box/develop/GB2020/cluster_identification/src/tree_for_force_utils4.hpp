namespace ParticleSimulator{

    template <class T>
    void MergeOrthotopeList(std::vector<T> & from_list,
                            std::vector<T> & to_list) {
    }
    
    template <>
    void MergeOrthotopeList<F64ort>(std::vector<F64ort> & from_list,
                                    std::vector<F64ort> & to_list) {
        while (!from_list.empty()) {
            // Extract vertex from `from_list`
            F64ort f = from_list.back();
            from_list.pop_back();
            // Find elements of dest which intersects with vertex above.
            std::vector<S32> delete_list;
            F64ort rslt = f;
            for (S32 i = 0; i < to_list.size(); i++) {
                if (f.overlapped(to_list[i])) {
                    delete_list.push_back(i);
                    rslt.merge(to_list[i]);
                }
            }
            // Delete
            for (S32 i = delete_list.size()-1; i >= 0; i--) {
                const S32 idx = delete_list[i];
                to_list.erase(to_list.begin() + idx);
            }
            // Add
            to_list.push_back(rslt);
        };
    
    }

    template<class T>
    void MergeOrthotopeList(std::vector<T> & from_list,
                            std::vector<std::vector<S32>> & from_lol, // lol := list of lists
                            std::vector<T> & to_list,
                            std::vector<std::vector<S32>> & to_lol) {
    }

    template <>
    void MergeOrthotopeList(std::vector<F64ort> & from_list,
                            std::vector<std::vector<S32>> & from_lol,
                            std::vector<F64ort> & to_list,
                            std::vector<std::vector<S32>> & to_lol) {
        while (!from_list.empty()) {
            // Extract vertex and rank list from `from_*list`
            F64ort f = from_list.back();
            from_list.pop_back();
            std::vector<S32> list_f = from_lol.back();
            from_lol.pop_back();
            // Find elements of dest which intersects with vertex above.
            std::vector<S32> delete_list;
            F64ort rslt = f;
            std::vector<S32> list_rslt = list_f;
            for (S32 i = 0; i < to_list.size(); i++) {
                if (f.overlapped(to_list[i])) {
                    delete_list.push_back(i);
                    rslt.merge(to_list[i]);
                    // Merge lists of ranks
                    for (S32 k = 0; k < to_lol[i].size(); k++) {
                        const S32 rnk = to_lol[i][k];
                        bool not_found {true};
                        for (S32 n = 0; n < list_rslt.size(); n++) {
                            if (list_rslt[n] == rnk) not_found = false;
                        }
                        if (not_found) list_rslt.push_back(rnk);
                    }
                }
            }
            // Delete
            for (S32 i = delete_list.size()-1; i >= 0; i--) {
                const S32 idx = delete_list[i];
                to_list.erase(to_list.begin() + idx);
                to_lol.erase(to_lol.begin() + idx);
            }
            // Add
            to_list.push_back(rslt);
            to_lol.push_back(list_rslt);
        };
    }


    template<class Ttc, class Tepj>
    std::vector<F64ort> GetClusterList(S32 adr_tc_level_partition[], 
                                       const S32 tc_size, 
                                       Ttc tc[],
                                       Tepj epj[],
                                       const S32 lev_max,
                                       const S32 n_leaf_limit,
                                       const F64 dt_crit){
        // Prepare buffers
        std::vector<F64ort> * info;
        info = new std::vector<F64ort>[tc_size];

        for(S32 i=lev_max; i>=0; --i){
            const S32 head = adr_tc_level_partition[i];
            const S32 next = adr_tc_level_partition[i+1];
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 j=head; j<next; j++){
                Ttc * tc_tmp = tc + j;
                const int n_tmp = tc_tmp->n_ptcl_;
                tc_tmp->mom_.init();
                if( n_tmp == 0) continue;
                else if(tc_tmp->isLeaf(n_leaf_limit)){
                    // Check if this leaf cell contains particles with small timesteps
                    // and calculate the size of a cuboid containing such particles if exit.
                    bool found = false;
                    F64ort rslt;
                    rslt.init();
                    const S32 adr = tc_tmp->adr_ptcl_;
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        const F64 dt = epj[k].getTimestep();
                        if (dt < dt_crit) {
                            const F64vec pos = epj[k].getPos();
                            const F64 h = epj[k].getRSearch();
                            rslt.merge(pos, h);
                            found = true;
                        }
                    }
                    if (found) info[j].push_back(rslt);
                }
                else{
                    // Merge clusters 
                    for (S32 k=0; k<N_CHILDREN; k++){
                        const S32 j_child = tc_tmp->adr_tc_ + k;
                        MergeOrthotopeList(info[j_child], info[j]);
                    }
                }
            }
        }
        // Store the result
        std::vector<F64ort> ret = info[0];

        // Free
        if (tc_size > 0) {
            for (S32 i = 0; i < tc_size; i++) 
                std::vector<F64ort>().swap(info[i]);
            delete [] info;
        }

        // Return the value
        return ret;
    }

} // END of namespace of ParticleSimulator
