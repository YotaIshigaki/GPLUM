#pragma once

#define RANK_LIST_SIZE 512

#ifdef USE_SIMPLEMAP
#include "simplemap.hpp"
#endif
class ExPair{
public:
    //PS::S64 id_in;
    //PS::S64 id_out;
    //PS::S64 id_in_rank_local;
    //PS::S64 id_out_rank_local;
    PS::S32 rank_in;
    PS::S32 id_local_in;
    PS::S32 rank_out;
    PS::S32 id_local_out;
    PS::S64 id_cluster;
    
    //PS::S32 * rank_list;
    PS::S32 rank_list[RANK_LIST_SIZE];
    
    static PS::S32 size;
    static PS::S32 rem;
    static PS::S32 n_bit;
    
    static void initialize() {
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        
        n_bit = 8 * sizeof(PS::S32);
        size = (PS::S32)std::ceil((PS::F64)n_proc/n_bit);
        rem  = n_bit*size - n_proc;
        assert( size <= RANK_LIST_SIZE );
    }
    static PS::S32 getSize() { return size+6; }
    
    ExPair(){
        //PS::S32 myrank = PS::Comm::getRank();
        //id_in = id_out = id_cluster = 0;
        //id_in_rank_local = id_out_rank_local = id_cluster = 0;
        rank_in = id_local_in = rank_out = id_local_out = -1;
        id_cluster = -1;
        
        //rank_list = new PS::S32[size];
        for ( PS::S32 i=0; i<size; i++ ) rank_list[i] = 0;
        //setFlag(myrank);
    }
    ExPair(PS::S32 myrank_in,
           PS::S32 id_loc_in,
           PS::S32 myrank_out,
           PS::S32 id_loc_out,
           PS::S64 id_cls){
        //PS::S32 myrank = PS::Comm::getRank();
        //id_in  = id_in0;
        //id_out = id_out0;
        //id_in_rank_local  = ((PS::S64)id_in_rank) * ((PS::S64)1<<32) + (PS::S64)id_in_local;
        //id_out_rank_local = ((PS::S64)id_out_rank)* ((PS::S64)1<<32) + (PS::S64)id_out_local;
        rank_in      = myrank_in;
        id_local_in  = id_loc_in;
        rank_out     = myrank_out;
        id_local_out = id_loc_out;
        id_cluster   = id_cls;
        
        //rank_list = new PS::S32[size];
        for ( PS::S32 i=0; i<size; i++ ) rank_list[i] = 0;
        //setFlag(myrank);
    }
    ExPair(const ExPair & ep){
        //id_in_rank_local  = ep.id_in_rank_local;
        //id_out_rank_local = ep.id_out_rank_local;
        rank_in      = ep.rank_in;
        id_local_in  = ep.id_local_in;
        rank_out     = ep.rank_out;
        id_local_out = ep.id_local_out;
        id_cluster   = ep.id_cluster;

        //rank_list = new PS::S32[size];
        for ( PS::S32 i=0; i<size; i++ ) rank_list[i] = ep.rank_list[i];
    }
    ExPair &operator=(const ExPair & ep){
        if ( this != &ep ){
            //id_in_rank_local  = ep.id_in_rank_local;
            //id_out_rank_local = ep.id_out_rank_local;
            rank_in      = ep.rank_in;
            id_local_in  = ep.id_local_in;
            rank_out     = ep.rank_out;
            id_local_out = ep.id_local_out;
            id_cluster   = ep.id_cluster;
        
            for ( PS::S32 i=0; i<size; i++ ) this->rank_list[i] = ep.rank_list[i];
        }
        return *this;
    }
    
    //~ExPair(){
    //    delete [] rank_list;
    //}

    PS::S32 getRank() const { return rank_in; }
    PS::S32 getIdLocal() const { return id_local_in; }
    PS::S32 getJRank() const { return rank_out; }
    PS::S32 getIdJLocal() const { return id_local_out; }
    
    std::pair<std::pair<PS::S32,PS::S32>, std::pair<PS::S32,PS::S32> > getPair() const {
        //PS::S32 id_in_local  = (PS::S32)(id_in_rank_local  % ((PS::S64)1<<32));
        //PS::S32 rank_in      = (PS::S32)(id_in_rank_local  / ((PS::S64)1<<32));
        //PS::S32 id_out_local = (PS::S32)(id_out_rank_local % ((PS::S64)1<<32));
        //PS::S32 rank_out     = (PS::S32)(id_out_rank_local / ((PS::S64)1<<32));
        return std::make_pair(std::make_pair(rank_in, id_local_in), std::make_pair(rank_out, id_local_out));
    }
    PS::S64 getIdCluster() const { return id_cluster; }
    PS::S64 setIdCluster(PS::S64 id_cluster0) { return id_cluster = id_cluster0; }

    PS::S32 input(std::vector<PS::S32> & inp,
                  const PS::S32 ofs){
        assert( ofs+size+6 <= inp.size() );
        //id_in_rank_local  = inp[1];
        //id_out_rank_local = inp[0];
        //id_cluster = inp[2];
        rank_in      = inp[ofs+2];
        id_local_in  = inp[ofs+3];
        rank_out     = inp[ofs+0];
        id_local_out = inp[ofs+1];
        id_cluster   = (PS::S64)inp[ofs+4] + (PS::S64)inp[ofs+5] * (1LL<<31);
        assert( inp[ofs+4] > -1 && inp[ofs+5] > -1);
        for ( PS::S32 i=0; i<size; i++ ) rank_list[i] = inp[ofs+i+6];
        return size+6;
    }
    PS::S32 output(std::vector<PS::S32> & outp,
                   const PS::S32 ofs){
        assert( ofs+size+6 <= outp.size() );
        //outp[0] = id_in_rank_local;
        //outp[1] = id_out_rank_local;
        //outp[2] = id_cluster;
        outp[ofs+0] = rank_in;
        outp[ofs+1] = id_local_in;
        outp[ofs+2] = rank_out;
        outp[ofs+3] = id_local_out;
        outp[ofs+4] = (PS::S32)(id_cluster % (1LL<<31));
        outp[ofs+5] = (PS::S32)(id_cluster / (1LL<<31));
        assert( id_cluster > -1 );
        for ( PS::S32 i=0; i<size; i++ ) outp[ofs+i+6] = rank_list[i];
        return size+6;
    }

    bool checkFlag(const PS::S32 i) const {
        PS::S32 n  = i / n_bit;
        PS::S32 ii = i % n_bit;
        return rank_list[n] & (1<<ii);
    }
    void setFlag(const PS::S32 i) {
        PS::S32 n  = i / n_bit;
        PS::S32 ii = i % n_bit;
        rank_list[n] |= (1<<ii);
    }
    void unsetFlag(const PS::S32 i) {
        PS::S32 n  = i / n_bit;
        PS::S32 ii = i % n_bit;
        rank_list[n] &= ~(1<<ii);
    }
    void resetFlag() {
        for ( PS::S32 i=0; i<size; i++ ) rank_list[i] = 0;
    }
    bool equalFlag(const ExPair & ep) const {
        bool check = true;
        for ( PS::S32 i=0; i<size; i++ ) check &= (rank_list[i]==ep.rank_list[i]);
        return check;
    }
    PS::S32 getMinFlag() const {
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        for (PS::S32 i=0; i<n_proc; i++) if ( checkFlag(i) ) return i;
        return n_proc;
    }

    void operator &= (const ExPair & ep) {
        for ( PS::S32 i=0; i<size; i++ ) this->rank_list[i] &= ep.rank_list[i];
	}
    void operator |= (const ExPair & ep) {
        for ( PS::S32 i=0; i<size; i++ ) this->rank_list[i] |= ep.rank_list[i];
	}

    bool exchange(const ExPair & ep) {
        bool check = (this->id_cluster != ep.id_cluster);
        this->id_cluster = std::min(this->id_cluster, ep.id_cluster);
        for ( PS::S32 i=0; i<size; i++ ) {
            check |= (this->rank_list[i] != ep.rank_list[i]);
            this->rank_list[i] |= ep.rank_list[i]; 
        }
        return check;
    }

    void dump(){
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        std::cout << "In: "    << rank_in  << "\t" << id_local_in
                  << "\tOut: " << rank_out << "\t" << id_local_out
                  << "\tCluster ID: "<< id_cluster
                  << "\tRnak List: ";
        for ( PS::S32 i=0; i<n_proc; i++ ) if ( checkFlag(i) ) std::cout << i << "\t";
        std::cout << std::endl;
    }
};

PS::S32 ExPair::size;
PS::S32 ExPair::rem;
PS::S32 ExPair::n_bit;


class NeighborId{
public:
    PS::S64 id;
    PS::S32 rank;
    PS::S32 id_local;

    NeighborId() {
        id       = -1;
        id_local = -1;
        rank     = -1;
    }
    NeighborId(PS::S32 myrank,
               PS::S32 id_loc) {
        id       = -1;
        rank     = myrank;
        id_local = id_loc;
    }
    NeighborId(PS::S64 id_0,
               PS::S32 myrank,
               PS::S32 id_loc) {
        id       = id_0;
        rank     = myrank;
        id_local = id_loc;
    }
    NeighborId(const NeighborId & ni) {
        id       = ni.id;
        rank     = ni.rank;
        id_local = ni.id_local;
    }
    NeighborId &operator=(const NeighborId & ni){
        if ( this != &ni ){
            id       = ni.id;
            rank     = ni.rank;
            id_local = ni.id_local;
        }
        return *this;
    }
};

class NeighborList{
public:
    //std::vector<std::vector<PS::S64> > n_list
    std::vector<std::vector<NeighborId> > n_list;
    std::vector<std::vector<NeighborId> > n_list_tmp;
    //#ifndef USE_SIMPLEMAP    
    //std::map<PS::S64, PS::S32> id_map;
    //#else
    //SimpleMapLib::Map<PS::S64, PS::S32> id_map;
    //#endif

    std::vector<std::vector<EPNgb> > ngb_send;
    std::vector<std::vector<EPNgb> > ngb_recv;

    std::vector<PS::S32> with_neighbor_list;
    std::vector<std::pair<PS::S32, PS::S32> > pair_list;                    // ptcl_id_loc, ngb_id_loc
    
    //std::vector<std::pair<PS::S64,PS::S64> > ex_list;
    std::vector<std::pair<PS::S32, std::pair<PS::S32, PS::S32> > > ex_list; // ptcl_id_loc, (ngb_rank, ngb_id_loc)
    std::vector<std::pair<PS::S32,PS::S32> > ex_adr_list;                   // ngb_rank, i_in_ex_data[ngb_rank]
    std::vector<PS::S32> connected_list;
    std::vector<std::vector<ExPair> > ex_data;
    //std::map<std::pair<PS::S32,PS::S32>, std::pair<PS::S32, PS::S32> > ex_data_map;

    std::vector<std::vector<PS::S32> > recv_list;
    std::vector<std::vector<PS::S32> > send_list;
    std::vector<PS::S32> recv_rank_list;
    std::vector<PS::S32> send_rank_list;

    std::vector<std::vector<PS::S32> > ex_data_send;
    std::vector<std::vector<PS::S32> > ex_data_recv;

    std::vector<NeighborId> & operator[](PS::S32 i){ return n_list[i]; }
    
    NeighborList() {
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        
        n_list.clear();
        n_list_tmp.clear();
        //id_map.clear();
        with_neighbor_list.clear();
        pair_list.clear();
        ex_list.clear();
        ex_adr_list.clear();
        connected_list.clear();
        //ex_data_map.clear();
        recv_rank_list.clear();
        send_rank_list.clear();
        ex_data_send.clear();
        ex_data_recv.clear();
        
        ex_data.resize(n_proc);
        recv_list.resize(n_proc);
        send_list.resize(n_proc);
        ngb_send.resize(n_proc);
        ngb_recv.resize(n_proc);
        ex_data_send.resize(n_proc);
        ex_data_recv.resize(n_proc);
        
#pragma omp parallel for
        for (PS::S32 i=0; i<n_proc; i++){
            ex_data[i].clear();
            recv_list[i].clear();
            send_list[i].clear();
            ngb_send[i].clear();
            ngb_recv[i].clear();
        }
      
        ExPair::initialize();
    }
    
    template <class Tpsys>
    void initializeList(Tpsys & pp) {
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        const PS::S32 n_loc  = pp.getNumberOfParticleLocal();

        n_list.clear();
        n_list_tmp.clear();
        //id_map.clear();
        with_neighbor_list.clear();
        pair_list.clear();
        ex_list.clear();
        ex_adr_list.clear();
        connected_list.clear();
        //ex_data_map.clear();
        recv_rank_list.clear();
        send_rank_list.clear();
            
#pragma omp parallel for
        for ( PS::S32 i=0; i<n_proc; i++ ){
            ex_data[i].clear();
            recv_list[i].clear();
            send_list[i].clear();
            ngb_send[i].clear();
            ngb_recv[i].clear();
        }
        
        n_list.resize(n_loc);
        n_list_tmp.resize(n_loc);
#pragma omp parallel for
        for(PS::S32 i=0; i<n_loc; i++) {
            n_list[i].clear();
            n_list_tmp[i].clear();
            n_list_tmp[i].resize(8);
        }
    }

    ExPair & getExData(std::pair<PS::S32, PS::S32> adr) {
        return ex_data[adr.first][adr.second];
    }
    ExPair & getExData(const PS::S32 i_id_loc,
                       const PS::S32 j_rank,
                       const PS::S32 j_id_loc) {
        assert( j_rank > -1 && j_rank < PS::Comm::getNumberOfProc() );
        assert( j_id_loc > -1 && i_id_loc > -1);
        const PS::U32 n_size = ex_data[j_rank].size();
        assert( n_size > 0 );

        for (PS::S32 i=0; i<n_size; i++)
            if ( ex_data[j_rank][i].getIdLocal() == i_id_loc
                 && ex_data[j_rank][i].getIdJLocal() == j_id_loc ) return ex_data[j_rank][i];
        assert( n_size == 0 );
        return ex_data[j_rank][n_size-1];
        
        //PS::S32 i = 0;
        //while ( ex_data[j_rank][i].getIdLocal() != i_id_loc
        //        || ex_data[j_rank][i].getIdJLocal() != j_id_loc ) i++;
        //assert ( i < n_size );
        //return ex_data[j_rank][i];
    }
    
    PS::S32 getNumberOfParticlesWithNeighbor() const { return with_neighbor_list.size(); }
    PS::S32 getNumberOfNeighborPairsLocal() const { return pair_list.size(); }
    
    PS::S32 getNumberOfRankSend() const { return send_rank_list.size(); }
    PS::S32 getNumberOfRankRecv() const { return recv_rank_list.size(); }
    PS::S32 getNumberOfRankConnected() const { return connected_list.size(); }
    PS::S32 getNumberOfPairConnected(const PS::S32 ii) const { return ex_data[connected_list.at(ii)].size(); }

    PS::U32 getNumberOfNeighbor(const PS::S32 i) const {return n_list[i].size(); }
    //PS::U32 getNumberOfTemporaryNeighbor(const PS::S32 i) const {return n_list_tmp[i].size(); }

    EPNgb & getNeighborInfo(PS::S32 rank,
                            PS::S32 id_loc) {
        assert( rank > -1 && rank < PS::Comm::getNumberOfProc() );
        assert( id_loc > -1 );
        PS::S32 n = ngb_send[rank].size();
        assert( n > 0 );
        
        for (PS::S32 i=0; i<n; i++) if ( ngb_send[rank][i].id_local == id_loc ) return ngb_send[rank][i];
        assert( n == 0 );
        return ngb_send[rank][0];
        //PS::S32 i = 0;
        //while ( ngb_recv[rank][i].id_local != id_loc ) i++;
        //assert( i < n );
    }

    template <class Tpsys, class Tptree>
    void makeTemporaryNeighborList(Tpsys & pp,
                                   Tptree & tree_grav) {
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();
        
#pragma omp parallel for
        for (PS::S32 i=0; i<n_loc; i++){
            PS::S32  neighbor = pp[i].neighbor.number;
            //n_list_tmp[i].resize(neighbor);
            assert( neighbor >= 0 );

            if ( neighbor > 0 && neighbor < 8 ) {
                for (PS::S32 j=0; j<neighbor; j++) {
                    PS::S32 rank   = pp[i].neighbor.getRank(j);
                    PS::S32 id_loc = pp[i].neighbor.getId(j);
                    assert( rank > -1 );
                    assert( id_loc > -1 );

                    //#pragma omp critical
                    //{
                        //n_list_tmp[i].push_back(NeighborId(rank, id_loc));
                    n_list_tmp[i][j] = NeighborId(rank, id_loc);
                    //}
                } 
                
            } else if ( neighbor >= 8 ) {
                PS::S32 n_ngb = 0;
                EPJ_t* next = NULL;
                n_ngb = tree_grav.getNeighborListOneParticle(pp[i], next);
                n_list_tmp[i].resize(n_ngb-1);

                PS::S32 jj = 0;
                for ( PS::S32 j=0; j<n_ngb; j++ ){
                    PS::S32 rank   = (next+j)->myrank;
                    PS::S32 id_loc = (next+j)->id_local;
                    if ( pp[i].myrank == rank && pp[i].id_local == id_loc ) continue;

#ifdef USE_INDIVIDUAL_CUTOFF
                    PS::F32    rsearchi = pp[i].r_search;
                    PS::F32    rsearchj = (next+j)->r_search;
                    PS::F32    rsearch  = std::max(rsearchi, rsearchj);
#else
                    PS::F32    rsearch  = FP_t::r_search;
#endif
                    PS::F32    rsearch2 = rsearch * rsearch * SAFTY_FACTOR2;
                    
                    PS::F64vec rij_64 = (next+j)->pos - pp[i].pos;
                    PS::F32vec rij    = (PS::F32vec)rij_64;
                    PS::F32    dr2 = rij * rij + (PS::F32)FP_t::eps2;

                    if ( dr2 < rsearch2 ) {
                        //#pragma omp critical
                        //{
                        //n_list_tmp[i].push_back(NeighborId(rank, id_loc));
                        n_list_tmp[i][jj] = NeighborId(rank, id_loc);
                        //}
                        jj ++ ;
                    }
                }
                assert( neighbor == jj );
            }
        }
    }

    template <class Tpsys>
    void exchangeNeighborInfo(Tpsys & pp) {
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();
        PS::S32 n_send[n_proc];
        PS::S32 n_recv[n_proc];     
        for (PS::S32 i=0; i<n_proc; i++) n_send[i] = 0;
        
#pragma omp parallel for reduction (+:n_send[:n_proc])
        for(PS::S32 i=0; i<n_loc; i++){
            PS::U32 n_size = pp[i].neighbor.number;
            for ( PS::S32 j=0; j<n_size; j++ ){
                PS::S32 rank   = n_list_tmp[i][j].rank;
                PS::S32 id_loc = n_list_tmp[i][j].id_local;
                assert( rank > -1 && rank < PS::Comm::getNumberOfProc() );
                assert( id_loc > -1 );
                if ( pp[i].myrank != rank ){
                    n_send[rank] ++;
#pragma omp critical
                    {
                        ngb_send[rank].push_back(EPNgb(id_loc));
                        //ngb_send[rank].push_back(EPNgb(pp[i]));
                        //EPNgb epn;
                        //epn.copyFromFP(pp[i]);
                        //ngb_send[rank].push_back(epn);
                    }
                } 
            }
        }
        
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Request req0[n_proc],  req1[n_proc];
        MPI_Status  stat0[n_proc], stat1[n_proc];
        for ( PS::S32 i=0; i<n_proc; i++ ) {
            PS::S32 rank = i;
            PS::S32 n_size = 1;
            PS::S32 TAG = 1024;
                    
            MPI_Isend(&n_send[i], n_size, PS::GetDataType(n_send[i]), rank, TAG, MPI_COMM_WORLD, &req0[i]);
            MPI_Irecv(&n_recv[i], n_size, PS::GetDataType(n_recv[i]), rank, TAG, MPI_COMM_WORLD, &req1[i]);
        }
        //MPI_Waitall(n_proc, req0, stat0);
        MPI_Waitall(n_proc, req1, stat1);
#else
        for ( PS::S32 i=0; i<n_proc; i++ ) assert ( n_send[i] == 0 );
#endif

        std::vector<PS::S32> send_rank_list;
        std::vector<PS::S32> recv_rank_list;
        send_rank_list.clear();
        recv_rank_list.clear();
        for(PS::S32 i=0; i<n_proc; i++){
            if ( n_send[i] ) send_rank_list.push_back(i);
            if ( n_recv[i] ) recv_rank_list.push_back(i);

            assert( ngb_send[i].size() == n_send[i] );
            ngb_recv[i].resize(n_recv[i]);
        }
            
        
        /*//#pragma omp parallel for
        for(PS::S32 i=0; i<n_proc; i++){
            if ( rank_c[i] > -1 ) {
                rank_c[i] = n_send;
                n_send ++;
                assert( i != PS::Comm::getRank() );
                
                ngb_recv[i].resize(ngb_send[i].size());
            }
            }*/

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_send_rank = send_rank_list.size();
        PS::S32 n_recv_rank = recv_rank_list.size();
        MPI_Request req2[n_send_rank],  req3[n_recv_rank];
        MPI_Status  stat2[n_send_rank], stat3[n_recv_rank];
        if ( n_send_rank ) {
            for ( PS::S32 i=0; i<n_send_rank; i++ ) {
                PS::S32 rank = send_rank_list.at(i);
                PS::S32 n_size = n_send[rank];
                PS::S32 TAG = 1025;
                    
                MPI_Isend(&ngb_send[rank][0], n_size, PS::GetDataType(ngb_send[rank][0]), rank, TAG, MPI_COMM_WORLD, &req2[i]);
            }
        }
        if ( n_recv_rank ) {
            for ( PS::S32 i=0; i<n_recv_rank; i++ ) {
                PS::S32 rank = recv_rank_list.at(i);
                PS::S32 n_size = n_recv[rank];
                PS::S32 TAG = 1025;
                    
                MPI_Irecv(&ngb_recv[rank][0], n_size, PS::GetDataType(ngb_recv[rank][0]), rank, TAG, MPI_COMM_WORLD, &req3[i]);
            }
        }
        //MPI_Waitall(n_send_rank, req2, stat2);
        MPI_Waitall(n_recv_rank, req3, stat3);
#endif

        
#pragma omp parallel for
        for(PS::S32 i=0; i<n_proc; i++){
            for(PS::S32 j=0; j<n_recv[i]; j++){
                ngb_recv[i][j].copyNeighbor(pp);
            }
        }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Request req4[n_send_rank],  req5[n_recv_rank];
        MPI_Status  stat4[n_send_rank], stat5[n_recv_rank];
        if ( n_send_rank ) {
            for ( PS::S32 i=0; i<n_send_rank; i++ ) {
                PS::S32 rank = send_rank_list.at(i);
                PS::S32 n_size = n_send[rank];
                PS::S32 TAG = 1026;
                
                MPI_Irecv(&ngb_send[rank][0], n_size, PS::GetDataType(ngb_send[rank][0]), rank, TAG, MPI_COMM_WORLD, &req4[i]);
            }
        }
        if ( n_recv_rank ) {
            for ( PS::S32 i=0; i<n_recv_rank; i++ ) {
                PS::S32 rank = recv_rank_list.at(i);
                PS::S32 n_size = n_recv[rank];
                PS::S32 TAG = 1026;
                    
                MPI_Isend(&ngb_recv[rank][0], n_size, PS::GetDataType(ngb_recv[rank][0]), rank, TAG, MPI_COMM_WORLD, &req5[i]);
            }
        }
        MPI_Waitall(n_send_rank, req4, stat4);
        //MPI_Waitall(n_recv_rank, req5, stat5);
#endif

        /*
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        if ( n_send ) {
            MPI_Request req0[n_send],  req1[n_send];
            MPI_Status  stat0[n_send], stat1[n_send];
            for ( PS::S32 i=0; i<n_proc; i++ ) if ( rank_c[i] > -1 ){
                    PS::S32 ii = rank_c[i];
                    PS::S32 rank = i;
                    PS::S32 n_size = ngb_send[i].size();
                    PS::S32 TAG = 1025;
                    assert( n_size > 0 );
                    
                    MPI_Isend(&ngb_recv[i][0], n_size, PS::GetDataType(ngb_recv[i][0]), rank, TAG, MPI_COMM_WORLD, &req0[ii]);
                    MPI_Irecv(&ngb_send[i][0], n_size, PS::GetDataType(ngb_send[i][0]), rank, TAG, MPI_COMM_WORLD, &req1[ii]);
                }
            //MPI_Waitall(n_send, req0, stat0);
            MPI_Waitall(n_send, req1, stat1);
        }
#else
        for ( PS::S32 i=0; i<n_proc; i++ ) assert ( ngb_send[i].size() == 0 );
#endif 
        */ 
    }

    template <class Tpsys>
    void addNeighbor(Tpsys & pp,
                     PS::S32 i,
                     PS::S64 j_id,
                     PS::S32 j_rank,
                     PS::S32 j_id_local) {
        assert( j_id > -1 );
        assert( j_rank > -1 );
        assert( j_id_local > -1 );
        n_list[i].push_back(NeighborId(j_id, j_rank, j_id_local));
        pp[i].neighbor.number ++;

        pp[i].id_cluster = std::min(pp[i].id_cluster, j_id);

        if ( j_rank != pp[i].myrank ) {
#pragma omp critical
            {
                ex_list.push_back(std::make_pair(i, std::make_pair(j_rank, j_id_local)));
                ex_adr_list.push_back(std::make_pair(j_rank, ex_data.at(j_rank).size()));
                //ex_data_map[std::make_pair(pp[i].id, j_id)] = std::make_pair(j_rank, ex_data.at(j_rank).size());
                ExPair ex_pair(pp[i].myrank, i, j_rank, j_id_local, pp[i].id_cluster);
                ex_pair.setFlag(pp[i].myrank);
                ex_pair.setFlag(j_rank);
                ex_data.at(j_rank).push_back(ex_pair);
            }
            pp[i].inDomain = false;
        } else {
            //if ( j_id_local < 0 ) j_id_local = id_map.at(j_id);
            if ( i<j_id_local ) {
#pragma omp critical
                {
                    pair_list.push_back(std::make_pair(i, j_id_local));
                }
            }
        }
    }
    
    template <class Tpsys>
    void checkNeighbor(Tpsys & pp) {
        const PS::S32 n_loc = n_list.size();
        bool check = true;
        PS::S32 nei_tot = 0;
        
        //for ( PS::S32 i=0; i<n_loc; i++ ) {
        //    if ( !pp[i].isDead )
                //assert ( id_map.at(pp[i].id) == i );
        //}

        for ( PS::S32 i=0; i<n_loc; i++ ) {
            PS::S32 n_ngb = n_list.at(i).size();
            //if ( pp[i].neighbor.number )
            //   std::cout << pp[i].id << "\t";
            nei_tot += n_ngb;
            
            for ( PS::S32 jj=0; jj<n_ngb; jj++ ) {
                PS::S64 j_id = n_list.at(i).at(jj).id;
                //if ( pp[i].neighbor.number )
                //    std::cout << j_id << " ";
                //auto itr = id_map.find(j_id);
                //if ( itr == id_map.end() ) continue;
                //#ifndef USE_SIMPLEMAP		
                //                PS::S32 j    = itr->second;
                //#else		
                //                PS::S32 j    = id_map.second(itr);
                //#endif
                PS::S32 j_rank = n_list.at(i).at(jj).rank;
                PS::S32 j      = n_list.at(i).at(jj).id_local;
                PS::S32 n_ngb_j = n_list.at(j).size();
                if ( j_rank != pp[i].myrank ) continue;

                PS::S32 n_p = 0;
                for ( PS::S32 kk=0; kk<n_ngb_j; kk++ ) {
                    PS::S64 k_id = n_list.at(j).at(kk).id;
                    /*                    auto itr1 = id_map.find(k_id);
                    if ( itr1 == id_map.end() ) continue;
#ifndef USE_SIMPLEMAP		    
		    auto ss = itr1->second;
#else		    
		    auto ss = id_map.second(itr1);
#endif		    
                    */
                    PS::S32 k_rank = n_list.at(j).at(kk).rank;
                    PS::S32 k = n_list.at(j).at(kk).id_local;
                    if ( k_rank != pp[j].myrank ) continue;

                    if ( k == i ) n_p ++ ;
                }
                if ( n_p != 1 && n_ngb > 0 ) {
                    std::cout << i << "\t" << pp[i].id << "\t" << pp[i].myrank << "\t" << j << "\t" << j_id << "\t" << j_rank << std::endl;
                    std::cout << "Neighbor of " << pp[i].id << ": ";
                    for (PS::S32 k=0; k<n_list.at(i).size(); k++) std::cout << n_list.at(i).at(k).id << "\t";
                    std::cout << std::endl;
                    std::cout << "Neighbor of " << j_id << ": ";
                    for (PS::S32 k=0; k<n_list.at(j).size(); k++) std::cout << n_list.at(j).at(k).id << "\t";
                    std::cout << std::endl;
                    check = check && false;
                }
            }
            //if ( pp[i].neighbor.number )
            //    std::cout << std::endl;
        }

        PS::S32 nei_tot_glb = PS::Comm::getSum(nei_tot);
        assert ( nei_tot_glb%2 == 0 );

        if ( !check ) {
            PS::Abort();
        }
    }

    void createConnectedRankList(){
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        connected_list.clear();

        for ( PS::S32 i=0; i<n_proc; i++ ) {
            if ( ex_data[i].size() ) {
                connected_list.push_back(i);
                assert( i != PS::Comm::getRank() );
            }
        }
    }

    void resizeExDataBuffer() {
        PS::S32 n_send = connected_list.size();
        
        //ex_data_send.resize(n_send);
        //ex_data_recv.resize(n_send);
        for ( PS::S32 i=0; i<n_send; i++ ) {
            PS::S32 n_size = ex_data[connected_list.at(i)].size() * ExPair::getSize();
            ex_data_send.at(i).resize(n_size);
            ex_data_recv.at(i).resize(n_size);
        }
    }

    /*    template <class Tpsys>
    void makeIdMap(Tpsys & pp){
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();
        //id_map.clear();
        //assert( (PS::S32)(n_list.size()) == n_loc );
#ifdef USE_SIMPLEMAP
	id_map.resize(n_loc);
#pragma omp parallel for schedule(static)
#endif	
        for(PS::S32 i=0; i<n_loc; i++){
            //assert( pp[i].neighbor.number == (PS::S32)(n_list[i].size()) );
            if ( !pp[i].isDead ) {
#ifndef USE_SIMPLEMAP
                id_map[pp[i].id] = i;
#else
		id_map.set(pp[i].id, i);
#endif		
            }else{
#ifdef USE_SIMPLEMAP
		id_map.set(-1, i);
#endif		
	    }
		
        }
#ifdef USE_SIMPLEMAP
	id_map.makemap();
#endif	
    }
    */

#if 1
    template <class Tpsys>
    void createNeighborCluster(Tpsys & pp){
        //const PS::S32 n_loc  = pp.getNumberOfParticleLocal();
        const PS::S32 n_wngb = with_neighbor_list.size();
        const PS::S32 n_pair = pair_list.size();
        
        bool check = true;
        while( check ){
            check = false;
            
#pragma omp parallel for reduction (||:check)
            for(PS::S32 ii=0; ii<n_pair; ii++){
                PS::S32 i  = pair_list.at(ii).first;
                PS::S32 j  = pair_list.at(ii).second;

                if ( pp[i].id_cluster != pp[j].id_cluster ) {
#pragma omp critical
                    {
                        pp[i].id_cluster = pp[j].id_cluster = std::min(pp[i].id_cluster, pp[j].id_cluster);
                    }
                    check = check || true;
                }
            }
        }

        if( ex_list.size() != 0 ){
            PS::S32 n_out = ex_list.size();
#pragma omp parallel for
            for(PS::S32 ii=0; ii<n_wngb; ii++){
                PS::S32 i = with_neighbor_list.at(ii);
                for(PS::S32 j=0; j<n_out; j++){
                    //PS::S32 i_out = id_map.at(ex_list.at(j).first);
                    PS::S32 i_out = ex_list[j].first;
                    PS::S32 id_cluster_out = pp[i_out].id_cluster;
                    if( pp[i].id_cluster == id_cluster_out ) pp[i].inDomain = false;
                }
            }
        }
        
    }
#else
    template <class Tpsys>
    void createNeighborCluster(Tpsys & pp){
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();

        PS::S64 j_id_cluster = 0;
        PS::S64 id_cluster[n_loc];
        bool check = true;
        while( check ){
            check = false;
#pragma omp parallel for
            for(PS::S32 i=0; i<n_loc; i++){
                PS::S64 j_id = 0;
                PS::S32 nei = 0;
                nei = pp[i].neighbor.number;
                id_cluster[i] = pp[i].id_cluster;
                
                if(nei == 0) continue;
                for(PS::S32 j=0; j<nei; j++){
                    /*auto itr = id_map.find(n_list[i].at(j));
                    if ( itr == id_map.end() ) continue;
#ifndef USE_SIMPLEMAP		
                    j_id = itr->second;
#else		    
                    j_id = id_map.second(itr);
                    #endif		    */
                    PS::S32 j_rank = n_list[i][j].rank;
                    PS::S32 j_id   = n_list[i][j].id_local;
                    if ( j_rank != pp[i].myrank ) continue;
                    
                    j_id_cluster = pp[j_id].id_cluster;
                    if( id_cluster[i] > j_id_cluster ) id_cluster[i] = j_id_cluster;
                }
            }
#pragma omp parallel for reduction (||:check)
            for(PS::S32 i=0; i<n_loc; i++){
                if ( pp[i].id_cluster != id_cluster[i] ) {
                    check = check || true;
                    pp[i].id_cluster = id_cluster[i];
                }
                assert( pp[i].id >= id_cluster[i] );
            }
        }

        if( ex_list.size() != 0 ){
            PS::S32 n_out = ex_list.size();
#pragma omp parallel for
            for(PS::S32 i=0; i<n_loc; i++){
                for(PS::S32 j=0; j<n_out; j++){
                    PS::S32 i_out = id_map.at(ex_list.at(j).first);
                    PS::S32 id_cluster_out = pp[i_out].id_cluster;
                    if( pp[i].id_cluster == id_cluster_out ) pp[i].inDomain = false;
                }
            }
        }
        
    }
#endif

    template <class Tpsys>
    void inputExData(Tpsys & pp){
        const PS::S32 n_out = ex_list.size();

#pragma omp parallel for
        for ( PS::S32 j=0; j<n_out; j++ ){
            //std::pair<PS::S64,PS::S64> pair = ex_list.at(j);
            PS::S32 ii = ex_list.at(j).first;
            std::pair<PS::S32,PS::S32> ex_adr = ex_adr_list.at(j);
            //assert( getExData(ex_adr).getId() == pair.first );
            getExData(ex_adr).setIdCluster(pp[ii].id_cluster);
        }
        
        for ( PS::S32 j=0; j<n_out; j++ ){
            //std::pair<PS::S32,PS::S32> pair = ex_list.at(j);
            std::pair<PS::S32,PS::S32> ex_adr = ex_adr_list.at(j);

            //assert( getExData(ex_adr).getId() == pair.first );
            //getExData(ex_adr).setIdCluster(pp[id_map.at(pair.first)].id_cluster);

            for ( PS::S32 k=0; k<n_out; k++ ){
                if ( k == j ) continue;
                //std::pair<PS::S32,PS::S32> pair2 = ex_list.at(k);
                std::pair<PS::S32,PS::S32> ex_adr2 = ex_adr_list.at(k);
                if ( getExData(ex_adr2).getIdCluster() == getExData(ex_adr).getIdCluster() ) {
                    getExData(ex_adr).exchange(getExData(ex_adr2));
                }
            }
        }
    }

    template <class Tpsys>
    bool exchangeExData(Tpsys & pp,
                        PS::S32 TAG){
        //PS::S32** & ex_data_send,
        //PS::S32** & ex_data_recv){
        //const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        const PS::S32 n_send = connected_list.size();
        //ex_data_send.resize(n_send);
        //ex_data_recv.resize(n_send);
        
        /*for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = connected_list.at(ii);
            PS::S32 n_size = ex_data[i].size() * ExPair::getSize();
            
            ex_data_send[ii].resize(n_size);
            ex_data_recv[ii].resize(n_size);
            }*/
#pragma omp parallel for
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = connected_list.at(ii);
            PS::S32 n_data = ex_data[i].size();
            
            PS::S32 jj = 0;
            for ( PS::S32 j=0; j<n_data; j++ ) {
                jj += ex_data[i][j].output(ex_data_send[ii], jj);
            }
        }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Request req0[n_send],  req1[n_send];
        MPI_Status  stat0[n_send], stat1[n_send];
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = connected_list.at(ii);
            PS::S32 n_size = ex_data[i].size() * ExPair::getSize();
            
            MPI_Isend(&ex_data_send[ii][0], n_size, PS::GetDataType(ex_data_send[ii][0]), i, TAG, MPI_COMM_WORLD, &req0[ii]);
            MPI_Irecv(&ex_data_recv[ii][0], n_size, PS::GetDataType(ex_data_recv[ii][0]), i, TAG, MPI_COMM_WORLD, &req1[ii]);
        }
        //MPI_Waitall(n_send, req0, stat0);
        MPI_Waitall(n_send, req1, stat1);
#else
        assert ( n_send == 0 );
#endif
        

        bool check = false;
        
#pragma omp parallel for reduction (||:check)
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = connected_list.at(ii);
            PS::S32 n_data = ex_data[i].size();

            PS::S32 jj = 0;
            for ( PS::S32 j=0; j<n_data; j++ ) {
                ExPair recv_pair;
                jj += recv_pair.input(ex_data_recv[ii], jj);

                PS::S32 i_id_loc = recv_pair.getIdLocal();
                PS::S32 j_rank   = recv_pair.getJRank();
                PS::S32 j_id_loc = recv_pair.getIdJLocal();
                ExPair & send_pair = getExData(i_id_loc, j_rank, j_id_loc);
                assert( send_pair.getRank() == recv_pair.getRank() );
                assert( send_pair.getIdLocal() == recv_pair.getIdLocal() );
                bool check_1 = send_pair.exchange(recv_pair);

                //std::pair<PS::S32,PS::S32> adr = ex_data_map.at(recv_pair.getPair());
                //assert ( adr.first == i );
                //assert ( recv_pair.getPair() == getExData(adr).getPair() );
                //bool check_1 = getExData(adr).exchange(recv_pair);
                check = check || check_1;
                //getExData(adr).show();
#pragma omp critical
                {
                    //PS::S32 i_loc = id_map.at(getExData(adr).getId());
                    PS::S32 i_loc = send_pair.getIdLocal();
                    pp[i_loc].id_cluster = std::min(pp[i_loc].id_cluster, send_pair.getIdCluster());
                }
            }
            
            //delete [] ex_data_send[ii];
            //delete [] ex_data_recv[ii];
        }
        //delete [] ex_data_send;
        //delete [] ex_data_recv;

        //PS::Comm::barrier();
        //bool check_glb = PS::Comm::synchronizeConditionalBranchOR(check);
        
        return check;
    }

    template <class Tpsys>
    void selectSendRecvParticle(Tpsys & pp){
        const PS::S32 myrank = PS::Comm::getRank();
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        const PS::S32 n_ptcl = ex_list.size();
        std::vector<PS::S64> ex_cluster;
        std::vector<std::pair<PS::S32,PS::S32> > ex_cluster_adr;
        ex_cluster.clear();
        ex_cluster_adr.clear();
        
        for ( PS::S32 ii=0; ii<n_ptcl; ii++ ) {
            //std::pair<PS::S32,PS::S32> pair = ex_list.at(ii);
            std::pair<PS::S32,PS::S32> adr =  ex_adr_list.at(ii);
            PS::S64 id_cluster = getExData(adr).id_cluster;

            PS::S32 n_l = ex_cluster.size();
            std::pair<PS::S32,PS::S32> adr2 = std::make_pair(-1,-1);
            for (PS::S32 j=0; j<n_l; j++){
                if ( id_cluster == ex_cluster.at(j) ){
                    adr2 = ex_cluster_adr.at(j); 
                    assert( getExData(adr).equalFlag(getExData(adr2)) );
                }
            }
            
            if ( adr2 == std::make_pair(-1,-1) ){
                ex_cluster.push_back(id_cluster);
                ex_cluster_adr.push_back(adr);
                
                PS::S32 min_rank = getExData(adr).getMinFlag();              
                if ( min_rank == myrank ) {
                    for ( PS::S32 j=0; j<n_proc; j++ ) {
                        if ( getExData(adr).checkFlag(j) ) {
                            if ( j == myrank ) continue;
                            recv_list[j].push_back(id_cluster);
                            assert ( j > myrank );
                        }
                    }
                } else {
                    assert ( min_rank < myrank );
                    send_list[min_rank].push_back(id_cluster);
                }
            }
        }

        for ( PS::S32 i=0; i<n_proc; i++ ) {
            if ( recv_list[i].size() ) recv_rank_list.push_back(i);
            if ( send_list[i].size() ) send_rank_list.push_back(i);
        }
    }

private:
    void operator =(const NeighborList& NL){}
    NeighborList(const NeighborList& NL) {}
};



template <class Tp>
class ExParticleSystem {
 public :
    PS::S32 n_send;
    PS::S32 n_recv;
    PS::S32 n_ex_ptcl_send_tot;
    //PS::S32 n_ex_nei_send_tot;
    PS::S32 n_ex_ptcl_recv_tot;
    //PS::S32 n_ex_nei_recv_tot;

    std::vector<Tp>         ex_ptcl_send;
    //std::vector<NeighborId> ex_nei_send;
    std::vector<Tp>         ex_ptcl_recv;
    //std::vector<NeighborId> ex_nei_recv;

    std::vector<std::vector<PS::S32> > ex_ptcl_send_list;
    //std::vector<NeighborId*> n_list;

    std::vector<PS::S32> n_ex_ptcl_send;
    //std::vector<PS::S32> n_ex_nei_send;
    std::vector<PS::S32> n_ex_ptcl_recv;
    //std::vector<PS::S32> n_ex_nei_recv;

    std::vector<PS::S32> adr_ex_ptcl_send;
    //std::vector<PS::S32> adr_ex_nei_send;
    std::vector<PS::S32> adr_ex_ptcl_recv;
    //std::vector<PS::S32> adr_ex_nei_recv;

    
    Tp & operator[](PS::S32 i){ return ex_ptcl_recv[i]; }
    PS::S32 getNumberOfParticleLocal() const { return n_ex_ptcl_recv_tot; }

    void initialize() {
        n_send = n_recv = 0;
        n_ex_ptcl_send_tot = n_ex_ptcl_recv_tot = 0;
        //n_ex_nei_send_tot  = n_ex_nei_recv_tot  = 0;
        
        ex_ptcl_send.clear();
        //ex_nei_send.clear();
        ex_ptcl_recv.clear();
        //ex_nei_recv.clear();
        
        ex_ptcl_send_list.clear();
    
        n_ex_ptcl_send.clear();
        //n_ex_nei_send.clear();
        n_ex_ptcl_recv.clear();
        //n_ex_nei_recv.clear();

        adr_ex_ptcl_send.clear();
        //adr_ex_nei_send.clear();
        adr_ex_ptcl_recv.clear();
        //adr_ex_nei_recv.clear();
    }
    

    void resize(PS::S32 n_send0,
                PS::S32 n_recv0){
        n_send = n_send0;
        n_ex_ptcl_send.resize(n_send);
        //n_ex_nei_send.resize(n_send);
        adr_ex_ptcl_send.resize(n_send);
        //adr_ex_nei_send.resize(n_send);

        ex_ptcl_send_list.resize(n_send);
#pragma omp parallel for
        for ( PS::S32 i=0; i<n_send; i++ ) ex_ptcl_send_list[i].clear();
        
        n_recv = n_recv0;
        n_ex_ptcl_recv.resize(n_recv);
        //n_ex_nei_recv.resize(n_recv);
        adr_ex_ptcl_recv.resize(n_recv);
        //adr_ex_nei_recv.resize(n_recv);
    }

    PS::S32 getNumberOfParticleSend() const { return n_ex_ptcl_send_tot; }
    PS::S32 getNumberOfParticleRecv() const { return n_ex_ptcl_recv_tot; }
    //PS::S32 getNumberOfNeighborSend() const { return n_ex_nei_send_tot; }
    //PS::S32 getNumberOfNeighborRecv() const { return n_ex_nei_recv_tot; }

    template <class Tpsys>
    void inputNumberOfExParticleSend(Tpsys & pp,
                                     NeighborList & NList){
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();
        
#pragma omp parallel for
        //for ( PS::S32 ii=0; ii<n_send; ii++ ) n_ex_ptcl_send[ii] = n_ex_nei_send[ii] = 0;
        for ( PS::S32 ii=0; ii<n_send; ii++ ) n_ex_ptcl_send[ii] = 0;

        if ( n_send ) {
#pragma omp parallel for
            for ( PS::S32 i=0; i<n_loc; i++) {
                if ( !pp[i].inDomain ) {
                    for ( PS::S32 jj=0; jj<n_send; jj++ ){
                        PS::S32 j = NList.send_rank_list[jj];
                        PS::S32 n_data = NList.send_list[j].size();
                        for ( PS::S32 k=0; k<n_data; k++ ) {
                            if ( NList.send_list[j][k] == pp[i].id_cluster ) {
#pragma omp critical       
                                {
                                    n_ex_ptcl_send[jj] ++;
                                    //n_ex_nei_send[jj] += pp[i].neighbor.number;
                                    //assert ( pp[i].neighbor.number == (PS::S32)(NList.n_list[i].size()) );
                                    ex_ptcl_send_list[jj].push_back(i);
                                }
                            }
                        }
                    }
                }
            }
        }
#pragma omp parallel for
        for ( PS::S32 ii=0; ii<n_send; ii++ ) assert( ex_ptcl_send_list[ii].size() );
    }

    void sendRecvNumberOfExParticle(NeighborList & NList,
                                    PS::S32 TAG = 0){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Request req0[n_send];
        MPI_Status  stat0[n_send];
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = NList.send_rank_list[ii];
            MPI_Isend(&n_ex_ptcl_send[ii], 1, PS::GetDataType(n_ex_ptcl_send[0]), i,
                      TAG,   MPI_COMM_WORLD, &req0[ii]);
            //MPI_Isend(&n_ex_nei_send[ii],  1, PS::GetDataType(n_ex_nei_send[0]),  i,
            //          TAG+1, MPI_COMM_WORLD, &req1[ii]);
        }
        MPI_Request req2[n_recv];
        MPI_Status  stat2[n_recv];
        for ( PS::S32 ii=0; ii<n_recv; ii++ ) {
            PS::S32 i = NList.recv_rank_list[ii];
            MPI_Irecv(&n_ex_ptcl_recv[ii], 1, PS::GetDataType(n_ex_ptcl_recv[0]), i,
                      TAG,   MPI_COMM_WORLD, &req2[ii]);
            //MPI_Irecv(&n_ex_nei_recv[ii],  1, PS::GetDataType(n_ex_nei_recv[0]),  i,
            //          TAG+1, MPI_COMM_WORLD, &req3[ii]);
        }
        //MPI_Waitall(n_send, req0, stat0);
        //MPI_Waitall(n_send, req1, stat1);
        MPI_Waitall(n_recv, req2, stat2);
        //MPI_Waitall(n_recv, req3, stat3);
#endif
    }
    
    void inputAdress(){
        //n_ex_ptcl_send_tot = n_ex_nei_send_tot = 0;
        n_ex_ptcl_send_tot = 0;
        for (PS::S32 i=0; i<n_send; i++){
            adr_ex_ptcl_send.at(i) = n_ex_ptcl_send_tot;
            //adr_ex_nei_send.at(i)  = n_ex_nei_send_tot;
            n_ex_ptcl_send_tot += n_ex_ptcl_send.at(i);
            //n_ex_nei_send_tot  += n_ex_nei_send.at(i);
        }
        
        //n_ex_ptcl_recv_tot = n_ex_nei_recv_tot = 0;
        n_ex_ptcl_recv_tot = 0;
        for (PS::S32 i=0; i<n_recv; i++){
            adr_ex_ptcl_recv.at(i) = n_ex_ptcl_recv_tot;
            // adr_ex_nei_recv.at(i)  = n_ex_nei_recv_tot;
            n_ex_ptcl_recv_tot += n_ex_ptcl_recv.at(i);
            //n_ex_nei_recv_tot  += n_ex_nei_recv.at(i);
        }

        ex_ptcl_send.resize(n_ex_ptcl_send_tot);
        //ex_nei_send.resize(n_ex_nei_send_tot);
        ex_ptcl_recv.resize(n_ex_ptcl_recv_tot);
        //ex_nei_recv.resize(n_ex_nei_recv_tot);
        //n_list.resize(n_ex_ptcl_recv_tot);
    }

    template <class Tpsys>
    void inputExParticleSend(Tpsys & pp,
                             NeighborList & NList){
#pragma omp parallel for
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 n_data   = n_ex_ptcl_send.at(ii);
            PS::S32 adr_ptcl = adr_ex_ptcl_send.at(ii);
            //PS::S32 adr_nei  = adr_ex_nei_send.at(ii);
            //PS::S32 n_nei = 0;
            for ( PS::S32 jj=0; jj<n_data; jj++ ) {
                PS::S32 j = ex_ptcl_send_list[ii].at(jj);

                pp[j].isSent = true;
                ex_ptcl_send.at(adr_ptcl + jj) = pp[j];
                assert( !pp[j].inDomain );
                
                //for ( PS::S32 k=0; k<pp[j].neighbor.number; k++ ) {
                //    ex_nei_send.at(adr_nei + n_nei) = NList.n_list[j][k];
                //    n_nei ++;
                //}
            }
            //assert ( n_ex_nei_send.at(ii) == n_nei );
        }
    }

    void sendRecvExParticle(NeighborList & NList,
                            PS::S32 TAG = 0){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Request req0[n_send];
        MPI_Status  stat0[n_send];
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = NList.send_rank_list[ii];       
            MPI_Isend(&ex_ptcl_send[adr_ex_ptcl_send[ii]], n_ex_ptcl_send[ii], PS::GetDataType(ex_ptcl_send[0]),
                      i, TAG+2, MPI_COMM_WORLD, &req0[ii]);
            //MPI_Isend(&ex_nei_send[adr_ex_nei_send[ii]],   n_ex_nei_send[ii],  PS::GetDataType(ex_nei_send[0]),
            //          i, TAG+3, MPI_COMM_WORLD, &req1[ii]);
        }
        MPI_Request req2[n_recv];
        MPI_Status  stat2[n_recv];
        for ( PS::S32 ii=0; ii<n_recv; ii++ ) {
            PS::S32 i = NList.recv_rank_list[ii];
            MPI_Irecv(&ex_ptcl_recv[adr_ex_ptcl_recv[ii]], n_ex_ptcl_recv[ii], PS::GetDataType(ex_ptcl_recv[0]),
                      i, TAG+2, MPI_COMM_WORLD, &req2[ii]);
            //MPI_Irecv(&ex_nei_recv[adr_ex_nei_recv[ii]],   n_ex_nei_recv[ii],  PS::GetDataType(ex_nei_recv[0]),
            //          i, TAG+3, MPI_COMM_WORLD, &req3[ii]);
        }
        //MPI_Waitall(n_send, req0, stat0);
        //MPI_Waitall(n_send, req1, stat1);
        MPI_Waitall(n_recv, req2, stat2);
        //MPI_Waitall(n_recv, req3, stat3);
#endif
    }

    /*
    void inputNeighborListOfExParticleRecv() { 
#pragma omp parallel for
        for ( PS::S32 ii=0; ii<n_recv; ii++ ) {
            PS::S32 n_data = n_ex_ptcl_recv.at(ii);
            PS::S32 adr_ptcl = adr_ex_ptcl_recv.at(ii);
            //PS::S32 n_nei    = adr_ex_nei_recv.at(ii);
            //for ( PS::S32 jj=0; jj<n_data; jj++ ) {
                //n_list.at(adr_ptcl + jj) = &(ex_nei_recv.at(n_nei));
                //n_nei += ex_ptcl_recv.at(adr_ptcl + jj).neighbor.number;
                //assert ( ex_ptcl_recv.at(adr_ptcl + jj).isSent );
            //}
            //if ( ii+1<n_recv ) assert ( adr_ex_nei_recv.at(ii+1) == n_nei );
        }
    }
    */

    void returnExParticle(NeighborList & NList,
                          PS::S32 TAG = 0){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Request req0[n_send];
        MPI_Status  stat0[n_send];
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = NList.send_rank_list[ii];       
            MPI_Irecv(&ex_ptcl_send[adr_ex_ptcl_send[ii]], n_ex_ptcl_send[ii], PS::GetDataType(ex_ptcl_send[0]),
                      i, TAG+4, MPI_COMM_WORLD, &req0[ii]);
            //MPI_Irecv(&ex_nei_send[adr_ex_nei_send[ii]],   n_ex_nei_send[ii],  PS::GetDataType(ex_nei_send[0]),
            //          i, TAG+5, MPI_COMM_WORLD, &req1[ii]);
        }
        MPI_Request req2[n_recv];
        MPI_Status  stat2[n_recv];
        for ( PS::S32 ii=0; ii<n_recv; ii++ ) {
            PS::S32 i = NList.recv_rank_list[ii];
            MPI_Isend(&ex_ptcl_recv[adr_ex_ptcl_recv[ii]], n_ex_ptcl_recv[ii], PS::GetDataType(ex_ptcl_recv[0]),
                      i, TAG+4, MPI_COMM_WORLD, &req2[ii]);
            //MPI_Isend(&ex_nei_recv[adr_ex_nei_recv[ii]],   n_ex_nei_recv[ii],  PS::GetDataType(ex_nei_recv[0]),
            //          i, TAG+5, MPI_COMM_WORLD, &req3[ii]);
        }
        MPI_Waitall(n_send, req0, stat0);
        //MPI_Waitall(n_send, req1, stat1);
        //MPI_Waitall(n_recv, req2, stat2);
        //MPI_Waitall(n_recv, req3, stat3);
#endif
    }

    template <class Tpsys>
    void outputExParticleSend(Tpsys & pp,
                              NeighborList & NList){
#pragma omp parallel for
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 n_data   = n_ex_ptcl_send.at(ii);
            PS::S32 adr_ptcl = adr_ex_ptcl_send.at(ii);
            for ( PS::S32 jj=0; jj<n_data; jj++ ) {
                PS::S32 j = ex_ptcl_send_list[ii].at(jj);
                PS::S32 id_pre = pp[j].id;
                
                pp[j] = ex_ptcl_send.at(adr_ptcl + jj);
                if (!pp[j].isDead) assert( pp[j].id == id_pre );
            }
        }
    }
    

};
