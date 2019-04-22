#pragma once

class ExPair{
public:
    PS::S32 id_in;
    PS::S32 id_out;
    PS::S32 id_cluster;
    
    PS::S32 * rank_list;
    
    static PS::S32 size;
    static PS::S32 rem;
    static PS::S32 n_bit;

    static void initialize() {
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        
        n_bit = 8 * sizeof(PS::S32);
        size = (PS::S32)std::ceil((PS::F64)n_proc/n_bit);
        rem  = n_bit*size - n_proc;
    }
    static PS::S32 getSize() { return size+3; }
    
    ExPair(){
        //PS::S32 myrank = PS::Comm::getRank();
        id_in = id_out = id_cluster = 0;
        
        rank_list = new PS::S32[size];
        for ( PS::S32 i=0; i<size; i++ ) rank_list[i] = 0;
        //setFlag(myrank);
    }
    ExPair(PS::S32 id_in0,
           PS::S32 id_out0,
           PS::S32 id_cluster0){
        //PS::S32 myrank = PS::Comm::getRank();
        id_in  = id_in0;
        id_out = id_out0;
        id_cluster = id_cluster0;

        rank_list = new PS::S32[size];
        for ( PS::S32 i=0; i<size; i++ ) rank_list[i] = 0;
        //setFlag(myrank);
    }
    ExPair(const ExPair & ep){
        id_in  = ep.id_in;
        id_out = ep.id_out;
        id_cluster = ep.id_cluster;

        rank_list = new PS::S32[size];
        for ( PS::S32 i=0; i<size; i++ ) rank_list[i] = ep.rank_list[i];
    }
    ExPair &operator=(const ExPair & ep){
        if ( this != &ep ){
            id_in  = ep.id_in;
            id_out = ep.id_out;
            id_cluster = ep.id_cluster;
        
            for ( PS::S32 i=0; i<size; i++ ) this->rank_list[i] = ep.rank_list[i];
        }
        return *this;
    }
    
    ~ExPair(){
        delete [] rank_list;
    }

    PS::S32 getId() const {
        return id_in;
    }
    std::pair<PS::S32,PS::S32> getPair() const {
        return std::make_pair(id_in, id_out);
    }
    PS::S32 getIdCluster() const {
        return id_cluster;
    }
    PS::S32 setIdCluster(PS::S32 id_cluster0) {
        return id_cluster = id_cluster0;
    }

    PS::S32 input(PS::S32 * inp){
        id_in  = inp[1];
        id_out = inp[0];
        id_cluster = inp[2];
        for ( PS::S32 i=0; i<size; i++ ) rank_list[i] = inp[i+3];
        return size+3;
    }
    PS::S32 output(PS::S32 * outp){
        outp[0] = id_in;
        outp[1] = id_out;
        outp[2] = id_cluster;
        for ( PS::S32 i=0; i<size; i++ ) outp[i+3] = rank_list[i];
        return size+3;
    }

    bool checkFlag(const PS::S32 i) const {
        PS::S32 n = i / n_bit;
        PS::S32 ii = i - n_bit * n;
        return rank_list[n] & (1<<ii);
    }
    void setFlag(const PS::S32 i) {
        PS::S32 n = i / n_bit;
        PS::S32 ii = i - n_bit * n;
        rank_list[n] |= (1<<ii);
    }
    void unsetFlag(const PS::S32 i) {
        PS::S32 n = i / n_bit;
        PS::S32 ii = i - n_bit * n;
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

    void show(){
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        std::cout << PS::Comm::getRank() << "\t" << id_in << "\t" << id_out << "\t" << id_cluster << "\t";
        for ( PS::S32 i=0; i<n_proc; i++ ) std::cout << (checkFlag(i));
        std::cout << std::endl;
    }
};

PS::S32 ExPair::size;
PS::S32 ExPair::rem;
PS::S32 ExPair::n_bit;


class NeighborList{
public:
    std::vector<std::vector<PS::S32> > n_list;
    std::map<PS::S32, PS::S32> id_map;
    
    std::vector<std::pair<PS::S32,PS::S32> > ex_list;
    std::vector<std::pair<PS::S32,PS::S32> > ex_adr_list;
    std::vector<PS::S32> connected_list;
    std::vector<std::vector<ExPair> > ex_data;
    std::map<std::pair<PS::S32,PS::S32>, std::pair<PS::S32, PS::S32> > ex_data_map;

    std::vector<std::vector<PS::S32> > recv_list;
    std::vector<std::vector<PS::S32> > send_list;
    std::vector<PS::S32> recv_rank_list;
    std::vector<PS::S32> send_rank_list;


    std::vector<PS::S32> & operator[](PS::S32 i){ return n_list[i]; }
    
    void initialize() {
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        
        n_list.clear();
        ex_data.resize(n_proc);
        recv_list.resize(n_proc);
        send_list.resize(n_proc);

#pragma omp parallel for
        for (PS::S32 i=0; i<n_proc; i++){
            ex_data[i].clear();
            recv_list[i].clear();
            send_list[i].clear();
        }
    }
    
    template <class Tpsys>
    void initializeList(Tpsys & pp) {
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        const PS::S32 n_loc  = pp.getNumberOfParticleLocal();

        n_list.clear();
        id_map.clear();
        ex_list.clear();
        ex_adr_list.clear();
        connected_list.clear();
        ex_data_map.clear();
        recv_rank_list.clear();
        send_rank_list.clear();
            
#pragma omp parallel for
        for ( PS::S32 i=0; i<n_proc; i++ ){
            ex_data[i].clear();
            recv_list[i].clear();
            send_list[i].clear();
        }
        
        n_list.resize(n_loc);
#pragma omp parallel for
        for(PS::S32 i=0; i<n_loc; i++) n_list.at(i).clear();
    }

    ExPair & getExData(std::pair<PS::S32, PS::S32> adr) {
        return ex_data[adr.first][adr.second];
    }
    
    PS::S32 getNumberOfRankSend() const { return send_rank_list.size(); }
    PS::S32 getNumberOfRankRecv() const { return recv_rank_list.size(); }
    PS::S32 getNumberOfRankConnected() const { return connected_list.size(); }
    PS::S32 getNumberOfPairConnected(const PS::S32 ii) const { return ex_data[connected_list.at(ii)].size(); }

    template <class Tpsys>
    void addNeighbor(Tpsys & pp,
                     PS::S32 i,
                     PS::S32 j_id,
                     PS::S32 j_rank) {
        n_list[i].push_back(j_id);
        pp[i].neighbor ++;

        pp[i].id_cluster = std::min(pp[i].id_cluster, j_id);

        if ( j_rank != pp[i].myrank ) {
#pragma omp critical
            {
                ex_list.push_back(std::make_pair(pp[i].id, j_id));
                ex_adr_list.push_back(std::make_pair(j_rank, ex_data.at(j_rank).size()));
                ex_data_map[std::make_pair(pp[i].id, j_id)] = std::make_pair(j_rank, ex_data.at(j_rank).size());
                ExPair ex_pair(pp[i].id, j_id, pp[i].id_cluster);
                ex_pair.setFlag(pp[i].myrank);
                ex_pair.setFlag(j_rank);
                ex_data.at(j_rank).push_back(ex_pair);
            }
            pp[i].inDomain = false;
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

    template <class Tpsys>
    void makeIdMap(Tpsys & pp){
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();
        assert( (PS::S32)(n_list.size()) == n_loc );
        
        for(PS::S32 i=0; i<n_loc; i++){
            assert( pp[i].neighbor == (PS::S32)(n_list[i].size()) );
            id_map[pp[i].id] = i;
        }
    }

    template <class Tpsys>
    void createNeighborCluster(Tpsys & pp){
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();

        PS::S32 j_id_cluster = 0;
        PS::S32 id_cluster[n_loc];
        bool check = true;
        while( check ){
            check = false;
#pragma omp parallel for
            for(PS::S32 i=0; i<n_loc; i++){
                PS::S32 j_id = 0;
                PS::S32 nei = 0;
                nei = pp[i].neighbor;
                id_cluster[i] = pp[i].id_cluster;
                
                if(nei == 0) continue;
                for(PS::S32 j=0; j<nei; j++){
                    auto itr = id_map.find(n_list[i].at(j));
                    if ( itr == id_map.end() ) continue;
                    j_id = itr->second;
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

    template <class Tpsys>
    void inputExData(Tpsys & pp){
        const PS::S32 n_out = ex_list.size();
        for ( PS::S32 j=0; j<n_out; j++ ){
            std::pair<PS::S32,PS::S32> pair = ex_list.at(j);
            std::pair<PS::S32,PS::S32> ex_adr = ex_adr_list.at(j);

            assert( getExData(ex_adr).getId() == pair.first );
            getExData(ex_adr).setIdCluster(pp[id_map.at(pair.first)].id_cluster);

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
                        PS::S32 TAG,
                        PS::S32** & ex_data_send,
                        PS::S32** & ex_data_recv){
        //const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        const PS::S32 n_send = connected_list.size();
        //PS::S32 ** ex_data_send = new PS::S32*[n_send];
        //PS::S32 ** ex_data_recv = new PS::S32*[n_send];
        
        //for ( PS::S32 ii=0; ii<n_send; ii++ ) {
        //    PS::S32 i = connected_list.at(ii);
        //    PS::S32 n_size = ex_data[i].size() * ExPair::getSize();
            
        //    ex_data_send[ii] = new PS::S32[n_size];
        //    ex_data_recv[ii] = new PS::S32[n_size];
        //}
#pragma omp parallel for
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = connected_list.at(ii);
            PS::S32 n_data = ex_data[i].size();
            
            PS::S32 jj = 0;
            for ( PS::S32 j=0; j<n_data; j++ ) {
                jj += ex_data[i][j].output(&ex_data_send[ii][jj]);
            }
        }

        MPI_Request req0[n_send],  req1[n_send];
        MPI_Status  stat0[n_send], stat1[n_send];
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = connected_list.at(ii);
            PS::S32 n_size = ex_data[i].size() * ExPair::getSize();
            
            MPI_Isend(&ex_data_send[ii][0], n_size, PS::GetDataType(*ex_data_send[ii]), i, TAG, MPI_COMM_WORLD, &req0[ii]);
            MPI_Irecv(&ex_data_recv[ii][0], n_size, PS::GetDataType(*ex_data_recv[ii]), i, TAG, MPI_COMM_WORLD, &req1[ii]);
        }
        MPI_Waitall(n_send, req0, stat0);
        MPI_Waitall(n_send, req1, stat1);
        

        bool check = false;
        
#pragma omp parallel for reduction (||:check)
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = connected_list.at(ii);
            PS::S32 n_data = ex_data[i].size();

            PS::S32 jj = 0;
            for ( PS::S32 j=0; j<n_data; j++ ) {
                ExPair recv_pair;
                jj += recv_pair.input(&ex_data_recv[ii][jj]);

                std::pair<PS::S32,PS::S32> adr = ex_data_map.at(recv_pair.getPair());
                assert ( adr.first == i );
                assert ( recv_pair.getPair() == getExData(adr).getPair() );
                bool check_1 = getExData(adr).exchange(recv_pair);
                check = check || check_1;
#pragma omp critical
                {
                    //getExData(adr).show();
                    PS::S32 i_loc = id_map.at(getExData(adr).getId());
                    pp[i_loc].id_cluster = std::min(pp[i_loc].id_cluster, getExData(adr).getIdCluster());
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
        std::vector<PS::S32> ex_cluster;
        std::vector<std::pair<PS::S32,PS::S32> > ex_cluster_adr;
        ex_cluster.clear();
        ex_cluster_adr.clear();
        
        for ( PS::S32 ii=0; ii<n_ptcl; ii++ ) {
            //std::pair<PS::S32,PS::S32> pair = ex_list.at(ii);
            std::pair<PS::S32,PS::S32> adr =  ex_adr_list.at(ii);
            PS::S32 id_cluster = getExData(adr).id_cluster;

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
    
};



template <class Tp>
class ExParticleSystem {
 public :
    PS::S32 n_send;
    PS::S32 n_recv;
    PS::S32 n_ex_ptcl_send_tot;
    PS::S32 n_ex_nei_send_tot;
    PS::S32 n_ex_ptcl_recv_tot;
    PS::S32 n_ex_nei_recv_tot;

    std::vector<Tp>      ex_ptcl_send;
    std::vector<PS::S32> ex_nei_send;
    std::vector<Tp>      ex_ptcl_recv;
    std::vector<PS::S32> ex_nei_recv;

    std::vector<std::vector<PS::S32> > ex_ptcl_send_list;
    std::vector<PS::S32*> n_list;

    std::vector<PS::S32> n_ex_ptcl_send;
    std::vector<PS::S32> n_ex_nei_send;
    std::vector<PS::S32> n_ex_ptcl_recv;
    std::vector<PS::S32> n_ex_nei_recv;

    std::vector<PS::S32> adr_ex_ptcl_send;
    std::vector<PS::S32> adr_ex_nei_send;
    std::vector<PS::S32> adr_ex_ptcl_recv;
    std::vector<PS::S32> adr_ex_nei_recv;

    
    Tp & operator[](PS::S32 i){ return ex_ptcl_recv[i]; }
    PS::S32 getNumberOfParticleLocal() const { return n_ex_ptcl_recv_tot; }

    void initialize() {
        n_send = n_recv = 0;
        n_ex_ptcl_send_tot = n_ex_ptcl_recv_tot = 0;
        n_ex_nei_send_tot  = n_ex_nei_recv_tot  = 0;
        
        ex_ptcl_send.clear();
        ex_nei_send.clear();
        ex_ptcl_recv.clear();
        ex_nei_recv.clear();
        
        ex_ptcl_send_list.clear();
    
        n_ex_ptcl_send.clear();
        n_ex_nei_send.clear();
        n_ex_ptcl_recv.clear();
        n_ex_nei_recv.clear();

        adr_ex_ptcl_send.clear();
        adr_ex_nei_send.clear();
        adr_ex_ptcl_recv.clear();
        adr_ex_nei_recv.clear();
    }
    

    void resize(PS::S32 n_send0,
                PS::S32 n_recv0){
        n_send = n_send0;
        n_ex_ptcl_send.resize(n_send);
        n_ex_nei_send.resize(n_send);
        adr_ex_ptcl_send.resize(n_send);
        adr_ex_nei_send.resize(n_send);

        ex_ptcl_send_list.resize(n_send);
#pragma omp parallel for
        for ( PS::S32 i=0; i<n_send; i++ ) ex_ptcl_send_list[i].clear();
        
        n_recv = n_recv0;
        n_ex_ptcl_recv.resize(n_recv);
        n_ex_nei_recv.resize(n_recv);
        adr_ex_ptcl_recv.resize(n_recv);
        adr_ex_nei_recv.resize(n_recv);
    }

    PS::S32 getNumberOfParticleSend() const { return n_ex_ptcl_send_tot; }
    PS::S32 getNumberOfParticleRecv() const { return n_ex_ptcl_recv_tot; }
    PS::S32 getNumberOfNeighborSend() const { return n_ex_nei_send_tot; }
    PS::S32 getNumberOfNeighborRecv() const { return n_ex_nei_recv_tot; }

    template <class Tpsys>
    void inputNumberOfExParticleSend(Tpsys & pp,
                                     NeighborList & NList){
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();
        
#pragma omp parallel for
        for ( PS::S32 ii=0; ii<n_send; ii++ ) n_ex_ptcl_send[ii] = n_ex_nei_send[ii] = 0;

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
                                    n_ex_nei_send[jj] += pp[i].neighbor;
                                    assert ( pp[i].neighbor == (PS::S32)(NList.n_list[i].size()) );
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
        MPI_Request req0[n_send],  req1[n_send];
        MPI_Status  stat0[n_send], stat1[n_send];
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = NList.send_rank_list[ii];
            MPI_Isend(&n_ex_ptcl_send[ii], 1, PS::GetDataType(n_ex_ptcl_send[0]), i,
                      TAG,   MPI_COMM_WORLD, &req0[ii]);
            MPI_Isend(&n_ex_nei_send[ii],  1, PS::GetDataType(n_ex_nei_send[0]),  i,
                      TAG+1, MPI_COMM_WORLD, &req1[ii]);
        }
        MPI_Request req2[n_recv],  req3[n_recv];
        MPI_Status  stat2[n_recv], stat3[n_recv];
        for ( PS::S32 ii=0; ii<n_recv; ii++ ) {
            PS::S32 i = NList.recv_rank_list[ii];
            MPI_Irecv(&n_ex_ptcl_recv[ii], 1, PS::GetDataType(n_ex_ptcl_recv[0]), i,
                      TAG,   MPI_COMM_WORLD, &req2[ii]);
            MPI_Irecv(&n_ex_nei_recv[ii],  1, PS::GetDataType(n_ex_nei_recv[0]),  i,
                      TAG+1, MPI_COMM_WORLD, &req3[ii]);
        }
        MPI_Waitall(n_send, req0, stat0);
        MPI_Waitall(n_send, req1, stat1);
        MPI_Waitall(n_recv, req2, stat2);
        MPI_Waitall(n_recv, req3, stat3);
    }
    
    void inputAdress(){
        n_ex_ptcl_send_tot = n_ex_nei_send_tot = 0;
        for (PS::S32 i=0; i<n_send; i++){
            adr_ex_ptcl_send.at(i) = n_ex_ptcl_send_tot;
            adr_ex_nei_send.at(i)  = n_ex_nei_send_tot;
            n_ex_ptcl_send_tot += n_ex_ptcl_send.at(i);
            n_ex_nei_send_tot  += n_ex_nei_send.at(i);
        }
        
        n_ex_ptcl_recv_tot = n_ex_nei_recv_tot = 0;
        for (PS::S32 i=0; i<n_recv; i++){
            adr_ex_ptcl_recv.at(i) = n_ex_ptcl_recv_tot;
            adr_ex_nei_recv.at(i)  = n_ex_nei_recv_tot;
            n_ex_ptcl_recv_tot += n_ex_ptcl_recv.at(i);
            n_ex_nei_recv_tot  += n_ex_nei_recv.at(i);
        }

        ex_ptcl_send.resize(n_ex_ptcl_send_tot);
        ex_nei_send.resize(n_ex_nei_send_tot);
        ex_ptcl_recv.resize(n_ex_ptcl_recv_tot);
        ex_nei_recv.resize(n_ex_nei_recv_tot);
        n_list.resize(n_ex_ptcl_recv_tot);
    }

    template <class Tpsys>
    void inputExParticleSend(Tpsys & pp,
                             NeighborList & NList){
#pragma omp parallel for
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 n_data   = n_ex_ptcl_send.at(ii);
            PS::S32 adr_ptcl = adr_ex_ptcl_send.at(ii);
            PS::S32 adr_nei  = adr_ex_nei_send.at(ii);
            PS::S32 n_nei = 0;
            for ( PS::S32 jj=0; jj<n_data; jj++ ) {
                PS::S32 j = ex_ptcl_send_list[ii].at(jj);

                pp[j].isSent = true;
                ex_ptcl_send.at(adr_ptcl + jj) = pp[j];
                assert( !pp[j].inDomain );
                
                for ( PS::S32 k=0; k<pp[j].neighbor; k++ ) {
                    ex_nei_send.at(adr_nei + n_nei) = NList.n_list[j].at(k);
                    n_nei ++;
                }
            }
            assert ( n_ex_nei_send.at(ii) == n_nei );
        }
    }

    void sendRecvExParticle(NeighborList & NList,
                            PS::S32 TAG = 0){
        MPI_Request req0[n_send],  req1[n_send];
        MPI_Status  stat0[n_send], stat1[n_send];
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = NList.send_rank_list[ii];       
            MPI_Isend(&ex_ptcl_send[adr_ex_ptcl_send[ii]], n_ex_ptcl_send[ii], PS::GetDataType(ex_ptcl_send[0]),
                      i, TAG+2, MPI_COMM_WORLD, &req0[ii]);
            MPI_Isend(&ex_nei_send[adr_ex_nei_send[ii]],   n_ex_nei_send[ii],  PS::GetDataType(ex_nei_send[0]),
                      i, TAG+3, MPI_COMM_WORLD, &req1[ii]);
        }
        MPI_Request req2[n_recv],  req3[n_recv];
        MPI_Status  stat2[n_recv], stat3[n_recv];
        for ( PS::S32 ii=0; ii<n_recv; ii++ ) {
            PS::S32 i = NList.recv_rank_list[ii];
            MPI_Irecv(&ex_ptcl_recv[adr_ex_ptcl_recv[ii]], n_ex_ptcl_recv[ii], PS::GetDataType(ex_ptcl_recv[0]),
                      i, TAG+2, MPI_COMM_WORLD, &req2[ii]);
            MPI_Irecv(&ex_nei_recv[adr_ex_nei_recv[ii]],   n_ex_nei_recv[ii],  PS::GetDataType(ex_nei_recv[0]),
                      i, TAG+3, MPI_COMM_WORLD, &req3[ii]);
        }
        MPI_Waitall(n_send, req0, stat0);
        MPI_Waitall(n_send, req1, stat1);
        MPI_Waitall(n_recv, req2, stat2);
        MPI_Waitall(n_recv, req3, stat3);
    }

    void inputNeighborListOfExParticleRecv() { 
#pragma omp parallel for
        for ( PS::S32 ii=0; ii<n_recv; ii++ ) {
            PS::S32 n_data = n_ex_ptcl_recv.at(ii);
            PS::S32 adr_ptcl = adr_ex_ptcl_recv.at(ii);
            PS::S32 n_nei    = adr_ex_nei_recv.at(ii);
            for ( PS::S32 jj=0; jj<n_data; jj++ ) {
                n_list.at(adr_ptcl + jj) = &(ex_nei_recv.at(n_nei));
                n_nei += ex_ptcl_recv.at(adr_ptcl + jj).neighbor;
                assert ( ex_ptcl_recv.at(adr_ptcl + jj).isSent );
            }
            if ( ii+1<n_recv ) assert ( adr_ex_nei_recv.at(ii+1) == n_nei );
        }
    }

    void returnExParticle(NeighborList & NList,
                          PS::S32 TAG = 0){
        MPI_Request req0[n_send],  req1[n_send];
        MPI_Status  stat0[n_send], stat1[n_send];
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 i = NList.send_rank_list[ii];       
            MPI_Irecv(&ex_ptcl_send[adr_ex_ptcl_send[ii]], n_ex_ptcl_send[ii], PS::GetDataType(ex_ptcl_send[0]),
                      i, TAG+4, MPI_COMM_WORLD, &req0[ii]);
            MPI_Irecv(&ex_nei_send[adr_ex_nei_send[ii]],   n_ex_nei_send[ii],  PS::GetDataType(ex_nei_send[0]),
                      i, TAG+5, MPI_COMM_WORLD, &req1[ii]);
        }
        MPI_Request req2[n_recv],  req3[n_recv];
        MPI_Status  stat2[n_recv], stat3[n_recv];
        for ( PS::S32 ii=0; ii<n_recv; ii++ ) {
            PS::S32 i = NList.recv_rank_list[ii];
            MPI_Isend(&ex_ptcl_recv[adr_ex_ptcl_recv[ii]], n_ex_ptcl_recv[ii], PS::GetDataType(ex_ptcl_recv[0]),
                      i, TAG+4, MPI_COMM_WORLD, &req2[ii]);
            MPI_Isend(&ex_nei_recv[adr_ex_nei_recv[ii]],   n_ex_nei_recv[ii],  PS::GetDataType(ex_nei_recv[0]),
                      i, TAG+5, MPI_COMM_WORLD, &req3[ii]);
        }
        MPI_Waitall(n_send, req0, stat0);
        MPI_Waitall(n_send, req1, stat1);
        MPI_Waitall(n_recv, req2, stat2);
        MPI_Waitall(n_recv, req3, stat3);
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
