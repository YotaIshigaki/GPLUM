//***************************************************************************************
//  This is the pair-ID data table class for calculate interactions.
//***************************************************************************************
#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>


//--- grobal object of parameter table
namespace IntraPair {

    template <class Tid>
    void check_invalid_connect(const Tid id_from, const Tid id_to){
        if(id_from == id_to){
            std::ostringstream oss;
            oss << "ID=" << id_from << " has connection to itself. invalid structure.";
            throw std::invalid_argument(oss.str());
        }
    }

    template <class Tptcl, class Tid>
    void check_nullptr_ptcl(const Tptcl *ptr_next,
                            const Tid   &id_from,
                            const Tid   &id_to   ){
        if(ptr_next == NULL){
            std::ostringstream oss;

            oss << "  next EPJ of id= "     << id_to   << " was not found in tree.\n";
            oss << "    searched from id= " << id_from << "\n";
            //throw std::out_of_range(oss.str());
        }
    }

   //! @brief dummy function for PS::TreeForForce::calcForceAll<>. use for communication without any force calculation.
    struct dummy_func {
        template <class TResult, class TEpi, class TEpj>
        void operator () (TEpi*    ptr_i, PS::S32 n_i,
                          TEpj*    ptr_j, PS::S32 n_j,
                          TResult* result             ){
            return;
        }
    };

    //--- generic function for construct table
    /**
    * @brief functor for construction the mask list for intermolecular interaction.
    * @tparam <Tid> data type of ID.
    * @tparam <_GetConnect> functor to get the next node list from particle.
    *                       MUST have:  <Container> operator () (const Tptcl &ptcl){}.
    *                                   <Container> is any container class with begin() and end() interface.
    *                                   such as std::vector<> or MD_EXT::fixed_vector<>.
    */
    template <class Tid, class _GetConnect>
    class IntraMaskMaker{
    private:
        std::vector<Tid>                     node_list_buff;
        std::vector<std::pair<Tid, PS::S32>> mask_order_buff;

        //
        //  NOTICE: PS::TreeForForce<>::getEpjFromId() is not thread safe.
        //     it is no meaning to modify accepting multi thread.

        //std::vector< std::vector<Tid> >                     node_list_buff;
        //std::vector< std::vector<std::pair<Tid, PS::S32>> > mask_order_buff;
        //
        //std::vector<Tid>& node_list(){
        //    #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
        //    const PS::S32 thread_id = omp_get_thread_num();
        //    #else
        //    const PS::S32 thread_id = 0;
        //    #endif
        //    return this->node_list_buff.at(thread_id);
        //}
        //std::vector<std::pair<Tid, PS::S32>>& mask_order(){
        //    #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
        //    const PS::S32 thread_id = omp_get_thread_num();
        //    #else
        //    const PS::S32 thread_id = 0;
        //    #endif
        //    return this->mask_order_buff.at(thread_id);
        //}
        std::vector<Tid>&                     node_list() { return this->node_list_buff;  }
        std::vector<std::pair<Tid, PS::S32>>& mask_order(){ return this->mask_order_buff; }

    public:
        void clear(){
            this->node_list().clear();
            this->mask_order().clear();
        }
        IntraMaskMaker(){
        //    #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
        //    const PS::S32 n_threads = omp_get_max_threads();
        //    #else
        //    const PS::S32 n_threads = 1;
        //    #endif
        //    this->node_list_buff.resize(n_threads);
        //    this->mask_order_buff.resize(n_threads);

            this->clear();
            this->node_list().reserve(4);
            this->mask_order().reserve(16);
        }
        ~IntraMaskMaker() = default;

        /**
        * @brief function for construction the mask list.
        * @tparam <Tptcl>   Particle class.
        *                   MUST have <Tid> getId(){} function.
        * @tparam <Ttree>   PS::TreeForForce<>
        * @tparam <Tcoef>   data type for mask.
        *                   MUST have <Tcoef> setId(const Tid &id){} function.
        * @tparam <Tresult> result data type.
        *                   MUST have void push_back(){} function.
        * @param[in]    ptcl target particle. the result must be used for this particle only.
        * @param[inout] tree PS::TreeForForce object. MUST NOT add const property.
        * @param[in]    mask_param the mask range = mask_param.size(),
        *                          1st level mask = mask_param[0],
        *                          2nd level mask = mask_param[1]...
        * @param[out]   result the intramolecular mask for ptcl.
        */
        template <class Tptcl, class Ttree, class Tcoef,
                  class Tresult>
        void operator () (const Tptcl              &ptcl,
                                Ttree              &tree,
                          const std::vector<Tcoef> &mask_param,
                                Tresult            &result     ){

            if(mask_param.size() < 1) return;
            if(this->node_list().size() == 0){
                this->node_list().push_back(ptcl.getId());
                result.clear();
            }

            const Tid     id_root   = this->node_list().front();
            const PS::S32 now_order = this->node_list().size() - 1;

            //--- search in connection lists gotten from get_func()
            const auto connect_list = _GetConnect()(ptcl);
            for(const Tid id_next : connect_list){

                check_invalid_connect(ptcl.getId(), id_next);

                //--- search in next node
                if(mask_param.size() - now_order > 1){
                    this->node_list().push_back(id_next);
                    const auto* ptr_next = tree.getEpjFromId(id_next);

                    check_nullptr_ptcl(ptr_next, id_root, id_next);

                    this->operator() (*ptr_next,
                                       tree,
                                       mask_param,
                                       result     );
                    this->node_list().pop_back();
                }

                //--- circulation check
                if( id_next == id_root) continue;

                //--- dupurication check
                auto itr = std::find_if(this->mask_order().begin(),
                                        this->mask_order().end(),
                                        [id_next](std::pair<Tid, PS::S32> &order){ return order.first == id_next; } );
                if(itr != this->mask_order().end()){
                    //--- record minimum order
                    itr->second = std::min(itr->second, now_order);
                } else {
                    //--- add mask_order list
                    this->mask_order().push_back(std::make_pair(id_next, now_order));
                }
            }

            //--- convert mask_order -> result
            if(this->node_list().size() == 1){
                result.clear();
                result.reserve(this->mask_order().size());
                for(const auto& mask_tgt : this->mask_order()){
                    auto param = mask_param.at(mask_tgt.second);
                    result.push_back(param.setId(mask_tgt.first));
                }

                //--- clear internal buffer
                this->clear();
            }
        }
    };

    /**
    * @brief functor for construction the angle pair list for intramolecular interaction.
    * @tparam <Tid> data type of ID.
    * @tparam <_GetConnect> functor to get the next node list from particle.
    *                       MUST have:  <Container> operator () (const Tptcl &ptcl){}.
    *                                   <Container> is any container class with begin() and end() interface.
    *                                   such as std::vector<> or MD_EXT::fixed_vector<>.
    */
    template <class Tid, class _GetConnect>
    class AngleListMaker{
    public:

        /**
        * @brief function for construction the angle pair list.
        * @tparam <Tptcl>   Particle class.
        *                   MUST have <Tid> getId(){} function.
        * @tparam <Ttree>   PS::TreeForForce<>
        * @tparam <Tresult> result data type.
        *                   container element is std::tuple<Tid, Tid, Tid>.
        *                   MUST have void push_back(){} function.
        * @param[in]    ptcl target particle. the result must be used for this particle only.
        * @param[inout] tree PS::TreeForForce object. MUST NOT add const property.
        * @param[out]   result the angle pair list for ptcl.
        */
        template <class Tptcl, class Ttree,
                  class Tresult            >
        void operator () (const Tptcl   &ptcl,
                                Ttree   &tree,
                                Tresult &result){
            result.clear();
            const Tid  id_root   = ptcl.getId();
            //--- I-J
            const auto connect_i = _GetConnect()(ptcl);
            for(size_t j=0; j<connect_i.size(); ++j){
                const Tid id_j = connect_i[j];

                //std::cout << "id_j = " << id_j << "  , id_root = " << id_root << std::endl;
                check_invalid_connect(id_root, id_j);

                //--- I-J-K
                const auto* ptr_j = tree.getEpjFromId(id_j);

                check_nullptr_ptcl(ptr_j, id_root, id_j);

                const auto connect_j = _GetConnect()(*ptr_j);
                for(size_t k=0; k<connect_j.size(); ++k){
                    const Tid id_k = connect_j[k];

                    check_invalid_connect(id_j, id_k);
                    if(id_k == id_root) continue;

                    result.push_back( std::make_tuple(id_root,
                                                      id_j,
                                                      id_k    ) );
                }

                //--- J-I-K
                for(size_t k=j+1; k<connect_i.size(); ++k){
                    const Tid id_k = connect_i[k];

                    check_invalid_connect(id_root, id_k);
                    if(id_k == id_j) continue;

                    result.push_back( std::make_tuple(id_j,
                                                      id_root,
                                                      id_k    )  );
                }
            }
        }
    };

    /**
    * @brief functor for construction the torsion pair list for intramolecular interaction.
    * @tparam <Tid> data type of ID.
    * @tparam <_GetConnect> functor to get the next node list from particle.
    *                       MUST have:  <Container> operator () (const Tptcl &ptcl){}.
    *                                   <Container> is any container class with begin() and end() interface.
    *                                   such as std::vector<> or MD_EXT::fixed_vector<>.
    */
    template <class Tid, class _GetConnect>
    class TorsionListMaker{
    public:

        /**
        * @brief function for construction the torsion pair list.
        * @tparam <Tptcl>   Particle class.
        *                   MUST have <Tid> getId(){} function.
        * @tparam <Ttree>   PS::TreeForForce<>
        * @tparam <Tresult> result data type.
        *                   container element is std::tuple<Tid, Tid, Tid, Tid>.
        *                   MUST have void push_back(){} function.
        * @param[in]    ptcl target particle. the result must be used for this particle only.
        * @param[inout] tree PS::TreeForForce object. MUST NOT add const property.
        * @param[out]   dihedral_result the dihedral torsion pair list for ptcl.
        * @param[out]   improper_result the improper torsion pair list for ptcl.
        */
        template <class Tptcl, class Ttree,
                  class Tresult            >
        void operator () (const Tptcl   &ptcl,
                                Ttree   &tree,
                                Tresult &dihedral_result,
                                Tresult &improper_result){
            dihedral_result.clear();
            improper_result.clear();
            const Tid  id_root   = ptcl.getId();
            //--- I-J
            const auto connect_i = _GetConnect()(ptcl);
            for(size_t j=0; j<connect_i.size(); ++j){
                const Tid id_j = connect_i[j];
                check_invalid_connect(id_root, id_j);

                //--- I-J-K
                auto* ptr_j = tree.getEpjFromId(id_j);
                check_nullptr_ptcl(ptr_j, id_root, id_j);

                const auto connect_j = _GetConnect()(*ptr_j);
                for(size_t k=0; k<connect_j.size(); ++k){
                    const Tid id_k = connect_j[k];
                    check_invalid_connect(id_j, id_k);
                    if(id_k == id_root) continue;

                    //--- dihedral, I-J-K-L
                    const auto* ptr_k = tree.getEpjFromId(id_k);
                    check_nullptr_ptcl(ptr_k, id_j, id_k);

                    const auto connect_k = _GetConnect()(*ptr_k);
                    for(size_t l=0; l<connect_k.size(); ++l){
                        const Tid id_l = connect_k[l];
                        check_invalid_connect(id_k, id_l);
                        if(id_l == id_j    ||
                           id_l == id_root   ) continue;

                        dihedral_result.push_back( std::make_tuple(id_root,
                                                                   id_j,
                                                                   id_k,
                                                                   id_l    ) );
                    }

                    //--- improper, I-J<KL
                    for(size_t l=k+1; l<connect_j.size(); ++l){
                        const Tid id_l = connect_j[l];
                        check_invalid_connect(id_j, id_l);
                        if(id_l == id_k   ||
                           id_l == id_root  ) continue;

                        improper_result.push_back( std::make_tuple(id_root,
                                                                   id_j,
                                                                   id_k,
                                                                   id_l    ) );  // root-j_k-l. j_k is axis.
                        improper_result.push_back( std::make_tuple(id_root,
                                                                   id_j,
                                                                   id_l,
                                                                   id_k    ) );  // root-j_l-k. j_l is axis.
                        improper_result.push_back( std::make_tuple(id_k,
                                                                   id_root,
                                                                   id_j,
                                                                   id_l    ) );  // k-root_j-l. root-j is axis.
                    }

                    //--- dihedral, L-I-J-K
                    for(size_t l=0; l<connect_i.size(); ++l){
                        const Tid id_l = connect_i[l];
                        check_invalid_connect(id_root, id_l);
                        if(id_l == id_j ||
                           id_l == id_k   ) continue;

                        dihedral_result.push_back( std::make_tuple(id_l,
                                                                   id_root,
                                                                   id_j,
                                                                   id_k    ) );
                    }
                }

                //--- improper, I center
                for(size_t k=j+1; k<connect_i.size(); ++k){
                    for(size_t l=k+1; l<connect_i.size(); ++l){
                        const Tid id_k = connect_i[k];
                        const Tid id_l = connect_i[l];

                        check_invalid_connect(id_root, id_k);
                        check_invalid_connect(id_root, id_l);
                        if(id_l == id_j) continue;
                        if(id_l == id_k) continue;

                        improper_result.push_back( std::make_tuple(id_j,
                                                                   id_root,
                                                                   id_k,
                                                                   id_l    ) );  // j-root_k-l. root_k is axis.
                        improper_result.push_back( std::make_tuple(id_j,
                                                                   id_root,
                                                                   id_l,
                                                                   id_k    ) );  // j-root_l-k. root_l is axis.
                        improper_result.push_back( std::make_tuple(id_k,
                                                                   id_root,
                                                                   id_j,
                                                                   id_l    ) );  // k-root_j-l. root-j is axis.
                    }
                }
            }
        }
    };

}
