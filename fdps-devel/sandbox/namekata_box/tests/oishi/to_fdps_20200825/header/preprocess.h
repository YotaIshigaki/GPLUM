#pragma once


void generate_ptcls(PS::ParticleSystem<FP> & sph_system){
  PS::S32 i = 0;
  // std::vector<FP> ptcl;
  // PS::F64vec min_pos(0.0);
  FileHeader header;

  char ip_path[] = "initial/initial_ptcls.bin";
  sph_system.readParticleBinary( ip_path );
  std::cout <<"readParticleBinary( "<< ip_path <<" ) is completed." <<std::endl;

  const PS::S32 nprocs = PS::Comm::getNumberOfProc();

  std::cout<<PS::Comm::getRank()<<std::endl;

  const PS::S32 nptotal = sph_system.getNumberOfParticleGlobal();
  std::cout << "nptotal:" << nptotal << "nprocs:" << nprocs <<std::endl;
  assert(nptotal % PS::Comm::getNumberOfProc() == 0);
  const PS::S32 nplocal = nptotal / PS::Comm::getNumberOfProc();
  std::cout <<"# of ptcls par process is..." <<nplocal <<std::endl;
  sph_system.setNumberOfParticleLocal(nplocal);
  
  // const PS::U32 i_head = nplocal*PS::Comm::getRank();
  // const PS::U32 i_tail = nplocal*(PS::Comm::getRank()+1);
  // for(PS::U32 i = 0 ; i < ptcl.size() ; ++i){
  //   if(i_head <= i && i < i_tail){
  //     const PS::U32 ii = i - nplocal * PS::Comm::getRank();
  //     sph_system[ii] = ptcl[i];
  //   }
  // }
  std::cout<<"setup..." <<PS::Comm::getRank()<<",  "
  <<sph_system.getNumberOfParticleLocal()<<std::endl;

}

