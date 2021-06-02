#pragma once
template <class ThisPtcl> void OutputFileWithTimeInterval(PS::ParticleSystem<ThisPtcl>& sph_system, const system_t& sysinfo, const PS::F64 end_time){
	static PS::F64 time = sysinfo.time;
	static PS::S64 step = sysinfo.step;

	if(sysinfo.time >= time){
		FileHeader header;
		header.time = sysinfo.time;
		header.Nbody = sph_system.getNumberOfParticleLocal();
		char filename[256];
		sprintf(filename, "result/%05d", step);
		sph_system.writeParticleAscii(filename, "%s_%05d_%05d.dat", header);
		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
			std::cout << "output " << filename << "." << std::endl;
			std::cout << "//================================" << std::endl;
		}
		time += end_time / PARAM::NUMBER_OF_SNAPSHOTS;
		++ step;
	}
}

template <class ThisPtcl> void InputFileWithTimeInterval(PS::ParticleSystem<ThisPtcl>& sph_system, system_t& sysinfo){
	FileHeader header;
	char filename[256];
	sprintf(filename, "result/%05d", sysinfo.step);
	sph_system.readParticleAscii(filename, "%s_%05d_%05d.dat", header);
	sysinfo.time = header.time;
	std::cout << header.time << std::endl;
}

template <class ThisPtcl> void DebugOutputFile(PS::ParticleSystem<ThisPtcl>& sph_system, const system_t& sysinfo){
	FileHeader header;
	header.time = sysinfo.time;
	header.Nbody = sph_system.getNumberOfParticleLocal();
	char filename[256];
	//sprintf(filename, "result/debug%05d", sysinfo.step);
	sprintf(filename, "result/debug%05d", 0);
	sph_system.writeParticleAscii(filename, "%s_%05d_%05d.dat", header);
	if(PS::Comm::getRank() == 0){
		std::cout << "DEBUG MODE... " << filename << "." << std::endl;
	}
}

