#include "header.h"

void OutputFileWithTimeInterval(PS::ParticleSystem<RealPtcl>& sph_system, const system_t& sysinfo){
	static double time = sysinfo.time;
	static int step = sysinfo.step;

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
		//DEBUG
		time += 500.0;
		//time += 0.0;
		++ step;
	}
}

void InputFileWithTimeInterval(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo){
	FileHeader header;
	char filename[256];
	sprintf(filename, "result/%05d", sysinfo.step);
	//sprintf(filename, "%05d", sysinfo.step);
	sph_system.readParticleAscii(filename, "%s_%05d_%05d.dat", header);
	sysinfo.time = header.time;
	sysinfo.end_time = 1.0e+5;
	std::cout << header.time << std::endl;
}

