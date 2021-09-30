typedef float real;

namespace Dens{
	struct EpiDev{
		real rx;
		real ry;
		real rz;
		real mass;
		real smth;
		int id_walk;
	};

	struct EpjDev{
		real rx;
		real ry;
		real rz;
		real mass;
		real smth;
	};
	struct ForceDev{
		real dens;
	};
}

namespace Drvt{
	struct EpiDev{
		real rx;
		real ry;
		real rz;
		real vx;
		real vy;
		real vz;
		real dens;
		real smth;
		int id_walk;
	};

	struct EpjDev{
		real rx;
		real ry;
		real rz;
		real vx;
		real vy;
		real vz;
		real mass;
		real smth;
	};
	struct ForceDev{
		real div_v;
		real rot_vx;
		real rot_vy;
		real rot_vz;
	};
}

namespace Hydr{
	struct EpiDev{
		real rx;
		real ry;
		real rz;
		real vx;
		real vy;
		real vz;
		real dens;
		real pres;
		real snds;
		real smth;
		real Bal;
		int  id_walk;
		real grad_smth;
	};

	struct EpjDev{
		real rx;
		real ry;
		real rz;
		real vx;
		real vy;
		real vz;
		real dens;
		real pres;
		real snds;
		real mass;
		real smth;
		real Bal;
		real grad_smth;
	};
	struct ForceDev{
		real ax;
		real ay;
		real az;
		real eng_dot;
		real dt;
	};	
}

namespace Grav{
	struct EpiDev{
		real rx;
		real ry;
		real rz;
		real eps2;
		int id_walk;
	};

	struct EpjDev{
		real rx;
		real ry;
		real rz;
		real mass;
	};
	struct ForceDev{
		real ax;
		real ay;
		real az;
		real pot;
	};	
}



