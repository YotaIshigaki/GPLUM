inline PS::F64 cutoff_poly_3rd(const PS::F64 rij,
                               const PS::F64 rout,
                               const PS::F64 rin){
    PS::F64 inv_dr = 1.0 / (rout-rin);
    PS::F64 x = (rij - rin)*inv_dr;
    PS::F64 K = 0.0;
    if(x <= 0.0)
        K = 0.0;
    else if(1.0 <= x)
        K = 1.0;
    else{
        PS::F64 x2 = x*x;
        PS::F64 x4 = x2*x2;
        PS::F64 x5 = x4*x;
        PS::F64 x6 = x4*x2;
        PS::F64 x7 = x6*x;
        K = -20.0*x7 + 70.0*x6 - 84.0*x5 + 35.0*x4;
    }
    return K;
}

inline PS::F64 cutoff_poly_3rd(const PS::F64 rij,
                               const PS::F64 rout,
                               const PS::F64 rin,
                               const PS::F64 inv_dr){
    PS::F64 x = (rij - rin)*inv_dr;
    PS::F64 K = 0.0;
    if(x <= 0.0)
        K = 0.0;
    else if(1.0 <= x)
        K = 1.0;
    else{
        PS::F64 x2 = x*x;
        PS::F64 x4 = x2*x2;
        PS::F64 x5 = x4*x;
        PS::F64 x6 = x4*x2;
        PS::F64 x7 = x6*x;
        K = -20.0*x7 + 70.0*x6 - 84.0*x5 + 35.0*x4;
    }
    return K;
}

inline void calc_acc_split(const PS::F64vec & posi,
                           PS::F64vec & acci_long,
                           PS::F64vec & acci_short,
                           PS::F64 & poti_tot,
                           const PS::F64vec & posj,
                           const PS::F64 massj,
                           const PS::F64 eps2,
                           const PS::F64 rcut_out,
                           const PS::F64 rcut_in,
                           const PS::F64 rsearch2){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2 = rij * rij;
    const PS::F64 r2_eps = r2 + eps2;
    //const PS::F64 r2_eps = rij*rij + eps2;
    const PS::F64 inv_dr = 1.0 / (rcut_out-rcut_in);
    //if(r2_eps <= rsearch2){
    if(r2 <= rsearch2){ // to be consisten with CalcForceEPEP()
        PS::F64 r_eps = sqrt(r2_eps);
        PS::F64 R = 1.0/r_eps;
        PS::F64 R2 = R*R;
        PS::F64 R3 = R2*R;
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
#else
        PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in, inv_dr);
#endif
        PS::F64vec F0 = -massj * R3 * rij;
        acci_short += F0 * (1.0-K);
        acci_long += F0 * K;
        poti_tot -= massj * R;
    }
}


inline PS::F64 cutoff_poly_3rd_dot(const PS::F64 &rij,
                                   const PS::F64 &rijvij,
                                   const PS::F64 &_rout,
                                   const PS::F64 &_rin){
    PS::F64 rout = _rout;
    PS::F64 rin = _rin;
    PS::F64 inv_dr = 1.0/(rout-rin);
    PS::F64 x = (rij - rin)*inv_dr;
    PS::F64 xdot = rijvij/rij*inv_dr;
    PS::F64 Kdot = 0.0;
    if(x <= 0.0)
        Kdot = 0.0;
    else if(1.0 <= x)
        Kdot = 0.0;
    else{
        PS::F64 x2 = x*x;
        PS::F64 x3 = x2*x;
        PS::F64 x4 = x2*x2;
        PS::F64 x5 = x4*x;
        PS::F64 x6 = x4*x2;
        Kdot = (-140.0*x6 + 420.0*x5 - 420.0*x4 + 140.0*x3) * xdot;
    }
    return Kdot;
}

#if 1
inline void CalcAcc0AndAcc1Cutoff(const PS::F64vec posi,
                                  const PS::F64vec veli,
                                  PS::F64vec & acci,
                                  PS::F64vec & jrki,
                                  const PS::F64vec posj, 
                                  const PS::F64vec velj, 
                                  const PS::F64 massj, 
                                  const PS::F64 eps2, 
                                  const PS::F64 rcut_out,
                                  const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2_eps = rij*rij + eps2;
    if(r2_eps <= rcut_out*rcut_out){
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
	//const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R*R;
        const PS::F64 R3 = R2*R;
        const PS::F64 A = (rijvij)*R2;
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -massj*R3*rij*(1.0-K);
        const PS::F64vec F1 = -massj*R3*vij*(1.0-K) - 3.0*A*F0 + massj*R3*rij*Kdot;
        acci += F0;
        jrki += F1;
    }
}
#endif

inline void CalcAcc0Acc1AndR2Cutoff(const PS::F64vec posi,
				    const PS::F64vec veli,
				    PS::F64vec & acci,
				    PS::F64vec & jrki,
				    PS::F64 & r2,
				    const PS::F64vec posj, 
				    const PS::F64vec velj,
				    const PS::F64 massj,
				    const PS::F64 eps2,
				    const PS::F64 rcut_out,
				    const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    r2 = rij*rij;
    const PS::F64 r2_eps = r2 + eps2;
    if(r2_eps <= rcut_out*rcut_out){
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
	//const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R*R;
        const PS::F64 R3 = R2*R;
        const PS::F64 A = (rijvij)*R2;
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -massj*R3*rij*(1.0-K);
        const PS::F64vec F1 = -massj*R3*vij*(1.0-K) - 3.0*A*F0 + massj*R3*rij*Kdot;
        acci += F0;
        jrki += F1;
    }
}

inline void CalcAcc0Acc1AndR2CutoffPair(const PS::F64 massi,
					const PS::F64vec posi,
					const PS::F64vec veli,
					PS::F64vec & acci,
					PS::F64vec & jrki,
					const PS::F64 massj,
					const PS::F64vec posj,
					const PS::F64vec velj,
					PS::F64vec & accj,
					PS::F64vec & jrkj,
					PS::F64 & r2,
					const PS::F64 eps2,
					const PS::F64 rcut_out,
					const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    r2 = rij*rij;
    const PS::F64 r2_eps = r2 + eps2;
    if(r2_eps <= rcut_out*rcut_out){
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
	//const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R*R;
        const PS::F64 R3 = R2*R;
        const PS::F64 A = (rijvij)*R2;
	//#ifdef FORDEBUG
#if 0
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -R3*rij*(1.0-K);
        const PS::F64vec F1 = -R3*vij*(1.0-K) - 3.0*A*F0 + R3*rij*Kdot;
        acci += massj*F0;
        jrki += massj*F1;
	accj -= massi*F0;
        jrkj -= massi*F1;
    }
}

inline void CalcAcc0Acc1AndPotCutoff(const PS::F64vec posi,
                                  const PS::F64vec veli,
                                  PS::F64vec & acci,
                                  PS::F64vec & jrki,
                                  PS::F64 & poti,
                                  const PS::F64vec posj, 
                                  const PS::F64vec velj, 
                                  const PS::F64 massj,
                                  const PS::F64 eps2,
                                  const PS::F64 rcut_out,
                                  const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2_eps = rij*rij + eps2;
#ifdef FORDEBUG
    if(1){
#else
    if(r2_eps <= rcut_out*rcut_out){
#endif
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
	//const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R * R;
        const PS::F64 R3 = R2 * R;
        const PS::F64 A = rijvij * R2;
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
        poti -= massj * R;
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -massj*R3*rij*(1.0-K);
        const PS::F64vec F1 = -massj*R3*vij*(1.0-K) - 3.0*A*F0 + massj*R3*rij*Kdot;
        acci += F0;
        jrki += F1;
    }
}

inline void CalcAcc0Acc1AndR2(const PS::F64vec posi,
			      const PS::F64vec veli,
			      PS::F64vec & acci,
			      PS::F64vec & jrki,
			      const PS::F64vec posj,
			      const PS::F64vec velj,
			      const PS::F64 massj,
			      const PS::F64 eps2){
    const PS::F64vec rij = posi - posj;
    const PS::F64vec vij = veli - velj;
    //const PS::F64 r2_eps = rij*rij + eps2;
    const PS::F64 r_inv = 1.0 / sqrt(rij * rij);
    const PS::F64 r2_inv = r_inv * r_inv;
    const PS::F64 r3_inv = r2_inv * r_inv;
    const PS::F64 m_r3 = massj * r3_inv;
    const PS::F64vec F0 = -m_r3*rij;
    const PS::F64vec F1 = -m_r3*vij - 3.0*rij*vij*r2_inv*F0;
    acci += F0;
    jrki += F1;
}

#if 0
inline PS::F64 calc_G(const PS::F64 y, const PS::F64 gamma){
    const PS::F64 y2 = y * y;
    const PS::F64 y3 = y2 * y;
    const PS::F64 y4 = y2 * y2;
    const PS::F64 y5 = y3 * y2;
    const PS::F64 y6 = y3 * y3;
    const PS::F64 y7 = y3 * y3;
    const PS::F64 gamma2 = gamma * gamma;
    const PS::F64 gamma3 = gamma2 * gamma;
    const PS::F64 c = 1.0 - gamma;
    const PS::F64 c7 = 1.0 / (c*c*c*c*c*c*c);
    PS::F64 G =
	(-10.0/3.0*y7
	 + 14.0*(gamma+1.0)*y6
	 - 21.0*(gamma2+3.0*gamma+1.0)*y5
	 + 35.0/3.0*(gamma3+9.0*gamma2+9.0*gamma+1.0)*y4
	 - 70.0*(gamma3+3.0*gamma2+gamma)*y3
	 + 210.0*(gamma3+gamma2)*y2
	 - 140.0*gamma3*y*log(y)
	 + (gamma7-7.0*gamma6+21.0*gamma5-35.0*gamma4)) * c7;
    return G;
}

inline PS::F64 calc_W(const PS::F64 y, const PS::F64 gamma){
    if(1.0 < y){
	return 1.0;
    }
    else if(gamma < y){
	calc_G
    }
    else{
	PS::F64 ret = 7.0*(gamma6 - 9.0*gamma5 + 45.0*gamma4)
    }
    const PS::F64 y2 = y * y;
    const PS::F64 y3 = y2 * y;
    const PS::F64 y4 = y2 * y2;
    const PS::F64 y5 = y3 * y2;
    const PS::F64 y6 = y3 * y3;
    const PS::F64 y7 = y3 * y3;
    const PS::F64 gamma2 = gamma * gamma;
    const PS::F64 gamma3 = gamma2 * gamma;
    const PS::F64 c = 1.0 - gamma;
    const PS::F64 c7 = 1.0 / (c*c*c*c*c*c*c);
    PS::F64 G =
	(-10.0/3.0*y7
	 + 14.0*(gamma+1.0)*y6
	 - 21.0*(gamma2+3.0*gamma+1.0)*y5
	 + 35.0/3.0*(gamma3+9.0*gamma2+9.0*gamma+1.0)*y4
	 - 70.0*(gamma3+3.0*gamma2+gamma)*y3
	 + 210.0*(gamma3+gamma2)*y2
	 - 140.0*gamma3*y*log(y)
	 + (gamma7-7.0*gamma6+21.0*gamma5-35.0*gamma4)) * c7;
    return G;
}


 
#endif
