#if KERNEL_TYPE == 0
static inline double __attribute__((always_inline)) SPHKernel(const double r, const double InvKerneli){ 

	double u = r*InvKerneli;
#if (DIMENSION == 1)
    const static double coef1d = 2.0/3.0;
	double coef = coef1d*InvKerneli;
#elif (DIMENSION == 2)   
    const static double coef2d = 10.0/(7.0*M_PI);
    double coef = coef2d*SQ(InvKerneli);
#elif (DIMENSION == 3)
	double coef = M_1_PI*CUBE(InvKerneli);
#endif

	if(u<1.e0){
		return (coef*(1.e0 - 1.5*SQ(u) + 0.75*CUBE(u)));
	} else if (u<2.e0){
		return (coef*(0.25*CUBE(2.e0-u)));
	} else {
	    return 0.e0;
    }
}

#define SuppressTensileInstability (ON)
static inline double __attribute__((always_inline)) dSPHKernel(const double r, const double InvKerneli){

	if(!(r>0.e0))
		return (0.e0);

	double u = r*InvKerneli;
#if (DIMENSION == 1)
    const static double coef1d = 2.0/3.0;
	double coef = coef1d*CUBE(InvKerneli);
#elif (DIMENSION == 2)   
    const static double coef2d = 10.0/(7.0*M_PI);
    double coef = coef2d*SQ(SQ(InvKerneli));
#elif (DIMENSION == 3)
	double coef = M_1_PI*CUBE(InvKerneli)*SQ(InvKerneli);
#endif
    /*
    */
	if(u<2.0/3.0){
		return (-coef/u);
    } else if(u<1.e0){
		return (-coef*(3.0 - 2.25*u));
	} else if (u<2.e0){
		return (-coef*(0.75*SQ(2.e0-u))/u);
	} else {
        return 0.e0;
    }
}

#elif KERNEL_TYPE == 1 // Wendland C2
#define WENDLANDC2_1D_COEF (1.620185)
#define WENDLANDC2_2D_COEF (1.897367)
#define WENDLANDC2_3D_COEF (1.936492)

static inline double __attribute__((always_inline)) SPHKernel(const double r, const double InvKerneli){ 

    double InvH = 0.5*InvKerneli; // 1.0/h->1.0/(2.0 h)
    double u = r*InvH;
    double uplus = fmax(0.0,(1.0-u));

#if (DIMENSION == 1)   
    static const double coef1d = 5.0/4.0;
    double coef = coef1d*InvH;
    return coef*CUBE(uplus)*(1.0+3.0*u);
#elif (DIMENSION == 2)   
    static const double coef2d = 7.0*M_1_PI;
    double coef = coef2d*SQ(InvH);
    return coef*SQ(SQ(uplus))*(1.0+4.0*u);
#elif (DIMENSION == 3)   
    static const double coef3d = 10.5*M_1_PI;
    double coef = coef3d*CUBE(InvH);
    return coef*SQ(SQ(uplus))*(1.0+4.0*u);
#else
#error DIMENSION is worng.
#endif
}

static inline double __attribute__((always_inline)) dSPHKernel(const double r, const double InvKerneli){

    double InvH = 0.5*InvKerneli;
    double u = r*InvH;
    double uplus = fmax(0.0,(1.0-u));

#if (DIMENSION == 1)   
    static const double coef1d = 5.0/4.0;
    double coef = coef1d*CUBE(InvH);
    return coef*(-12.0)*SQ(uplus);
#elif (DIMENSION == 2)   
    static const double coef2d = 7.0*M_1_PI; 
    double coef = coef2d*SQ(SQ(InvH));
    return coef*(-20.0)*CUBE(uplus);
#elif (DIMENSION == 3)   
    static const double coef3d = 10.5*M_1_PI; 
    double coef = coef3d*SQ(InvH)*CUBE(InvH);
    return coef*(-20.0)*CUBE(uplus);
#else
#error DIMENSION is worng.
#endif
}

#elif KERNEL_TYPE == 2 // Wendland C4
#define WENDLANDC4_1D_COEF (1.936492)
#define WENDLANDC4_2D_COEF (2.171239)
#define WENDLANDC4_3D_COEF (2.207940)

static inline double __attribute__((always_inline)) SPHKernel(const double r, const double InvKerneli){ 

    double InvH = 0.5*InvKerneli;
    double u = r*InvH;
    double uplus = fmax(0.0,(1.0-u));

#if (DIMENSION == 1)   
    static const double coef1d = 1.5;
    double coef = coef1d*InvH;
    return coef*CUBE(uplus)*SQ(uplus)*(1.0+5.0*u+8.0*SQ(u));
#elif (DIMENSION == 2)   
    static const double coef2d = 9.0*M_1_PI;
    double coef = coef2d*SQ(InvH);
    return coef*SQ(CUBE(uplus))*(1.0+6.0*u+11.6666666666667*SQ(u)); //35/3 = 11.6666666666667
#elif (DIMENSION == 3)   
    static const double coef3d = 15.46875*M_1_PI; // 495.0/32.0 = 15.46875
    double coef = coef3d*CUBE(InvH);
    return coef*SQ(CUBE(uplus))*(1.0+6.0*u+11.6666666666667*SQ(u)); //35/3 = 11.6666666666667
#else
#error DIMENSION is worng.
#endif
}

static inline double __attribute__((always_inline)) dSPHKernel(const double r, const double InvKerneli){

    double InvH = 0.5*InvKerneli;
    double u = r*InvH;
    double uplus = fmax(0.0,(1.0-u));

#if (DIMENSION == 1)   
    static const double coef1d = 1.5;
    double coef = coef1d*CUBE(InvH);
    //return coef*(-14.0)*SQ(SQ(uplus))*(1.0+4.0*u);
    return coef*SQ(SQ(uplus))*(-14.0-56.0*u);
#elif (DIMENSION == 2)   
    static const double coef2d = 9.0*M_1_PI;
    double coef = coef2d*SQ(SQ(InvH));
    //return coef*(-6.666666667)*SQ(uplus)*CUBE(uplus)*(1.0+14.0*u);
    return coef*SQ(uplus)*CUBE(uplus)*(-18.66666666666667-93.33333333333333*u);
#elif (DIMENSION == 3)   
    static const double coef3d = 15.46875*M_1_PI; // 495.0/32.0 = 15.46875
    double coef = coef3d*CUBE(InvH)*SQ(InvH);
    //return coef*(-6.666666667)*SQ(uplus)*CUBE(uplus)*(1.0+14.0*u);
    return coef*SQ(uplus)*CUBE(uplus)*(-18.66666666666667-93.33333333333333*u);
#else
#error DIMENSION is worng.
#endif
}

#elif KERNEL_TYPE == 3 // Wendland C6

#define WENDLANDC6_1D_COEF (2.207940)
#define WENDLANDC6_2D_COEF (2.415230)
#define WENDLANDC6_3D_COEF (2.449490)

static inline double __attribute__((always_inline)) SPHKernel(const double r, const double InvKerneli){ 

    double InvH = 0.5*InvKerneli;
    double u = r*InvH;
    double uplus = fmax(0.0,(1.0-u));

#if (DIMENSION == 1)   
    static const double coef1d = 55.0/32.0;
    double coef = coef1d*InvH;
    return coef*CUBE(SQ(uplus))*(uplus)*(1.0+7.0*u+19.0*SQ(u)+21.0*CUBE(u));
#elif (DIMENSION == 2)   
    static const double coef2d = (78.0/7.0)*M_1_PI;
    double coef = coef2d*SQ(InvH);
    return coef*SQ(SQ(SQ(uplus)))*(1.0+8.0*u+25.0*SQ(u)+32.0*CUBE(u));
#elif (DIMENSION == 3)   
    static const double coef3d = (1365.0/64.0)*M_1_PI; 
    double coef = coef3d*CUBE(InvH);
    return coef*SQ(SQ(SQ(uplus)))*(1.0+8.0*u+25.0*SQ(u)+32.0*CUBE(u));
#else
#error DIMENSION is worng.
#endif
}

static inline double __attribute__((always_inline)) dSPHKernel(const double r, const double InvKerneli){

    double InvH = 0.5*InvKerneli;
    double u = r*InvH;
    double uplus = fmax(0.0,(1.0-u));

#if (DIMENSION == 1)   
    static const double coef1d = 55.0/32.0;
    double coef = coef1d*CUBE(InvH);
    return coef*SQ(CUBE(uplus))*(-18.0-108.0*u-210.0*SQ(u));
#elif (DIMENSION == 2)   
    static const double coef2d = (78.0/7.0)*M_1_PI;
    double coef = coef2d*SQ(SQ(InvH));
    return coef*SQ(SQ(uplus))*CUBE(uplus)*(-22.0-154.0*u-352.0*SQ(u));
#elif (DIMENSION == 3)   
    static const double coef3d = (1365.0/64.0)*M_1_PI;
    double coef = coef3d*CUBE(InvH)*SQ(InvH);
    return coef*SQ(SQ(uplus))*CUBE(uplus)*(-22.0-154.0*u-352.0*SQ(u));
#else
#error DIMENSION is worng.
#endif
}




#else
#error Kernel type undefined. 
#endif


