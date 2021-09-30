/*
 
 transfer_eisenstein.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 */

#include "transfer_function.hh"

// forward declaration of WDM class
class transfer_eisenstein_wdm_plugin;


struct eisenstein_transfer
{
  //Cosmology m_Cosmology;
  double  m_h0;
  double	omhh,		/* Omega_matter*h^2 */
  obhh,		/* Omega_baryon*h^2 */
  theta_cmb,	/* Tcmb in units of 2.7 K */
  z_equality,	/* Redshift of matter-radiation equality, really 1+z */
  k_equality,	/* Scale of equality, in Mpc^-1 */
  z_drag,		/* Redshift of drag epoch */
  R_drag,		/* Photon-baryon ratio at drag epoch */
  R_equality,	/* Photon-baryon ratio at equality epoch */
  sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
  k_silk,		/* Silk damping scale, in Mpc^-1 */
  alpha_c,	/* CDM suppression */
  beta_c,		/* CDM log shift */
  alpha_b,	/* Baryon suppression */
  beta_b,		/* Baryon envelope shift */
  beta_node,	/* Sound horizon shift */
  k_peak,		/* Fit to wavenumber of first peak, in Mpc^-1 */
  sound_horizon_fit,	/* Fit to sound horizon, in Mpc */
  alpha_gamma;	/* Gamma suppression in approximate TF */
  
  //! private member function: sets internal quantities for Eisenstein & Hu fitting
  void TFset_parameters(double omega0hh, double f_baryon, double Tcmb)
  /* Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
  /* Input: omega0hh -- The density of CDM and baryons, in units of critical dens,
   multiplied by the square of the Hubble constant, in units
   of 100 km/s/Mpc */
  /* 	  f_baryon -- The fraction of baryons to CDM */
  /*        Tcmb -- The temperature of the CMB in Kelvin.  Tcmb<=0 forces use
   of the COBE value of  2.728 K. */
  /* Output: Nothing, but set many global variables used in TFfit_onek().
   You can access them yourself, if you want. */
  /* Note: Units are always Mpc, never h^-1 Mpc. */
  {
    double z_drag_b1, z_drag_b2;
    double alpha_c_a1, alpha_c_a2, beta_c_b1, beta_c_b2, alpha_b_G, y;
    
    if (f_baryon<=0.0 || omega0hh<=0.0) {
      fprintf(stderr, "TFset_parameters(): Illegal input.\n");
      exit(1);
    }
    omhh = omega0hh;
    obhh = omhh*f_baryon;
    if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
    theta_cmb = Tcmb/2.7;
    
    z_equality = 2.50e4*omhh/POW4(theta_cmb);  /* Really 1+z */
    k_equality = 0.0746*omhh/SQR(theta_cmb);
    
    z_drag_b1 = 0.313*pow((double)omhh,-0.419)*(1+0.607*pow((double)omhh,0.674));
    z_drag_b2 = 0.238*pow((double)omhh,0.223);
    z_drag = 1291*pow(omhh,0.251)/(1+0.659*pow((double)omhh,0.828))*
    (1+z_drag_b1*pow((double)obhh,(double)z_drag_b2));
    
    R_drag = 31.5*obhh/POW4(theta_cmb)*(1000/(1+z_drag));
    R_equality = 31.5*obhh/POW4(theta_cmb)*(1000/z_equality);
    
    sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
    log((sqrt(1+R_drag)+sqrt(R_drag+R_equality))/(1+sqrt(R_equality)));
    
    k_silk = 1.6*pow((double)obhh,0.52)*pow((double)omhh,0.73)*(1+pow((double)10.4*omhh,-0.95));
    
    alpha_c_a1 = pow((double)46.9*omhh,0.670)*(1+pow(32.1*omhh,-0.532));
    alpha_c_a2 = pow((double)12.0*omhh,0.424)*(1+pow(45.0*omhh,-0.582));
    alpha_c = pow(alpha_c_a1,-f_baryon)*
    pow(alpha_c_a2,-CUBE(f_baryon));
    
    beta_c_b1 = 0.944/(1+pow(458*omhh,-0.708));
    beta_c_b2 = pow(0.395*omhh, -0.0266);
    beta_c = 1.0/(1+beta_c_b1*(pow(1-f_baryon, beta_c_b2)-1));
    
    y = z_equality/(1+z_drag);
    alpha_b_G = y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
    alpha_b = 2.07*k_equality*sound_horizon*pow(1+R_drag,-0.75)*alpha_b_G;
    
    beta_node = 8.41*pow(omhh, 0.435);
    beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*sqrt(pow(17.2*omhh,2.0)+1);
    
    k_peak = 2.5*3.14159*(1+0.217*omhh)/sound_horizon;
    sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1+10.0*pow(obhh,0.75));
    
    alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*
    SQR(f_baryon);
    
    return;
  }
  
  //! private member function: computes transfer function for mode k (k in Mpc)
  inline double TFfit_onek(double k, double *tf_baryon, double *tf_cdm)
  /* Input: k -- Wavenumber at which to calculate transfer function, in Mpc^-1.
   *tf_baryon, *tf_cdm -- Input value not used; replaced on output if
   the input was not NULL. */
  /* Output: Returns the value of the full transfer function fitting formula.
   This is the form given in Section 3 of Eisenstein & Hu (1997).
   *tf_baryon -- The baryonic contribution to the full fit.
   *tf_cdm -- The CDM contribution to the full fit. */
  /* Notes: Units are Mpc, not h^-1 Mpc. */
  {
    double T_c_ln_beta, T_c_ln_nobeta, T_c_C_alpha, T_c_C_noalpha;
    double q, xx, xx_tilde;//, q_eff;
    double T_c_f, T_c, s_tilde, T_b_T0, T_b, f_baryon, T_full;
    //double T_0_L0, T_0_C0, T_0, gamma_eff;
    //double T_nowiggles_L0, T_nowiggles_C0, T_nowiggles;
    
    k = fabs(k);	/* Just define negative k as positive */
    if (k==0.0) {
      if (tf_baryon!=NULL) *tf_baryon = 1.0;
      if (tf_cdm!=NULL) *tf_cdm = 1.0;
      return 1.0;
    }
    
    q = k/13.41/k_equality;
    xx = k*sound_horizon;
    
    T_c_ln_beta = log(2.718282+1.8*beta_c*q);
    T_c_ln_nobeta = log(2.718282+1.8*q);
    T_c_C_alpha = 14.2/alpha_c + 386.0/(1+69.9*pow(q,1.08));
    T_c_C_noalpha = 14.2 + 386.0/(1+69.9*pow(q,1.08));
    
    T_c_f = 1.0/(1.0+POW4(xx/5.4));
    T_c = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*SQR(q)) +
    (1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*SQR(q));
    
    s_tilde = sound_horizon*pow(1.+CUBE(beta_node/xx),-1./3.);
    xx_tilde = k*s_tilde;
    
    T_b_T0 = T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*SQR(q));
    T_b = sin(xx_tilde)/(xx_tilde)*(T_b_T0/(1.+SQR(xx/5.2))+
                                    alpha_b/(1.+CUBE(beta_b/xx))*exp(-pow(k/k_silk,1.4)));
    
    f_baryon = obhh/omhh;
    T_full = f_baryon*T_b + (1-f_baryon)*T_c;
    
    /* Now to store these transfer functions */
    if (tf_baryon!=NULL) *tf_baryon = T_b;
    if (tf_cdm!=NULL) *tf_cdm = T_c;
    return T_full;
  }
  
  double fb_, fc_;
  
  eisenstein_transfer()
  {  }
  
  void set_parameters( const cosmology& cosmo, double Tcmb )
  {
    m_h0 = cosmo.H0*0.01;
    TFset_parameters( (cosmo.Omega_m)*cosmo.H0*cosmo.H0*(0.01*0.01),
                     cosmo.Omega_b/cosmo.Omega_m, Tcmb);
    
    fb_ = cosmo.Omega_b/(cosmo.Omega_m);
    fc_ = (cosmo.Omega_m-cosmo.Omega_b)/(cosmo.Omega_m) ;
  }
  
  inline double at_k( double k )
  {
    double tfb, tfcdm;
    TFfit_onek( k*m_h0, &tfb, &tfcdm );
    return fb_*tfb+fc_*tfcdm;
  }
};


//! Implementation of abstract base class TransferFunction for the Eisenstein & Hu transfer function 
/*!
 This class implements the analytical fit to the matter transfer
 function by Eisenstein & Hu (1999). In fact it is their code.
 */
class transfer_eisenstein_plugin : public transfer_function_plugin
{
protected:
	using transfer_function_plugin::cosmo_;
   eisenstein_transfer etf_;
	
public:
	//! Constructor for Eisenstein & Hu fitting for transfer function
	/*!
	 \param aCosm structure of type Cosmology carrying the cosmological parameters
	 \param Tcmb mean temperature of the CMB fluctuations (defaults to
	 Tcmb = 2.726 if not specified)
	 */
	transfer_eisenstein_plugin( config_file &cf )//Cosmology aCosm, double Tcmb = 2.726 )
    :  transfer_function_plugin(cf)
	{
		double Tcmb = pcf_->getValueSafe("cosmology","Tcmb",2.726);
		
        etf_.set_parameters( cosmo_, Tcmb );
		
		tf_distinct_ = false;
		tf_withvel_  = false;
	}
	
	//! Computes the transfer function for k in Mpc/h by calling TFfit_onek
	inline double compute( double k, tf_type type ){
        return etf_.at_k( k );
	}
	
	inline double get_kmin( void ){
		return 1e-4;
	}
	
	inline double get_kmax( void ){
		return 1.e4;
	}
	
};


#include <map>
class transfer_eisenstein_wdm_plugin : public transfer_function_plugin
{
protected:
	real_t m_WDMalpha, m_h0;
	double omegam_, wdmm_, wdmgx_, wdmnu_, H0_, omegab_;
    std::string type_;
    std::map< std::string, int > typemap_;
  
    eisenstein_transfer etf_;
  
    enum wdmtyp { wdm_bode, wdm_viel, wdm_bode_wrong=99};
    
public:
	transfer_eisenstein_wdm_plugin( config_file &cf )
	: transfer_function_plugin(cf), m_h0( cosmo_.H0*0.01 )
	{
        double Tcmb = pcf_->getValueSafe("cosmology","Tcmb",2.726);
        etf_.set_parameters( cosmo_, Tcmb );
      
        typemap_.insert( std::pair<std::string,int>( "BODE", wdm_bode ) );
        typemap_.insert( std::pair<std::string,int>( "VIEL", wdm_viel ) ); // add the other types
        typemap_.insert( std::pair<std::string,int>( "BODE_WRONG", wdm_bode_wrong ) ); // add the other types
        
		omegam_ = cf.getValue<double>("cosmology","Omega_m");
		omegab_ = cf.getValue<double>("cosmology","Omega_b");
		wdmm_   = cf.getValue<double>("cosmology","WDMmass");
      
		
        H0_     = cf.getValue<double>("cosmology","H0");
        type_   = cf.getValueSafe<std::string>("cosmology","WDMtftype","BODE");
        
        //type_ = std::string( toupper( type_.c_str() ) );
        
        if( typemap_.find( type_ ) == typemap_.end() )
            throw std::runtime_error("unknown transfer function fit for WDM");
        
        m_WDMalpha = 1.0;
        
        switch( typemap_[type_] )
        {
            //... parameterisation from Bode et al. (2001), ApJ, 556, 93
            case wdm_bode:
                    wdmnu_  = cf.getValueSafe<double>("cosmology","WDMnu",1.0);
                    wdmgx_  = cf.getValueSafe<double>("cosmology","WDMg_x",1.5);
                    m_WDMalpha = 0.05 * pow( omegam_/0.4,0.15)
                        *pow(H0_*0.01/0.65,1.3)*pow(wdmm_,-1.15)
                        *pow(1.5/wdmgx_,0.29);
                
                break;
            
            //... parameterisation from Viel et al. (2005), Phys Rev D, 71
            case wdm_viel:
                    wdmnu_  = cf.getValueSafe<double>("cosmology","WDMnu",1.12);
                    m_WDMalpha = 0.049 * pow( omegam_/0.25,0.11)
                        *pow(H0_*0.01/0.7,1.22)*pow(wdmm_,-1.11);
                break;
            
            
            //.... below is for historical reasons due to the buggy parameterisation
            //.... in early versions of MUSIC, but apart from H instead of h, Bode et al.
            case wdm_bode_wrong:
                    wdmnu_  = cf.getValueSafe<double>("cosmology","WDMnu",1.0);
                    wdmgx_  = cf.getValueSafe<double>("cosmology","WDMg_x",1.5);
                    m_WDMalpha = 0.05 * pow( omegam_/0.4,0.15)
                        *pow(H0_/0.65,1.3)*pow(wdmm_,-1.15)
                        *pow(1.5/wdmgx_,0.29);
                break;
                
            default:
                    wdmnu_  = cf.getValueSafe<double>("cosmology","WDMnu",1.0);
                    wdmgx_  = cf.getValueSafe<double>("cosmology","WDMg_x",1.5);
                    m_WDMalpha = 0.05 * pow( omegam_/0.4,0.15)
                        *pow(H0_*0.01/0.65,1.3)*pow(wdmm_,-1.15)
                        *pow(1.5/wdmgx_,0.29);
                break;
        }
        std::cerr << "WDM alpha = " << m_WDMalpha << std::endl;
	}
	
	inline double compute( double k, tf_type type )
	{
		return etf_.at_k( k )*pow(1.0+pow(m_WDMalpha*k,2.0*wdmnu_),-5.0/wdmnu_);
	}
  
    inline double get_kmin( void ){
      return 1e-4;
    }
    
    inline double get_kmax( void ){
      return 1.e4;
    }
	
};

// CDM Bino type WIMP small-scale damped spectrum from Green, Hofmann & Schwarz (2004)
class transfer_eisenstein_cdmbino_plugin : public transfer_function_plugin
{
protected:
	real_t m_h0;
	double omegam_, H0_, omegab_, mcdm_, Tkd_, kfs_, kd_;
    eisenstein_transfer etf_;
    
public:
	transfer_eisenstein_cdmbino_plugin( config_file &cf )
	: transfer_function_plugin(cf), m_h0( cosmo_.H0*0.01 )
	{
        double Tcmb = pcf_->getValueSafe("cosmology","Tcmb",2.726);
        etf_.set_parameters( cosmo_, Tcmb );
        
		omegam_ = cf.getValue<double>("cosmology","Omega_m");
		omegab_ = cf.getValue<double>("cosmology","Omega_b");
		H0_     = cf.getValue<double>("cosmology","H0");
        
        mcdm_   = cf.getValueSafe<double>("cosmology","CDM_mass", 100.0); // bino particle mass in GeV
        Tkd_    = cf.getValueSafe<double>("cosmology","CDM_Tkd", 33.0); // temperature at which CDM particle kinetically decouples (in MeV)
        
        kfs_    = 1.7e6 / m_h0 * sqrt( mcdm_ / 100. * Tkd_ / 30. ) / (1.0 + log( Tkd_ / 30. ) / 19.2 );
        kd_     = 3.8e7 / m_h0 * sqrt( mcdm_ / 100. * Tkd_ / 30. );
        
        LOGINFO(" bino CDM: k_fs = %g, k_d = %g", kfs_, kd_ );
		
	}
	
	inline double compute( double k, tf_type type )
	{
        double kkfs = k/kfs_;
        double kkfs2 = kkfs*kkfs;
        double kkd2 = (k/kd_)*(k/kd_);
        
        // in principle the Green et al. (2004) works only up to k/k_fs < 1
        // the fit crosses zero at (k/k_fs)**2 = 3/2, we just zero it there...
        if( kkfs2 < 1.5 )
            return etf_.at_k( k ) * (1.0-2.0/3.0*kkfs2) * exp( -kkfs2 - kkd2 );
        else
            return 0.0;
	}
    
    inline double get_kmin( void ){
        return 1e-4;
    }
    
    inline double get_kmax( void ){
        return 1.e8;
    }
	
};



namespace{
	transfer_function_plugin_creator_concrete< transfer_eisenstein_plugin > creator("eisenstein");
	transfer_function_plugin_creator_concrete< transfer_eisenstein_wdm_plugin > creator2("eisenstein_wdm");
    transfer_function_plugin_creator_concrete< transfer_eisenstein_cdmbino_plugin > creator3("eisenstein_cdmbino");
}

