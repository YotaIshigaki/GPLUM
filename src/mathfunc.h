#pragma once

PS::F64 pow2(PS::S32 x)
{
    if ( x == 0 ){
        return 1.;
    } else if ( x > 0 ) {
        return 2.*pow2(x-1);
    } else {
        return 0.5*pow2(x+1);
    }
}

inline PS::F64 calcErf(PS::F64 x)
{
    const PS::F64 p  = 0.47047;
    const PS::F64 a1 = 0.3480242;
    const PS::F64 a2 = -0.0958798;
    const PS::F64 a3 = 0.7478556;
    PS::F64 t = 1./(1.+p*x);

    return 1. - ( a1*t + a2*t*t + a3*t*t*t )*exp(-x*x);
}

inline PS::F64 getGaussian(PS::F64 sigma)
{
    PS::F64 R = drand48();

    PS::F64 X;
    if ( R < 0.2 ) {
        X = R / M_2_SQRTPI;
    } else {
        PS::F64 X0 = R;
        for ( PS::S32 i=0; i<10; i++ ){
            X = X0 + ( R-calcErf(X0) )*exp(-X0*X0);
            X0 = X;
        }
    }

    return sqrt(2.*sigma*sigma) * X;
}

PS::F64 getvalue(std::string value,
                 PS::F64 MKS_UNIT,
                 PS::F64 CGS_UNIT)
{
    PS::F64 result = 1.;
    std::vector<std::string> numer, denom;
    numer.clear();
    denom.clear();

    if ( value.size() > 4 ){
        std::string last3 = value.substr(value.size()-3, 3);
        if ( last3 == "MKS" || last3 == "mks" ) {
            result /= MKS_UNIT;
            value.erase(value.end()-3, value.end());
        }
        if ( last3 == "CGS" || last3 == "cgs" ) {
            result /= CGS_UNIT;
            value.erase(value.end()-3, value.end());
        }
    }

    PS::S32 begin = 0;
    PS::F64 x = 0;
    bool isNumer = true;
    PS::S32 size = value.size();
    for ( PS::S32 i=0; i<size+1; i++ ){
        if ( i == size || value[i] == '*' || value[i] == '/' ){
            std::string x_str = value.substr(begin, i-begin);
            
            if ( x_str.size() > 2 && x_str.substr(0, 2) == "2^" ) {
                x = pow2( std::atoi( x_str.substr(2, x_str.size()-2).c_str() ) );
            } else {
                x = std::atof( x_str.c_str() );
            }
            begin = i+1;

            if ( isNumer ) {
                result *= x;
            } else {
                result /= x;
            }
            
            if ( value[i] == '*' ) isNumer = true;
            if ( value[i] == '/' ) isNumer = false;
        }
    }
    
    return result;
}
