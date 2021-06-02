#pragma once


inline PS::F64 W_3d (const PS::F64 s){
  // const PS::F64 s = sqrt( dr * dr )/ smth_len;
  // PS::F64 value = 10.0 /( 7.0 * pi * smth_len * smth_len );//for 2D
  PS::F64 ret = 0.0;
  if (s >= zero && s <= one){
    ret = val_W_3d * (3.0 * s * s * (s - 2.0) + 4.0);
  }
  else if (s > one && s <= two){
    ret = val_W_3d *( 2.0 - s )*( 2.0 - s )*( 2.0 - s );// =a(2.0-s)^3
  }
  return ret;
}

inline PS::F64vec gradW_3d (const PS::F64vec dr, const PS::F64 s){
  // const PS::F64 dis = sqrt(dr * dr);
  // const PS::F64 s = dis / smth_len;
  // PS::F64 value=45.0/(14.0*pi*smth_len*smth_len*smth_len*smth_len);
  PS::F64vec ret(0.0);
  if (s >= zero && s <= one){
    ret = dr * (val_gradW_3d * (3.0 * s - 4.0) );//=a(s-4/3)r
  }
  else if (s > one && s <= two){
    ret = dr * (val_gradW_3d * ( -s + 4.0 - 4.0 / s));//
  }
  return ret;
}
