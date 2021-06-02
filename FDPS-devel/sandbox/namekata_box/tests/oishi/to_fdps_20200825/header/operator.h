#pragma once


inline PS::F64mat product_02_3d( const PS::F64 s, const PS::F64mat m3 ){
  PS::F64mat ret(s * m3.xx, s * m3.yy, s * m3.zz, s * m3.xy, s * m3.xz, s * m3.yz);
  return ret;
}
// PS::F64 dotpro_vec3( const PS::F64vec v3 ){
//   return v3.x*v3.x + v3.y*v3.y + v3.z*v3.z ;
// }

// inline PS::F64vec dot_product_21_2d(const PS::F64mat m, const PS::F64vec v){
//   PS::F64vec ret(m.xx * v.x + m.xy * v.y, m.xy * v.x + m.yy * v.y);
//   return ret;
// }

inline PS::F64vec dot_product_21_3d(const PS::F64mat m, const PS::F64vec v){
  PS::F64vec ret(m.xx * v.x + m.xy * v.y + m.xz * v.z,
                 m.xy * v.x + m.yy * v.y + m.yz * v.z,
                 m.xz * v.x + m.yz * v.y + m.zz * v.z);
  return ret;
}

inline PS::F64 J2_3d(const PS::F64mat m){
  return 0.5 * (m.xx * m.xx + m.yy * m.yy + m.zz * m.zz) 
         + m.xy * m.xy + m.xz * m.xz + m.yz * m.yz;
}

