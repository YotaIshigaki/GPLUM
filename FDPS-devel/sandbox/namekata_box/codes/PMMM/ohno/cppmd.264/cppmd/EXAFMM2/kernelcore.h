#ifdef SPARC_SIMD
#include <emmintrin.h>
#endif

void getCoef(real *C, const real *dist, real &invR2, const real &invR) {
  C[0] = invR;
  invR2 = -invR2;
  real x = dist[0], y = dist[1], z = dist[2];
  real invR3 = invR * invR2;
  C[1] = x * invR3;
  C[2] = y * invR3;
  C[3] = z * invR3;

  real invR5 = 3 * invR3 * invR2;
  real t = x * invR5;
  C[4] = x * t + invR3;
  C[5] = y * t;
  C[6] = z * t;
  t = y * invR5;
  C[7] = y * t + invR3;
  C[8] = z * t;
  C[9] = z * z * invR5 + invR3;

  real invR7 = 5 * invR5 * invR2;
  t = x * x * invR7;
  C[10] = x * (t + 3 * invR5);
  C[11] = y * (t +     invR5);
  C[12] = z * (t +     invR5);
  t = y * y * invR7;
  C[13] = x * (t +     invR5);
  C[16] = y * (t + 3 * invR5);
  C[17] = z * (t +     invR5);
  t = z * z * invR7;
  C[15] = x * (t +     invR5);
  C[18] = y * (t +     invR5);
  C[19] = z * (t + 3 * invR5);
  C[14] = x * y * z * invR7;

  real invR9 = 7 * invR7 * invR2;
  t = x * x * invR9;
  C[20] = x * x * (t + 6 * invR7) + 3 * invR5;
  C[21] = x * y * (t + 3 * invR7);
  C[22] = x * z * (t + 3 * invR7);
  C[23] = y * y * (t +     invR7) + x * x * invR7 + invR5;
  C[24] = y * z * (t +     invR7);
  C[25] = z * z * (t +     invR7) + x * x * invR7 + invR5;
  t = y * y * invR9;
  C[26] = x * y * (t + 3 * invR7);
  C[27] = x * z * (t +     invR7);
  C[30] = y * y * (t + 6 * invR7) + 3 * invR5;
  C[31] = y * z * (t + 3 * invR7);
  C[32] = z * z * (t +     invR7) + y * y * invR7 + invR5;
  t = z * z * invR9;
  C[28] = x * y * (t +     invR7);
  C[29] = x * z * (t + 3 * invR7);
  C[33] = y * z * (t + 3 * invR7);
  C[34] = z * z * (t + 6 * invR7) + 3 * invR5;

  real invR11 = 9 * invR9 * invR2;
  t = x * x * invR11;
  C[35] = x * x * x * (t + 10 * invR9) + 15 * x * invR7;
  C[36] = x * x * y * (t +  6 * invR9) +  3 * y * invR7;
  C[37] = x * x * z * (t +  6 * invR9) +  3 * z * invR7;
  C[38] = x * y * y * (t +  3 * invR9) + x * x * x * invR9 + 3 * x * invR7;
  C[39] = x * y * z * (t +  3 * invR9);
  C[40] = x * z * z * (t +  3 * invR9) + x * x * x * invR9 + 3 * x * invR7;
  C[41] = y * y * y * (t +      invR9) + 3 * x * x * y * invR9 + 3 * y * invR7;
  C[42] = y * y * z * (t +      invR9) + x * x * z * invR9 + z * invR7;
  C[43] = y * z * z * (t +      invR9) + x * x * y * invR9 + y * invR7;
  C[44] = z * z * z * (t +      invR9) + 3 * x * x * z * invR9 + 3 * z * invR7;
  t = y * y * invR11;
  C[45] = x * y * y * (t +  6 * invR9) +  3 * x * invR7;
  C[46] = x * y * z * (t +  3 * invR9);
  C[47] = x * z * z * (t +      invR9) + x * y * y * invR9 + x * invR7;
  C[50] = y * y * y * (t + 10 * invR9) + 15 * y * invR7;
  C[51] = y * y * z * (t +  6 * invR9) + 3 * z * invR7;
  C[52] = y * z * z * (t +  3 * invR9) + y * y * y * invR9 + 3 * y * invR7;
  C[53] = z * z * z * (t +      invR9) + 3 * y * y * z * invR9 + 3 * z * invR7;
  t = z * z * invR11;
  C[48] = x * y * z * (t +  3 * invR9);
  C[49] = x * z * z * (t +  6 * invR9) +  3 * x * invR7;
  C[54] = y * z * z * (t +  6 * invR9) +  3 * y * invR7;
  C[55] = z * z * z * (t + 10 * invR9) + 15 * z * invR7;

  real invR13 = 11 * invR11 * invR2;
  t = x * x * invR13;
  C[56] = x * x * x * x * (t + 15 * invR11) + 45 * x * x * invR9 + 15 * invR7;
  C[57] = x * x * x * y * (t + 10 * invR11) + 15 * x * y * invR9;
  C[58] = x * x * x * z * (t + 10 * invR11) + 15 * x * z * invR9;
  C[59] = x * x * y * y * (t +  6 * invR11) + x * x * x * x * invR11 + (6 * x * x + 3 * y * y) * invR9 + 3 * invR7;
  C[60] = x * x * y * z * (t +  6 * invR11) + 3 * y * z * invR9;
  C[61] = x * x * z * z * (t +  6 * invR11) + x * x * x * x * invR11 + (6 * x * x + 3 * z * z) * invR9 + 3 * invR7;
  C[62] = x * y * y * y * (t +  3 * invR11) + 3 * x * x * x * y * invR11 + 9 * x * y * invR9;
  C[63] = x * y * y * z * (t +  3 * invR11) + x * x * x * z * invR11 + 3 * x * z * invR9;
  C[64] = x * y * z * z * (t +  3 * invR11) + x * x * x * y * invR11 + 3 * x * y * invR9;
  C[65] = x * z * z * z * (t +  3 * invR11) + 3 * x * x * x * z * invR11 + 9 * x * z * invR9;
  C[66] = y * y * y * y * (t +      invR11) + 6 * x * x * y * y * invR11 + (3 * x * x + 6 * y * y) * invR9 + 3 * invR7;
  C[67] = y * y * y * z * (t +      invR11) + 3 * x * x * y * z * invR11 + 3 * y * z * invR9;
  C[68] = y * y * z * z * (t +      invR11) + (x * x * y * y + x * x * z * z) * invR11 + (x * x + y * y + z * z) * invR9 + invR7;
  C[69] = y * z * z * z * (t +      invR11) + 3 * x * x * y * z * invR11 + 3 * y * z * invR9;
  C[70] = z * z * z * z * (t +      invR11) + 6 * x * x * z * z * invR11 + (3 * x * x + 6 * z * z) * invR9 + 3 * invR7;
  t = y * y * invR13;
  C[71] = x * y * y * y * (t + 10 * invR11) + 15 * x * y * invR9;
  C[72] = x * y * y * z * (t +  6 * invR11) + 3 * x * z * invR9;
  C[73] = x * y * z * z * (t +  3 * invR11) + x * y * y * y * invR11 + 3 * x * y * invR9;
  C[74] = x * z * z * z * (t +      invR11) + 3 * x * y * y * z * invR11 + 3 * x * z * invR9;
  C[77] = y * y * y * y * (t + 15 * invR11) + 45 * y * y * invR9 + 15 * invR7;
  C[78] = y * y * y * z * (t + 10 * invR11) + 15 * y * z * invR9;
  C[79] = y * y * z * z * (t +  6 * invR11) + y * y * y * y * invR11 + (6 * y * y + 3 * z * z) * invR9 + 3 * invR7;
  C[80] = y * z * z * z * (t +  3 * invR11) + 3 * y * y * y * z * invR11 + 9 * y * z * invR9;
  C[81] = z * z * z * z * (t +      invR11) + 6 * y * y * z * z * invR11 + (3 * y * y + 6 * z * z) * invR9 + 3 * invR7;
  t = z * z * invR13;
  C[75] = x * y * z * z * (t +  6 * invR11) + 3 * x * y * invR9;
  C[76] = x * z * z * z * (t + 10 * invR11) + 15 * x * z * invR9;
  C[82] = y * z * z * z * (t + 10 * invR11) + 15 * y * z * invR9;
  C[83] = z * z * z * z * (t + 15 * invR11) + 45 * z * z * invR9 + 15 * invR7;
}

void M2LSum(real *L, const real *C, const real*M) {
  for( int l=0; l<LTERM; l++ ) L[l] += M[0] * C[l];

  L[0] += M[1]*C[1]+M[2]*C[2]+M[3]*C[3];
  L[1] += M[1]*C[4]+M[2]*C[5]+M[3]*C[6];
  L[2] += M[1]*C[5]+M[2]*C[7]+M[3]*C[8];
  L[3] += M[1]*C[6]+M[2]*C[8]+M[3]*C[9];

  for( int l=4; l<10; l++ ) L[0] += M[l] * C[l];
  L[1] += M[4]*C[10]+M[5]*C[11]+M[6]*C[12]+M[7]*C[13]+M[8]*C[14]+M[9]*C[15];
  L[2] += M[4]*C[11]+M[5]*C[13]+M[6]*C[14]+M[7]*C[16]+M[8]*C[17]+M[9]*C[18];
  L[3] += M[4]*C[12]+M[5]*C[14]+M[6]*C[15]+M[7]*C[17]+M[8]*C[18]+M[9]*C[19];
  L[4] += M[1]*C[10]+M[2]*C[11]+M[3]*C[12];
  L[5] += M[1]*C[11]+M[2]*C[13]+M[3]*C[14];
  L[6] += M[1]*C[12]+M[2]*C[14]+M[3]*C[15];
  L[7] += M[1]*C[13]+M[2]*C[16]+M[3]*C[17];
  L[8] += M[1]*C[14]+M[2]*C[17]+M[3]*C[18];
  L[9] += M[1]*C[15]+M[2]*C[18]+M[3]*C[19];

  for( int l=10; l<20; l++ ) L[0] += M[l] * C[l];
  L[1] += M[10]*C[20]+M[11]*C[21]+M[12]*C[22]+M[13]*C[23]+M[14]*C[24]+M[15]*C[25]+M[16]*C[26]+M[17]*C[27]+M[18]*C[28]+M[19]*C[29];
  L[2] += M[10]*C[21]+M[11]*C[23]+M[12]*C[24]+M[13]*C[26]+M[14]*C[27]+M[15]*C[28]+M[16]*C[30]+M[17]*C[31]+M[18]*C[32]+M[19]*C[33];
  L[3] += M[10]*C[22]+M[11]*C[24]+M[12]*C[25]+M[13]*C[27]+M[14]*C[28]+M[15]*C[29]+M[16]*C[31]+M[17]*C[32]+M[18]*C[33]+M[19]*C[34];
  L[4] += M[4]*C[20]+M[5]*C[21]+M[6]*C[22]+M[7]*C[23]+M[8]*C[24]+M[9]*C[25];
  L[5] += M[4]*C[21]+M[5]*C[23]+M[6]*C[24]+M[7]*C[26]+M[8]*C[27]+M[9]*C[28];
  L[6] += M[4]*C[22]+M[5]*C[24]+M[6]*C[25]+M[7]*C[27]+M[8]*C[28]+M[9]*C[29];
  L[7] += M[4]*C[23]+M[5]*C[26]+M[6]*C[27]+M[7]*C[30]+M[8]*C[31]+M[9]*C[32];
  L[8] += M[4]*C[24]+M[5]*C[27]+M[6]*C[28]+M[7]*C[31]+M[8]*C[32]+M[9]*C[33];
  L[9] += M[4]*C[25]+M[5]*C[28]+M[6]*C[29]+M[7]*C[32]+M[8]*C[33]+M[9]*C[34];
  L[10] += M[1]*C[20]+M[2]*C[21]+M[3]*C[22];
  L[11] += M[1]*C[21]+M[2]*C[23]+M[3]*C[24];
  L[12] += M[1]*C[22]+M[2]*C[24]+M[3]*C[25];
  L[13] += M[1]*C[23]+M[2]*C[26]+M[3]*C[27];
  L[14] += M[1]*C[24]+M[2]*C[27]+M[3]*C[28];
  L[15] += M[1]*C[25]+M[2]*C[28]+M[3]*C[29];
  L[16] += M[1]*C[26]+M[2]*C[30]+M[3]*C[31];
  L[17] += M[1]*C[27]+M[2]*C[31]+M[3]*C[32];
  L[18] += M[1]*C[28]+M[2]*C[32]+M[3]*C[33];
  L[19] += M[1]*C[29]+M[2]*C[33]+M[3]*C[34];

 for( int l=20; l<35; l++ ) L[0] += M[l] * C[l];
  L[1] += M[20]*C[35]+M[21]*C[36]+M[22]*C[37]+M[23]*C[38]+M[24]*C[39]+M[25]*C[40]+M[26]*C[41]+M[27]*C[42]+M[28]*C[43]+M[29]*C[44]+M[30]*C[45]+M[31]*C[46]+M[32]*C[47]+M[33]*C[48]+M[34]*C[49];
  L[2] += M[20]*C[36]+M[21]*C[38]+M[22]*C[39]+M[23]*C[41]+M[24]*C[42]+M[25]*C[43]+M[26]*C[45]+M[27]*C[46]+M[28]*C[47]+M[29]*C[48]+M[30]*C[50]+M[31]*C[51]+M[32]*C[52]+M[33]*C[53]+M[34]*C[54];
  L[3] += M[20]*C[37]+M[21]*C[39]+M[22]*C[40]+M[23]*C[42]+M[24]*C[43]+M[25]*C[44]+M[26]*C[46]+M[27]*C[47]+M[28]*C[48]+M[29]*C[49]+M[30]*C[51]+M[31]*C[52]+M[32]*C[53]+M[33]*C[54]+M[34]*C[55];
  L[4] += M[10]*C[35]+M[11]*C[36]+M[12]*C[37]+M[13]*C[38]+M[14]*C[39]+M[15]*C[40]+M[16]*C[41]+M[17]*C[42]+M[18]*C[43]+M[19]*C[44];
  L[5] += M[10]*C[36]+M[11]*C[38]+M[12]*C[39]+M[13]*C[41]+M[14]*C[42]+M[15]*C[43]+M[16]*C[45]+M[17]*C[46]+M[18]*C[47]+M[19]*C[48];
  L[6] += M[10]*C[37]+M[11]*C[39]+M[12]*C[40]+M[13]*C[42]+M[14]*C[43]+M[15]*C[44]+M[16]*C[46]+M[17]*C[47]+M[18]*C[48]+M[19]*C[49];
  L[7] += M[10]*C[38]+M[11]*C[41]+M[12]*C[42]+M[13]*C[45]+M[14]*C[46]+M[15]*C[47]+M[16]*C[50]+M[17]*C[51]+M[18]*C[52]+M[19]*C[53];
  L[8] += M[10]*C[39]+M[11]*C[42]+M[12]*C[43]+M[13]*C[46]+M[14]*C[47]+M[15]*C[48]+M[16]*C[51]+M[17]*C[52]+M[18]*C[53]+M[19]*C[54];
  L[9] += M[10]*C[40]+M[11]*C[43]+M[12]*C[44]+M[13]*C[47]+M[14]*C[48]+M[15]*C[49]+M[16]*C[52]+M[17]*C[53]+M[18]*C[54]+M[19]*C[55];
  L[10] += M[4]*C[35]+M[5]*C[36]+M[6]*C[37]+M[7]*C[38]+M[8]*C[39]+M[9]*C[40];
  L[11] += M[4]*C[36]+M[5]*C[38]+M[6]*C[39]+M[7]*C[41]+M[8]*C[42]+M[9]*C[43];
  L[12] += M[4]*C[37]+M[5]*C[39]+M[6]*C[40]+M[7]*C[42]+M[8]*C[43]+M[9]*C[44];
  L[13] += M[4]*C[38]+M[5]*C[41]+M[6]*C[42]+M[7]*C[45]+M[8]*C[46]+M[9]*C[47];
  L[14] += M[4]*C[39]+M[5]*C[42]+M[6]*C[43]+M[7]*C[46]+M[8]*C[47]+M[9]*C[48];
  L[15] += M[4]*C[40]+M[5]*C[43]+M[6]*C[44]+M[7]*C[47]+M[8]*C[48]+M[9]*C[49];
  L[16] += M[4]*C[41]+M[5]*C[45]+M[6]*C[46]+M[7]*C[50]+M[8]*C[51]+M[9]*C[52];
  L[17] += M[4]*C[42]+M[5]*C[46]+M[6]*C[47]+M[7]*C[51]+M[8]*C[52]+M[9]*C[53];
  L[18] += M[4]*C[43]+M[5]*C[47]+M[6]*C[48]+M[7]*C[52]+M[8]*C[53]+M[9]*C[54];
  L[19] += M[4]*C[44]+M[5]*C[48]+M[6]*C[49]+M[7]*C[53]+M[8]*C[54]+M[9]*C[55];
  L[20] += M[1]*C[35]+M[2]*C[36]+M[3]*C[37];
  L[21] += M[1]*C[36]+M[2]*C[38]+M[3]*C[39];
  L[22] += M[1]*C[37]+M[2]*C[39]+M[3]*C[40];
  L[23] += M[1]*C[38]+M[2]*C[41]+M[3]*C[42];
  L[24] += M[1]*C[39]+M[2]*C[42]+M[3]*C[43];
  L[25] += M[1]*C[40]+M[2]*C[43]+M[3]*C[44];
  L[26] += M[1]*C[41]+M[2]*C[45]+M[3]*C[46];
  L[27] += M[1]*C[42]+M[2]*C[46]+M[3]*C[47];
  L[28] += M[1]*C[43]+M[2]*C[47]+M[3]*C[48];
  L[29] += M[1]*C[44]+M[2]*C[48]+M[3]*C[49];
  L[30] += M[1]*C[45]+M[2]*C[50]+M[3]*C[51];
  L[31] += M[1]*C[46]+M[2]*C[51]+M[3]*C[52];
  L[32] += M[1]*C[47]+M[2]*C[52]+M[3]*C[53];
  L[33] += M[1]*C[48]+M[2]*C[53]+M[3]*C[54];
  L[34] += M[1]*C[49]+M[2]*C[54]+M[3]*C[55];

  for( int l=35; l<56; l++ ) L[0] += M[l] * C[l];
  L[1] += M[35]*C[56]+M[36]*C[57]+M[37]*C[58]+M[38]*C[59]+M[39]*C[60]+M[40]*C[61]+M[41]*C[62]+M[42]*C[63]+M[43]*C[64]+M[44]*C[65]+M[45]*C[66]+M[46]*C[67]+M[47]*C[68]+M[48]*C[69]+M[49]*C[70]+M[50]*C[71]+M[51]*C[72]+M[52]*C[73]+M[53]*C[74]+M[54]*C[75]+M[55]*C[76];
  L[2] += M[35]*C[57]+M[36]*C[59]+M[37]*C[60]+M[38]*C[62]+M[39]*C[63]+M[40]*C[64]+M[41]*C[66]+M[42]*C[67]+M[43]*C[68]+M[44]*C[69]+M[45]*C[71]+M[46]*C[72]+M[47]*C[73]+M[48]*C[74]+M[49]*C[75]+M[50]*C[77]+M[51]*C[78]+M[52]*C[79]+M[53]*C[80]+M[54]*C[81]+M[55]*C[82];
  L[3] += M[35]*C[58]+M[36]*C[60]+M[37]*C[61]+M[38]*C[63]+M[39]*C[64]+M[40]*C[65]+M[41]*C[67]+M[42]*C[68]+M[43]*C[69]+M[44]*C[70]+M[45]*C[72]+M[46]*C[73]+M[47]*C[74]+M[48]*C[75]+M[49]*C[76]+M[50]*C[78]+M[51]*C[79]+M[52]*C[80]+M[53]*C[81]+M[54]*C[82]+M[55]*C[83];
  L[4] += M[20]*C[56]+M[21]*C[57]+M[22]*C[58]+M[23]*C[59]+M[24]*C[60]+M[25]*C[61]+M[26]*C[62]+M[27]*C[63]+M[28]*C[64]+M[29]*C[65]+M[30]*C[66]+M[31]*C[67]+M[32]*C[68]+M[33]*C[69]+M[34]*C[70];
  L[5] += M[20]*C[57]+M[21]*C[59]+M[22]*C[60]+M[23]*C[62]+M[24]*C[63]+M[25]*C[64]+M[26]*C[66]+M[27]*C[67]+M[28]*C[68]+M[29]*C[69]+M[30]*C[71]+M[31]*C[72]+M[32]*C[73]+M[33]*C[74]+M[34]*C[75];
  L[6] += M[20]*C[58]+M[21]*C[60]+M[22]*C[61]+M[23]*C[63]+M[24]*C[64]+M[25]*C[65]+M[26]*C[67]+M[27]*C[68]+M[28]*C[69]+M[29]*C[70]+M[30]*C[72]+M[31]*C[73]+M[32]*C[74]+M[33]*C[75]+M[34]*C[76];
  L[7] += M[20]*C[59]+M[21]*C[62]+M[22]*C[63]+M[23]*C[66]+M[24]*C[67]+M[25]*C[68]+M[26]*C[71]+M[27]*C[72]+M[28]*C[73]+M[29]*C[74]+M[30]*C[77]+M[31]*C[78]+M[32]*C[79]+M[33]*C[80]+M[34]*C[81];
  L[8] += M[20]*C[60]+M[21]*C[63]+M[22]*C[64]+M[23]*C[67]+M[24]*C[68]+M[25]*C[69]+M[26]*C[72]+M[27]*C[73]+M[28]*C[74]+M[29]*C[75]+M[30]*C[78]+M[31]*C[79]+M[32]*C[80]+M[33]*C[81]+M[34]*C[82];
  L[9] += M[20]*C[61]+M[21]*C[64]+M[22]*C[65]+M[23]*C[68]+M[24]*C[69]+M[25]*C[70]+M[26]*C[73]+M[27]*C[74]+M[28]*C[75]+M[29]*C[76]+M[30]*C[79]+M[31]*C[80]+M[32]*C[81]+M[33]*C[82]+M[34]*C[83];
  L[10] += M[10]*C[56]+M[11]*C[57]+M[12]*C[58]+M[13]*C[59]+M[14]*C[60]+M[15]*C[61]+M[16]*C[62]+M[17]*C[63]+M[18]*C[64]+M[19]*C[65];
  L[11] += M[10]*C[57]+M[11]*C[59]+M[12]*C[60]+M[13]*C[62]+M[14]*C[63]+M[15]*C[64]+M[16]*C[66]+M[17]*C[67]+M[18]*C[68]+M[19]*C[69];
  L[12] += M[10]*C[58]+M[11]*C[60]+M[12]*C[61]+M[13]*C[63]+M[14]*C[64]+M[15]*C[65]+M[16]*C[67]+M[17]*C[68]+M[18]*C[69]+M[19]*C[70];
  L[13] += M[10]*C[59]+M[11]*C[62]+M[12]*C[63]+M[13]*C[66]+M[14]*C[67]+M[15]*C[68]+M[16]*C[71]+M[17]*C[72]+M[18]*C[73]+M[19]*C[74];
  L[14] += M[10]*C[60]+M[11]*C[63]+M[12]*C[64]+M[13]*C[67]+M[14]*C[68]+M[15]*C[69]+M[16]*C[72]+M[17]*C[73]+M[18]*C[74]+M[19]*C[75];
  L[15] += M[10]*C[61]+M[11]*C[64]+M[12]*C[65]+M[13]*C[68]+M[14]*C[69]+M[15]*C[70]+M[16]*C[73]+M[17]*C[74]+M[18]*C[75]+M[19]*C[76];
  L[16] += M[10]*C[62]+M[11]*C[66]+M[12]*C[67]+M[13]*C[71]+M[14]*C[72]+M[15]*C[73]+M[16]*C[77]+M[17]*C[78]+M[18]*C[79]+M[19]*C[80];
  L[17] += M[10]*C[63]+M[11]*C[67]+M[12]*C[68]+M[13]*C[72]+M[14]*C[73]+M[15]*C[74]+M[16]*C[78]+M[17]*C[79]+M[18]*C[80]+M[19]*C[81];
  L[18] += M[10]*C[64]+M[11]*C[68]+M[12]*C[69]+M[13]*C[73]+M[14]*C[74]+M[15]*C[75]+M[16]*C[79]+M[17]*C[80]+M[18]*C[81]+M[19]*C[82];
  L[19] += M[10]*C[65]+M[11]*C[69]+M[12]*C[70]+M[13]*C[74]+M[14]*C[75]+M[15]*C[76]+M[16]*C[80]+M[17]*C[81]+M[18]*C[82]+M[19]*C[83];
  L[20] += M[4]*C[56]+M[5]*C[57]+M[6]*C[58]+M[7]*C[59]+M[8]*C[60]+M[9]*C[61];
  L[21] += M[4]*C[57]+M[5]*C[59]+M[6]*C[60]+M[7]*C[62]+M[8]*C[63]+M[9]*C[64];
  L[22] += M[4]*C[58]+M[5]*C[60]+M[6]*C[61]+M[7]*C[63]+M[8]*C[64]+M[9]*C[65];
  L[23] += M[4]*C[59]+M[5]*C[62]+M[6]*C[63]+M[7]*C[66]+M[8]*C[67]+M[9]*C[68];
  L[24] += M[4]*C[60]+M[5]*C[63]+M[6]*C[64]+M[7]*C[67]+M[8]*C[68]+M[9]*C[69];
  L[25] += M[4]*C[61]+M[5]*C[64]+M[6]*C[65]+M[7]*C[68]+M[8]*C[69]+M[9]*C[70];
  L[26] += M[4]*C[62]+M[5]*C[66]+M[6]*C[67]+M[7]*C[71]+M[8]*C[72]+M[9]*C[73];
  L[27] += M[4]*C[63]+M[5]*C[67]+M[6]*C[68]+M[7]*C[72]+M[8]*C[73]+M[9]*C[74];
  L[28] += M[4]*C[64]+M[5]*C[68]+M[6]*C[69]+M[7]*C[73]+M[8]*C[74]+M[9]*C[75];
  L[29] += M[4]*C[65]+M[5]*C[69]+M[6]*C[70]+M[7]*C[74]+M[8]*C[75]+M[9]*C[76];
  L[30] += M[4]*C[66]+M[5]*C[71]+M[6]*C[72]+M[7]*C[77]+M[8]*C[78]+M[9]*C[79];
  L[31] += M[4]*C[67]+M[5]*C[72]+M[6]*C[73]+M[7]*C[78]+M[8]*C[79]+M[9]*C[80];
  L[32] += M[4]*C[68]+M[5]*C[73]+M[6]*C[74]+M[7]*C[79]+M[8]*C[80]+M[9]*C[81];
  L[33] += M[4]*C[69]+M[5]*C[74]+M[6]*C[75]+M[7]*C[80]+M[8]*C[81]+M[9]*C[82];
  L[34] += M[4]*C[70]+M[5]*C[75]+M[6]*C[76]+M[7]*C[81]+M[8]*C[82]+M[9]*C[83];
  L[35] += M[1]*C[56]+M[2]*C[57]+M[3]*C[58];
  L[36] += M[1]*C[57]+M[2]*C[59]+M[3]*C[60];
  L[37] += M[1]*C[58]+M[2]*C[60]+M[3]*C[61];
  L[38] += M[1]*C[59]+M[2]*C[62]+M[3]*C[63];
  L[39] += M[1]*C[60]+M[2]*C[63]+M[3]*C[64];
  L[40] += M[1]*C[61]+M[2]*C[64]+M[3]*C[65];
  L[41] += M[1]*C[62]+M[2]*C[66]+M[3]*C[67];
  L[42] += M[1]*C[63]+M[2]*C[67]+M[3]*C[68];
  L[43] += M[1]*C[64]+M[2]*C[68]+M[3]*C[69];
  L[44] += M[1]*C[65]+M[2]*C[69]+M[3]*C[70];
  L[45] += M[1]*C[66]+M[2]*C[71]+M[3]*C[72];
  L[46] += M[1]*C[67]+M[2]*C[72]+M[3]*C[73];
  L[47] += M[1]*C[68]+M[2]*C[73]+M[3]*C[74];
  L[48] += M[1]*C[69]+M[2]*C[74]+M[3]*C[75];
  L[49] += M[1]*C[70]+M[2]*C[75]+M[3]*C[76];
  L[50] += M[1]*C[71]+M[2]*C[77]+M[3]*C[78];
  L[51] += M[1]*C[72]+M[2]*C[78]+M[3]*C[79];
  L[52] += M[1]*C[73]+M[2]*C[79]+M[3]*C[80];
  L[53] += M[1]*C[74]+M[2]*C[80]+M[3]*C[81];
  L[54] += M[1]*C[75]+M[2]*C[81]+M[3]*C[82];
  L[55] += M[1]*C[76]+M[2]*C[82]+M[3]*C[83];
}

void powerM(real *C, const real *dist) {
  C[1] = C[0] * dist[0];
  C[2] = C[0] * dist[1];
  C[3] = C[0] * dist[2];

  C[4] = C[1] * dist[0] / 2;
  C[5] = C[2] * dist[0];
  C[6] = C[3] * dist[0];
  C[7] = C[2] * dist[1] / 2;
  C[8] = C[3] * dist[1];
  C[9] = C[3] * dist[2] / 2;
  
  C[10] = C[4] * dist[0] / 3;
  C[11] = C[5] * dist[0] / 2;
  C[12] = C[6] * dist[0] / 2;
  C[13] = C[7] * dist[0];
  C[14] = C[8] * dist[0];
  C[15] = C[9] * dist[0];
  C[16] = C[7] * dist[1] / 3;
  C[17] = C[8] * dist[1] / 2;
  C[18] = C[9] * dist[1];
  C[19] = C[9] * dist[2] / 3;
  
  C[20] = C[10] * dist[0] / 4;
  C[21] = C[11] * dist[0] / 3;
  C[22] = C[12] * dist[0] / 3;
  C[23] = C[13] * dist[0] / 2;
  C[24] = C[14] * dist[0] / 2;
  C[25] = C[15] * dist[0] / 2;
  C[26] = C[16] * dist[0];
  C[27] = C[17] * dist[0];
  C[28] = C[18] * dist[0];
  C[29] = C[19] * dist[0];
  C[30] = C[16] * dist[1] / 4;
  C[31] = C[17] * dist[1] / 3;
  C[32] = C[18] * dist[1] / 2;
  C[33] = C[19] * dist[1];
  C[34] = C[19] * dist[2] / 4;

  C[35] = C[20] * dist[0] / 5;
  C[36] = C[21] * dist[0] / 4;
  C[37] = C[22] * dist[0] / 4;
  C[38] = C[23] * dist[0] / 3;
  C[39] = C[24] * dist[0] / 3;
  C[40] = C[25] * dist[0] / 3;
  C[41] = C[26] * dist[0] / 2;
  C[42] = C[27] * dist[0] / 2;
  C[43] = C[28] * dist[0] / 2;
  C[44] = C[29] * dist[0] / 2;
  C[45] = C[30] * dist[0];
  C[46] = C[31] * dist[0];
  C[47] = C[32] * dist[0];
  C[48] = C[33] * dist[0];
  C[49] = C[34] * dist[0];
  C[50] = C[30] * dist[1] / 5;
  C[51] = C[31] * dist[1] / 4;
  C[52] = C[32] * dist[1] / 3;
  C[53] = C[33] * dist[1] / 2;
  C[54] = C[34] * dist[1];
  C[55] = C[34] * dist[2] / 5;
}

void powerL(real *C, const real *dist) {
  C[1] = C[0] * dist[0];
  C[2] = C[0] * dist[1];
  C[3] = C[0] * dist[2];

  C[4] = C[1] * dist[0] / 2;
  C[5] = C[2] * dist[0];
  C[6] = C[3] * dist[0];
  C[7] = C[2] * dist[1] / 2;
  C[8] = C[3] * dist[1];
  C[9] = C[3] * dist[2] / 2;
  
  C[10] = C[4] * dist[0] / 3;
  C[11] = C[5] * dist[0] / 2;
  C[12] = C[6] * dist[0] / 2;
  C[13] = C[7] * dist[0];
  C[14] = C[8] * dist[0];
  C[15] = C[9] * dist[0];
  C[16] = C[7] * dist[1] / 3;
  C[17] = C[8] * dist[1] / 2;
  C[18] = C[9] * dist[1];
  C[19] = C[9] * dist[2] / 3;
  
  C[20] = C[10] * dist[0] / 4;
  C[21] = C[11] * dist[0] / 3;
  C[22] = C[12] * dist[0] / 3;
  C[23] = C[13] * dist[0] / 2;
  C[24] = C[14] * dist[0] / 2;
  C[25] = C[15] * dist[0] / 2;
  C[26] = C[16] * dist[0];
  C[27] = C[17] * dist[0];
  C[28] = C[18] * dist[0];
  C[29] = C[19] * dist[0];
  C[30] = C[16] * dist[1] / 4;
  C[31] = C[17] * dist[1] / 3;
  C[32] = C[18] * dist[1] / 2;
  C[33] = C[19] * dist[1];
  C[34] = C[19] * dist[2] / 4;

  C[35] = C[20] * dist[0] / 5;
  C[36] = C[21] * dist[0] / 4;
  C[37] = C[22] * dist[0] / 4;
  C[38] = C[23] * dist[0] / 3;
  C[39] = C[24] * dist[0] / 3;
  C[40] = C[25] * dist[0] / 3;
  C[41] = C[26] * dist[0] / 2;
  C[42] = C[27] * dist[0] / 2;
  C[43] = C[28] * dist[0] / 2;
  C[44] = C[29] * dist[0] / 2;
  C[45] = C[30] * dist[0];
  C[46] = C[31] * dist[0];
  C[47] = C[32] * dist[0];
  C[48] = C[33] * dist[0];
  C[49] = C[34] * dist[0];
  C[50] = C[30] * dist[1] / 5;
  C[51] = C[31] * dist[1] / 4;
  C[52] = C[32] * dist[1] / 3;
  C[53] = C[33] * dist[1] / 2;
  C[54] = C[34] * dist[1];
  C[55] = C[34] * dist[2] / 5;

  C[56] = C[35] * dist[0] / 6;
  C[57] = C[36] * dist[0] / 5;
  C[58] = C[37] * dist[0] / 5;
  C[59] = C[38] * dist[0] / 4;
  C[60] = C[39] * dist[0] / 4;
  C[61] = C[40] * dist[0] / 4;
  C[62] = C[41] * dist[0] / 3;
  C[63] = C[42] * dist[0] / 3;
  C[64] = C[43] * dist[0] / 3;
  C[65] = C[44] * dist[0] / 3;
  C[66] = C[45] * dist[0] / 2;
  C[67] = C[46] * dist[0] / 2;
  C[68] = C[47] * dist[0] / 2;
  C[69] = C[48] * dist[0] / 2;
  C[70] = C[49] * dist[0] / 2;
  C[71] = C[50] * dist[0];
  C[72] = C[51] * dist[0];
  C[73] = C[52] * dist[0];
  C[74] = C[53] * dist[0];
  C[75] = C[54] * dist[0];
  C[76] = C[55] * dist[0];
  C[77] = C[50] * dist[1] / 6;
  C[78] = C[51] * dist[1] / 5;
  C[79] = C[52] * dist[1] / 4;
  C[80] = C[53] * dist[1] / 3;
  C[81] = C[54] * dist[1] / 2;
  C[82] = C[55] * dist[1];
  C[83] = C[55] * dist[2] / 6;
}

void M2MSum(real *MI, const real *C, const real *MJ) {
  for( int i=1; i<MTERM; i++ ) MI[i] += MJ[i];

  MI[4] += C[1]*MJ[1];
  MI[5] += C[1]*MJ[2]+C[2]*MJ[1];
  MI[6] += C[1]*MJ[3]+C[3]*MJ[1];
  MI[7] += C[2]*MJ[2];
  MI[8] += C[2]*MJ[3]+C[3]*MJ[2];
  MI[9] += C[3]*MJ[3];

  MI[10] += C[1]*MJ[4]+C[4]*MJ[1];
  MI[11] += C[1]*MJ[5]+C[2]*MJ[4]+C[4]*MJ[2]+C[5]*MJ[1];
  MI[12] += C[1]*MJ[6]+C[3]*MJ[4]+C[4]*MJ[3]+C[6]*MJ[1];
  MI[13] += C[1]*MJ[7]+C[2]*MJ[5]+C[5]*MJ[2]+C[7]*MJ[1];
  MI[14] += C[1]*MJ[8]+C[2]*MJ[6]+C[3]*MJ[5]+C[5]*MJ[3]+C[6]*MJ[2]+C[8]*MJ[1];
  MI[15] += C[1]*MJ[9]+C[3]*MJ[6]+C[6]*MJ[3]+C[9]*MJ[1];
  MI[16] += C[2]*MJ[7]+C[7]*MJ[2];
  MI[17] += C[2]*MJ[8]+C[3]*MJ[7]+C[7]*MJ[3]+C[8]*MJ[2];
  MI[18] += C[2]*MJ[9]+C[3]*MJ[8]+C[8]*MJ[3]+C[9]*MJ[2];
  MI[19] += C[3]*MJ[9]+C[9]*MJ[3];

  MI[20] += C[1]*MJ[10]+C[4]*MJ[4]+C[10]*MJ[1];
  MI[21] += C[1]*MJ[11]+C[2]*MJ[10]+C[4]*MJ[5]+C[5]*MJ[4]+C[10]*MJ[2]+C[11]*MJ[1];
  MI[22] += C[1]*MJ[12]+C[3]*MJ[10]+C[4]*MJ[6]+C[6]*MJ[4]+C[10]*MJ[3]+C[12]*MJ[1];
  MI[23] += C[1]*MJ[13]+C[2]*MJ[11]+C[4]*MJ[7]+C[5]*MJ[5]+C[7]*MJ[4]+C[11]*MJ[2]+C[13]*MJ[1];
  MI[24] += C[1]*MJ[14]+C[2]*MJ[12]+C[3]*MJ[11]+C[4]*MJ[8]+C[5]*MJ[6]+C[6]*MJ[5]+C[8]*MJ[4]+C[11]*MJ[3]+C[12]*MJ[2]+C[14]*MJ[1];
  MI[25] += C[1]*MJ[15]+C[3]*MJ[12]+C[4]*MJ[9]+C[6]*MJ[6]+C[9]*MJ[4]+C[12]*MJ[3]+C[15]*MJ[1];
  MI[26] += C[1]*MJ[16]+C[2]*MJ[13]+C[5]*MJ[7]+C[7]*MJ[5]+C[13]*MJ[2]+C[16]*MJ[1];
  MI[27] += C[1]*MJ[17]+C[2]*MJ[14]+C[3]*MJ[13]+C[5]*MJ[8]+C[6]*MJ[7]+C[7]*MJ[6]+C[8]*MJ[5]+C[13]*MJ[3]+C[14]*MJ[2]+C[17]*MJ[1];
  MI[28] += C[1]*MJ[18]+C[2]*MJ[15]+C[3]*MJ[14]+C[5]*MJ[9]+C[6]*MJ[8]+C[8]*MJ[6]+C[9]*MJ[5]+C[14]*MJ[3]+C[15]*MJ[2]+C[18]*MJ[1];
  MI[29] += C[1]*MJ[19]+C[3]*MJ[15]+C[6]*MJ[9]+C[9]*MJ[6]+C[15]*MJ[3]+C[19]*MJ[1];
  MI[30] += C[2]*MJ[16]+C[7]*MJ[7]+C[16]*MJ[2];
  MI[31] += C[2]*MJ[17]+C[3]*MJ[16]+C[7]*MJ[8]+C[8]*MJ[7]+C[16]*MJ[3]+C[17]*MJ[2];
  MI[32] += C[2]*MJ[18]+C[3]*MJ[17]+C[7]*MJ[9]+C[8]*MJ[8]+C[9]*MJ[7]+C[17]*MJ[3]+C[18]*MJ[2];
  MI[33] += C[2]*MJ[19]+C[3]*MJ[18]+C[8]*MJ[9]+C[9]*MJ[8]+C[18]*MJ[3]+C[19]*MJ[2];
  MI[34] += C[3]*MJ[19]+C[9]*MJ[9]+C[19]*MJ[3];

  MI[35] += C[1]*MJ[20]+C[4]*MJ[10]+C[10]*MJ[4]+C[20]*MJ[1];
  MI[36] += C[1]*MJ[21]+C[2]*MJ[20]+C[4]*MJ[11]+C[5]*MJ[10]+C[10]*MJ[5]+C[11]*MJ[4]+C[20]*MJ[2]+C[21]*MJ[1];
  MI[37] += C[1]*MJ[22]+C[3]*MJ[20]+C[4]*MJ[12]+C[6]*MJ[10]+C[10]*MJ[6]+C[12]*MJ[4]+C[20]*MJ[3]+C[22]*MJ[1];
  MI[38] += C[1]*MJ[23]+C[2]*MJ[21]+C[4]*MJ[13]+C[5]*MJ[11]+C[7]*MJ[10]+C[10]*MJ[7]+C[11]*MJ[5]+C[13]*MJ[4]+C[21]*MJ[2]+C[23]*MJ[1];
  MI[39] += C[1]*MJ[24]+C[2]*MJ[22]+C[3]*MJ[21]+C[4]*MJ[14]+C[6]*MJ[12]+C[8]*MJ[11]+C[10]*MJ[10]+C[11]*MJ[8]+C[12]*MJ[6]+C[14]*MJ[4]+C[21]*MJ[3]+C[22]*MJ[2]+C[24]*MJ[1];
  MI[40] += C[1]*MJ[25]+C[3]*MJ[22]+C[4]*MJ[15]+C[6]*MJ[12]+C[9]*MJ[10]+C[10]*MJ[9]+C[12]*MJ[6]+C[15]*MJ[4]+C[22]*MJ[3]+C[25]*MJ[1];
  MI[41] += C[1]*MJ[26]+C[2]*MJ[23]+C[4]*MJ[16]+C[5]*MJ[13]+C[7]*MJ[11]+C[11]*MJ[7]+C[13]*MJ[5]+C[16]*MJ[4]+C[23]*MJ[2]+C[26]*MJ[1];
  MI[42] += C[1]*MJ[27]+C[2]*MJ[24]+C[3]*MJ[23]+C[4]*MJ[17]+C[5]*MJ[14]+C[6]*MJ[13]+C[7]*MJ[12]+C[8]*MJ[11]+C[11]*MJ[8]+C[12]*MJ[7]+C[13]*MJ[6]+C[14]*MJ[5]+C[17]*MJ[4]+C[23]*MJ[3]+C[24]*MJ[2]+C[27]*MJ[1];
  MI[43] += C[1]*MJ[28]+C[2]*MJ[25]+C[3]*MJ[24]+C[4]*MJ[18]+C[5]*MJ[15]+C[6]*MJ[14]+C[8]*MJ[12]+C[9]*MJ[11]+C[11]*MJ[9]+C[12]*MJ[8]+C[14]*MJ[6]+C[15]*MJ[5]+C[18]*MJ[4]+C[24]*MJ[3]+C[25]*MJ[2]+C[28]*MJ[1];
  MI[44] += C[1]*MJ[29]+C[3]*MJ[25]+C[4]*MJ[19]+C[6]*MJ[15]+C[9]*MJ[12]+C[12]*MJ[9]+C[15]*MJ[6]+C[19]*MJ[4]+C[25]*MJ[3]+C[29]*MJ[1];
  MI[45] += C[1]*MJ[30]+C[2]*MJ[26]+C[5]*MJ[16]+C[7]*MJ[13]+C[13]*MJ[7]+C[16]*MJ[5]+C[26]*MJ[2]+C[30]*MJ[1];
  MI[46] += C[1]*MJ[31]+C[2]*MJ[27]+C[3]*MJ[26]+C[5]*MJ[17]+C[7]*MJ[16]+C[8]*MJ[14]+C[13]*MJ[13]+C[14]*MJ[8]+C[16]*MJ[7]+C[17]*MJ[5]+C[26]*MJ[3]+C[27]*MJ[2]+C[31]*MJ[1];
  MI[47] += C[1]*MJ[32]+C[2]*MJ[28]+C[3]*MJ[27]+C[5]*MJ[18]+C[6]*MJ[17]+C[7]*MJ[15]+C[8]*MJ[14]+C[9]*MJ[13]+C[13]*MJ[9]+C[14]*MJ[8]+C[15]*MJ[7]+C[17]*MJ[6]+C[18]*MJ[5]+C[27]*MJ[3]+C[28]*MJ[2]+C[32]*MJ[1];
  MI[48] += C[1]*MJ[33]+C[2]*MJ[29]+C[3]*MJ[28]+C[5]*MJ[19]+C[8]*MJ[18]+C[9]*MJ[15]+C[14]*MJ[14]+C[15]*MJ[9]+C[18]*MJ[8]+C[19]*MJ[5]+C[28]*MJ[3]+C[29]*MJ[2]+C[33]*MJ[1];
  MI[49] += C[1]*MJ[34]+C[3]*MJ[29]+C[6]*MJ[19]+C[9]*MJ[15]+C[15]*MJ[9]+C[19]*MJ[6]+C[29]*MJ[3]+C[34]*MJ[1];
  MI[50] += C[2]*MJ[30]+C[7]*MJ[16]+C[16]*MJ[7]+C[30]*MJ[2];
  MI[51] += C[2]*MJ[31]+C[3]*MJ[30]+C[7]*MJ[17]+C[8]*MJ[16]+C[16]*MJ[8]+C[17]*MJ[7]+C[30]*MJ[3]+C[31]*MJ[2];
  MI[52] += C[2]*MJ[32]+C[3]*MJ[31]+C[7]*MJ[18]+C[8]*MJ[17]+C[9]*MJ[16]+C[16]*MJ[9]+C[17]*MJ[8]+C[18]*MJ[7]+C[31]*MJ[3]+C[32]*MJ[2];
  MI[53] += C[2]*MJ[33]+C[3]*MJ[32]+C[7]*MJ[19]+C[8]*MJ[18]+C[9]*MJ[17]+C[17]*MJ[9]+C[18]*MJ[8]+C[19]*MJ[7]+C[32]*MJ[3]+C[33]*MJ[2];
  MI[54] += C[2]*MJ[34]+C[3]*MJ[33]+C[8]*MJ[19]+C[9]*MJ[18]+C[18]*MJ[9]+C[19]*MJ[8]+C[33]*MJ[3]+C[34]*MJ[2];
  MI[55] += C[3]*MJ[34]+C[9]*MJ[19]+C[19]*MJ[9]+C[34]*MJ[3];
}

void L2LSum(real *LI, const real *C, const real *LJ) {
  LI[1] += C[1]*LJ[4]+C[2]*LJ[5]+C[3]*LJ[6];
  LI[2] += C[1]*LJ[5]+C[2]*LJ[7]+C[3]*LJ[8];
  LI[3] += C[1]*LJ[6]+C[2]*LJ[8]+C[3]*LJ[9];

  LI[1] += C[4]*LJ[10]+C[5]*LJ[11]+C[6]*LJ[12]+C[7]*LJ[13]+C[8]*LJ[14]+C[9]*LJ[15];
  LI[2] += C[4]*LJ[11]+C[5]*LJ[13]+C[6]*LJ[14]+C[7]*LJ[16]+C[8]*LJ[17]+C[9]*LJ[18];
  LI[3] += C[4]*LJ[12]+C[5]*LJ[14]+C[6]*LJ[15]+C[7]*LJ[17]+C[8]*LJ[18]+C[9]*LJ[19];
  LI[4] += C[1]*LJ[10]+C[2]*LJ[11]+C[3]*LJ[12];
  LI[5] += C[1]*LJ[11]+C[2]*LJ[13]+C[3]*LJ[14];
  LI[6] += C[1]*LJ[12]+C[2]*LJ[14]+C[3]*LJ[15];
  LI[7] += C[1]*LJ[13]+C[2]*LJ[16]+C[3]*LJ[17];
  LI[8] += C[1]*LJ[14]+C[2]*LJ[17]+C[3]*LJ[18];
  LI[9] += C[1]*LJ[15]+C[2]*LJ[18]+C[3]*LJ[19];

  LI[1] += C[10]*LJ[20]+C[11]*LJ[21]+C[12]*LJ[22]+C[13]*LJ[23]+C[14]*LJ[24]+C[15]*LJ[25]+C[16]*LJ[26]+C[17]*LJ[27]+C[18]*LJ[28]+C[19]*LJ[29];
  LI[2] += C[10]*LJ[21]+C[11]*LJ[23]+C[12]*LJ[24]+C[13]*LJ[26]+C[14]*LJ[27]+C[15]*LJ[28]+C[16]*LJ[30]+C[17]*LJ[31]+C[18]*LJ[32]+C[19]*LJ[33];
  LI[3] += C[10]*LJ[22]+C[11]*LJ[24]+C[12]*LJ[25]+C[13]*LJ[27]+C[14]*LJ[28]+C[15]*LJ[29]+C[16]*LJ[31]+C[17]*LJ[32]+C[18]*LJ[33]+C[19]*LJ[34];
  LI[4] += C[4]*LJ[20]+C[5]*LJ[21]+C[6]*LJ[22]+C[7]*LJ[23]+C[8]*LJ[24]+C[9]*LJ[25];
  LI[5] += C[4]*LJ[21]+C[5]*LJ[23]+C[6]*LJ[24]+C[7]*LJ[26]+C[8]*LJ[27]+C[9]*LJ[28];
  LI[6] += C[4]*LJ[22]+C[5]*LJ[24]+C[6]*LJ[25]+C[7]*LJ[27]+C[8]*LJ[28]+C[9]*LJ[29];
  LI[7] += C[4]*LJ[23]+C[5]*LJ[26]+C[6]*LJ[27]+C[7]*LJ[30]+C[8]*LJ[31]+C[9]*LJ[32];
  LI[8] += C[4]*LJ[24]+C[5]*LJ[27]+C[6]*LJ[28]+C[7]*LJ[31]+C[8]*LJ[32]+C[9]*LJ[33];
  LI[9] += C[4]*LJ[25]+C[5]*LJ[28]+C[6]*LJ[29]+C[7]*LJ[32]+C[8]*LJ[33]+C[9]*LJ[34];
  LI[10] += C[1]*LJ[20]+C[2]*LJ[21]+C[3]*LJ[22];
  LI[11] += C[1]*LJ[21]+C[2]*LJ[23]+C[3]*LJ[24];
  LI[12] += C[1]*LJ[22]+C[2]*LJ[24]+C[3]*LJ[25];
  LI[13] += C[1]*LJ[23]+C[2]*LJ[26]+C[3]*LJ[27];
  LI[14] += C[1]*LJ[24]+C[2]*LJ[27]+C[3]*LJ[28];
  LI[15] += C[1]*LJ[25]+C[2]*LJ[28]+C[3]*LJ[29];
  LI[16] += C[1]*LJ[26]+C[2]*LJ[30]+C[3]*LJ[31];
  LI[17] += C[1]*LJ[27]+C[2]*LJ[31]+C[3]*LJ[32];
  LI[18] += C[1]*LJ[28]+C[2]*LJ[32]+C[3]*LJ[33];
  LI[19] += C[1]*LJ[29]+C[2]*LJ[33]+C[3]*LJ[34]; // 282

  LI[1] += C[20]*LJ[35]+C[21]*LJ[36]+C[22]*LJ[37]+C[23]*LJ[38]+C[24]*LJ[39]+C[25]*LJ[40]+C[26]*LJ[41]+C[27]*LJ[42]+C[28]*LJ[43]+C[29]*LJ[44]+C[30]*LJ[45]+C[31]*LJ[46]+C[32]*LJ[47]+C[33]*LJ[48]+C[34]*LJ[49];
  LI[2] += C[20]*LJ[36]+C[21]*LJ[38]+C[22]*LJ[39]+C[23]*LJ[41]+C[24]*LJ[42]+C[25]*LJ[43]+C[26]*LJ[45]+C[27]*LJ[46]+C[28]*LJ[47]+C[29]*LJ[48]+C[30]*LJ[50]+C[31]*LJ[51]+C[32]*LJ[52]+C[33]*LJ[53]+C[34]*LJ[54];
  LI[3] += C[20]*LJ[37]+C[21]*LJ[39]+C[22]*LJ[40]+C[23]*LJ[42]+C[24]*LJ[43]+C[25]*LJ[44]+C[26]*LJ[46]+C[27]*LJ[47]+C[28]*LJ[48]+C[29]*LJ[49]+C[30]*LJ[51]+C[31]*LJ[52]+C[32]*LJ[53]+C[33]*LJ[54]+C[34]*LJ[55]; // 90
  LI[4] += C[10]*LJ[35]+C[11]*LJ[36]+C[12]*LJ[37]+C[13]*LJ[38]+C[14]*LJ[39]+C[15]*LJ[40]+C[16]*LJ[41]+C[17]*LJ[42]+C[18]*LJ[43]+C[19]*LJ[44];
  LI[5] += C[10]*LJ[36]+C[11]*LJ[38]+C[12]*LJ[39]+C[13]*LJ[41]+C[14]*LJ[42]+C[15]*LJ[43]+C[16]*LJ[45]+C[17]*LJ[46]+C[18]*LJ[47]+C[19]*LJ[48];
  LI[6] += C[10]*LJ[37]+C[11]*LJ[39]+C[12]*LJ[40]+C[13]*LJ[42]+C[14]*LJ[43]+C[15]*LJ[44]+C[16]*LJ[46]+C[17]*LJ[47]+C[18]*LJ[48]+C[19]*LJ[49];
  LI[7] += C[10]*LJ[38]+C[11]*LJ[41]+C[12]*LJ[42]+C[13]*LJ[45]+C[14]*LJ[46]+C[15]*LJ[47]+C[16]*LJ[50]+C[17]*LJ[51]+C[18]*LJ[52]+C[19]*LJ[53];
  LI[8] += C[10]*LJ[39]+C[11]*LJ[42]+C[12]*LJ[43]+C[13]*LJ[46]+C[14]*LJ[47]+C[15]*LJ[48]+C[16]*LJ[51]+C[17]*LJ[52]+C[18]*LJ[53]+C[19]*LJ[54];
  LI[9] += C[10]*LJ[40]+C[11]*LJ[43]+C[12]*LJ[44]+C[13]*LJ[47]+C[14]*LJ[48]+C[15]*LJ[49]+C[16]*LJ[52]+C[17]*LJ[53]+C[18]*LJ[54]+C[19]*LJ[55];
  LI[10] += C[4]*LJ[35]+C[5]*LJ[36]+C[6]*LJ[37]+C[7]*LJ[38]+C[8]*LJ[39]+C[9]*LJ[40];
  LI[11] += C[4]*LJ[36]+C[5]*LJ[38]+C[6]*LJ[39]+C[7]*LJ[41]+C[8]*LJ[42]+C[9]*LJ[43];
  LI[12] += C[4]*LJ[37]+C[5]*LJ[39]+C[6]*LJ[40]+C[7]*LJ[42]+C[8]*LJ[43]+C[9]*LJ[44];
  LI[13] += C[4]*LJ[38]+C[5]*LJ[41]+C[6]*LJ[42]+C[7]*LJ[45]+C[8]*LJ[46]+C[9]*LJ[47];
  LI[14] += C[4]*LJ[39]+C[5]*LJ[42]+C[6]*LJ[43]+C[7]*LJ[46]+C[8]*LJ[47]+C[9]*LJ[48];
  LI[15] += C[4]*LJ[40]+C[5]*LJ[43]+C[6]*LJ[44]+C[7]*LJ[47]+C[8]*LJ[48]+C[9]*LJ[49];
  LI[16] += C[4]*LJ[41]+C[5]*LJ[45]+C[6]*LJ[46]+C[7]*LJ[50]+C[8]*LJ[51]+C[9]*LJ[52];
  LI[17] += C[4]*LJ[42]+C[5]*LJ[46]+C[6]*LJ[47]+C[7]*LJ[51]+C[8]*LJ[52]+C[9]*LJ[53];
  LI[18] += C[4]*LJ[43]+C[5]*LJ[47]+C[6]*LJ[48]+C[7]*LJ[52]+C[8]*LJ[53]+C[9]*LJ[54];
  LI[19] += C[4]*LJ[44]+C[5]*LJ[48]+C[6]*LJ[49]+C[7]*LJ[53]+C[8]*LJ[54]+C[9]*LJ[55];
  LI[20] += C[1]*LJ[35]+C[2]*LJ[36]+C[3]*LJ[37];
  LI[21] += C[1]*LJ[36]+C[2]*LJ[38]+C[3]*LJ[39];
  LI[22] += C[1]*LJ[37]+C[2]*LJ[39]+C[3]*LJ[40];
  LI[23] += C[1]*LJ[38]+C[2]*LJ[41]+C[3]*LJ[42];
  LI[24] += C[1]*LJ[39]+C[2]*LJ[42]+C[3]*LJ[43];
  LI[25] += C[1]*LJ[40]+C[2]*LJ[43]+C[3]*LJ[44];
  LI[26] += C[1]*LJ[41]+C[2]*LJ[45]+C[3]*LJ[46];
  LI[27] += C[1]*LJ[42]+C[2]*LJ[46]+C[3]*LJ[47];
  LI[28] += C[1]*LJ[43]+C[2]*LJ[47]+C[3]*LJ[48];
  LI[29] += C[1]*LJ[44]+C[2]*LJ[48]+C[3]*LJ[49];
  LI[30] += C[1]*LJ[45]+C[2]*LJ[50]+C[3]*LJ[51];
  LI[31] += C[1]*LJ[46]+C[2]*LJ[51]+C[3]*LJ[52];
  LI[32] += C[1]*LJ[47]+C[2]*LJ[52]+C[3]*LJ[53];
  LI[33] += C[1]*LJ[48]+C[2]*LJ[53]+C[3]*LJ[54];
  LI[34] += C[1]*LJ[49]+C[2]*LJ[54]+C[3]*LJ[55];

  LI[1] += C[35]*LJ[56]+C[36]*LJ[57]+C[37]*LJ[58]+C[38]*LJ[59]+C[39]*LJ[60]+C[40]*LJ[61]+C[41]*LJ[62]+C[42]*LJ[63]+C[43]*LJ[64]+C[44]*LJ[65]+C[45]*LJ[66]+C[46]*LJ[67]+C[47]*LJ[68]+C[48]*LJ[69]+C[49]*LJ[70]+C[50]*LJ[71]+C[51]*LJ[72]+C[52]*LJ[73]+C[53]*LJ[74]+C[54]*LJ[75]+C[55]*LJ[76];
  LI[2] += C[35]*LJ[57]+C[36]*LJ[59]+C[37]*LJ[60]+C[38]*LJ[62]+C[39]*LJ[63]+C[40]*LJ[64]+C[41]*LJ[66]+C[42]*LJ[67]+C[43]*LJ[68]+C[44]*LJ[69]+C[45]*LJ[71]+C[46]*LJ[72]+C[47]*LJ[73]+C[48]*LJ[74]+C[49]*LJ[75]+C[50]*LJ[77]+C[51]*LJ[78]+C[52]*LJ[79]+C[53]*LJ[80]+C[54]*LJ[81]+C[55]*LJ[82];
  LI[3] += C[35]*LJ[58]+C[36]*LJ[60]+C[37]*LJ[61]+C[38]*LJ[63]+C[39]*LJ[64]+C[40]*LJ[65]+C[41]*LJ[67]+C[42]*LJ[68]+C[43]*LJ[69]+C[44]*LJ[70]+C[45]*LJ[72]+C[46]*LJ[73]+C[47]*LJ[74]+C[48]*LJ[75]+C[49]*LJ[76]+C[50]*LJ[78]+C[51]*LJ[79]+C[52]*LJ[80]+C[53]*LJ[81]+C[54]*LJ[82]+C[55]*LJ[83];
  LI[4] += C[20]*LJ[56]+C[21]*LJ[57]+C[22]*LJ[58]+C[23]*LJ[59]+C[24]*LJ[60]+C[25]*LJ[61]+C[26]*LJ[62]+C[27]*LJ[63]+C[28]*LJ[64]+C[29]*LJ[65]+C[30]*LJ[66]+C[31]*LJ[67]+C[32]*LJ[68]+C[33]*LJ[69]+C[34]*LJ[70];
  LI[5] += C[20]*LJ[57]+C[21]*LJ[59]+C[22]*LJ[60]+C[23]*LJ[62]+C[24]*LJ[63]+C[25]*LJ[64]+C[26]*LJ[66]+C[27]*LJ[67]+C[28]*LJ[68]+C[29]*LJ[69]+C[30]*LJ[71]+C[31]*LJ[72]+C[32]*LJ[73]+C[33]*LJ[74]+C[34]*LJ[75];
  LI[6] += C[20]*LJ[58]+C[21]*LJ[60]+C[22]*LJ[61]+C[23]*LJ[63]+C[24]*LJ[64]+C[25]*LJ[65]+C[26]*LJ[67]+C[27]*LJ[68]+C[28]*LJ[69]+C[29]*LJ[70]+C[30]*LJ[72]+C[31]*LJ[73]+C[32]*LJ[74]+C[33]*LJ[75]+C[34]*LJ[76];
  LI[7] += C[20]*LJ[59]+C[21]*LJ[62]+C[22]*LJ[63]+C[23]*LJ[66]+C[24]*LJ[67]+C[25]*LJ[68]+C[26]*LJ[71]+C[27]*LJ[72]+C[28]*LJ[73]+C[29]*LJ[74]+C[30]*LJ[77]+C[31]*LJ[78]+C[32]*LJ[79]+C[33]*LJ[80]+C[34]*LJ[81];
  LI[8] += C[20]*LJ[60]+C[21]*LJ[63]+C[22]*LJ[64]+C[23]*LJ[67]+C[24]*LJ[68]+C[25]*LJ[69]+C[26]*LJ[72]+C[27]*LJ[73]+C[28]*LJ[74]+C[29]*LJ[75]+C[30]*LJ[78]+C[31]*LJ[79]+C[32]*LJ[80]+C[33]*LJ[81]+C[34]*LJ[82];
  LI[9] += C[20]*LJ[61]+C[21]*LJ[64]+C[22]*LJ[65]+C[23]*LJ[68]+C[24]*LJ[69]+C[25]*LJ[70]+C[26]*LJ[73]+C[27]*LJ[74]+C[28]*LJ[75]+C[29]*LJ[76]+C[30]*LJ[79]+C[31]*LJ[80]+C[32]*LJ[81]+C[33]*LJ[82]+C[34]*LJ[83];
  LI[10] += C[10]*LJ[56]+C[11]*LJ[57]+C[12]*LJ[58]+C[13]*LJ[59]+C[14]*LJ[60]+C[15]*LJ[61]+C[16]*LJ[62]+C[17]*LJ[63]+C[18]*LJ[64]+C[19]*LJ[65];
  LI[11] += C[10]*LJ[57]+C[11]*LJ[59]+C[12]*LJ[60]+C[13]*LJ[62]+C[14]*LJ[63]+C[15]*LJ[64]+C[16]*LJ[66]+C[17]*LJ[67]+C[18]*LJ[68]+C[19]*LJ[69];
  LI[12] += C[10]*LJ[58]+C[11]*LJ[60]+C[12]*LJ[61]+C[13]*LJ[63]+C[14]*LJ[64]+C[15]*LJ[65]+C[16]*LJ[67]+C[17]*LJ[68]+C[18]*LJ[69]+C[19]*LJ[70];
  LI[13] += C[10]*LJ[59]+C[11]*LJ[62]+C[12]*LJ[63]+C[13]*LJ[66]+C[14]*LJ[67]+C[15]*LJ[68]+C[16]*LJ[71]+C[17]*LJ[72]+C[18]*LJ[73]+C[19]*LJ[74];
  LI[14] += C[10]*LJ[60]+C[11]*LJ[63]+C[12]*LJ[64]+C[13]*LJ[67]+C[14]*LJ[68]+C[15]*LJ[69]+C[16]*LJ[72]+C[17]*LJ[73]+C[18]*LJ[74]+C[19]*LJ[75];
  LI[15] += C[10]*LJ[61]+C[11]*LJ[64]+C[12]*LJ[65]+C[13]*LJ[68]+C[14]*LJ[69]+C[15]*LJ[70]+C[16]*LJ[73]+C[17]*LJ[74]+C[18]*LJ[75]+C[19]*LJ[76];
  LI[16] += C[10]*LJ[62]+C[11]*LJ[66]+C[12]*LJ[67]+C[13]*LJ[71]+C[14]*LJ[72]+C[15]*LJ[73]+C[16]*LJ[77]+C[17]*LJ[78]+C[18]*LJ[79]+C[19]*LJ[80];
  LI[17] += C[10]*LJ[63]+C[11]*LJ[67]+C[12]*LJ[68]+C[13]*LJ[72]+C[14]*LJ[73]+C[15]*LJ[74]+C[16]*LJ[78]+C[17]*LJ[79]+C[18]*LJ[80]+C[19]*LJ[81];
  LI[18] += C[10]*LJ[64]+C[11]*LJ[68]+C[12]*LJ[69]+C[13]*LJ[73]+C[14]*LJ[74]+C[15]*LJ[75]+C[16]*LJ[79]+C[17]*LJ[80]+C[18]*LJ[81]+C[19]*LJ[82];
  LI[19] += C[10]*LJ[65]+C[11]*LJ[69]+C[12]*LJ[70]+C[13]*LJ[74]+C[14]*LJ[75]+C[15]*LJ[76]+C[16]*LJ[80]+C[17]*LJ[81]+C[18]*LJ[82]+C[19]*LJ[83];
  LI[20] += C[4]*LJ[56]+C[5]*LJ[57]+C[6]*LJ[58]+C[7]*LJ[59]+C[8]*LJ[60]+C[9]*LJ[61];
  LI[21] += C[4]*LJ[57]+C[5]*LJ[59]+C[6]*LJ[60]+C[7]*LJ[62]+C[8]*LJ[63]+C[9]*LJ[64];
  LI[22] += C[4]*LJ[58]+C[5]*LJ[60]+C[6]*LJ[61]+C[7]*LJ[63]+C[8]*LJ[64]+C[9]*LJ[65];
  LI[23] += C[4]*LJ[59]+C[5]*LJ[62]+C[6]*LJ[63]+C[7]*LJ[66]+C[8]*LJ[67]+C[9]*LJ[68];
  LI[24] += C[4]*LJ[60]+C[5]*LJ[63]+C[6]*LJ[64]+C[7]*LJ[67]+C[8]*LJ[68]+C[9]*LJ[69];
  LI[25] += C[4]*LJ[61]+C[5]*LJ[64]+C[6]*LJ[65]+C[7]*LJ[68]+C[8]*LJ[69]+C[9]*LJ[70];
  LI[26] += C[4]*LJ[62]+C[5]*LJ[66]+C[6]*LJ[67]+C[7]*LJ[71]+C[8]*LJ[72]+C[9]*LJ[73];
  LI[27] += C[4]*LJ[63]+C[5]*LJ[67]+C[6]*LJ[68]+C[7]*LJ[72]+C[8]*LJ[73]+C[9]*LJ[74];
  LI[28] += C[4]*LJ[64]+C[5]*LJ[68]+C[6]*LJ[69]+C[7]*LJ[73]+C[8]*LJ[74]+C[9]*LJ[75];
  LI[29] += C[4]*LJ[65]+C[5]*LJ[69]+C[6]*LJ[70]+C[7]*LJ[74]+C[8]*LJ[75]+C[9]*LJ[76];
  LI[30] += C[4]*LJ[66]+C[5]*LJ[71]+C[6]*LJ[72]+C[7]*LJ[77]+C[8]*LJ[78]+C[9]*LJ[79];
  LI[31] += C[4]*LJ[67]+C[5]*LJ[72]+C[6]*LJ[73]+C[7]*LJ[78]+C[8]*LJ[79]+C[9]*LJ[80];
  LI[32] += C[4]*LJ[68]+C[5]*LJ[73]+C[6]*LJ[74]+C[7]*LJ[79]+C[8]*LJ[80]+C[9]*LJ[81];
  LI[33] += C[4]*LJ[69]+C[5]*LJ[74]+C[6]*LJ[75]+C[7]*LJ[80]+C[8]*LJ[81]+C[9]*LJ[82];
  LI[34] += C[4]*LJ[70]+C[5]*LJ[75]+C[6]*LJ[76]+C[7]*LJ[81]+C[8]*LJ[82]+C[9]*LJ[83];
  LI[35] += C[1]*LJ[56]+C[2]*LJ[57]+C[3]*LJ[58];
  LI[36] += C[1]*LJ[57]+C[2]*LJ[59]+C[3]*LJ[60];
  LI[37] += C[1]*LJ[58]+C[2]*LJ[60]+C[3]*LJ[61];
  LI[38] += C[1]*LJ[59]+C[2]*LJ[62]+C[3]*LJ[63];
  LI[39] += C[1]*LJ[60]+C[2]*LJ[63]+C[3]*LJ[64];
  LI[40] += C[1]*LJ[61]+C[2]*LJ[64]+C[3]*LJ[65];
  LI[41] += C[1]*LJ[62]+C[2]*LJ[66]+C[3]*LJ[67];
  LI[42] += C[1]*LJ[63]+C[2]*LJ[67]+C[3]*LJ[68];
  LI[43] += C[1]*LJ[64]+C[2]*LJ[68]+C[3]*LJ[69];
  LI[44] += C[1]*LJ[65]+C[2]*LJ[69]+C[3]*LJ[70];
  LI[45] += C[1]*LJ[66]+C[2]*LJ[71]+C[3]*LJ[72];
  LI[46] += C[1]*LJ[67]+C[2]*LJ[72]+C[3]*LJ[73];
  LI[47] += C[1]*LJ[68]+C[2]*LJ[73]+C[3]*LJ[74];
  LI[48] += C[1]*LJ[69]+C[2]*LJ[74]+C[3]*LJ[75];
  LI[49] += C[1]*LJ[70]+C[2]*LJ[75]+C[3]*LJ[76];
  LI[50] += C[1]*LJ[71]+C[2]*LJ[77]+C[3]*LJ[78];
  LI[51] += C[1]*LJ[72]+C[2]*LJ[78]+C[3]*LJ[79];
  LI[52] += C[1]*LJ[73]+C[2]*LJ[79]+C[3]*LJ[80];
  LI[53] += C[1]*LJ[74]+C[2]*LJ[80]+C[3]*LJ[81];
  LI[54] += C[1]*LJ[75]+C[2]*LJ[81]+C[3]*LJ[82];
  LI[55] += C[1]*LJ[76]+C[2]*LJ[82]+C[3]*LJ[83];
}

void L2PSum(real *TRG, const real *C, const real *L) {
  TRG[1] += C[1]*L[4]+C[2]*L[5]+C[3]*L[6];
  TRG[2] += C[1]*L[5]+C[2]*L[7]+C[3]*L[8];
  TRG[3] += C[1]*L[6]+C[2]*L[8]+C[3]*L[9];

  TRG[1] += C[4]*L[10]+C[5]*L[11]+C[6]*L[12]+C[7]*L[13]+C[8]*L[14]+C[9]*L[15];
  TRG[2] += C[4]*L[11]+C[5]*L[13]+C[6]*L[14]+C[7]*L[16]+C[8]*L[17]+C[9]*L[18];
  TRG[3] += C[4]*L[12]+C[5]*L[14]+C[6]*L[15]+C[7]*L[17]+C[8]*L[18]+C[9]*L[19];

  TRG[1] += C[10]*L[20]+C[11]*L[21]+C[12]*L[22]+C[13]*L[23]+C[14]*L[24]+C[15]*L[25]+C[16]*L[26]+C[17]*L[27]+C[18]*L[28]+C[19]*L[29];
  TRG[2] += C[10]*L[21]+C[11]*L[23]+C[12]*L[24]+C[13]*L[26]+C[14]*L[27]+C[15]*L[28]+C[16]*L[30]+C[17]*L[31]+C[18]*L[32]+C[19]*L[33];
  TRG[3] += C[10]*L[22]+C[11]*L[24]+C[12]*L[25]+C[13]*L[27]+C[14]*L[28]+C[15]*L[29]+C[16]*L[31]+C[17]*L[32]+C[18]*L[33]+C[19]*L[34];

  TRG[1] += C[20]*L[35]+C[21]*L[36]+C[22]*L[37]+C[23]*L[38]+C[24]*L[39]+C[25]*L[40]+C[26]*L[41]+C[27]*L[42]+C[28]*L[43]+C[29]*L[44]+C[30]*L[45]+C[31]*L[46]+C[32]*L[47]+C[33]*L[48]+C[34]*L[49];
  TRG[2] += C[20]*L[36]+C[21]*L[38]+C[22]*L[39]+C[23]*L[41]+C[24]*L[42]+C[25]*L[43]+C[26]*L[45]+C[27]*L[46]+C[28]*L[47]+C[29]*L[48]+C[30]*L[50]+C[31]*L[51]+C[32]*L[52]+C[33]*L[53]+C[34]*L[54];
  TRG[3] += C[20]*L[37]+C[21]*L[39]+C[22]*L[40]+C[23]*L[42]+C[24]*L[43]+C[25]*L[44]+C[26]*L[46]+C[27]*L[47]+C[28]*L[48]+C[29]*L[49]+C[30]*L[51]+C[31]*L[52]+C[32]*L[53]+C[33]*L[54]+C[34]*L[55];

  TRG[1] += C[35]*L[56]+C[36]*L[57]+C[37]*L[58]+C[38]*L[59]+C[39]*L[60]+C[40]*L[61]+C[41]*L[62]+C[42]*L[63]+C[43]*L[64]+C[44]*L[65]+C[45]*L[66]+C[46]*L[67]+C[47]*L[68]+C[48]*L[69]+C[49]*L[70]+C[50]*L[71]+C[51]*L[72]+C[52]*L[73]+C[53]*L[74]+C[54]*L[75]+C[55]*L[76];
  TRG[2] += C[35]*L[57]+C[36]*L[59]+C[37]*L[60]+C[38]*L[62]+C[39]*L[63]+C[40]*L[64]+C[41]*L[66]+C[42]*L[67]+C[43]*L[68]+C[44]*L[69]+C[45]*L[71]+C[46]*L[72]+C[47]*L[73]+C[48]*L[74]+C[49]*L[75]+C[50]*L[77]+C[51]*L[78]+C[52]*L[79]+C[53]*L[80]+C[54]*L[81]+C[55]*L[82];
  TRG[3] += C[35]*L[58]+C[36]*L[60]+C[37]*L[61]+C[38]*L[63]+C[39]*L[64]+C[40]*L[65]+C[41]*L[67]+C[42]*L[68]+C[43]*L[69]+C[44]*L[70]+C[45]*L[72]+C[46]*L[73]+C[47]*L[74]+C[48]*L[75]+C[49]*L[76]+C[50]*L[78]+C[51]*L[79]+C[52]*L[80]+C[53]*L[81]+C[54]*L[82]+C[55]*L[83];
}
