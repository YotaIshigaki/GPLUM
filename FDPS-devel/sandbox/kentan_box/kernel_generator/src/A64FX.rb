require "common.rb"
require "loop_fission.rb"

def get_declare_type_a64fx(type)
  case type
  when "F64" then
    decl = "svfloat64_t"
  when "F32" then
    decl = "svfloat32_t"
  when "S64" then
    decl = "svint64_t"
  when "S32" then
    decl = "svint32_t"
  when "U64" then
    decl = "svuint64_t"
  when "U32" then
    decl = "svuint32_t"
  when "F64vec1" then 
    decl = "svfloat64_t"
  when "F32vec1" then 
    decl = "svfloat32_t"
  when "F64vec2" then 
    decl = "svfloat64x2_t"
  when "F32vec2" then 
    decl = "svfloat32x2_t"
  when "F64vec3" then 
    decl = "svfloat64x3_t"
  when "F32vec3" then 
    decl = "svfloat32x3_t"
  when "F64vec4" then 
    decl = "svfloat64x4_t"
  when "F32vec4" then 
    decl = "svfloat32x4_t"
  when "F64vec" then 
    decl = "svfloat64x3_t"
  when "F32vec" then 
    decl = "svfloat32x3_t"
  else
    #abort "error: unsupported vector type of #{@type} for A64FX"
    decl = "UNSUPPORTED_VECTOR_TYPE"
  end
  decl
end

def get_declare_scalar_type(type,h=$varhash)
  case type
  when "F64" then
    decl = "float64_t"
  when "F32" then
    decl = "float32_t"
  when "S64" then
    decl = "int64_t"
  when "S32" then
    decl = "int32_t"
  when "U64" then
    decl = "uint64_t"
  when "U32" then
    decl = "uint32_t"
  when "F64vec1" then 
    decl = "float64_t"
  when "F32vec1" then 
    decl = "float32_t"
  when "F64vec2" then 
    decl = "float64x2_t"
  when "F32vec2" then 
    decl = "float32x2_t"
  when "F64vec3" then 
    decl = "float64x3_t"
  when "F32vec3" then 
    decl = "float32x3_t"
  when "F64vec4" then 
    decl = "float64x4_t"
  when "F32vec4" then 
    decl = "float32x4_t"
  when "F64vec" then 
    decl = "float64x3_t"
  when "F32vec" then 
    decl = "float32x3_t"
  else
    abort "error: unsupported scalar type of #{@type} for A64FX (get_declare_scalar_type)"
  end
  decl
end

def get_type_suffix(type)
  case type
  when "F64" then
    suffix = "f64"
  when "F32" then
    suffix = "f32"
  when "S64" then
    suffix = "s64"
  when "S32" then
    suffix = "s32"
  when "U64" then
    suffix = "u64"
  when "U32" then
    suffix = "u32"
  when /F64vec/ then
    suffix = "f64"
  when /F32vec/ then
    suffix = "f32"
  else
    abort "error: unsupported scalar type of #{@type} for A64FX (get_type_suffix)"
  end
  suffix
end

def get_pointer_cast(type)
  case type
  when "F64" then
    cast = "(double*)"
  when "F32" then
    cast = "(float*)"
  when "S64" then
    cast = "(long long*)"
  when "S32" then
    cast = "(int*)"
  when "U64" then
    cast = "(unsigned long long*)"
  when "U32" then
    cast = "(unsigned int*)"
  when "F64vec" then
    cast = "(double*)"
  when "F32vec" then
    cast = "(float*)"
  else
    abort "error: unsupported scalar type of #{@type} for A64FX"
  end
  cast
end


class Kernelprogram
  def reserved_func_def_a64fx(conversion_type)
    code = ""
    code += "svfloat32_t rsqrt(svbool_t pg,svfloat32_t op){\n"
    code += "svfloat32_t rinv = svrsqrte_f32(op);\n"
    code += "svfloat32_t h = svmul_f32_z(pg,op,rinv);\n"
    code += "h = svmsb_n_f32_z(pg,h,rinv,1.f);\n"
    code += "svfloat32_t poly = svmad_n_f32_z(pg,h,svdup_f32(0.375f),0.5f);\n"
    code += "poly = svmul_f32_z(pg,poly,h);\n"
    code += "rinv = svmad_f32_z(pg,rinv,poly,rinv);\n"
    code += "return rinv;\n"
    code += "}\n"

    code += "svfloat32_t sqrt(svbool_t pg,svfloat32_t op){ return svsqrt_f32_z(pg,op); }\n"
    # http://math-koshimizu.hatenablog.jp/entry/2017/07/28/083000
    code += "svfloat32_t inv(svbool_t pg,svfloat32_t op){\n"
    code += "svfloat32_t x1 = svrecpe_f32(op);\n"
    code += "svfloat32_t x2 = svmsb_n_f32_z(pg,op,x1,2.f);\n"
    code += "x2 = svmul_f32_z(pg,x2,x1);\n"
    code += "svfloat32_t ret = svmsb_n_f32_z(pg,op,x2,2.f);\n"
    code += "ret = svmul_f32_z(pg,ret,x2);\n"
    code += "return ret;\n"
    code += "}\n"
    code += "svfloat32_t max(svbool_t pg,svfloat32_t a,svfloat32_t b){ return svmax_f32_z(pg,a,b);}\n"
    code += "svfloat32_t min(svbool_t pg,svfloat32_t a,svfloat32_t b){ return svmin_f32_z(pg,a,b);}\n"  

    # code += "void transpose4x4(svfloat32x4_t& v){\n"
    # code += "const unsigned int tmp[16] = { 0, 2, 1, 3, 4, 6, 5, 7, 8,10, 9,11,12,14,13,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "v.v0 = svtbl_f32(v.v0,index);\n" # { 0, 2, 1, 3, 4, 6, 5, 7, 8,10, 9,11,12,14,13,15}
    # code += "v.v1 = svtbl_f32(v.v1,index);\n" #
    # code += "v.v2 = svtbl_f32(v.v2,index);\n" #
    # code += "v.v3 = svtbl_f32(v.v3,index);\n" #
    # code += "svfloat64_t xy0 = svreinterpret_f64_f32(svtrn1_f32(v.v0,v.v1));\n" # { 0,16, 1,17, 4,20, 5,21, 8,24, 9,25,12,28,13,29}
    # code += "svfloat64_t xy1 = svreinterpret_f64_f32(svtrn2_f32(v.v0,v.v1));\n" # { 2,18, 3,19, 6,22, 7,23,10,26,11,27,14,30,15,31}
    # code += "svfloat64_t zw0 = svreinterpret_f64_f32(svtrn1_f32(v.v2,v.v3));\n" # {32,48,33,49,36,52,37,53,40,56,41,57,44,60,45,61}
    # code += "svfloat64_t zw1 = svreinterpret_f64_f32(svtrn2_f32(v.v2,v.v3));\n" # {34,50,35,51,38,54,39,55,42,58,43,59,46,62,47,63}
    # code += "v.v0 = svreinterpret_f32_f64(svtrn1_f64(xy0,zw0));\n" # { 0,16,32,48, 4,20,36,52, 8,24,40,56,12,28,44,60}
    # code += "v.v1 = svreinterpret_f32_f64(svtrn2_f64(xy0,zw0));\n" # { 1,17,33,49, 5,21,37,53, 9,25,41,57,13,29,45,61}
    # code += "v.v2 = svreinterpret_f32_f64(svtrn1_f64(xy1,zw1));\n" # { 2,18,34,50, 6,22,38,54,10,26,42,58,14,30,46,62}
    # code += "v.v3 = svreinterpret_f32_f64(svtrn2_f64(xy1,zw1));\n" # { 3,19,35,51, 7,23,39,55,11,27,43,59,15,31,47,63}
    # code += "}\n"

    # code += "void gather8(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "const unsigned int tmp[16] = {  0,  8,  1,  9,  2, 10,  3, 11, 4, 12,  5, 13,  6, 14,  7, 15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "svfloat64_t a = svreinterpret_f64_f32(svtbl_f32(v0,index));\n"
    # code += "svfloat64_t b = svreinterpret_f64_f32(svtbl_f32(v1,index));\n"
    # code += "svfloat64_t c = svreinterpret_f64_f32(svtbl_f32(v2,index));\n"
    # code += "svfloat64_t d = svreinterpret_f64_f32(svtbl_f32(v3,index));\n"
    # code += "svfloat64_t e = svreinterpret_f64_f32(svtbl_f32(v4,index));\n"
    # code += "svfloat64_t f = svreinterpret_f64_f32(svtbl_f32(v5,index));\n"
    # code += "svfloat64_t g = svreinterpret_f64_f32(svtbl_f32(v6,index));\n"
    # code += "svfloat64_t h = svreinterpret_f64_f32(svtbl_f32(v7,index));\n"

    # code += "svfloat64_t ae0 = svzip1_f64(a,e);\n" # {  0,  8, 64, 72,  1,  9, 65, 73,  2, 10, 66, 74,  3, 11, 67, 75} is composed of v0,v4
    # code += "svfloat64_t ae1 = svzip2_f64(a,e);\n" # {  4, 12, 68, 76,  5, 13, 69, 77,  6, 14, 70, 78,  7, 15, 71, 79} is composed of v0,v4
    # code += "svfloat64_t bf0 = svzip1_f64(b,f);\n"
    # code += "svfloat64_t bf1 = svzip2_f64(b,f);\n" 
    # code += "svfloat64_t cg0 = svzip1_f64(c,g);\n" # { 32, 40, 96,104, 33, 41, 97,105, 34, 42, 98,106, 35, 43, 99,107} is composed of v2,v6
    # code += "svfloat64_t cg1 = svzip2_f64(c,g);\n" # { 36, 44,100,108, 37, 45,101,109, 38, 46,102,110, 39, 47,103,111} is composed of v2,v6
    # code += "svfloat64_t dh0 = svzip1_f64(d,h);\n"
    # code += "svfloat64_t dh1 = svzip2_f64(d,h);\n" 
    
    # code += "svfloat64_t aceg0 = svzip1_f64(ae0,cg0);\n" # {  0,  8, 32, 40, 64, 72, 96,104,  1,  9, 33, 41, 65, 73, 97,105} is composed of v0,v2,v4,v6
    # code += "svfloat64_t aceg1 = svzip2_f64(ae0,cg0);\n" # {  2, 10, 34, 42, 66, 74, 98,106,  3, 11, 35, 43, 67, 75, 99,107} is composed of v0,v2,v4,v6
    # code += "svfloat64_t aceg2 = svzip1_f64(ae1,cg1);\n"
    # code += "svfloat64_t aceg3 = svzip2_f64(ae1,cg1);\n"
    # code += "svfloat64_t bdfh0 = svzip1_f64(bf0,dh0);\n" # { 16, 24, 48, 56, 80, 88,112,120, 17, 25, 49, 57, 81, 89,113,121} is composed of v1,v3,v5,v7
    # code += "svfloat64_t bdfh1 = svzip2_f64(bf0,dh0);\n" # { 18, 26, 50, 58, 82, 90,114,122, 19, 27, 51, 59, 83, 91,115,123} is composed of v1,v3,v5,v7
    # code += "svfloat64_t bdfh2 = svzip1_f64(bf1,dh1);\n"
    # code += "svfloat64_t bdfh3 = svzip2_f64(bf1,dh1);\n"
    
    # code += "v0 = svreinterpret_f32_f64(svzip1_f64(aceg0,bdfh0));\n" # {  0,  8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96,104,112,120}
    # code += "v1 = svreinterpret_f32_f64(svzip2_f64(aceg0,bdfh0));\n" # {  1,  9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, 97,105,113,121}
    # code += "v2 = svreinterpret_f32_f64(svzip1_f64(aceg1,bdfh1));\n" # {  2, 10, 18, 26, 34, 42, 50, 58, 66, 74, 82, 90, 98,106,113,122}
    # code += "v3 = svreinterpret_f32_f64(svzip2_f64(aceg1,bdfh1));\n" # {  3, 11, 19, 27}
    # code += "v4 = svreinterpret_f32_f64(svzip1_f64(aceg2,bdfh2));\n" # {  4, 12, 20, 28}
    # code += "v5 = svreinterpret_f32_f64(svzip2_f64(aceg2,bdfh2));\n" # {  5, 13, 21, 29}
    # code += "v6 = svreinterpret_f32_f64(svzip1_f64(aceg3,bdfh3));\n" # {  6, 14, 22, 30}
    # code += "v7 = svreinterpret_f32_f64(svzip2_f64(aceg3,bdfh3));\n" # {  7, 15, 23, 31}
    # code += "}\n"

    # code += "void gather5(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,10,11,12,5,6,7,8,9,13,14,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code ++ "v0 = svtbl_f32(v0,index);\n"
    # code ++ "v1 = svtbl_f32(v1,index);\n"
    # code ++ "v2 = svtbl_f32(v2,index);\n"
    # code ++ "v3 = svtbl_f32(v3,index);\n"
    # code ++ "v4 = svtbl_f32(v4,index);\n"
    # code ++ "v5 = svtbl_f32(v5,index);\n"
    # code ++ "v6 = svtbl_f32(v6,index);\n"
    # code ++ "v7 = svtbl_f32(v7,index);\n"
    # code += "gather8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "}\n"     
    # code += "void gather6(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,5,12,13,6,7,8,9,10,11,14,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code ++ "v0 = svtbl_f32(v0,index);\n"
    # code ++ "v1 = svtbl_f32(v1,index);\n"
    # code ++ "v2 = svtbl_f32(v2,index);\n"
    # code ++ "v3 = svtbl_f32(v3,index);\n"
    # code ++ "v4 = svtbl_f32(v4,index);\n"
    # code ++ "v5 = svtbl_f32(v5,index);\n"
    # code ++ "v6 = svtbl_f32(v6,index);\n"
    # code ++ "v7 = svtbl_f32(v7,index);\n"
    # code += "gather8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "}\n"     
    # code += "void gather7(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,5,6,14,7,8,9,10,11,12,13,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code ++ "v0 = svtbl_f32(v0,index);\n"
    # code ++ "v1 = svtbl_f32(v1,index);\n"
    # code ++ "v2 = svtbl_f32(v2,index);\n"
    # code ++ "v3 = svtbl_f32(v3,index);\n"
    # code ++ "v4 = svtbl_f32(v4,index);\n"
    # code ++ "v5 = svtbl_f32(v5,index);\n"
    # code ++ "v6 = svtbl_f32(v6,index);\n"
    # code ++ "v7 = svtbl_f32(v7,index);\n"
    # code += "gather8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "}\n"     
    
    # code += "void scatter8(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "svfloat32_t ae0 = svzip1_f32(v0,v4);\n" # {  0, 64,  1, 65,  2, 66,  3, 67,  4, 68,  5, 69,  6, 70,  7, 71} is composed of v0,v4
    # code += "svfloat32_t ae1 = svzip2_f32(v0,v4);\n" # {  8, 72,  9, 73, 10, 74, 11, 75, 12, 76, 13, 77, 14, 78, 15, 79} is composed of v0,v4
    # code += "svfloat32_t bf0 = svzip1_f32(v1,v5);\n"
    # code += "svfloat32_t bf1 = svzip2_f32(v1,v5);\n"
    # code += "svfloat32_t cg0 = svzip1_f32(v2,v6);\n"
    # code += "svfloat32_t cg1 = svzip2_f32(v2,v6);\n"
    # code += "svfloat32_t dh0 = svzip1_f32(v3,v7);\n"
    # code += "svfloat32_t dh1 = svzip2_f32(v3,v7);\n"
    # code += "svfloat32_t aceg0 = svzip1_f32(ae0,cg0);\n" # {  0, 32, 64, 96,  1, 33, 65, 97,  2, 34, 66, 98,  3, 35, 67, 99} is composed of v0,v2,v4,v6
    # code += "svfloat32_t aceg1 = svzip2_f32(ae0,cg0);\n" # {  4, 36, 68,100,  5, 37, 69,101,  6, 38, 70,102,  7, 39, 71,103} is composed of v0,v2,v4,v6
    # code += "svfloat32_t aceg2 = svzip1_f32(ae1,cg1);\n"
    # code += "svfloat32_t aceg3 = svzip2_f32(ae1,cg1);\n"
    # code += "svfloat32_t bdfh0 = svzip1_f32(bf0,dh0);\n" # { 16, 48, 80,112, 17, 49, 81,113, 18, 50, 82,114, 19, 51, 83,115} is composed of v1,v3,v5,v7
    # code += "svfloat32_t bdfh1 = svzip2_f32(bf0,dh0);\n" # { 20, 52, 84,116, 21, 53, 85,117, 22, 54, 86,118, 23, 55, 87,119} is composed of v1,v3,v5,v7
    # code += "svfloat32_t bdfh2 = svzip1_f32(bf1,dh1);\n"
    # code += "svfloat32_t bdfh3 = svzip2_f32(bf1,dh1);\n"
    # code += "v0 = svzip1_f32(aceg0,bdfh0);\n" # {  0, 16, 32, 48, 64, 80, 96,112,  1, 17, 33, 49, 65, 81, 97,113}
    # code += "v1 = svzip2_f32(aceg0,bdfh0);\n" # {  2, 18, 34, 50, 66, 82, 98,114,  3, 19, 35, 51, 67, 83, 99,115} 
    # code += "v2 = svzip1_f32(aceg1,bdfh1);\n" # {  4, 20, 36, 52, 68, 84,100,116,  5, 21, 37, 53, 69, 85,101,117}
    # code += "v3 = svzip2_f32(aceg1,bdfh1);\n" # {  6, 22, 38, 54, 70, 86,102,118,  7, 23, 39, 55, 71, 87,103,119}
    # code += "v4 = svzip1_f32(aceg2,bdfh2);\n" # {  8, 24, 40, 56, 72, 88,104,120,  9, 25, 41, 57, 73, 89,105,121}
    # code += "v5 = svzip2_f32(aceg2,bdfh2);\n" # { 10, 26, 42, 58, 74, 90,106,122, 11, 27, 43, 59, 75, 91,107,123}
    # code += "v6 = svzip1_f32(aceg3,bdfh3);\n" # { 12, 28, 44, 60, 76, 92,108,124, 13, 29, 45, 61, 77, 93,109,125}
    # code += "v7 = svzip2_f32(aceg3,bdfh3);\n" # { 14, 30, 46, 62, 78, 94,110,126, 15, 31, 47, 63, 79, 94,111,127
    # code += "}\n"
    
    # code += "void scatter5(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "scatter8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,8,9,10,11,12,5,6,7,13,14,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "v0 = svtbl_f32(v0,index);\n"
    # code += "v1 = svtbl_f32(v1,index);\n"
    # code += "v2 = svtbl_f32(v2,index);\n"
    # code += "v3 = svtbl_f32(v3,index);\n"
    # code += "v4 = svtbl_f32(v4,index);\n"
    # code += "v5 = svtbl_f32(v5,index);\n"
    # code += "v6 = svtbl_f32(v6,index);\n"
    # code += "v7 = svtbl_f32(v7,index);\n"
    # code += "}\n"

    # code += "void scatter6(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "scatter8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,5,8,9,10,11,12,13,6,7,14,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "v0 = svtbl_f32(v0,index);\n"
    # code += "v1 = svtbl_f32(v1,index);\n"
    # code += "v2 = svtbl_f32(v2,index);\n"
    # code += "v3 = svtbl_f32(v3,index);\n"
    # code += "v4 = svtbl_f32(v4,index);\n"
    # code += "v5 = svtbl_f32(v5,index);\n"
    # code += "v6 = svtbl_f32(v6,index);\n"
    # code += "v7 = svtbl_f32(v7,index);\n"
    # code += "}\n"

    # code += "void scatter7(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "scatter8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,5,6,8,9,10,11,12,13,14,7,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "v0 = svtbl_f32(v0,index);\n"
    # code += "v1 = svtbl_f32(v1,index);\n"
    # code += "v2 = svtbl_f32(v2,index);\n"
    # code += "v3 = svtbl_f32(v3,index);\n"
    # code += "v4 = svtbl_f32(v4,index);\n"
    # code += "v5 = svtbl_f32(v5,index);\n"
    # code += "v6 = svtbl_f32(v6,index);\n"
    # code += "}\n"

    # code += "void transpose8x8(svfloat32x4_t& v0,svfloat32x4_t v1){\n"
    # code += "const unsigned int tmp[16] = {0, 8, 4,12, 1, 9, 5,13, 2,10, 6,14, 3,11, 7,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "v0.v0 = svtbl_f32(v0.v0,index);\n" # {0, 8, 4,12, 1, 9, 5,13, 2,10, 6,14, 3,11, 7,15}
    # code += "v0.v1 = svtbl_f32(v0.v1,index);\n"
    # code += "v0.v2 = svtbl_f32(v0.v2,index);\n"
    # code += "v0.v3 = svtbl_f32(v0.v3,index);\n"
    # code += "v1.v0 = svtbl_f32(v1.v0,index);\n" # {0, 8, 4,12, 1, 9, 5,13, 2,10, 6,14, 3,11, 7,15}
    # code += "v1.v1 = svtbl_f32(v1.v1,index);\n"
    # code += "v1.v2 = svtbl_f32(v1.v2,index);\n"
    # code += "v1.v3 = svtbl_f32(v1.v3,index);\n"
    # code += "svfloat64_t ab0 = svreinterpret_f64_f32(svtrn1_f32(v0.v0,v0.v1));\n" # ab{0, 4, 1, 5, 2, 6, 3, 7}
    # code += "svfloat64_t ab1 = svreinterpret_f64_f32(svtrn2_f32(v0.v0,v0.v1));\n" # ab{8,12, 9,13,10,14,11,15}
    # code += "svfloat64_t cd0 = svreinterpret_f64_f32(svtrn1_f32(v0.v2,v0.v3));\n" # cd{0, 4, 1, 5, 2, 6, 3, 7}
    # code += "svfloat64_t cd1 = svreinterpret_f64_f32(svtrn2_f32(v0.v2,v0.v3));\n" # cd{8,12, 9,13,10,14,11,15}
    # code += "svfloat64_t ef0 = svreinterpret_f64_f32(svtrn1_f32(v1.v0,v1.v1));\n" # ef{0, 4, 1, 5, 2, 6, 3, 7}
    # code += "svfloat64_t ef1 = svreinterpret_f64_f32(svtrn2_f32(v1.v0,v1.v1));\n" # ef{8,12, 9,13,10,14,11,15}
    # code += "svfloat64_t gh0 = svreinterpret_f64_f32(svtrn1_f32(v1.v2,v1.v3));\n" # gh{0, 4, 1, 5, 2, 6, 3, 7}
    # code += "svfloat64_t gh1 = svreinterpret_f64_f32(svtrn2_f32(v1.v2,v1.v3));\n" # gh{8,12, 9,13,10,14,11,15}
    # code += "v0.v0 = svreinterpret_f32_f64(svtrn1_f64(ab0,cd0));\n" # abcd{ 0, 1, 2, 3}
    # code += "v0.v1 = svreinterpret_f32_f64(svtrn2_f64(ab0,cd0));\n" # abcd{ 4, 5, 6, 7}
    # code += "v0.v2 = svreinterpret_f32_f64(svtrn1_f64(ab1,cd1));\n" # abcd{ 8, 9,10,11}
    # code += "v0.v3 = svreinterpret_f32_f64(svtrn2_f64(ab1,cd1));\n" # abcd{12,13,14,15}
    # code += "v1.v0 = svreinterpret_f32_f64(svtrn1_f64(ef0,gh0));\n" # efgh{ 0, 1, 2, 3}
    # code += "v1.v1 = svreinterpret_f32_f64(svtrn2_f64(ef0,gh0));\n" # efgh{ 4, 5, 6, 7}
    # code += "v1.v2 = svreinterpret_f32_f64(svtrn1_f64(ef1,gh1));\n" # efgh{ 8, 9,10,11}
    # code += "v1.v3 = svreinterpret_f32_f64(svtrn2_f64(ef1,gh1));\n" # efgh{12,13,14,15}
    # code += "v0.v0 = svreinterpret_f32_f64(svtrn1_f64(ab0,cd0));\n" # abcdefgh{ 0, 1}
    # code += "v0.v1 = svreinterpret_f32_f64(svtrn2_f64(ab0,cd0));\n" # abcdefgh{ 2, 3}
    # code += "v0.v2 = svreinterpret_f32_f64(svtrn1_f64(ab1,cd1));\n" # abcdefgh{ 4, 5}
    # code += "v0.v3 = svreinterpret_f32_f64(svtrn2_f64(ab1,cd1));\n" # abcdefgh{ 6, 7}
    # code += "v1.v0 = svreinterpret_f32_f64(svtrn1_f64(ef0,gh0));\n" # abcdefgh{ 8, 9}
    # code += "v1.v1 = svreinterpret_f32_f64(svtrn2_f64(ef0,gh0));\n" # abcdefgh{10,11}
    # code += "v1.v2 = svreinterpret_f32_f64(svtrn1_f64(ef1,gh1));\n" # abcdefgh{12,13}
    # code += "v1.v3 = svreinterpret_f32_f64(svtrn2_f64(ef1,gh1));\n" # abcdefgh{14,15}
    # code += ""
    # code += "}\n"

    # code += "void transpose16x16(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7,svfloat32_t& v8,svfloat32_t& v9,svfloat32_t& v10,svfloat32_t& v11,svfloat32_t& v12,svfloat32_t& v13,svfloat32_t& v14,svfloat32_t& v15){\n"
    # code += "svfloat32_t ai0 = svzip1_f32(v0,v8);\n" # {  0, 64,  1, 65,  2, 66,  3, 67,  4, 68,  5, 69,  6, 70,  7, 71} is composed of v0,v8
    # code += "svfloat32_t ai1 = svzip2_f32(v0,v8);\n" # {  8, 72,  9, 73, 10, 74, 11, 75, 12, 76, 13, 77, 14, 78, 15, 79} is composed of v0,v8
    # code += "svfloat32_t bj0 = svzip1_f32(v1,v9);\n"
    # code += "svfloat32_t bj1 = svzip2_f32(v1,v9);\n"
    # code += "svfloat32_t ck0 = svzip1_f32(v2,v10);\n"
    # code += "svfloat32_t ck1 = svzip2_f32(v2,v10);\n"
    # code += "svfloat32_t dl0 = svzip1_f32(v3,v11);\n"
    # code += "svfloat32_t dl1 = svzip2_f32(v3,v11);\n"
    # code += "svfloat32_t em0 = svzip1_f32(v4,v12);\n"
    # code += "svfloat32_t em1 = svzip2_f32(v4,v12);\n"
    # code += "svfloat32_t fn0 = svzip1_f32(v5,v13);\n"
    # code += "svfloat32_t fn1 = svzip2_f32(v5,v13);\n"
    # code += "svfloat32_t go0 = svzip1_f32(v6,v14);\n"
    # code += "svfloat32_t go1 = svzip2_f32(v6,v14);\n"
    # code += "svfloat32_t hp0 = svzip1_f32(v7,v15);\n"
    # code += "svfloat32_t hp1 = svzip2_f32(v7,v15);\n"
    
    # code += "svfloat32_t aeim0 = svzip1_f32(ai0,em0);\n" # {  0, 64,128,192,  1, 65,129,193,  2, 68,130,194,  3, 67,131,195} is composed of v0,v4,v8,v12
    # code += "svfloat32_t aeim1 = svzip2_f32(ai0,em0);\n"
    # code += "svfloat32_t aeim2 = svzip1_f32(ai1,em1);\n"
    # code += "svfloat32_t aeim3 = svzip2_f32(ai1,em1);\n"
    # code += "svfloat32_t bfjn0 = svzip1_f32(bj0,fn0);\n"
    # code += "svfloat32_t bfjn1 = svzip2_f32(bj0,fn0);\n"
    # code += "svfloat32_t bfjn2 = svzip1_f32(bj1,fn1);\n"
    # code += "svfloat32_t bfjn3 = svzip2_f32(bj1,fn1);\n"
    # code += "svfloat32_t cgko0 = svzip1_f32(ck0,go0);\n"
    # code += "svfloat32_t cgko1 = svzip2_f32(ck0,go0);\n"
    # code += "svfloat32_t cgko2 = svzip1_f32(ck1,go1);\n"
    # code += "svfloat32_t cgko3 = svzip2_f32(ck1,go1);\n"
    # code += "svfloat32_t dhlp0 = svzip1_f32(dl0,hp0);\n"
    # code += "svfloat32_t dhlp1 = svzip2_f32(dl0,hp0);\n"
    # code += "svfloat32_t dhlp2 = svzip1_f32(dl1,hp1);\n"
    # code += "svfloat32_t dhlp3 = svzip2_f32(dl1,hp1);\n"
    
    # code += "svfloat32_t acegikmo0 = svzip1_f32(aeim0,cgko0);\n" # {  0, 32, 64, 96,128,160,192,224,  1, 33, 65, 97,129,161,193,225} is composed of v0,v2,v4,v6,v8,v10,v12,v14
    # code += "svfloat32_t acegikmo1 = svzip2_f32(aeim0,cgko0);\n" # {  2, 34, 68, 98,130,162,194,226,  3, 35, 67, 99,131,163,195,227} is composed of v0,v2,v4,v6,v8,v10,v12,v14
    # code += "svfloat32_t acegikmo2 = svzip1_f32(aeim1,cgko1);\n"
    # code += "svfloat32_t acegikmo3 = svzip2_f32(aeim1,cgko1);\n"
    # code += "svfloat32_t acegikmo4 = svzip1_f32(aeim2,cgko2);\n"
    # code += "svfloat32_t acegikmo5 = svzip2_f32(aeim2,cgko2);\n"
    # code += "svfloat32_t acegikmo6 = svzip1_f32(aeim3,cgko3);\n"
    # code += "svfloat32_t acegikmo7 = svzip2_f32(aeim3,cgko3);\n"
    # code += "svfloat32_t bdfhjlnp0 = svzip1_f32(bfjn0,dhlp0);\n" # { 16, 48, 80,112,144,176,208,240, 17, 49, 81,113,145,177,209,241} is composed of v1,v3,v5,v7,v9,v11,v13,v15
    # code += "svfloat32_t bdfhjlnp1 = svzip2_f32(bfjn0,dhlp0);\n" # { 18, 50, 82,114,146,178,210,242, 19, 51, 83,115,147,179,211,243} is composed of v1,v3,v5,v7,v9,v11,v13,v15
    # code += "svfloat32_t bdfhjlnp2 = svzip1_f32(bfjn1,dhlp1);\n"
    # code += "svfloat32_t bdfhjlnp3 = svzip2_f32(bfjn1,dhlp1);\n"
    # code += "svfloat32_t bdfhjlnp4 = svzip1_f32(bfjn2,dhlp2);\n"
    # code += "svfloat32_t bdfhjlnp5 = svzip2_f32(bfjn2,dhlp2);\n"
    # code += "svfloat32_t bdfhjlnp6 = svzip1_f32(bfjn3,dhlp3);\n"
    # code += "svfloat32_t bdfhjlnp7 = svzip2_f32(bfjn3,dhlp3);\n"
    # code += "v0  = svzip1_f32(acegikmo0,bdfhjlnp0);\n" # {  0, 16, 32, 48, 64, 80, 96,112,128,144,160,176,192,208,224,240}
    # code += "v1  = svzip2_f32(acegikmo0,bdfhjlnp0);\n" # {  1, 17, 33, 49, 65, 81, 97,113,129,145,161,177,193,209,225,241}
    # code += "v2  = svzip1_f32(acegikmo1,bdfhjlnp1);\n" # {  2, 18, 34, 50, 66, 82, 98,114,130,146,162,178,194,210,226,242}
    # code += "v3  = svzip2_f32(acegikmo1,bdfhjlnp1);\n" # {  3, 19, 35, 51, 67, 83, 99,115,131,147,163,179,195,211,227,243}
    # code += "v4  = svzip1_f32(acegikmo2,bdfhjlnp2);\n" # {  4, 20, 36, 52, 68, 84,100,116,132,148,164,180,196,212,228,244}
    # code += "v5  = svzip2_f32(acegikmo2,bdfhjlnp2);\n" # {  5, 21, 37, 53, 69, 85,101,117,133,149,165,181,197,213,229,245}
    # code += "v6  = svzip1_f32(acegikmo3,bdfhjlnp3);\n" # {  6, 22, 38, 54, 70, 86,102,118,134,150,166,182,198,214,230,246}
    # code += "v7  = svzip2_f32(acegikmo3,bdfhjlnp3);\n" # {  7, 23, 39, 55, 71, 87,103,119,135,151,167,183,199,215,231,247}
    # code += "v8  = svzip1_f32(acegikmo4,bdfhjlnp4);\n" # {  8, 24, 40, 56, 72, 88,104,120,136,152,168,184,200,216,232,248}
    # code += "v9  = svzip2_f32(acegikmo4,bdfhjlnp4);\n" # {  9, 25, 41, 57, 73, 89,105,121,137,153,169,185,201,217,233,249}
    # code += "v10 = svzip1_f32(acegikmo5,bdfhjlnp5);\n" # { 10, 26, 42, 58, 74, 90,106,122,138,154,170,186,202,218,234,250} 
    # code += "v11 = svzip2_f32(acegikmo5,bdfhjlnp5);\n" # { 11, 27, 43, 59, 75, 91,107,123,139,155,171,187,203,219,235,251}
    # code += "v12 = svzip1_f32(acegikmo6,bdfhjlnp6);\n" # { 12, 28, 44, 60, 76, 92,108,124,140,156,172,188,204,220,236,252}
    # code += "v13 = svzip2_f32(acegikmo6,bdfhjlnp6);\n" # { 13, 29, 45, 61, 77, 93,109,125,141,157,173,189,205,221,237,253}
    # code += "v14 = svzip1_f32(acegikmo7,bdfhjlnp7);\n" # { 14, 30, 46, 62, 78, 94,110,126,142,158,174,190,206,222,238,254}
    # code += "v15 = svzip2_f32(acegikmo7,bdfhjlnp7);\n" # { 15, 31, 47, 63, 79, 95,111,127,143,159,175,191,207,223,239,255}
    # code += "}\n"

    # code += "void gather16(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7,svfloat32_t& v8,svfloat32_t& v9,svfloat32_t& v10,svfloat32_t& v11,svfloat32_t& v12,svfloat32_t& v13,svfloat32_t& v14,svfloat32_t& v15){ transpose16x16(v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15); };\n"
    # code += "void scatter16(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7,svfloat32_t& v8,svfloat32_t& v9,svfloat32_t& v10,svfloat32_t& v11,svfloat32_t& v12,svfloat32_t& v13,svfloat32_t& v14,svfloat32_t& v15){ transpose16x16(v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15); };\n"
    code
  end

  def kernel_body_a64fx_f32(conversion_type)
    code = ""
    epi_tab = []
    xyzw = ["v0","v1","v2","v3"]
    epi_count = 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "EPI"
        name = v[0]
        type = v[1][1]
        if type == "F32vec"
          epi_tab += [name+".v0","#{epi_count}"]
          epi_count += 1
          epi_tab += [name+".v1","#{epi_count}"]
          epi_count += 1
          epi_tab += [name+".v2","#{epi_count}"]
          epi_count += 1
        elsif type == "F32" || type == "U32" || type == "S32"
          epi_tab += [name,"#{epi_count}"]
          epi_count += 1
        else
          abort "FP64 exist in #{conversion_type} optimization"
        end
      end
    }
    epi_tab = Hash[*epi_tab]
    abort "# of EPI class member > 16 is not supported" if epi_count > 16

    force_count = 0
    force_tab = []
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "FORCE"
        name = v[0]
        type = v[1][1]
        if type == "F32vec"
          force_tab += [name+".v0","#{force_count}"]
          force_count += 1
          force_tab += [name+".v1","#{force_count}"]
          force_count += 1
          force_tab += [name+".v2","#{force_count}"]
          force_count += 1
        elsif type == "F32" || type == "U32" || type == "S32"
          force_tab += [name,"#{force_count}"]
          force_count += 1
        else
          abort "Only F32[vec] is supported in #{conversion_type} optimization"
        end
      end
    }
    force_tab = Hash[*force_tab]
    abort "# of FORCE class member > 16 is not supported" if force_count > 16
    code += "for(int i=0;i<((ni+15)/16)*16;i+=16){\n"
    # declare i particle
    code += "#{$current_predicate} = svwhilelt_b32_s32(i,ni);\n"
    abort "pg at the top of i loop must be pg0" if $pg_count > 0

    if epi_count <= 4 then
      code += "svfloat32x#{  epi_count}_t r = svld#{  epi_count}_f32(#{$current_predicate},(float*)&epi[i]);\n"
      for i in 0...epi_count do
        code += "svfloat32_t& r#{i} = r.v#{i};\n"
      end
    else
      if epi_count <= 8 then
        offset = 2
      else
        offset = 1
      end
      for i in 0...(16/offset) do
        code += "svfloat32_t r#{i};\n"
      end
      for i in 0...(16/offset) do
        code += "r#{i} = svld1_f32(#{$current_predicate},(float*)&epi[i+#{offset*i}]);\n"
      end
      code += "gather#{16/offset}("
      for i in 0...(16/offset) do
        code += "," if i > 0
        code += "r#{i}"
      end
      code += ");\n"
    end

    if force_count <= 4 then
      code += "svfloat32x#{force_count}_t f = svld#{force_count}_f32(#{$current_predicate},(float*)&force[i]);\n"
      for i in 0...force_count do
        code += "svfloat32_t& f#{i} = f.v#{i};\n"
      end
    else
      if force_count <= 8 then
        offset = 2
      else
        offset = 1
      end
      for i in 0...(16/offset) do
        code += "svfloat32_t f#{i};\n"
      end
      for i in 0...(16/offset) do
        code += "f#{i} = svld1_f32(#{$current_predicate},(float*)&force[i+#{offset*i}]);\n"
      end
      code += "gather#{16/offset}("
      for i in 0...(16/offset) do
        code += "," if i > 0
        code += "f#{i}"
      end
      code += ");\n"
    end
    $varhash.each{ |v|
      iotype = v[1][0]
      if iotype == "ACCUM"
        name = v[0]
        type = v[1][1]
        case type
        when "F32"
          code += "svfloat32_t " + name + ";\n"
        when "F32vec"
          code += "svfloat32x3_t " + name + ";\n"
        else
          abort "error: unsupported type #{type} for ACCUM variable"
        end
      end
    }
    accum_init,ss,accum_finalize = split_accum_statement(@statements)
    accum_init.each{|s|
      line = s.convert_to_code(conversion_type) + "\n"
      epi_tab.each_with_index{ |v,i|
        orig = v[0]
        new  = "r" + v[1]
        line.gsub!(orig,new)
      }
      force_tab.each_with_index{ |v,i|
        orig = v[0]
        new  = "f" + v[1]
        line.gsub!(orig,new)
      }
      code += line
    }
    code += "#pragma loop loop_fission\n"
    code += "#pragma loop loop_fission_stripmining L1\n"
    code += "for(int j=0;j<nj;j++){\n"
    abort "pg at the top of j loop must be pg0" if $pg_count > 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "EPJ"
        name = v[0]
        type = v[1][1]
        fdpsname = v[1][2]
        case type
        when "F32vec" then
          code += "svfloat32x3_t " + name + ";\n"
          code += name + ".v0 = svdup_n_f32(epj[j]." + fdpsname + ".x);\n"
          code += name + ".v1 = svdup_n_f32(epj[j]." + fdpsname + ".y);\n"
          code += name + ".v2 = svdup_n_f32(epj[j]." + fdpsname + ".z);\n"
        when "F32" then
          code += "svfloat32_t " + name + " = svdup_n_f32(epj[j]." + fdpsname + ");\n"
        end
      end
    }
    ss.each{|s|
      line = s.convert_to_code(conversion_type)+"\n"
      if isStatement(s)
        epi_tab.each_with_index{ |v,i|
          orig = v[0]
          new  = "r" + v[1]
          line.gsub!(orig,new)
        }
        force_tab.each_with_index{ |v,i|
          orig = v[0]
          new  = "f" + v[1]
          line.gsub!(orig,new)
        }
      end
      code += line
    }
    code += "}\n" # j loop
    abort "pg at the end of j loop must be pg0" if $pg_count >  0
    accum_finalize.each{|s|
      line = s.convert_to_code(conversion_type) + "\n"
      epi_tab.each_with_index{ |v,i|
        orig = v[0]
        new  = "r" + v[1]
        line.gsub!(orig,new)
      }
      force_tab.each_with_index{ |v,i|
        orig = v[0]
        new  = "f" + v[1]
        line.gsub!(orig,new)
      }
      code += line
    }
    if force_count <= 4 then
      code += "svst#{force_count}_f32(#{$current_predicate},(float*)&force[i],f);\n"
    else
      if force_count <= 8 then
        offset = 2
      else
        offset = 1
      end
      #abort "number of FORCE class member more than 4 is not supported"
      code += "scatter#{16/offset}("
      for i in 0...(16/offset) do
        code += "," if i > 0
        code += "f#{i}"
      end
      code += ");\n"
      for i in 0...(16/offset) do
        code += "svst1_f32(svptrue_b32(),(float*)&force[i+#{offset*i}], f#{i});\n"
      end
    end
    code += "}\n" # i loop
    code += "}\n" # fuction body
    code
  end

  def generate_load_store_code_i_vars(io,fvars,conversion_type,h=$varhash)
    lcode = ""
    scode = ""
    # count bytes of EPI and FORCE member variable
    tot = 0
    is_uniform = true
    prev_type = nil
    max_byte_size = 0
    h.each{ |v|
      iotype = v[1][0]
      type   = v[1][1]
      modifier = v[1][3]
      if iotype == io
        prev_type = type.delete("vec") if tot == 0
        byte = byte_count(type)
        tot += byte_count(type) if modifier == nil
        byte = byte / 3 if type =~ /vec/
        max_byte_size = byte if byte > max_byte_size
        is_uniform = false if prev_type != nil && type.delete("vec") != prev_type
      end
    }
    warn "max byte size of #{io} member is #{max_byte_size}, #{tot} byte in total"
    # declare and load EPI, FORCE and local variables
    if is_uniform && (tot/byte_count(prev_type)) <= 4
      nelem = tot/byte_count(prev_type)
      vname = "__fkg_tmp_pos"   if io == "EPI"
      val  = "epi"   if io == "EPI"
      vname = "__fkg_tmp_force" if io == "FORCE"
      val  = "force" if io == "FORCE"
      tmp_type = "#{prev_type}vec#{nelem}"
      lcode += "#{get_declare_type(tmp_type,conversion_type)} #{vname} = svld#{nelem}_#{get_type_suffix(prev_type)}(#{$current_predicate},#{get_pointer_cast(prev_type)}&#{val}[i]);\n"

      count = 0
      h.each_with_index{ |v,n|
        name     = v[0]
        iotype   = v[1][0]
        type     = v[1][1]
        fdpsname = v[1][2]
        modifier = v[1][3]
        if iotype == io
          if fvars.index(name)
            case modifier
            when "local"
              lcode += "#{get_declare_type(type,conversion_type)} #{name};\n"
              if type =~ /vec/
                lcode += "#{name} = svld3_#{get_type_suffix(type)}(#{$current_predicate},#{get_pointer_cast(type)}(#{name}_tmp+i));\n"
              else
                lcode += "#{name} = svld1_#{get_type_suffix(type)}(#{$current_predicate},#{get_pointer_cast(type)}(#{name}_tmp+i));\n"
              end
            else
              if type =~ /vec/
                lcode += "#{get_declare_type(type,conversion_type)} #{name};\n"
                lcode += "#{name}.v0 = #{vname}.v#{count};\n"
                lcode += "#{name}.v1 = #{vname}.v#{count+1};\n"
                lcode += "#{name}.v2 = #{vname}.v#{count+2};\n"
                scode += "#{vname}.v#{count  } = #{name}.v0;\n" if io == "FORCE"
                scode += "#{vname}.v#{count+1} = #{name}.v1;\n" if io == "FORCE"
                scode += "#{vname}.v#{count+2} = #{name}.v2;\n" if io == "FORCE"
              else
                lcode += "#{get_declare_type(type,conversion_type)} #{name} = #{vname}.v#{count};\n"
                scode += "#{vname}.v#{count} = #{name};\n" if io == "FORCE"
              end
            end
          end
          if modifier != "local"
            count += 1
            count += 2 if type=~ /vec/
          end
        end
      }
      #scode += "svst#{nelem}_#{get_type_suffix(prev_type)}(#{$current_predicate},#{get_pointer_cast(prev_type)}&force[i],#{vname});\n"
      if io == "FORCE"
        for i in 0...nelem
          scode += "#{get_declare_scalar_type(prev_type)} __fkg_tmp_force_store#{i}[#{512/byte_count(prev_type)}];\n"
          scode += "svst1_#{get_type_suffix(prev_type)}(#{$current_predicate},__fkg_tmp_force_store#{i},__fkg_tmp_force.v#{i});\n"
        end
        scode += "for(int j=0;j<std::min(16,ni-i);j++){\n"
        count = 0
        h.each{ |v|
          name     = v[0]
          iotype   = v[1][0]
          type     = v[1][1]
          fdpsname = v[1][2]
          modifier = v[1][3]
          if iotype == "FORCE"
            if type =~ /vec/
              scode += "force[i+j].#{fdpsname}.x = __fkg_tmp_force_store#{count}[j];\n"
              count += 1
              scode += "force[i+j].#{fdpsname}.y = __fkg_tmp_force_store#{count}[j];\n"
              count += 1
              scode += "force[i+j].#{fdpsname}.z = __fkg_tmp_force_store#{count}[j];\n"
              count += 1
            else
              scode += "force[i+j].#{fdpsname} = __fkg_tmp_force_store#{count}[j];\n"
              count += 1
            end
          end
        }
        scode += "}\n"
      end
    else
      count = 0
      h.each_with_index{ |v,n|
      name     = v[0]
        iotype   = v[1][0]
        type     = v[1][1]
        fdpsname = v[1][2]
        modifier = v[1][3]
        if iotype == io
          if fvars.index(name)
            lcode += "#{get_declare_type(type,conversion_type)} #{name};\n"
            case modifier
            when "local"
              if type =~ /vec/
                lcode += "#{name} = svld3_#{get_type_suffix(type)}(#{$current_predicate},#{get_pointer_cast(type)}(#{name}_tmp+i));\n"
              else
                lcode += "#{name} = svld1_#{get_type_suffix(type)}(#{$current_predicate},#{get_pointer_cast(type)}(#{name}_tmp+i));\n"
              end
            else
              byte = byte_count(type)
              byte = byte / 3 if type =~ /vec/
              interval = ((tot+max_byte_size-1)/max_byte_size)*max_byte_size / byte
              offset = count / byte
              val = "epi"   if io == "EPI"
              val = "force" if io == "FORCE"
              if type =~ /vec/
                ["v0","v1","v2"].each_with_index{ |dim,dcount|
                  index_name = "index_#{name}#{dim}"
                  lcode += "uint#{byte}_t #{index_name}[#{512/byte}] = {"
                  for i in 0...(512/byte)
                    lcode += "," if i > 0
                    lcode += "#{offset + dcount + i*interval}"
                  end
                  lcode += "};\n"
                  vindex_name = "vindex_#{name}#{dim}"
                  index_type = "U#{byte}"
                  lcode += "svuint#{byte}_t #{vindex_name} = svld1_#{get_type_suffix(index_type)}(svptrue_b32(),#{index_name});\n"
                  lcode += "#{name}.#{dim} = svld1_gather_u32index_#{get_type_suffix(type)}(#{$current_predicate},#{get_pointer_cast(type)}&#{val}[i],#{vindex_name});\n"
                  scode += "svst1_scatter_u#{byte}index_#{get_type_suffix(type)}(#{$current_predicate},#{get_pointer_cast(type)}&#{val}[i],#{vindex_name},#{name}.#{dim});\n"
                }
              else
                index_name = "index_#{name}"
                index_type = "U#{byte}"
                lcode += "uint#{byte}_t #{index_name}[#{512/byte}] = {"
                for i in 0...(512/byte)
                  lcode += "," if i > 0
                  lcode += "#{offset + i*interval}"
                end
                lcode += "};\n"
                vindex_name = "vindex_#{name}"
                lcode += "svuint#{byte}_t #{vindex_name} = svld1_#{get_type_suffix(index_type)}(svptrue_b32(),#{index_name});\n"
                lcode += "#{name} = svld1_gather_u#{byte}index_#{get_type_suffix(type)}(#{$current_predicate},#{get_pointer_cast(type)}&#{val}[i],#{vindex_name});\n"
                scode += "svst1_scatter_u#{byte}index_#{get_type_suffix(type)}(#{$current_predicate},#{get_pointer_cast(type)}&#{val}[i],#{vindex_name},#{name});\n"
              end # isVector
            end # local or not
          end # if fvars.index(name)
          count += byte_count(type)
        end # if iotype == io
      }
    end # if !uniform
    [lcode,scode]
  end


  def j_loop_a64fx_mp(statements,j,fvars,conversion_type,h=$varhash)
    code = ""
    h.each{ |v|
      name     = v[0]
      iotype   = v[1][0]
      type     = v[1][1]
      fdpsname = v[1][2]
      modifier = v[1][3]
      if iotype == "EPJ" && fvars.index(name)
        code += "#{get_declare_type(type,conversion_type)} #{name};\n"
        case modifier
        when "local"
          if type =~ /vec/
            code += "#{name}.v0 = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}].x);\n"
            code += "#{name}.v1 = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}].y);\n"
            code += "#{name}.v2 = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}].z);\n"
          else
            code += "#{name} = svdup_#{get_type_suffix(type)}(#{name}_tmp[j]);\n"
          end
        else
          if type =~ /vec/
            code += "#{name}.v0 = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname}.x);\n"
            code += "#{name}.v1 = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname}.y);\n"
            code += "#{name}.v2 = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname}.z);\n"
          else
            code += "#{name} = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname});\n"
          end
        end
      end
    }
    statements.each{ |s|
      code += s.convert_to_code(conversion_type) + "\n" if !isStatement(s) || h[get_name(s)][3] != "local"
    }
    code
  end

  def j_loop_fissioned_a64fx_mp(statements,j,fvars,conversion_type,h=$varhash)
    loop_fission_vars = find_loop_fission_load_store_vars(statements)
    fission_count = 0
    code = ""

    tmpvars = []
    loop_fission_vars.each{ |vs|
      vs[0].each{ |v|
        tmpvars += [v]
      }
    }
    tmpvars.uniq.each{ |v|
      iotype = $varhash[v][0]
      type = $varhash[v][1]
      if type =~ /vec/
        type = type.delete("vec")
        code += "PS::#{type} #{v}_tmp_x[#{512/byte_count(type)}*#{$strip_mining}];\n"
        code += "PS::#{type} #{v}_tmp_y[#{512/byte_count(type)}*#{$strip_mining}];\n"
        code += "PS::#{type} #{v}_tmp_z[#{512/byte_count(type)}*#{$strip_mining}];\n"
      else
        code += "PS::#{type} #{v}_tmp[#{512/byte_count(type)}*#{$strip_mining}];\n"
      end
    }

    software_pipelining(statements,loop_fission_vars,$swpl_stage) if $swpl_stage != nil

    code += "for(int jj=0;jj<#{$strip_mining};jj++){\n"
    code += "const int jjj = j + jj;\n"
    loop_fission_vars[fission_count][1].each{ |v|
      name     = v
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      modifier = h[v][3]
      if iotype == "EPJ" && fvars.index(name)
        code += "#{get_declare_type(type,conversion_type)} #{name};\n"
        case modifier
        when "local"
          if type =~ /vec/
            code += "#{name}.v0 = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}].x);\n"
            code += "#{name}.v1 = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}].y);\n"
            code += "#{name}.v2 = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}].z);\n"
          else
            code += "#{name} = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}]);\n"
          end
        else
          if type =~ /vec/
            code += "#{name}.v0 = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname}.x);\n"
            code += "#{name}.v1 = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname}.y);\n"
            code += "#{name}.v2 = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname}.z);\n"
          else
            code += "#{name} = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname});\n"
          end
        end
      end
    }
    fission_count += 1

    statements.each{ |s|
      if s.class == Pragma
        loop_fission_vars[fission_count][0].each{ |v|
          type = h[v][1]
          if type =~ /vec/
            code += "svst1_#{get_type_suffix(type)}(#{$current_predicate},#{v}_tmp_x+16*jj,#{v}.v0);\n"
            code += "svst1_#{get_type_suffix(type)}(#{$current_predicate},#{v}_tmp_y+16*jj,#{v}.v1);\n"
            code += "svst1_#{get_type_suffix(type)}(#{$current_predicate},#{v}_tmp_z+16*jj,#{v}.v2);\n"
          else
            code += "svst1_#{get_type_suffix(type)}(#{$current_predicate},#{v}_tmp+16*jj,#{v});\n"
          end
        }
        code += "}\n"
        code += "for(int jj=0;jj<#{$strip_mining};jj++){\n"
        code += "const int jjj = j + jj;\n"
        loop_fission_vars[fission_count][1].each{ |v|
          #p v
          iotype = h[v][0]
          type   = h[v][1]
          if iotype == "declared"
            if type =~ /vec/
              #code += "#{get_daclare_type(type)} #{v};\n"
              code += "#{v}.v0 = svld1_#{get_type_suffix(type)}(#{$current_predicate},#{v}_tmp_x+16*jj);\n"
              code += "#{v}.v1 = svld1_#{get_type_suffix(type)}(#{$current_predicate},#{v}_tmp_y+16*jj);\n"
              code += "#{v}.v2 = svld1_#{get_type_suffix(type)}(#{$current_predicate},#{v}_tmp_z+16*jj);\n"
            else
              #code += "#{get_declare_type(type)} #{v};\n"
              code += "#{v} = svld1_#{get_type_suffix(type)}(#{$current_predicate},#{v}_tmp+16*jj);\n"
            end
          elsif iotype == "EPJ"
            name     = v
            iotype   = h[v][0]
            type     = h[v][1]
            fdpsname = h[v][2]
            modifier = h[v][3]
            if iotype == "EPJ" && fvars.index(name)
              code += "#{get_declare_type(type,conversion_type)} #{name};\n"
              case modifier
              when "local"
                if type =~ /vec/
                  code += "#{name}.v0 = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}].x);\n"
                  code += "#{name}.v1 = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}].y);\n"
                  code += "#{name}.v2 = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}].z);\n"
                else
                  code += "#{name} = svdup_#{get_type_suffix(type)}(#{name}_tmp[#{j}]);\n"
                end
              else
                if type =~ /vec/
                  code += "#{name}.v0 = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname}.x);\n"
                  code += "#{name}.v1 = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname}.y);\n"
                  code += "#{name}.v2 = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname}.z);\n"
                else
                  code += "#{name} = svdup_#{get_type_suffix(type)}(epj[#{j}].#{fdpsname});\n"
                end
              end
            end
          end
        }
        fission_count += 1
      else
        code += s.convert_to_code(conversion_type) + "\n" if !isStatement(s) || h[get_name(s)][3] != "local"
      end
    }
    code += "}\n"
    code
  end

  def kernel_body_a64fx_mp(conversion_type,h=$varhash)
    code = ""
    force_related_vars = generate_force_related_map(@statements)
    accum_init,ss,accum_finalize = split_accum_statement(@statements)
    code += "for(int i=0;i<((ni+15)/16)*16;i+=16){\n"
    code += "#{$current_predicate} = svwhilelt_b32_s32(i,ni);\n"
    load_code,store_code = generate_load_store_code_i_vars("EPI",force_related_vars,conversion_type)
    code += load_code
    load_code,store_code = generate_load_store_code_i_vars("FORCE",force_related_vars,conversion_type)
    code += load_code
    accum_init.each{ |s|
      code += s.convert_to_code(conversion_type)
    }
    code += "int j=0;\n"
    #code += "#pragma loop loop_fission\n"
    #code += "#pragma loop loop_fission_stripmining L1\n"
    if $strip_mining != nil
      code += "int j_end = (nj/#{$strip_mining})*(#{$strip_mining});\n"
      code += "for(;j<j_end;j+=#{$strip_mining}){\n"
      code += j_loop_fissioned_a64fx_mp(ss,"jjj",force_related_vars,conversion_type,h)
      code += "}\n"
    end
    code += "for(;j<nj;j++){\n"
    code += j_loop_a64fx_mp(ss,"j",force_related_vars,conversion_type,h)
    code += "}\n"
    accum_finalize.each{ |s|
      code += s.convert_to_code(conversion_type)
    }
    code += store_code
    code += "}\n" # i loop
    code += "}\n" # end of operator()
    code
  end
end

class Funcdeclaration
end

class Statement
  def declare_variable_a64fx(h = $varhash)
    declpart=""
    name = get_name(@name)
    tail = get_tail(@name)
    if h[name][0] == nil
      if tail == nil
        case @type
        when "F64" then
          declpart = "svfloat64_t #{name};\n"
        when "F32" then
          declpart = "svfloat32_t #{name};\n"
        when "S64" then
          declpart = "svint64_t #{name};\n"
        when "S32" then
          declpart = "svint32_t #{name};\n"
        when "U64" then
          declpart = "svuint64_t #{name};\n"
        when "U32" then
          declpart = "svuint32_t #{name};\n"
        else
          abort "error: unsupported scalar type of #{@type} for A64FX in declaration"
        end
      else
        case @type
        when "F64" then 
          declpart = "svfloat64x3_t #{name};\n"
        when "F32" then 
          declpart = "svfloat32x3_t #{name};\n"
        else
          abort "error: unsupported vector type of #{@type} for A64FX in declaration"
        end
      end
      h[name][0] = "declared"
    end
    declpart
  end
end

class FuncCall
  def convert_to_code_a64fx(conversion_type)
    retval = @name
    retval += "(#{$current_predicate},"
    if @ops.length > 0
      @ops.each_with_index{ |op,i| 
        retval += "," if i > 0
        retval += op.convert_to_code(conversion_type)
      }
    end
    retval += ")"
    retval
  end
end

class IfElseState
  def convert_to_code_a64fx(conversion_type)
    ret = ""
    if $if_count == nil
      $if_count = -1
      $if_exp = [[]]
    end
    case @operator
    when :if
      $pg_count += 1
      $if_count += 1
      $if_exp[$if_count] = [@expression]
      pgtmp = "pg#{$pg_count}"
      #ret += "svbool_t " + pgtmp + ";\n" if $pg_count > $max_pg_count
      $max_pg_count = $pg_count if $pg_count > $max_pg_count
      ret += pgtmp + "=" + @expression.convert_to_code(conversion_type) + ";"
      $current_predicate = pgtmp
    when :elsif
      pgtmp  = "pg#{$pg_count}" 
      pgprev = "pg#{$pg_count-1}"
      # svand(pg,A,svand(pg,B,C)
      tail = $if_exp[$if_count].length - 1
      line = "svnot_b_z(#{pgprev},"
      $if_exp[$if_count].each_with_index{ |exp,i|
        line += "svand_b_z(#{pgprev}," if i < tail
        line += exp.convert_to_code(conversion_type)
        line += "," if i < tail
        line += ")" if i == tail
      }
      ret += "#{pgtmp} = svand_b_z(#{pgprev},#{line}," + @expression.convert_to_code(conversion_type) + ");"
      $if_exp[$if_count] += [@expression]
    when :else
      pgtmp = "pg#{$pg_count}"
      pgprev = "pg#{$pg_count-1}"
      tail = $if_exp[$if_count].length - 1
      line = "svnot_b_z(#{pgprev},"
      $if_exp[$if_count].each_with_index{ |exp,i|
        line += "svand_b_z(#{pgprev}," if i < tail
        line += exp.convert_to_code(conversion_type)
        line += "," if i < tail
        line += ")" if i == tail
      }
      ret += "#{pgtmp} = #{line};"
    when :endif
      $if_count -= 1
      $pg_count -= 1
      $current_predicate = "pg#{$pg_count}"
    else
      abort "error: undefined if operator #{@operator}"
    end
    ret
  end
end

class Expression
  def convert_to_code_a64fx(conversion_type)
    if @operator != :array && @operator != :func
      case self.get_type
      when "F64"
        type = "f64"
      when "F32"
        type = "f32"
      when "S64"
        type = "s64"
      when "S32"
        type = "s32"
      when "U64"
        type = "u64"
      when "U32"
        type = "u32"
      when "B64"
        type = "b64"
      when "B32"
        type = "b32"
      when "F64vec"
        type = "f64x3"
      when "F32vec"
        type = "f32x3"
      else
        #abort "error: unsupported type #{@type} and operator #{@operator}"
      end
    end

    case @operator
    when :uminus then
      retval="svneg_#{type}_z(" + $current_predicate + "," + @lop.convert_to_code(conversion_type) + ")"
    when :plus  then
      retval="svadd_#{type}_z(" + $current_predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :minus then
      retval="svsub_#{type}_z(" + $current_predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :mult  then
      retval="svmul_#{type}_z(" + $current_predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :div   then
      retval="svdiv_#{type}_z(" + $current_predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :lt    then
      retval="svcmplt_#{type}(" + $current_predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :le    then
      retval="svcmple_#{type}(" + $current_predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :gt    then
      retval="svcmpgt_#{type}(" + $current_predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :ge    then
      retval="svcmpge_#{type}(" + $current_predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :eq    then
      retval="svcmpeq_#{type}(" + $current_predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :neq   then
      retval="svcmpne_#{type}(" + $current_predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :dot   then
      if @lop =~ /\d+/ && @rop =~ /\d+(f|h)?/
        case rop
        when /f/
          retval="svdup_n_f32(#{@lop}.#{@rop})"
        when /h/
          retval="svdup_n_f16(#{@lop}.#{@rop.gsub!("h","")})"
        else
          retval="svdup_n_f64(#{@lop}.#{@rop})"
        end
      elsif @rop == "x" || @rop == "y" || @rop == "z" || @rop == "w"
        retval=@lop.convert_to_code(conversion_type)+"."
        retval += ["v0","v1","v2","v3"][["x","y","z","w"].index(@rop)]
        #@rop.convert_to_code(conversion_type)
      elsif @rop == "v0" || @rop == "v1" || @rop == "v2" || @rop == "v3"
        retval=@lop.convert_to_code(conversion_type)+"." + @rop
        #@rop.convert_to_code(conversion_type)
      else
        #retval = "#{@lop.convert_to_code(conversion_type)}.#{@rop.convert_to_code(conversion_type)}"
        retval = "svdup_n_#{type}(#{@lop.convert_to_code(conversion_type)}.#{@rop.convert_to_code(conversion_type)})"
      end
    when :func  then
      retval=@lop+"("+@rop.convert_to_code(conversion_type)+")"
    when :array
      retval = @lop+"["+@rop.convert_to_code(conversion_type)+"]"
    else
      abort "error: unsupported operator #{@operator}"
    end
    retval
  end
end

class String
  def convert_to_code_a64fx(conversion_type,h=$varhash)
    name = get_name(self)
    #abort "error: undefined reference to #{name} in convert_to_code_a64fx of String"
    s = self
    if h[name] != nil
      iotype = h[name][0]
      if iotype == "MEMBER"
        type = self.get_type
        case type
        when "F64"
          s = "svdup_n_f64("+self+")"
        when "F32"
          s = "svdup_n_f32("+self+")"
        when "S64"
          s = "svdup_n_s64("+self+")"
        when "S32"
          s = "svdup_n_s32("+self+")"
        when "U64"
          s = "svdup_n_u64("+self+")"
        when "U32"
          s = "svdup_n_u32("+self+")"
        end
      end
    end
    s
  end
end


def count_class_member(io,h=$varhash)
  tot = 0
  is_uniform = true
  prev_type = nil
  max_byte_size = 0
  h.each{ |v|
    iotype = v[1][0]
    type   = v[1][1]
    modifier = v[1][3]
    if iotype == io
      prev_type = type.delete("vec") if tot == 0
      byte = byte_count(type)
      tot += byte_count(type) if modifier == nil
      byte = byte / 3 if type =~ /vec/
      max_byte_size = byte if byte > max_byte_size
      is_uniform = false if prev_type != nil && type.delete("vec") != prev_type
    end
  }
  warn "size of #{io} member is #{max_byte_size}, # of elemennts are #{tot/max_byte_size}" if is_uniform
  [tot,max_byte_size,is_uniform]
end

def aos2soa_a64fx(fvars,h=$varhash)
  ret = []
  # count bytes of EPI and FORCE member variable
  epi_tot, epi_max_byte_size, is_epi_uniform = count_class_member("EPI")
  force_tot, force_max_byte_size, is_force_uniform = count_class_member("FORCE")

  tot = [epi_tot,force_tot]
  max_byte_size = [epi_max_byte_size,force_max_byte_size]
  is_uniform = [is_epi_uniform,is_force_uniform]
  nelems = [epi_tot/epi_max_byte_size, force_tot/force_max_byte_size]
  iotypes = ["EPI","FORCE"]

  # declare and load EPI, FORCE and local variables
  ["EPI","FORCE"].each{ |io|
    nelem  = nelems[iotypes.index(io)]
    offset = (tot[iotypes.index(io)] + max_byte_size[iotypes.index(io)] - 1) / max_byte_size[iotypes.index(io)]
    if is_uniform[iotypes.index(io)] && nelem <= 4
      #if false # structure pattern is not supported yet
      vname = ["__fkg_tmp_epi","__fkg_tmp_force"][iotypes.index(io)]
      type = "F#{max_byte_size[iotypes.index(io)]}vec#{nelem}"
      ret += [Declaration.new([type,vname])]
      ret += [StructureLoad.new([vname,PointerOf.new([type,Expression.new([:array,get_iotype_array(io),"i",type])]),nelem,type])]

      fmem = []
      h.each{ |v|
        if v[1][0] == io
          fmem += [v[0]] 
          ret += [Declaration.new([v[1][1],v[0]])]
        end
      }
      count = 0
      fvars.each{ |v|
        if h[v][0] == io
          if h[v][1] =~ /vec/
            ret += [Statement.new([Expression.new([:dot,v,"x"]),Expression.new([:dot,vname,"v#{count+0}",type]),type])]
            ret += [Statement.new([Expression.new([:dot,v,"y"]),Expression.new([:dot,vname,"v#{count+1}",type]),type])]
            ret += [Statement.new([Expression.new([:dot,v,"z"]),Expression.new([:dot,vname,"v#{count+2}",type]),type])]
            count += 3
          else
            ret += [Statement.new([v,Expression.new([:dot,vname,"v#{count}",type]),type])]
            count += 1
          end
        end
      }
    else # is not uniform
      fvars.each{ |v|
        name = v
        iotype   = h[v][0]
        type     = h[v][1]
        fdpsname = h[v][2]
        modifier = h[v][3]
        if iotype == io
          ret += [Declaration.new([type,name])]
          if modifier == "local"
            if type =~ /vec/
              tmp_type = type.delete("vec")
              ret += [LoadState.new([Expression.new([:dot,name,"x"]),PointerOf.new([type,Expression.new([:array,"#{name}_tmp_x","i",tmp_type])]),tmp_type])]
              ret += [LoadState.new([Expression.new([:dot,name,"y"]),PointerOf.new([type,Expression.new([:array,"#{name}_tmp_y","i",tmp_type])]),tmp_type])]
              ret += [LoadState.new([Expression.new([:dot,name,"z"]),PointerOf.new([type,Expression.new([:array,"#{name}_tmp_z","i",tmp_type])]),tmp_type])]
            else
              ret += [LoadState.new([name,PointerOf.new([type,Expression.new([:array,"#{name}_tmp","i",type])]),type])]
            end
          else
            tmp_type = type.delete("vec")
            offset_tmp = offset * max_byte_size[iotypes.index(io)] / byte_count(tmp_type)
            if type =~ /vec/
              ret += [GatherLoad.new([Expression.new([:dot,name,"x",tmp_type]),PointerOf.new([type,Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type]),"x",tmp_type])]),"#{offset_tmp}",tmp_type])]
              ret += [GatherLoad.new([Expression.new([:dot,name,"y",tmp_type]),PointerOf.new([type,Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type]),"y",tmp_type])]),"#{offset_tmp}",tmp_type])]
              ret += [GatherLoad.new([Expression.new([:dot,name,"z",tmp_type]),PointerOf.new([type,Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type]),"z",tmp_type])]),"#{offset_tmp}",tmp_type])]
            else
              ret += [GatherLoad.new([name,PointerOf.new([type,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type])]),"#{offset_tmp}",type])]
            end
          end
        end
      }
    end
  }
  ret
end

def soa2aos_a64fx(fvars,h=$varhash)
  ret = []
  # count bytes of EPI and FORCE member variable
  tot, max_byte_size, is_uniform = count_class_member("FORCE")

  nelem = tot/max_byte_size
  io = "FORCE"

  # declare and load EPI, FORCE and local variables
  offset = (tot + max_byte_size - 1) / max_byte_size
  if false #is_uniform && nelem <= 4 # structure store does not work with FCCpx
    vname = "__fkg_tmp_force"
    type = "F#{max_byte_size}"
    ret += [StructureStore.new([PointerOf.new([type,Expression.new([:array,get_iotype_array(io),"i",type])]),vname,nelem,type])]
  else # is not uniform
    h.each{ |v|
      name = v[0]
      iotype   = v[1][0]
      type     = v[1][1]
      fdpsname = v[1][2]
      modifier = v[1][3]
      if iotype == io
        tmp_type = type.delete("vec")
        tmp_offset = offset * max_byte_size / byte_count(tmp_type)
        if type =~ /vec/
          ret += [ScatterStore.new([PointerOf.new([type,Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type]),"x",tmp_type])]),Expression.new([:dot,name,"x",tmp_type]),"#{tmp_offset}",tmp_type])]
          ret += [ScatterStore.new([PointerOf.new([type,Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type]),"y",tmp_type])]),Expression.new([:dot,name,"y",tmp_type]),"#{tmp_offset}",tmp_type])]
          ret += [ScatterStore.new([PointerOf.new([type,Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type]),"z",tmp_type])]),Expression.new([:dot,name,"z",tmp_type]),"#{tmp_offset}",tmp_type])]
        else
          ret += [ScatterStore.new([PointerOf.new([type,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type])]),name,"#{tmp_offset}",type])]
        end
      end
    }
  end

  ret
end
