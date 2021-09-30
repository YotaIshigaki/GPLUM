#
# convert_f90_struct_to_c.rb
#
def print_header
end
def print_trailer
end


$type_conversion_table =<<-EOF
integer(kind=c_long_long) long long
real(kind=c_double)   double
real(kind=c_float)    float
type(fdps_f64vec) Cvec_Float64
type(fdps_f32vec) Cvec_Float32
type(fdps_f64mat) Cmat_Float64
type(fdps_f32mat) Cmat_Float32
EOF


$type_conversion_hash = $type_conversion_table.split("\n").map{|s|
  a=s.split
  [a[0], a[1..(a.length-1)].join(" ") ]}.to_h


def print_c_member(s)
  a = s.split
  return if a.length==0
  ss = "     "
  if a[0][0] != "!"
    # member variable line
    varname = a[2]
    typename = $type_conversion_hash[a[0]]
    raise "unparsable line: #{s}"    if a[1] != "::" || typename == nil
    ss += " "+typename + " "+ varname + ";"
    a = a[3..(a.length)]
  end
  if a.length >0
    a[0]= "//"+a[0]
    ss += " "+a.join(" ")
  end
  print ss,"\n"
end
        
      
      
    
print_header

while s=gets
  if s =~ /!\*\*\*\* PS::(\S+)/
    print "// "+ s
    struct_name = $1
    s = gets
    lines=[]
    instruct=true
    while instruct
      s=gets
      print "// "+ s
      if s=~/^(\s*)end type/
        instruct =false
      else
        lines.push s.chomp
      end
    end
    cname = struct_name
    struct_name= struct_name.downcase
    print "typedef   struct #{struct_name}{\n"
    lines.each{|s| print_c_member(s)}
    print "}#{cname},*P#{cname};\n"
  end
end

print_trailer

