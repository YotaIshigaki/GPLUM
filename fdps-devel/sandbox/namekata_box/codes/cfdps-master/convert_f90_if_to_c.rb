#
# convert_f90_if_to_c.rb
#

$type_conversion_table =<<-EOF
void          void
int           int
bool          bool
PS::S64       long long
PS::S32       int
PS::F32       float
PS::F64       double
double        double
PS::TimeProfile TimeProfile
PS::F32vec Cvec_Float32
char          char
size_t        size_t
full_particle  Full_particle
float         float
EOF


$type_conversion_hash = $type_conversion_table.split("\n").map{|s|
  a= s.split
  [a[0], a[1..(a.length-1)].join(" ")]
}.to_h

$functions_to_skip=<<-EOF
get_psys_cptr
EOF
def convert_functype(s)
  ss = $type_conversion_hash[s]
  if ss == nil
    raise "unknown functype #{s}"
  end
  ss
end

def preprocess_single_line_def(s)
  a=[]
  if s =~/^([a-z]\w*) fdps_(\S*)\((.*)\)/
    a.push $1
    a.push $2
    a += $3.split(",")
  end
  a
end

def preprocess_funarg(a)
  b=[]
  s= a.shift
#  print "Funarg[0]=",s,"\n"
  while s
    if s =~/([a-z]\w*) \(\*/
#      print "Funarg complex type"
      ss = s
      s= a.shift
      while s !~ /\)/ 
        ss += s
        s= a.shift
      end
      ss += s
      s=ss
    end
    b.push s.chomp
    s = a.shift
  end
  b
end
def preprocess_multi_line_def(a)
  b=[]
  s = a.shift
  if s =~/^([a-z]\w*) fdps_(\S*)\((.*)/
    b.push $1
    b.push $2
    a.unshift $3
  else
    raise "unparsable line: #{s}"
  end
  a=preprocess_funarg(a)
  b+a
  
end
def convert_to_prototype(a)
  s=a.join()
  s.gsub!(/{/,";")
  s.split("\n")
end

class String
  def tocap
    s=self
    s=self.capitalize if self[0] =~/[a-z]/
    s
  end
end      
    
def convert_function_pointers(s)
  if s.index("(*pfunc_comp)")
    s =~ /pfunc_comp\)\(const struct (\S+) (.*)\n(\s*) const struct (\S+)/
    args = [$1, $4].map{|s| s.tocap}
    s = "bool  (*pfunc_comp)(#{args[0]}*, #{args[1]}*)"
  elsif s.index("(*pfunc_ep")
    if  s =~/(pfunc_ep\S+)\)\(struct (\S+) .*\n\s* int ,\n\s* struct (\S+) .*\n\s*int.*\n\s* struct (\S+)/
    args=[$1,$2,$3,$4]
    elsif s =~/(pfunc_ep\S+)\)\(struct (\S+) .*\n\s* int ,\n\s* PS::(\S+) .*\n\s*int.*\n\s* struct (\S+)/
      args=[$1,$2,$3,$4]
    end
    funname = args.shift
    args = args.map{|s| s.tocap}
    s = <<EOF
         void (*#{funname})(#{args[0]}*, int,
                         #{args[1]}*, int,
                         #{args[2]}*)
EOF
  end
  s.chomp
end

def convert_to_c_arg_def(s)
  result=""
  s.gsub!(/\*(\s+)/," *")
  if s =~ /\n/
    result= convert_function_pointers(s)
  else
    a=s.split(")")[0].split(",")[0].split
    a.shift if a[0] == "const"
    varname = a.pop
    a.pop if a[a.length-1]== "const"
    typename = $type_conversion_hash[a.join(" ")]
    if typename  == nil
      a.shift if a[0] == "enum"
      typename = a.join(" ")
    end
#    if varname[0]=="*"
#     varname = varname[1,varname.length-1]
#      typename = "Pointer(#{typename})"
#    end
    result= "#{typename} #{varname}"
  end
  result
end

def convert_one_func(a)
  funtype = convert_functype(a[0])
  funname = a[1]
  if $functions_to_skip.index(funname)
    print "// function #{funname} skipped --- should be written manually\n"
    return
  end
  funargs = a[2..(a.length-1)].map{|s|  convert_to_c_arg_def(s)}
  print "#{funtype} fdps_#{funname}("
  print funargs[0] if funargs.length >0
  funargs.shift
  funargs.each{|s|
    print ",\n           ", s
  }
  print ");\n"
  
end

def print_header
end
def print_trailer
end
print_header
while s=gets
  if s =~  /^([a-z]\w*) fdps_(\S*)\(/
#    p s
    a=[]
    a.push s
    while ! ( a[a.length-1] =~ /\{/ )
      a.push gets
    end
#    print "func body:\n"
    #    print  a.join
    print "\n"
    a.each{|s|print "// "+s}
    print "\n"
      
    if a.length==1
      a =preprocess_single_line_def(a[0])
    else
      a =preprocess_multi_line_def(a)
    end
    convert_one_func(a)
#    print convert_to_prototype(a).join("\n")    
  end
end
print_trailer
