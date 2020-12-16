require "kernelparser.rb"

class Kernelprogram
  attr_accessor :nsdeclarations, :iodeclarations, :statements
  def initialize(x)
    @nsdeclarations, @iodeclarations, @statements =*x
  end
end
class Nsdeclaration
  attr_accessor :type,:name
  def initialize(x)
    @type,@name=x
  end
end
class Iodeclaration
  attr_accessor :iotype,:type,:name,:fdpsname
  def initialize(x)
    @iotype,@type,@name,@fdpsname=x
  end
  
end
class Statement
  attr_accessor :type, :name, :expression
  def initialize(x)
    @type, @name, @expression=x
  end
  def get_type
    @type = @name.get_type if @name.get_type
    @type = @expression.get_type if @type==nil
    @type
  end
  def add_to_varhash
    $varhash[@name] = [nil, @type, nil] if $varhash[@name] ==nil
  end
  def convert_to_code(conversion_type)
#    p self
#    p @name
#    p $varhash[@name]
    declpart=""
    declpart= @type+" "if $varhash[@name][0]==nil 
    declpart+@name.convert_to_code(conversion_type) +
      " = " + @expression.convert_to_code(conversion_type)
  end
end

class Expression
  attr_accessor :operator, :lop, :rop,:type
  def initialize(x)
    @type=nil
    @operator, @lop, @rop=x
  end
end


class Kernelprogram
  def generate_hash(kerneltype)
    #    print "print kernel\n"
    #    p self
    $nshash=process_nsdecl(@nsdeclarations)
    p $nshash
    $varhash=process_iodecl(@iodeclarations)
    p $varhash
    @statements.each{|s|
      s.get_type
      s.add_to_varhash
    }
    p self
  end
  def print_statements(convertion_type)
    @statements.each{|s|  print s.convert_to_code(convertion_type)+"\n"}
  end
end
class String
  def get_type
    propaties=$varhash[self]
    if propaties
      propaties[1]
    else
      nil
    end
  end
  def convert_to_code(conversion_type)
    s=self
#    print "convert to code #{s}\n"
#    p $varhash[self]
    if conversion_type == "noconversion"
      if $varhash[self]
        t=$varhash[self][0]
        s += "[i]" if t== "iin"|| t== "iout" 
        s += "[j]" if t== "jin"
      end
    end
#    print "result=", s, "\n"
    s 
  end
end
class Expression

  def derive_type(operator, lt, rt)
    type=nil
    #    print "delive type ", operator," ", lt," ", rt, "\n"
    if [:plus, :minus, :mult, :div, :uminus].index(operator)
      if operator == :mult && lt.index("vec")  && rt.index("vec")
        if lt == "f64vec" || rt == "f64vec"
          type="f64"
        else
          type="f32"
        end
      else
        type="i"
        type= "f" if lt.index("f") ||  rt.index("f")
        if lt.index("64") ||  rt.index("64")  
          type += "64"
        else
          type += "32"
        end
        type += "vec" if lt.index("vec") ||  rt.index("vec")  
      end
    else
      type=rt
    end
    #    print "final type=", type, "\n"
    type
  end     
  def  get_type
    if @type == nil
      @type = derive_type(@operator, @lop.get_type, @rop.get_type)
    end
    @type
  end

  def convert_to_code(conversion_type)
    retval=""
    if type=[:plus,:minus,:mult,:div].index(@operator)
      opcode="+-*/"[type]
      retval="("+@lop.convert_to_code(conversion_type)+opcode+
             @rop.convert_to_code(conversion_type)+")"
    elsif @operator== :func
      retval= @lop+"("+ @rop.convert_to_code(conversion_type)+")"
    end
    retval
  end
end

def process_nsdecl(n)
  a=[]
  n.each{|x| a += [x.type, x.name]}
  Hash[*a]
end
def process_iodecl(ios)
  a=[]
  ios.each{|x|
    a +=  [x.name, [x.iotype, x.type, x.fdpsname]]
  }
  p a   
  Hash[*a]
end

parser=KernelParser.new
s=<<EOF
epj jparticle
epi iparticle
fi  force
iin f64vec xi r
jin f64vec xj r
jin f64    mj m
iin f64 eps2
iout f64vec f
iout f64 phi
rij = xi-xj
r2 = rij*rij+eps2
rinv = rsqrtinv(r2)
mrinv = mj*rinv
phi = phi-mrinv
f = f-mrinv*rinv*rinv* rij
EOF

program=parser.parse(s)
p program
p program.class
program.generate_hash("noconversion")
program.print_statements("noconversion")
