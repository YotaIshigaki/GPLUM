# class definition
class Kernelprogram
  attr_accessor :iodeclarations, :functions, :statements
  def initialize(x)
    #@iodeclarations, @functions, @statements =*x
    @iodeclarations = x.shift
    if x[0][0].class == Function
      @functions = x.shift
    else
      @functions = []
    end
    @statements = x.shift
  end
end

class Loop
  attr_accessor :index,:loop_beg,:loop_end,:interval,:statements
  def initialize(x)
    @index,@loop_beg,@loop_end,@interval,@statements = x
  end
  def convert_to_code(conversion_type)
    ret = ""
    ret += "for("
    ret += "#{@index} = #{@loop_beg}" if @loop_beg  != nil
    ret += ";"
    ret += "#{@index} < #{@loop_end};"
    if @interval == 1
      ret += "++#{@index}"
    else
      ret += "#{@index} += #{@interval}"
    end
    ret += "){\n"
    statements.each{|s|
      ret += s.convert_to_code(conversion_type) + "\n"
    }
    ret += "}\n"
  end
end

class ILoop
  attr_accessor :statements
  def initialize(x)
    @statements = x
  end
  def convert_to_code(conversion_type)
    ret = ""
    case conversion_type
    when "reference"
      ret += "for(i=0;i<ni;i++){\n"
    when /A64FX/
      ret += "for(i=0;i<((ni+15)/16)*16;i+=16){\n"
      ret += "#{$current_predicate} = svwhilelt_b32(i,ni);\n"
    else
      abort "unsupporeted conversion_type (#{conversion_type}) in ILoop conversion"
    end
    statements.each{|s|
      ret += s.convert_to_code(conversion_type) + "\n"
    }
    ret += "}\n"
  end
end

class JLoop
  attr_accessor :statements
  def initialize(x)
    @statements = x
  end
  def convert_to_code(conversion_type)
    ret = ""
    case conversion_type
    when "reference"
      ret += "for(j=0;j<nj;j++){\n"
    when /A64FX/
    else
      abort "unsupporeted conversion_type (#{conversion_type}) in JLoop conversion"
    end
    statements.each{|s|
      ret += s.convert_to_code(conversion_type) + "\n"
    }
    ret += "}\n"
  end
end

class ConditionalBranch
  attr_accessor :conditions,:bodies
  def initialize(x)
    @conditions,@bodies = x
    abort "# of bodies is not correspond to # of conditions" if @conditions.length != @bodies.length
  end
  def push_body(x)
    @bodies.last.push(x)
  end
  def push_condition(x)
    @conditions.push(x)
    @bodies.push([])
    #p self
  end

  def get_related_variable
    ret = []
    @conditions.zip(@bodies){ |c,b|
      ret.push(["",c.expression.get_related_variable]) if c.expression != nil
      b.each{|s|
        if isStatement(s)
          exp = s.expression.get_related_variable
          exp = [get_name(s)] if exp == []
          ret.push([get_name(s),exp]) if exp != []
        elsif s.class == ConditionalBranch
          ret += s.get_related_variable
        end
      }
    }
    ret
  end

  def get_cond_related_variable
    ret = []
    @conditions.zip(@bodies){ |c,b|
      ret += c.expression.get_related_variable if c.expression != nil
      b.each{|s|
        ret += s.get_cond_related_variable if s.class == ConditionalBranch
      }
    }
    ret
  end

  def declare_temporal_var(h=$varhash)
    ret = []
    @bodies.each{ |b|
      b.each{ |s|
        ret += s.declare_temporal_var(h) if s.class == Statement || s.class == ConditionalBranch
      }
    }
    ret
  end

  def calc_max_predicate_count(count = 1)
    ret = count + 1
    @bodies.each{ |b|
      b.each{ |s|
        if s.class == ConditionalBranch
          tmp = s.calc_max_predicate_count(ret)
          ret = tmp if tmp > ret
        end
      }
    }
    ret
  end

  def convert_to_code(conversion_type)
    ret = ""
    conditions.zip(bodies){ |c,b|
      #p [c,b]
      ret += c.convert_to_code(conversion_type) + "\n"
      b.each{ |s|
        ret += s.convert_to_code(conversion_type) + "\n"
      }
    }
    ret += IfElseState.new([:endif,nil]).convert_to_code(conversion_type)
  end
end

class Iodeclaration
  attr_accessor :iotype,:type,:name,:fdpsname,:modifier
  def initialize(x)
    #p x
    @iotype,@type,@name,@fdpsname,@modifier=x
  end
end

class Declaration
  attr_accessor :type, :name
  def initialize(x)
    #p x
    @type,@name = x
  end
  def append(x)
    @name += [x]
  end
  def get_type(h = $varhash)
    @type = @name.get_type(h) if @type==nil
    @type
  end
  def convert_to_code(conversion_type)
    #p self
    "#{get_declare_type(@type,conversion_type)} #{@name.convert_to_code(conversion_type)};\n"
  end
end

class NonSimdDecl < Declaration
  def convert_to_code(conversion_type)
    #p self
    "#{get_declare_type(@type,"reference")} #{@name.convert_to_code("reference")};\n"
  end
end

class Function
  attr_accessor :decl, :statements, :retval
  def initialize(x)
    @decl, @statements, @retval = x
  end
end

class Funcdeclaration
  attr_accessor :name, :vars, :attr
  def initialize(x)
    #p x
    @name,@vars,@attr = x
  end
end

class ReturnState
  attr_accessor :ret
  def initialize(x)
    #p x
    @ret = x
  end

  def convert_to_code(conversion_type)
    retval = "return " + ret.convert_to_code(conversion_type) + ";\n"
    retval
  end
end

class Statement
  attr_accessor :name, :expression, :type
  def initialize(x)
    #p x
    @name, @expression, @type=x
  end
  def get_type(h = $varhash)
    @type = @expression.get_type(h) if @type==nil
    @type
  end
  def add_to_varhash(h=$varhash)
    #p $varhash
    name = get_name(@name)
    if h[name] == nil
      h[name] = [nil, @type, nil]
    end
  end

  def expand_inner_prod(orig)
    #p "expand_inner_prod:"
    exp = orig
    if !isLeaf(orig)
      lt = orig.lop.get_type
      if !(lt =~ /vec/)
        orig.lop = expand_inner_prod(orig.lop)
      end
      if orig.rop != nil
        rt = orig.rop.get_type
        if !(rt =~ /vec/)
          orig.rop = expand_inner_prod(orig.rop)
        end
        if lt =~ /vec/ && rt =~ /vec/
          if orig.operator == :mult
            type = lt.delete("vec")
            xx  = Expression.new([:mult,Expression.new([:dot,orig.lop,"x",type]),Expression.new([:dot,orig.rop,"x",type]),type])
            yy  = Expression.new([:mult,Expression.new([:dot,orig.lop,"y",type]),Expression.new([:dot,orig.rop,"y",type]),type])
            zz  = Expression.new([:mult,Expression.new([:dot,orig.lop,"z",type]),Expression.new([:dot,orig.rop,"z",type]),type])
            exp = Expression.new([:plus,xx,Expression.new([:plus,yy,zz,type]),type])
          else
            exp = orig
          end
        else
          exp = orig
        end
      end
    end
    exp
  end

  def expand_tree
    ret = expand_inner_prod(@expression)
    @expression = ret
  end

  def declare_variable(conversion_type,h=$varhash)
    declpart=""
    name = get_name(@name)
    tail = get_tail(@name)
    #p self
    #p h[name]
    if h[name][0]==nil
      case conversion_type
      when "reference" then
        if tail == nil
          declpart= "PS::" + @type + " " + name + ";\n"
        else
          declpart= "PS::" + @type + "vec " + name + ";\n"
        end
        h[name][0] = "declared"
      when /A64FX/
        declpart += self.declare_variable_a64fx
      end
    end
    declpart
  end

  def replace_name(orig,replaced)
    abort "name must be Strinng class to replace" if orig.class != String || replaced.class != String
    if @name.class == String && @name == orig
      @name =  replaced
    end
  end

  def declare_temporal_var(h = $varhash)
    ret = []
    name = get_name(self)
    tail = get_tail(self)
    if h[name][0] == nil && h[name][0] != "declared"
      type = @type
      type += "vec" if tail != nil
      ret += [Declaration.new([type,name])]
      h[name][0] = "declared"
    end
    ret
  end

  def convert_to_code(conversion_type)
    #p self
    #p @name
    #p $varhash[@name]
    ret=""
    ret += @name.convert_to_code(conversion_type) + " = " + @expression.convert_to_code(conversion_type) + ";"
    ret
  end
end

class NonSimdState < Statement
  def convert_to_code(conversion_type)
    @name.convert_to_code("reference") + " = " + @expression.convert_to_code("reference") + ";"
  end
end

class FuncCall
  attr_accessor :name, :ops, :type
  def initialize(x)
    @name, @ops = x
    @type = nil
  end
  def get_type(h = $varhash)
    ret_type = nil
    if $reserved_function.index(@name)
      #ret_type = h[@ops[0]][1]
      ret_type = @ops[0].get_type(h)
    else
      if $funchash[@name] != nil
        func = $funchash[@name]
        #p func
        abort "error: number of operands(= #{@ops.length}) for function #{@name} is wrong, #{func.decl.vars.length}" if @ops.length != func.decl.vars.length

        a = []
        func.decl.vars.zip(@ops).each{|tmp, orig|
          type = nil
          if orig.class == Expression
            type = [nil,orig.get_type,nil]
          else
            type = h[orig]
          end
          abort "error: undefined reference to #{orig} in get_type of FuncCall" if type == nil
          a += [tmp, type]
        }
        tmphash = Hash[*a]
        #p tmphash
        if func.statements != nil
          func.statements.each{ |stmt|
            #p stmt
            type = stmt.get_type(tmphash)
            if tmphash[stmt.name] == nil
              tmphash[stmt.name] = [nil,type,nil]
            end
          }
        end
        ret_type = func.retval.ret.get_type(tmphash)
        abort "error: function returns vector type is not allowed" if ret_type =~ /vec/
      else
        warn "error: undefined reference to function #{@name}"
      end
    end
    abort "error: function returns vector variable is not supported" if ret_type =~ /vec/
    abort "error: type inference failed for #{self}" if ret_type == nil
    ret_type
  end

  def find_function
    [self]
  end

  
  def expand_tree
  end

  def get_related_variable
    ret = []
    @ops.each{ |op|
      ret += op.get_related_variable
    }
    ret.sort!
    ret.uniq!
    ret
  end
  
  def isJRelated(list)
    ret = false
    @ops.each{ |op|
      ret |= op.isJRelated(list)
    }
    ret
  end
  
  def convert_to_code(conversion_type)
    case conversion_type
    when "reference" then
      retval = @name
      # retval += "<PS::" + get_type
      # @ops.each{ |op|
      #   retval += ",PS::" + op.get_type
      # }
      # retval += ">"
      retval += "("
      if @ops.length > 0
        @ops.each_with_index{ |op,i|
          retval += "," if i > 0
          retval += op.convert_to_code(conversion_type)
        }
      end
      retval += ")"
    when /A64FX/
      retval = self.convert_to_code_a64fx(conversion_type)
    else
      abort "error: unsupported conversion_type #{conversion_type}"
    end
    retval
  end
end

class Pragma
  attr_accessor :name, :option
  def initialize(x)
    @name,@option = x
  end

  def convert_to_code(conversion_type)
    ret = "#pragma #{@name}"
    if @option != nil then
      @option.each{ |x|
        ret += " #{x}"
      }
    end
    ret
  end
end

class IfElseState
  attr_accessor :operator, :expression
  def initialize(x)
    @operator, @expression = x
  end

  def convert_to_code(conversion_type)
    ret = ""
    case conversion_type
    when "reference"
      case @operator
      when :if
        ret = "if(" + @expression.convert_to_code(conversion_type) + "){"
      when :elsif
        ret = "}else if(" + @expression.convert_to_code(conversion_type) + "){"
      when :else
        ret = "}else{"        
      when :endif
        ret = "}"
      else
        abort "error: undefined if operator #{@operator}"
      end
    when /A64FX/
      ret = self.convert_to_code_a64fx(conversion_type)
    else
      abort "error: unsupported conversion_type #{conversion_type}"
    end
    ret
  end
end

class StoreState
  attr_accessor :dest,:src,:type
  def initialize(x)
    @dest, @src, @type = x
  end
  def get_type(h = $varhash)
    @type = @dest.get_type if @type == nil
    @type
  end
  def convert_to_code(conversion_type)
    case conversion_type
    when "reference"
      ret = @dest.convert_to_code(conversion_type) + "=" + @src.convert_to_code(conversion_type) + ";"
    when /A64FX/
      ret = "svst1_#{get_type_suffix(type)}(#{$current_predicate},#{@dest.convert_to_code(conversion_type)},#{@src.convert_to_code(conversion_type)});"
    end
  end
end

class LoadState
  attr_accessor :dest,:src,:type
  def initialize(x)
    @dest, @src, @type = x
  end
  def get_type(h = $varhash)
    @type = @dest.get_type if @type == nil
    @type
  end
  def convert_to_code(conversion_type)
    case conversion_type
    when "reference"
      ret = @dest.convert_to_code(conversion_type) + "=" + @src.convert_to_code(conversion_type) + ";"
    when /A64FX/
      ret = "#{@dest.convert_to_code(conversion_type)} = svld1_#{get_type_suffix(type)}(#{$current_predicate},#{@src.convert_to_code(conversion_type)});"
    end
  end
end

class PointerOf
  attr_accessor :type,:exp
  def initialize(x)
    @type,@exp = x
  end
  def get_pointer_type(type)
    case type
    when /F64/
      "double*"
    when /F32/
      "float*"
    when /S64/
      "long*"
    when /S32/
      "int*"
    when /U64/
      "unsigned long long*"
    when /U32/
      "unsigned int*"
    end
  end

  def convert_to_code(conversion_type)
    #"(#{get_pointer_type(@type)})&#{@exp.convert_to_code(conversion_type)}"
    case conversion_type
    when "reference"
      "#{@exp.convert_to_code("reference")}"
    when /A64FX/
      "(#{get_pointer_type(@type)})&#{@exp.convert_to_code("reference")}"
    end
  end
end

class GatherLoad
  attr_accessor :dest,:src,:offset,:type
  def initialize(x)
    @dest,@src,@offset,@type = x
  end
  def convert_to_code(conversion_type)
    $gather_load_count = 0 if $gather_load_count == nil
    ret = ""
    case conversion_type
    when /A64FX/
      index_name = "index_gather_load#{$gather_load_count}"
      vindex_name = "v" + index_name
      nsimd = get_simd_width(conversion_type)
      index = "0"
      for i in 1...nsimd
        index +=  ",#{i*offset.to_i}"
      end
      ret += "static uint32_t #{index_name}[#{nsimd}] = {#{index}};\n"
      ret += "svuint32_t #{vindex_name} = svld1_u32(#{$current_predicate},#{index_name});\n"
      ret += "#{@dest.convert_to_code(conversion_type)} = svld1_gather_u32index_#{get_type_suffix(@type)}(#{$current_predicate},#{@src.convert_to_code(conversion_type)},#{vindex_name});"
    else
      abort "unsupported conversion type for GatherLoad"
    end
    $gather_load_count += 1
    ret
  end
end

class ScatterStore
  attr_accessor :dest,:src,:offset,:type
  def initialize(x)
    @dest,@src,@offset,@type = x
  end
  def convert_to_code(conversion_type)
    $scatter_store_count = 0 if $scatter_store_count == nil
    ret = ""
    case conversion_type
    when /A64FX/
      index_name = "index_scatter_store#{$scatter_store_count}"
      vindex_name = "v" + index_name
      nsimd = get_simd_width(conversion_type)
      index = "0"
      for i in 1...nsimd
        index +=  ",#{i*offset.to_i}"
      end
      ret += "static uint32_t #{index_name}[#{nsimd}] = {#{index}};\n"
      ret += "svuint32_t #{vindex_name} = svld1_u32(#{$current_predicate},#{index_name});\n"
      ret += "svst1_scatter_u32index_#{get_type_suffix(@type)}(#{$current_predicate},#{@dest.convert_to_code(conversion_type)},#{vindex_name},#{@src.convert_to_code(conversion_type)});"
    else
      abort "unsupported conversion type for GatherLoad"
    end
    $scatter_store_count += 1
    ret
  end
end

class StructureLoad
  attr_accessor :dest,:src,:nelem,:type
  def initialize(x)
    @dest,@src,@nelem,@type = x
    abort "nelem must be 1 to 4" if nelem > 4 || nelem <= 0
  end
  def convert_to_code(conversion_type)
    "#{@dest.convert_to_code(conversion_type)} = svld#{nelem}_#{get_type_suffix(@type)}(#{$current_predicate},#{@src.convert_to_code(conversion_type)});"
  end
end

class StructureStore
  attr_accessor :dest,:src,:nelem,:type
  def initialize(x)
    @dest,@src,@nelem,@type = x
    abort "nelem must be 1 to 4" if nelem > 4 || nelem <= 0
  end
  def convert_to_code(conversion_type)
    "svst#{nelem}_#{get_type_suffix(@type)}(#{$current_predicate},#{@dest.convert_to_code(conversion_type)},#{@src.convert_to_code(conversion_type)});"
  end
end

class Duplicate
  attr_accessor :name,:expression,:type
  def initialize(x)
    @name,@expression,@type = x
  end
  def get_type(h = $varhash)
    @type = @expression.get_type(h) if @type == nil
    @type
  end
  def convert_to_code(conversion_type)
    case conversion_type
    when "reference"
      "#{@name.convert_to_code(conversion_type)} = #{@expression.convert_to_code(conversion_type)};"
    when /A64FX/
      "#{@name.convert_to_code(conversion_type)} = svdup_n_#{get_type_suffix(@type)}(#{@expression.convert_to_code(conversion_type)});"
    end
  end
end

class MADD
  attr_accessor :operator, :aop, :bop, :cop, :type
  def initialize(x)
    @operator, @aop, @bop, @cop, @type = x
  end

  def derive_type(op,aop,bop,cop,h=$varhash)
    type=nil
    at = aop.get_type(h)
    bt = bop.get_type(h)
    ct = cop.get_type(h)
    if at == "F64" || bt == "F64" || ct == "F64"
      type="F64"
    elsif at == "F32" || bt == "F32" || ct == "F32"      
      type="F32"
    else
      type="F16"
    end
    type
  end
  
  def get_type(h = $varhash)
    if @type == nil
      @type = derive_type(@operator,@aop,@bop,@cop,h)
    end
    @type
  end

  def get_related_variable
    ret = []
    ret += @aop.get_related_variable
    ret += @bop.get_related_variable
    ret += @cop.get_related_variable
    ret.sort!
    ret.uniq!
    ret
  end
  def isJRelated(list)
    ret = false
    ret |= @aop.isJRelated(list)
    ret |= @bop.isJRelated(list)
    ret |= @cop.isJRelated(list)
    ret
  end
  
  def convert_to_code(conversion_type)
    retval=""
    case conversion_type
    when "reference" then
      case @operator
      when :madd then
        retval = "(#{@cop.convert_to_code(conversion_type)} + #{@aop.convert_to_code(conversion_type)}*#{@bop.convert_to_code(conversion_type)})"
      when :msub then
        retval = "(#{@cop.convert_to_code(conversion_type)} - #{@aop.convert_to_code(conversion_type)}*#{@bop.convert_to_code(conversion_type)})"
      when :nmadd then
        retval = "(-(#{@cop.convert_to_code(conversion_type)} + #{@aop.convert_to_code(conversion_type)}*#{@bop.convert_to_code(conversion_type)}))"
      when :nmsub then
        retval = "(#{@aop.convert_to_code(conversion_type)}*#{@bop.convert_to_code(conversion_type)}-#{@cop.convert_to_code(conversion_type)})"
      else
        abort "error: unsupported operator for MADD class"
      end
    when /A64FX/
      self.get_type
      case @operator
      when :madd
        retval += "svmad_#{get_type_suffix(@type)}_z("
      when :msub
        retval += "svmsb_#{get_type_suffix(@type)}_z("
      when :nmadd
        retval += "svnmad_#{get_type_suffix(@type)}_z("
      when :nmsub
        retval += "svnmsb_#{get_type_suffix(@type)}_z("
      end
      retval += $current_predicate + ","
      retval += @aop.convert_to_code(conversion_type) + ","
      retval += @bop.convert_to_code(conversion_type) + ","
      retval += @cop.convert_to_code(conversion_type) + ")"
    end
    retval
  end
end

class Merge
  attr_accessor :op1,:op2,:type
  def initialize(x)
    @op1,@op2,@type = x
  end
  def get_related_variable
    [@op1,@op2]
  end
  def convert_to_code(conversion_type)
    ret = ""
    case conversion_type
    when /reference/
      ret = "#{@op1.convert_to_code(conversion_type)}"
    when /A64FX/
      ret = "svsel_#{get_type_suffix(@type)}(#{$current_predicate},#{@op1.convert_to_code(conversion_type)},#{@op2.convert_to_code(conversion_type)});"
    end
    ret
  end
end

class Expression
  attr_accessor :operator, :lop, :rop, :type
  def initialize(x)
    @operator, @lop, @rop, @type=x
  end

  def derive_type (operator,lop,rop,h=$varhash)
    type=nil
    lt = lop.get_type(h)
    #print "derive type ", operator," ", lop," ", rop, "\n"
    #p self
    if [:plus, :minus, :mult, :div,:and,:or,:eq,:neq,:lt,:le,:gt,:ge].index(operator)
      rt = rop.get_type(h)
      abort "type is not derived for #{lop} = #{lt}, #{rop} = #{rt} in derive_type" if lt == nil || rt == nil
      #p lop,rop
      #print "#{lt},#{rt}\n"
      #p h
      if operator == :mult && lt.index("vec")  && rt.index("vec")
        if lt == "F64vec" || rt == "F64vec"
          type="F64"
        elsif lt == "F32vec" || rt == "F32vec"
          type="F32"
        else
          type="F16"
        end
      else
        type="U"
        type="S" if lt.index("S") || rt.index("S")
        type="F" if lt.index("F") || rt.index("F")
        if lt.index("64") ||  rt.index("64")  
          type += "64"
        elsif lt.index("32") || rt.index("32")
          type += "32"
        else
          p lop,rop,h
          nil[0]
          abort "U16 is currently not supported"
          type += "16"
        end
        type += "vec" if lt.index("vec") ||  rt.index("vec")  
      end
    elsif operator == :uminus || operator == :array
      type = lt
    elsif operator == :dot
      if rop == "x" || rop == "y" || rop == "z"
        if lt == "F64vec"
          type = "F64"
        elsif lt == "F32vec"
          type = "F32"
        elsif lt == "F16vec"
          type = "F16"
        else
          warn "error: #{get_name(lop.name)} is not vector type!"
          abort
        end
      elsif lop =~ /^\d+$/ && rop =~ /^\d+(f|h)?$/
        if rop =~ /f/ then
          type = "F32"
        elsif rop =~ /h/ then
          type = "F16"
        else
          type = "F64"
        end
      else
        abort "error: unsupported vector member \"#{rop.name}\""
      end
    else
      p self
      abort "error: undefined reference to #{self} in derive_type"
      type=rt
    end
    #print "final type=", type, "\n"
    type
  end

  def  get_type(h = $varhash)
    if @operator == :dot && @lop =~ /^\d+$/
      if @rop =~ /^\d+h$/
        @type = "F16"
      elsif @rop =~ /^\d+f$/
        @type = "F32"
      elsif @rop =~ /^\d+$/
        @type = "F64"
      else
        p self
        abort "unsupported floatinng point number! (must be (\d+).(\d+[f,h]))"
      end
    end
    if @type == nil
      @type = derive_type(@operator, @lop, @rop, h)
    end
    abort "get_type failed" if @type == nil
    @type
  end

  def find_function
    ret = []
    if @lop.class == FuncCall
      ret += [@lop]
    elsif !isLeaf(@lop)
      ret += @lop.find_function
    end
    if @rop != nil
      if @rop.class == FuncCall
        ret += [@rop]
      elsif !isLeaf(@rop)
        ret += @rop.find_function
      end
    end
    ret
  end

  def split_vector_expression(dim)
    ret = self
    if self.class != String
      
    else
      if @type =~ /vec/
        
      end
    end
  end

  def get_related_variable
    ret = []
    case @operator
    when :dot
      ret += [@lop] if @lop =~ /[a-z_].*/
    when :array
      ret += [@lop]
    else
      ret += @lop.get_related_variable
      ret += @rop.get_related_variable if @rop != nil
    end
    ret.sort!
    ret.uniq!
    ret
  end
  def isJRelated(list)
    ret = false
    if isLeaf(self)
      abort "only . and [] operator can be Leaf for Expression class" if @operator != :dot && @operator != :array
      if @lop =~ /\d+/
      else
        ret |= @lop.isJRelated(list)
      end
    else
      ret |= @lop.isJRelated(list)
      ret |= @rop.isJRelated(list) if @rop != nil
    end
    ret
  end

  def get_opcode(op)
    case op
    when :plus
      "+"
    when :minus
      "-"
    when :mult
      "*"
    when :div
      "/"
    when :lt
      "<"
    when :le
      "<="
    when :gt
      ">"
    when :ge
      ">="
    when :eq
      "=="
    when :neq
      "!="
    when :eq
      "=="
    when :and
      "&&"
    when :or
      "||"
    else
      abort "error: undefined opcode #{op}"
    end
  end

  def replace_recursive(orig,replaced)
    #p self
    @lop = @lop.replace_recursive(orig,replaced)
    @rop = @rop.replace_recursive(orig,replaced) if @rop != nil
    self
  end
  def replace_fdpsname_recursive(h=$varhash)
    #p self
    ret = self.dup
    if @operator == :array
      name = ret.lop.dup
      iotype = h[name][0]
      type   = h[name][1]
      fdpsname = h[name][2]
      abort "array expression must be used for EPI, EPJ, or FORCE variable" if fdpsname == nil || iotype == nil
      ret = Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),ret.rop]),fdpsname,type])
    elsif @operator == :dot
      ret.lop = ret.lop.replace_fdpsname_recursive(h)
    else
      ret.lop = ret.lop.replace_fdpsname_recursive(h)
      ret.rop = ret.rop.replace_fdpsname_recursive(h) if ret.rop != nil
    end
    ret
  end
  
  def convert_to_code(conversion_type)
    #p self
    retval=""
    case conversion_type
    when "reference"
      if type=[:plus,:minus,:mult,:div,:eq,:neq,:lt,:le,:gt,:ge,:and,:or].index(@operator)
        opcode=get_opcode(@operator)
        retval="("+@lop.convert_to_code(conversion_type)+opcode+
               @rop.convert_to_code(conversion_type)+")"
      elsif @operator == :uminus
        retval="-(" + @lop.convert_to_code(conversion_type) + ")"
      elsif @operator == :dot
        retval=@lop.convert_to_code(conversion_type)+"."+
               @rop.convert_to_code(conversion_type)
      elsif @operator == :array
        retval=lop.convert_to_code(conversion_type) + "[" + @rop.convert_to_code(conversion_type) + "]"
      elsif @operator == :func
        retval= @lop+"("+ @rop.convert_to_code(conversion_type)+")"
      else
        abort "error: unsupported operator #{@operator} in expression code conversion"
      end
    when /A64FX/
      retval = self.convert_to_code_a64fx(conversion_type)
    end
    retval
  end
end

class NonSimdExp < Expression
  def convert_to_code(conversion_type)
    super("reference")
  end
end
