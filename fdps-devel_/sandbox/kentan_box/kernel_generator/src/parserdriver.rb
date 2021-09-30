require "common.rb"
require "intermediate_exp_class.rb"
require "disassemble_statement.rb"
require "kernelparser.rb"
require "expand_function.rb"
require "reduce_madd.rb"
require "loop_fission.rb"
require "software_pipelining.rb"
require "A64FX.rb"

$dimension = 3
$reserved_function = ["rsqrt","sqrt","inv","max","min","madd","msub","nmadd","nmsub"]
$tmp_val_count = 0

def accumulate_related_variable(orig,vars,h)
  ret = []
  if vars != nil
    ret = vars
    vars.each{ |v|
      if  h[v] != nil
        if !($varhash[v][0] =~ /(EPI|EPJ|FORCE)/)
          ret += accumulate_related_variable(v,h[v][0],h) if v != orig && h[v][1] == true
        end
        h[v][1] = false
      end      
    }
  end
  ret
end

def generate_related_map(fs,ss,h=$varhash)
  tmp = []
  ss.each{ |s|
    #p get_name(s),h[get_name(s)] if isStatement(s)
    if isStatement(s) && h[get_name(s)][3] != "local"
      exp = s.expression.get_related_variable
      tmp += [[get_name(s),exp]] if exp != []
    elsif s.class == ConditionalBranch
      tmp += s.get_related_variable
    end
  }
  tmp.sort!
  #tmp.sort_by{ |v| v[0]}
  tmp.uniq!
  tmp2 = Hash.new()
  tmp.each{ |v|
    if tmp2[v[0]] == nil
      tmp2[v[0]] = [v[1],true]
    else
      tmp2[v[0]][0] += v[1]
    end
  }

  tmp3 = []
  fs.each{ |f|
    tmp3 += accumulate_related_variable(f,tmp2[f][0],tmp2)
  }
  ss.each{ |s|
    tmp3 += s.expression.get_related_variable if s.class == IfElseState && s.expression != nil
    tmp3 += s.get_cond_related_variable if s.class == ConditionalBranch
  }
  ret = tmp3.sort.uniq

  ret
end

def generate_force_related_map(ss,h=$varhash)
  fs = []
  h.each{ |v|
    iotype = v[1][0]
    fs += [v[0]] if iotype == "FORCE"
  }
  ret = generate_related_map(fs,ss,h)
  ret
end

class Kernelprogram
  def generate_hash(kerneltype)
    #    print "print kernel\n"
    #p self
    #p $funchash
    $varhash=process_iodecl(@iodeclarations)
    $funchash=process_funcdecl(@functions)
    #p $varhash
    @statements.each{|s|
      #p s
      if isStatement(s)
        if s.name.class == Expression
          if s.name.operator == :dot
            if ["x","y","z"].index(s.name.rop)
              type = s.get_type
              s.type = type + "vec"
              s.add_to_varhash
              s.type = type
            else
              abort  "left value must be vector or scalar variable"
            end
          else
            abort "only vector expression is allowed for left value"
          end
        else
          s.get_type
          s.add_to_varhash
        end
      end
    }
    #reserved variables
    ["ni","nj","i","j","jj","jjj"].each{ |x|
      $varhash[x] = [nil,"S32",nil]
    }
    #p $funchash
  end
  def print_function(conversion_type)
    @functions.each{ |f|
      print f.convert_to_code(conversion_type)
    }
  end
  def print_statements(conversion_type)
    @statements.each{|s|
      #p s
      print s.convert_to_code(conversion_type)+"\n"
    }
  end
 
  def split_coefficients(orig)
    #p "split_coefficients",orig
    ret = []
    exp = orig.dup
    if !isLeaf(orig)
      if orig.operator == :mult
        lt = orig.lop.get_type
        rt = orig.rop.get_type
        if ((lt =~ /vec/) && !(rt =~ /vec/))
          #p "split candidate lop: #{orig.rop} with #{orig.lop}"
          if !isLeaf(orig.rop)
            tmp_name = add_new_tmpvar(rt)
            tmp = Statement.new([tmp_name,Expression.new([orig.rop.operator,orig.rop.lop,orig.rop.rop])])
            tmp.type = rt
            tmp.expression.type = orig.rop.get_type
            tmp.add_to_varhash
            orig = Expression.new([orig.operator,orig.lop,tmp_name,orig.type])
            ret += split_coefficients(orig.lop)
            ret.append(tmp)
          end
        elsif (!(lt =~ /vec/) && (rt =~ /vec/))
          #p "split candidate lop: #{orig.lop} with #{orig.rop}"
          if !isLeaf(orig.lop)
            tmp_name = add_new_tmpvar(lt) #"__fkg_tmp#{$tmp_val_count}"
            tmp = Statement.new([tmp_name,Expression.new([orig.lop.operator,orig.lop.lop,orig.lop.rop])])
            tmp.type = lt
            tmp.expression.type = orig.lop.get_type
            tmp.add_to_varhash
            #orig.lop = tmp_name
            orig = Expression.new([orig.operator,tmp_name,orig.rop,orig.type])
            ret += split_coefficients(orig.rop)
            ret.append(tmp)
          end
        else
          ret += split_coefficients(orig.lop)
          ret += split_coefficients(orig.rop) if orig.rop != nil
        end
      else
        ret += split_coefficients(orig.lop)
        ret += split_coefficients(orig.rop) if orig.rop != nil
      end
    else
      #if orig.class == FuncCall
      #tmp_name = "__fkg_tmp#{$tmp_val_count}"
      #tmp = Statement.new([tmp_name,orig])
      #tmp.type = orig.get_type
      #tmp.expression.type = orig.get_type
      #tmp.add_to_varhash
      #p orig.get_type
      #p tmp
      #orig = tmp_name
      #ret.append(tmp)
      #$tmp_val_count += 1
      #end
    end
    #p "return value:", ret, orig
    ret
  end

  def vector_to_scalar(exp,dim)
    #p "vector_to_scalar:"
    ret = exp.dup
    if !isLeaf(ret)
      #p ret
      ret.type = ret.type.delete("vec") if isVector(ret)
      ret.lop = vector_to_scalar(ret.lop,dim)
      ret.rop = vector_to_scalar(ret.rop,dim)
    else
      if isVector(ret)
        ret = Expression.new([:dot,ret,dim])
        ret.type = exp.get_type.delete("vec")
        #ret.type = ret.type.delete("vec")
      end
    end
    ret
  end
  
  def split_vector_expression(orig,dim = $dimension)
    ret = []
    if isVector(orig)
      if !(orig.expression.class == FuncCall)
        ["x","y","z"].each{ |d|
          val = Expression.new([:dot,orig.name,d])
          exp = vector_to_scalar(orig.expression,d)
          #p exp
          type = orig.type.delete("vec")
          tmp = Statement.new([val,exp])
          tmp.type = type
          ret.append(tmp)
        }
      else
        ret.append(orig)
      end
    else
      ret = [orig]
    end
    ret
  end
  
  def expand_vector_statement(orig)
    ret = []
    val = orig.name
    exp = orig.expression
    ret += split_coefficients(exp)
    ret += split_vector_expression(orig)
    ret
  end

  def expand_tree
    @statements.each{ |s|
      if isStatement(s)
        s.expand_tree
      end
    }
    new_s = []
    @statements.each{ |orig|
      if isStatement(orig)
        expand_vector_statement(orig).each{ |s|
          s.reduce_madd_recursive
          s.reduce_negate_recursive
          new_s.append(s)
        }
      else
        new_s.append(orig)
      end
    }
    @statements = new_s
    #p "expanded statements:",@statements
  end

  def calc_max_predicate_count(ss)
    ret = 1
    count = 1
    ss.each{ |s|
      if s.class == IfElseState
        if s.operator == :if
          count += 1
        elsif s.operator == :endif
          count -= 1
        end
        ret = count if count > ret
      elsif s.class == ConditionalBranch
        tmp = s.calc_max_predicate_count
        ret = tmp if ret < tmp
      end
    }
    ret
  end
  
  def kernel_class_def(conversion_type)
    code = ""
    code += "struct #{$kernel_name}{\n"
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += "PS::" + type + " " + name + ";\n"
      end
    }
   
    code += "#{$kernel_name}("
    member_count = 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += "," if member_count > 0
        code += "PS::" + type + " " + name
        member_count = member_count+1
      end
    }
    code += ")"
    member_count = 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        code += ":" if member_count == 0
        code += "," if member_count > 0
        code += name + "(" + name +")"
        member_count = member_count+1
      end
    }
    code += "{}\n"
    code += "void operator()(const #{$epi_name}* __restrict__ epi,const int ni,const #{$epj_name}* __restrict__ epj,const int nj,#{$force_name}* __restrict__ force){\n"
    if conversion_type =~ /A64FX/
      $pg_count = 0
      $current_predicate = "pg#{$pg_count}"
      $max_pg_count = calc_max_predicate_count(@statements)
      for i in 0...$max_pg_count
        code += "svbool_t pg#{i};\n"
      end
    end

    code
  end

  def reserved_func_def(conversion_type)
    code = ""
    case conversion_type
    when "reference" then
      # code += "template<typename Tret,typename Top>\n"
      # code += "Tret rsqrt(Top op){ return (Tret)1.0/std::sqrt(op); }\n"
      # code += "template<typename Tret,typename Top>\n"
      # code += "Tret sqrt(Top op){ return std::sqrt(op); }\n"
      # code += "template<typename Tret,typename Top>\n"
      # code += "Tret inv(Top op){ return 1.0/op; }\n"
      # code += "template<typename Tret,typename Ta,typename Tb>\n"
      # code += "Tret max(Ta a,Tb b){ return std::max(a,b);}\n"
      # code += "template<typename Tret,typename Ta,typename Tb>\n"
      # code += "Tret min(Ta a,Tb b){ return std::min(a,b);}\n"
      code += "PS::F64 rsqrt(PS::F64 op){ return 1.0/std::sqrt(op); }\n"
      code += "PS::F64 sqrt(PS::F64 op){ return std::sqrt(op); }\n"
      code += "PS::F64 inv(PS::F64 op){ return 1.0/op; }\n"
      code += "PS::F64 max(PS::F64 a,PS::F64 b){ return std::max(a,b);}\n"
      code += "PS::F64 min(PS::F64 a,PS::F64 b){ return std::min(a,b);}\n"

      code += "PS::F32 rsqrt(PS::F32 op){ return 1.f/std::sqrt(op); }\n"
      code += "PS::F32 sqrt(PS::F32 op){ return std::sqrt(op); }\n"
      code += "PS::F32 inv(PS::F32 op){ return 1.f/op; }\n"

      code += "PS::S64 max(PS::S64 a,PS::S64 b){ return std::max(a,b);}\n"
      code += "PS::S64 min(PS::S64 a,PS::S64 b){ return std::min(a,b);}\n"
      code += "PS::S32 max(PS::S32 a,PS::S32 b){ return std::max(a,b);}\n"
      code += "PS::S32 min(PS::S32 a,PS::S32 b){ return std::min(a,b);}\n"
    when /A64FX/ then
      code += self.reserved_func_def_a64fx(conversion_type)
    end
    code
  end

  def find_j_variable(exp)
  end
  

  def split_accum_statement(ss,h = $varhash)
    j_list = []
    h.each{ |v|
      iotype = v[1][0]
      if iotype == "EPJ"
        name = v[0]
        j_list += [name]
      end
    }

    init = []
    body = []
    accum = []
    ss.each{|s|
      if isStatement(s)
        name = get_name(s.name)
        if h[name][3] != "local"
          if s.expression.isJRelated(j_list)
            j_list += [name] if j_list.index(name) == nil
            body += [s]
          else
            v = h[name]
            iotype = v[0]
            case iotype
            when "FORCE"
              accum += [s]
            when nil
              init += [s]
            else
              body += [s]
            end
          end
        end
      else
        body += [s]
      end
    }
    [init,body,accum]
  end
  
  def kernel_body(conversion_type,h=$varhash,istart=0)
    code = ""
    @statements.each{ |s|
      #p s
      code += s.declare_variable(conversion_type) if isStatement(s)
    }
    h.each{ |v|
      modify = v[1][3]
      if modify == "local"
        name   = v[0]
        iotype = v[1][0]
        type   = v[1][1]
        code += "PS::#{type} #{name}_tmp[ni];\n" if iotype == "EPI"
        code += "PS::#{type} #{name}_tmp[nj];\n" if iotype == "EPJ"
      end
    }
    @statements.each{ |s|
      name = get_name(s.name) if isStatement(s)
      if name != nil && h[name][3] == "local"
        tail = get_tail(s.name)
        iotype = h[name][0]
        exp = s.expression
        code += "for(int i=0;i<ni;i++){\n" if iotype == "EPI" || iotype == "FORCE"
        code += "for(int i=0;i<nj;i++){\n" if iotype == "EPJ"
        if tail != nil
          code += "#{name}_tmp[i].#{tail} = "
        else
          code += "#{name}_tmp[i] = "
        end
        tmp = exp.convert_to_code("reference")
        h.each{ |v|
          name = v[0]
          iotype = v[1][0]
          fdpsname = v[1][2]
          case iotype
          when "EPI"
            tmp.gsub!(name,"epi[i].#{fdpsname}")
          when "EPJ"
            tmp.gsub!(name,"epj[i].#{fdpsname}")
          when "FORCE"
            tmp.gsub!(name,"force[i].#{fdpsname}")
          end
        }
        code += tmp + ";\n"
        code += "}\n"
      end
    }
    case conversion_type
    when "reference" then
      fvars = generate_force_related_map(@statements)
      code += "for(int i=#{istart};i<ni;i++){\n"
      $varhash.each{|v|
        iotype = v[1][0]
        if iotype == "EPI" || iotype == "FORCE"
          name     = v[0]
          type     = v[1][1]
          modifier = v[1][3]
          if modifier == "local"
            code += "PS::#{type} #{name} = #{name}_tmp[i];\n"
          else
            if fvars.index(name)
              fdpsname = "epi[i]."   + v[1][2] if iotype == "EPI"
              fdpsname = "force[i]." + v[1][2] if iotype == "FORCE"
              code += "PS::" + type + " " + name + " = " + fdpsname + ";\n"
            end
          end
        elsif iotype == "ACCUM"
          name = v[0]
          type = v[1][1]
          code += "PS::" + type + " " + name + ";\n"
        end
      }

      accum_init,ss,accum_finalize = split_accum_statement(@statements)
      accum_init.each{|s|
        code += s.convert_to_code(conversion_type)+"\n"
      }
      code += "for(int j=0;j<nj;j++){\n"
      $varhash.each{|v|
        name     = v[0]
        iotype   = v[1][0]
        type     = v[1][1]
        fdpsname = v[1][2]
        modifier = v[1][3]
        if iotype == "EPJ"
          if modifier == "local"
            code += "PS::" + type + " " + name + " = " + name + "_tmp[j];\n"
          elsif fvars.index(name)
            code += "PS::" + type + " " + name + " = " + "epj[j]." + fdpsname + ";\n"
          end
        end
      }
      ss.each{|s|
        code += s.convert_to_code(conversion_type)+"\n"
      }
      code += "}\n" # j loop
      accum_finalize.each{|s|
        code += s.convert_to_code(conversion_type) + "\n"
      }
      $varhash.each{|v|
        iotype = v[1][0]
        if iotype == "FORCE"
          name = v[0]
          type = v[1][1]
          fdpsname = "epi[i]."   + v[1][2] if iotype == "EPI"
          fdpsname = "force[i]." + v[1][2] if iotype == "FORCE"
          code += fdpsname + " = " + name + ";\n"
        end
      }
      code += "}\n" # i loop
      code += "}\n" # end of operator()

      #p @functions
    when "A64FX_F32_I16_J1" then
      code += self.kernel_body_a64fx_f32(conversion_type)
    when "A64FX_MP" then
      code += self.kernel_body_a64fx_mp(conversion_type)
    else
      warn "error: unsupported conversion_type #{conversion_type}"
    end
    code
  end
  
  def generate_optimized_code(conversion_type)
    code = ""
    if conversion_type =~ /A64FX/
      code += "#include <arm_sve.h>\n"
    end
    code += "#include \"user_defined_class.h\"\n"
    code += kernel_class_def(conversion_type)
    code += kernel_body(conversion_type)
    #@functions.each{ |f|
    #  code += f.convert_to_code(conversion_type)
    #}
    code += reserved_func_def(conversion_type)
    code += "};\n"
    print code
  end

  def aos2soa(fvars,conversion_type,h=$varhash)
    ret = []
    case conversion_type
    when "reference"
      h.each{|v|
        iotype = v[1][0]
        if iotype == "EPI" || iotype == "FORCE"
          name     = v[0]
          type     = v[1][1]
          modifier = v[1][3]
          ret += [Declaration.new([type,name])]
          if modifier == "local"
            if type =~ /vec/
              ret += [Statement.new([Expression.new([:dot,name,"x"]),Expression.new([:array,"#{name}_tmp_x",get_io_index(iotype),type]),type])]
              ret += [Statement.new([Expression.new([:dot,name,"y"]),Expression.new([:array,"#{name}_tmp_y",get_io_index(iotype),type]),type])]
              ret += [Statement.new([Expression.new([:dot,name,"z"]),Expression.new([:array,"#{name}_tmp_z",get_io_index(iotype),type]),type])]
            else
              ret += [Statement.new([name,Expression.new([:array,"#{name}_tmp",get_io_index(iotype),type]),type])]
            end
          elsif fvars.index(name)
            fdpsname = v[1][2]
            ret += [Statement.new([name,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype),type]),fdpsname,type]),type])]
          end
        end
      }
    when /A64FX/
      ret += aos2soa_a64fx(fvars,h)
    end
    ret
  end

  def soa2aos(fvars,conversion_type,h=$varhash)
    ret = []
    case conversion_type
    when "reference"
      h.each{|v|
        iotype = v[1][0]
        if iotype == "FORCE"
          name = v[0]
          type = v[1][1]
          fdpsname = v[1][2]
          if type =~ /vec/
            ret += [StoreState.new([Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype),type]),fdpsname,type]),"x"]),Expression.new([:dot,name,"x"]),type])]
            ret += [StoreState.new([Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype),type]),fdpsname,type]),"y"]),Expression.new([:dot,name,"y"]),type])]
            ret += [StoreState.new([Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype),type]),fdpsname,type]),"z"]),Expression.new([:dot,name,"z"]),type])]
          else
            ret += [StoreState.new([Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype),type]),fdpsname,type]),name,type])]
          end
        end
      }
    when /A64FX/
      ret += soa2aos_a64fx(fvars,h)
    end
    ret
  end

  def generate_iloop_begin(conversion_type)
    case conversion_type
    when "reference"
      ret = Loop.new(["i","0","ni",1,[]])
    when  /A64FX/
      nsimd = get_simd_width(conversion_type)
      ret = Loop.new(["i","0","((ni+#{nsimd-1})/#{nsimd})*#{nsimd}","#{nsimd}",["#{$current_predicate}=svwhilelt_b32_s32(i,ni);"]])
    end
    ret
  end
  
  def kernel_body2(conversion_type,h=$varhash,istart=0)
    code = ""
    ret = []
    ret.push(NonSimdDecl.new(["S32","i"]))
    ret.push(NonSimdDecl.new(["S32","j"]))
    # declare "local" variables
    h.each{ |v|
      modify = v[1][3]
      if modify == "local"
        name = v[0]
        iotype = v[1][0]
        type   = v[1][1]
        array_size = ["ni","nj"][["EPI","EPJ"].index(iotype)]
        if type =~ /vec/
          type = type.delete("vec")
          ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp_x",array_size,type])])]
          ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp_y",array_size,type])])]
          ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp_z",array_size,type])])]
        else
          ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp",array_size,type])])]
        end
      end
    }
    # calc or copy local variables from EPI, EPJ, or FORCE
    @statements.each{ |s|
      name = get_name(s) if s.class == Statement
      if name != nil && h[name][3] == "local"
        tail = get_tail(s.name)
        iotype = h[name][0]
        type   = h[name][1]
        exp = s.expression
        index = get_io_index(iotype)
        loop_tmp = Loop.new([index,"0","n#{index}",1,[]])
        if tail != nil
          new_name = Expression.new([:array,"#{name}_tmp_#{tail}",index,type])
        else
          new_name = Expression.new([:array,"#{name}_tmp",index,type]) #"#{name}_tmp[i]"
        end
        new_exp = exp.replace_fdpsname_recursive(h)
        loop_tmp.statements += [Statement.new([new_name,new_exp])]
        ret += [loop_tmp]
      end
    }

    ret.each{ |s|
      code += s.convert_to_code("reference")
    }
    ret = []
    # declare temporal variables
    @statements.each{ |s|
      if s.class == Statement || s.class == ConditionalBranch
        ret += s.declare_temporal_var
      end
    }

    nsimd = get_simd_width(conversion_type)
    # building i loop
    accum_init,ss,accum_finalize = split_accum_statement(@statements)
    fvars = generate_force_related_map(ss)
    # iloop = Loop.new(["i=#{istart}","i<((ni+#{nsimd-1})/#{nsimd})*#{nsimd}","i+=#{nsimd}",[]])
    iloop = generate_iloop_begin(conversion_type)

    # load EPI and FORCE variable
    iloop.statements += aos2soa(fvars,conversion_type)

    accum_init.each{|s|
      iloop.statements += [s]
    }

    iloop.statements += [NonSimdState.new(["j","0"])]
    if $strip_mining != nil
      warn "strip mining is applied"
      loop_fission_vars = find_loop_fission_load_store_vars(ss)
      fission_count = 0
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
          iloop.statements += [NonSimdDecl.new([type,"#{v}_tmp_x[#{nsimd}*#{$strip_mining}]"])]
          iloop.statements += [NonSimdDecl.new([type,"#{v}_tmp_y[#{nsimd}*#{$strip_mining}]"])]
          iloop.statements += [NonSimdDecl.new([type,"#{v}_tmp_z[#{nsimd}*#{$strip_mining}]"])]
        else
          iloop.statements += [NonSimdDecl.new([type,"#{v}_tmp[#{nsimd}*#{$strip_mining}]"])]
        end
      }

      jloop = Loop.new(["j",nil,"(nj/#{$strip_mining})*#{$strip_mining}","#{$strip_mining}",[NonSimdDecl.new(["S32","jj"])]])
      jjloop = Loop.new(["jj","0","#{$strip_mining}",1,[]])
      jjj = NonSimdExp.new([:plus,"j","jj"])
      first_loop = true
      ss.each{ |s|
        $unroll_stage = s.option[0].to_i if s.class == Pragma && s.name == "unroll"

        if (s.class == Pragma && s.name == "statement" && s.option == ["loop_fission_point"])|| first_loop
          if !first_loop
            loop_fission_vars[fission_count][0].each{ |v|
              type = h[v][1]
              if type =~ /vec/
                jjloop.statements += [StoreState.new([PointerOf.new([type,Expression.new([:array,"#{v}_tmp_x",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),Expression.new([:dot,"#{v}","x"]),type])]
                jjloop.statements += [StoreState.new([PointerOf.new([type,Expression.new([:array,"#{v}_tmp_y",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),Expression.new([:dot,"#{v}","y"]),type])]
                jjloop.statements += [StoreState.new([PointerOf.new([type,Expression.new([:array,"#{v}_tmp_z",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),Expression.new([:dot,"#{v}","z"]),type])]
              else
                jjloop.statements += [StoreState.new([PointerOf.new([type,Expression.new([:array,"#{v}_tmp",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),"#{v}",type])]
              end
            }
            #jjloop = software_pipelining(jjloop) if $swpl_stage > 1
            jjloop = loop_unroll(jjloop,$unroll_stage) if $unroll_stage > 1
            jloop.statements += [jjloop.dup]
          end
          first_loop = false
          jjloop = Loop.new(["jj","0","#{$strip_mining}",1,[]])

          loop_fission_vars[fission_count][1].each{ |v|
            iotype = h[v][0]
            type   = h[v][1]
            if iotype == "declared"
              if type =~ /vec/
                jjloop.statements += [LoadState.new([Expression.new([:dot,"#{v}","x"]),PointerOf.new([type,Expression.new([:array,"#{v}_tmp_x",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),type])]
                jjloop.statements += [LoadState.new([Expression.new([:dot,"#{v}","y"]),PointerOf.new([type,Expression.new([:array,"#{v}_tmp_y",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),type])]
                jjloop.statements += [LoadState.new([Expression.new([:dot,"#{v}","z"]),PointerOf.new([type,Expression.new([:array,"#{v}_tmp_z",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),type])]
              else
                jjloop.statements += [LoadState.new(["#{v}",PointerOf.new([type,Expression.new([:array,"#{v}_tmp",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),type])]
              end
            elsif iotype == "EPJ"
              name     = v
              fdpsname = h[v][2]
              modifier = h[v][3]

              jjloop.statements += [Declaration.new([type,name])]
              case modifier
              when "local"
                if type =~ /vec/
                  jjloop.statements += [Duplicate.new([Expression.new([:dot,name,"x"]),Expression.new([:array,"#{name}_tmp_x",jjj,type]),type])]
                  jjloop.statements += [Duplicate.new([Expression.new([:dot,name,"y"]),Expression.new([:array,"#{name}_tmp_y",jjj,type]),type])]
                  jjloop.statements += [Duplicate.new([Expression.new([:dot,name,"z"]),Expression.new([:array,"#{name}_tmp_z",jjj,type]),type])] 
                else
                  jjloop.statements += [Duplicate.new([name,Expression.new([:array,"#{name}_tmp",jjj,type]),type])]
                end
              else
                if type =~ /vec/
                  stype = type.delete("vec")
                  jjloop.statements += [Statement.new([Expression.new([:dot,name,"x"]),Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),jjj]),fdpsname,type]),"x"]),stype])]
                  jjloop.statements += [Statement.new([Expression.new([:dot,name,"y"]),Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),jjj]),fdpsname,type]),"y"]),stype])]
                  jjloop.statements += [Statement.new([Expression.new([:dot,name,"z"]),Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),jjj]),fdpsname,type]),"z"]),stype])]
                else
                  jjloop.statements += [Statement.new([name,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),jjj]),fdpsname,type]),type])]
                end
              end
            elsif fvars.index(name)
              fdpsname = v[1][2]
              jjloop.statements += [Duplicate.new([name,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype)]),fdpsname,type]),type])]
            end
          }
          fission_count += 1
        end # s.class != Pragma
        jjloop.statements += [s] if (!isStatement(s) || h[get_name(s)][3] != "local") && s.class != Pragma
      } # ss.each
      #jjloop = software_pipelining(jjloop) if $swpl_stage > 1
      jjloop = loop_unroll(jjloop,$unroll_stage) if $unroll_stage > 1
      jloop.statements += [jjloop.dup]

      iloop.statements += [jloop]
    end # strip_mining

    # tail loop
    jloop = Loop.new(["j",nil,"nj",1,[]])
    fvars.each{|v|
      iotype = h[v][0]
      if iotype == "EPJ"
        name     = v
        type     = h[v][1]
        modifier = h[v][3]
        jloop.statements += [Declaration.new([type,name])]
        if modifier == "local"
          if type =~ /vec/
            jloop.statements += [Duplicate.new([Expression.new([:dot,name,"x"]),Expression.new([:array,"#{name}_tmp_x",get_io_index(iotype),type]),type])]
            jloop.statements += [Duplicate.new([Expression.new([:dot,name,"y"]),Expression.new([:array,"#{name}_tmp_y",get_io_index(iotype),type]),type])]
            jloop.statements += [Duplicate.new([Expression.new([:dot,name,"z"]),Expression.new([:array,"#{name}_tmp_z",get_io_index(iotype),type]),type])]
          else
            jloop.statements += [Duplicate.new([name,Expression.new([:array,"#{name}_tmp",get_io_index(iotype),type]),type])]
          end
        elsif fvars.index(name)
          fdpsname = h[v][2]
          jloop.statements += [Statement.new([name,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype)]),fdpsname,type]),type])]
        end
      end
    }
    ss.each{|s|
      jloop.statements += [s] if s.class != Pragma
    }
    iloop.statements += [jloop]

    accum_finalize.each{|s|
      iloop.statements += [s]
    }
    iloop.statements += soa2aos(fvars,conversion_type)

    ret += [iloop]
    ret.each{ |s|
      code += s.convert_to_code(conversion_type)
    }
    #abort "kernel_body2 is testing"
    code
  end
  def generate_optimized_code2(conversion_type)
    code = ""
    if conversion_type =~ /A64FX/
      code += "#include <arm_sve.h>\n"
    end
    code += "#include \"user_defined_class.h\"\n"
    code += kernel_class_def(conversion_type)
    code += kernel_body2(conversion_type)
    code += "}\n"
    code += reserved_func_def(conversion_type)
    code += "};\n"
    print code
  end

  def make_conditional_branch_block
    @statements = make_conditional_branch_block_recursive(@statements)
  end
  
  def make_conditional_branch_block_recursive(ss)
    new_s = []
    nest_level = 0
    cbb = ConditionalBranch.new([[],[]])
    isLastBranch = true
    computed_vars = []
    ss.each{ |s|
      if s.class == IfElseState
        nest_level += 1 if s.operator == :if
        nest_level -= 1 if s.operator == :endif
        abort "nest level is less than 0" if nest_level < 0
      end
      if nest_level == 0
        if isStatement(s)
          computed_vars.push(get_name(s))
        end
        if s.class == IfElseState
          abort "operator #{s.operator} of #{s} must be :endif" if s.operator != :endif
          new_b = []
          cbb.bodies.each{ |b|
            new_b.push(make_conditional_branch_block_recursive(b))
          }
          cbb.bodies = new_b
          new_s.push(cbb)
          cbb.bodies.each{ |bss|
            bss.each{ |bs|
              if isStatement(bs) && bs.expression.class != Merge
                name = get_name(bs)
                tail = get_tail(bs)
                tmp = computed_vars.find(name)
                if tmp != nil
                  abort "merging after conditional branch is not supported for vector variable" if bs.type =~  /vec/
                  tmp_name = add_new_tmpvar(bs.type)
                  bs.replace_name(name,tmp_name)
                  bss.push(Statement.new([name,Merge.new([tmp_name,name,bs.type]),bs.type]))
                end
              end
            }
          }
          
          cbb = ConditionalBranch.new([[],[]])
        else
          new_s.push(s)
        end
      elsif nest_level == 1
        if s.class == IfElseState && s.operator != :endif
          cbb.push_condition(s)
          #p cbb
        else
          cbb.push_body(s)
        end
      else
        cbb.push_body(s)
      end
    }
    abort "ConditionalBranch is not terminated" if nest_level > 0
    new_s
  end
end

class String
  def get_type(h = $varhash)
    propaties=h[self]
    if propaties
      propaties[1]
    else
      p h[self]
      p self
      nil[0]
      abort "undefined reference to #{self} at get_type of String class"
    end
  end

  def get_related_variable
    [self]
  end
  def isJRelated(list)
    if list.index(self)
      true
    else
      false
    end
  end

  def replace_recursive(orig,replaced)
    if self == orig
      replaced
    else
      self
    end
  end
  def replace_fdpsname_recursive(h=$varhash)
    name = self.dup
    ret = name
    return ret if h[name] == nil
    iotype = h[name][0]
    type = h[name][1]
    fdpsname = h[name][2]
    if iotype != nil && fdpsname != nil
      op = "i" if iotype == "EPI" || iotype == "FORCE"
      op = "j" if iotype == "EPJ"
      ret = Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),op]),fdpsname,type])
    end
    ret
  end
  
  def convert_to_code(conversion_type,h=$varhash)
    s=self
    #print "convert to code #{s}\n"
    #p $varhash[self]
    case conversion_type
    when "reference" then
      # nothing to do
    when /A64FX/
      s = self.convert_to_code_a64fx(h)
    end
    #    print "result=", s, "\n"
    s 
  end
end


def process_iodecl(ios)
  a=[]
  ios.each{|x|
    a +=  [x.name, [x.iotype, x.type, x.fdpsname, x.modifier]]
  }
  #p a
  Hash[*a]
end

def process_funcdecl(func)
  a=[]
  #p func
  func.each{|x|
    #p x
    decl = x.decl
    stmt = x.statements
    ret  = x.retval
    #p decl
    #p stmt
    #p ret
    a += [decl.name, x]
  }
  # reserved name function

  Hash[*a]
end

parser=KernelParser.new
$kernel_name="Kernel"
$epi_name="EPI"
$epj_name="EPJ"
$force_name="Force"
$conversion_type = "reference"
$swpl_stage = 1
$unroll_stage = 1
while true
  opt = ARGV.shift
  break if opt == nil
  case opt
  when "-i"
    filename = ARGV.shift
    warn "input file: #{filename}\n"
  when "--kernel-name"
    $kernel_name = ARGV.shift
    warn "kernel name: #{$kernel_name}\n"
  when "--epi-name"
    $epi_name = ARGV.shift
    warn "epi name: #{$epi_name}\n"
  when "--epj-name"
    $epj_name = ARGV.shift
    warn "epj name: #{$epj_name}\n"
  when "--force-name"
    $force_name = ARGV.shift
    warn "force name: #{$force_name}\n"
  when "--conversion-type"
    $conversion_type = ARGV.shift
    warn "conversion type: #{$conversion_type}\n"
  when "--strip-mining"
    $strip_mining = ARGV.shift
    warn "strip mining size: #{$strip_mining}\n"
  when "--software-pipelining"
    $swpl_stage = ARGV.shift.to_i
    abort "software pipelining is not available"
    warn "software pipelining stage: #{$swpl_stage}\n"
  when "--unroll"
    $unroll_stage = ARGV.shift.to_i
    warn "software pipelining stage: #{$unroll_stage}\n"
  else
    warn "error: unsupported option #{opt}"
    abort
  end
end

src = ""
File.foreach(filename){|line|
  src = src + line
}
program=parser.parse(src)
program.generate_hash("noconversion")
program.expand_function
program.expand_tree
program.make_conditional_branch_block
program.disassemble_statement
program.generate_optimized_code2($conversion_type)
#program.generate_optimized_code($conversion_type)

__END__
