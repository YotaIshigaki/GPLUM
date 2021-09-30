def isVector(val)
  if val.get_type =~ /vec/
    true
  else
    false
  end
end

def isVariable(val)
  if val.class == String
    true
  else
    false
  end
end

def isLeaf(val)
  ret = false
  if val.class == String || val.class == FuncCall
    ret = true
  else
    if val.class == Expression
      if val.operator == :dot || val.operator == :array
        ret = true
      end
    end
  end
  ret
end

def isStatement(val)
  if val.class == Statement
    true
  else
    false
  end
end

def get_name(x)
  ret = nil
  if x.class == Expression
    if x.operator == :dot || x.operator == :array
      ret = get_name(x.lop)
    else
      abort "get_name for not :dot epxpression"
    end
  elsif x.class ==  Statement
    if x.name.class == Expression
      ret = get_name(x.name)
    else
      ret = x.name
    end
  elsif x.class == Duplicate
    ret = get_name(x.name)
  elsif x.class == Declaration
    ret = get_name(x.name)
  elsif x.class == LoadState
    ret = get_name(x.dest)
  elsif x.class == StoreState
    ret = get_name(x.src)
  elsif x.class == PointerOf
    ret = get_name(x.exp)
  elsif x.class == String
    ret = x
  else
    abort "get_name is not allowed to use for #{x.class}"
  end
  ret
end

def get_tail(x)
  ret = nil
  if x.class == Expression
    if x.operator == :dot
      ret = x.rop
    else
      abort "get_name for not :dot epxpression"
    end
  elsif x.class == Statement
    if x.name.class == Expression
      if x.name.operator == :dot
        ret = x.name.lop
      else
        ret = x.name
      end
    else
      ret = nil
    end
  elsif x.class == String
    ret = nil
  else
    abort "get_name is not allowed to use for #{x.class}"
  end
  ret
end

def get_declare_type(type,conversion_type)
  case conversion_type
  when "reference"
    case type
    when "S32"
      decl = "PS::S32"
    when "S64"
      decl = "PS::S64"
    when "U32"
      decl = "PS::U32"
    when "U64"
      decl = "PS::U64"
    when "F32"
      decl = "PS::F32"
    when "F64"
      decl = "PS::F64"
    when "F32vec"
      decl = "PS::F32vec"
    when "F64vec"
      decl = "PS::F64vec"
    else
      abort "error: unsupported type #{type} for get_declare_type"
    end
  when /A64FX/
    decl = get_declare_type_a64fx(type)
  else
    abort "error: unsupported conversion_type #{conversion_type}"
  end
  decl
end

def get_iotype_array(iotype)
  case iotype
  when "EPI"
    "epi"
  when "EPJ"
    "epj"
  when "FORCE"
    "force"
  else
    abort "error: undefined iotype #{iotype}"
  end
end
def get_io_index(iotype)
  case iotype
  when "EPI"
    "i"
  when "EPJ"
    "j"
  when "FORCE"
    "i"
  else
    abort "error: undefined iotype #{iotype}"
  end
end

def get_simd_width(conversion_type)
  case conversion_type
  when "reference"
    1
  when /A64FX/
    16
  end
end


def add_new_tmpvar(type,h=$varhash)
  ret = "__fkg_tmp#{$tmp_val_count}"
  abort "type of #{ret} is nil at add_new_tmpvar!" if type == nil
  h[ret] = [nil,type,nil]
  $tmp_val_count += 1
  ret
end
