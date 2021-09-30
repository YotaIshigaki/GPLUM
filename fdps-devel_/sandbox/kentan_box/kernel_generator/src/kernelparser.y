# coding: utf-8
class KernelParser
  prechigh
  left '.'
  nonassoc UMINUS
  left '*' '/'
  left '+' '-'
  left '<' '>' '<=' '>='
  left '&' '|'
  left '==' '!='
  left '+=' '-='
    
  preclow
  start innerkernel

  rule

  innerkernel: iodeclarations functions statements {result=Kernelprogram.new(val)}
             | iodeclarations statements           {result=Kernelprogram.new(val)}
  iodeclarations: iodeclaration
                | iodeclarations iodeclaration {result = val[0]+val[1]}
  iodeclaration: iotype type varname ':' fdpsname EOL          {result = [Iodeclaration.new([val[0],val[1],val[2],val[4],nil])]}
               | iotype "static" type varname ':' fdpsname EOL {result = [Iodeclaration.new([val[0],val[2],val[3],val[5],"static"])]}
               | iotype "local" type varname EOL { result = [Iodeclaration.new([val[0],val[2],val[3],nil,"local"])]}
               | declaration

  functions : function
            | function functions {result = val[0]+val[1]}
  function : funcdeclaration statements ret_state end_state {result = [Function.new([val[0],val[1],val[2]])]}
           | funcdeclaration ret_state end_state            {result = [Function.new([val[0],[],val[1]])]}
  ret_state : "return" expression EOL {result = ReturnState.new(val[1])}
  end_state : "end" EOL
  funcdeclaration : "function" TEXT '(' operands ')' EOL {result = Funcdeclaration.new([val[1],val[3]])}
                  | "function" func_attr TEXT '(' operands ')' EOL {result = Funcdeclaration.new([val[2],val[4],val[1]])}
  func_attr : "__ppa__"

  operands : operand
           | operands ',' operand {result = val[0] + val[2]}
  operand : varname {result = [val[0]]}

  declaration: type varname EOL {result = [Iodeclaration.new(["MEMBER",val[0],val[1],nil,nil])]}

  iotype : EPI
         | EPJ
         | FORCE

  type : F64VEC
       | F64
       | F32VEC
       | F32
       | F16VEC
       | F16
       | S64
       | U64
       | S32
       | U32
       | S16
       | U16

  varname : TEXT

  fdpsname : TEXT

  statements : statement
             | statement statements {result = val[0]+val[1]}

  statement : var "="  expression EOL {result = [Statement.new([val[0],val[2]])]}
            | var '+=' expression EOL {result = [Statement.new([val[0],Expression.new([:plus, val[0], val[2]])])]}
            | var '-=' expression EOL {result = [Statement.new([val[0],Expression.new([:minus,val[0], val[2]])])]}
            | var '*=' expression EOL {result = [Statement.new([val[0],Expression.new([:mult,val[0], val[2]])])]}
            | var '/=' expression EOL {result = [Statement.new([val[0],Expression.new([:div,val[0], val[2]])])]}
            | pragma EOL {result = [val[0]]}
            | "if" expression EOL    {result = [IfElseState.new([:if,val[1]])]}
            | "else" EOL             {result = [IfElseState.new([:else,nil])]}
            | "elsif" expression EOL {result = [IfElseState.new([:elsif,val[1]])]}
            | "endif" EOL            {result = [IfElseState.new([:endif,nil])]}

  pragma: "#pragma" TEXT {result = Pragma.new([val[1],nil])}
        | "#pragma" TEXT options {result = Pragma.new([val[1],val[2]])}
  options: option { result = val[0]}
         | options option { result = val[0] + val[1] }
  option: TEXT { result = [val[0]]}

  expression: binary

  binary: binary '+'  binary {result = Expression.new([:plus,  val[0], val[2]])}
        | binary '-'  binary {result = Expression.new([:minus, val[0], val[2]])}
        | binary '*'  binary {result = Expression.new([:mult,  val[0], val[2]])}
        | binary '/'  binary {result = Expression.new([:div,   val[0], val[2]])}
        | binary '==' binary {result = Expression.new([:eq,    val[0], val[2]])}
        | binary '!=' binary {result = Expression.new([:neq,   val[0], val[2]])}
        | binary '&'  binary {result = Expression.new([:and,   val[0], val[2]])}
        | binary '|'  binary {result = Expression.new([:or,    val[0], val[2]])}
        | binary '>'  binary {result = Expression.new([:gt,    val[0], val[2]])}
        | binary '<'  binary {result = Expression.new([:lt,    val[0], val[2]])}
        | binary '>=' binary {result = Expression.new([:ge,    val[0], val[2]])}
        | binary '<=' binary {result = Expression.new([:le,    val[0], val[2]])}
        | '(' binary ')'     {result = val[1]}
        | funccall           {result = val[0]}
	| '-' binary =UMINUS {result = Expression.new([:uminus,val[1], nil])}
        | var

  funccall: TEXT '(' args ')' {result = FuncCall.new([val[0], val[2]])}
  args: arg
      | args ',' arg {result = val[0] + val[2]}
  arg: binary {result = [val[0]]}

  var: TEXT                {result = val[0]}
     | TEXT '.' TEXT       {result = Expression.new([:dot,val[0],val[2]])}
     | TEXT '[' TEXT ']' {result = Expression.new([:array,val[0],val[2]])}
end

---- inner
def vartype(string)
[[:F64VEC,:F64,:F32VEC,:F32,:F16VEC,:F16,:S64,:U64,:S32,:U32,:S16,:U16][["F64vec","F64","F32vec","F32","F16vec","F16","S64","U64","S32","U32","S16","U16"].index(string)],string]
end
def iotype(string)
  [[:EPI,:EPJ,:FORCE][["EPI","EPJ","FORCE"].index(string)], string]
end
def dim(string)
  [[:x,:y,:z,:xx,:yy,:zz,:xy,:xz,:yz,:yx,:zx,:zy][["x","y","z","xx","yy","zz","xy","xz","yz","yx","zx","zy"].index(string)],string]
end

def parse(program)
  @q=[]
  program.each_line{|str|
#    p str.chomp
    a=str.chomp.split(/(\s|:|\=\=|\!\=|\+\=|-\=|\*\=|\&|\||\+|-|\*|\/|\=|\(|\)|,|\.|>\=|<\=|>|<|\[|\])/).select{|s| s=~/\S+/}
    symbols = /(\=\=|\!\=|\&|\||\+|-|\*|\/|\(|\)|,|\.|>\=|<\=|>|<|\[|\])/
    if a == []
      next
    end
    if a[0] =~ /(EPI|EPJ|FORCE|MEMBER|ACCUM)/
      io = a[0]
      @q << iotype(a.shift) # iotype
    end
    if a[0] == "static"
      @q << ["static",a.shift]
    elsif a[0] == "local"
      @q << ["local",a.shift]
    end

    if a[0] =~ /(F(64|32|16)|F(64|32|16)vec|(U|S)(64|32|16))/
      @q << vartype(a.shift) # type
      @q << [:TEXT, a.shift] # varname
      if io != nil && a[0] == ":"
        @q << [":",a.shift] # :
        @q << [:TEXT, a.shift] # fdpsname
      end
    elsif a[0] == "function"
    #      print "function decl"
      @q << ["function", a.shift]
      @q << ["__ppa__",a.shift] if a[0] == "__ppa__"
      @q << [:TEXT, a.shift]
      @q << ["(", a.shift]
      a.each{|x|
        if x =~ /(F(64|32|16)|F(64|32|16)vec|(U|S)(64|32|16))/
          @q << vartype(x)
        elsif x  == ','
          @q << [',', x]
        elsif x  == ')'
          @q << [')', x]
        else
          @q << [:TEXT, x]
        end
      }
    elsif a[0] == "return"
      @q << ["return", a.shift]
      a.each{ |x|
        if x =~ symbols
          @q << [x, x]
        else
          @q << [:TEXT, x]
        end
      }
    elsif a[0] == "end"
      @q << ["end", a.shift]
    elsif a[0] == "#pragma"
      #p a
      @q << ["#pragma",a.shift]
      a.each{|x|
	@q << [:TEXT,x]
      }
    elsif a[0] == "if"
      @q << ["if","if"]
      a.shift
      a.each{|x|
	if x =~ symbols
	  @q << [x,x]
	else
	  @q << [:TEXT,x]
	end
      }
    elsif a[0] == "elsif"
      @q << ["elsif","elsif"]
      a.shift
      a.each{|x|
	if x =~ symbols
	  @q << [x,x]
	else
	  @q << [:TEXT,x]
	end
      }
    elsif a[0] == "else"
      @q << ["else","else"]
    elsif a[0] == "endif"
      @q << ["endif","endif"]
    elsif a[1] =~ /(\=|\+\=|-\=)/
      #print "statement \n"
      @q << [:TEXT,a.shift]
      #p a
      @q << [a[0], a[0]]
      a.shift
      a.each{|x|
        if x =~ symbols
          @q << [x,x]
        else
          @q << [:TEXT,x]
        end
      }
    elsif a[3] =~ /(\=|\+\=|-\=)/
      #print "statement \n"
      @q << [:TEXT,a.shift] # var
      @q << [a[0], a[0]]    #.
      a.shift
      @q << [:TEXT,a.shift]
      @q << [a[0], a[0]]    # =|+=|-=|*=|/=
      a.shift
      a.each{|x|
        if x =~ symbols
          @q << [x,x]
        else
          @q << [:TEXT,x]
        end
      }
    else
      abort "error: unsupported DSL description \"#{str.chomp}\""
    end
    @q << [:EOL,:EOL]
  }
  @q << nil
  #p @q
  do_parse
end

def next_token
  @q.shift
end

