# coding: utf-8
class KernelParser
  prechigh
  nonassoc UMINUS
  left '*' '/'
  left '+' '-'
  preclow
  start innerkernel
  rule
  innerkernel: nsdeclarations iodeclarations statements {result=Kernelprogram.new(val)}
  nsdeclarations:nsdeclaration 
                |nsdeclaration nsdeclarations {result = val[0]+val[1]}
  iodeclarations:iodeclaration 
                 |iodeclaration iodeclarations {result = val[0]+val[1]}
   statements: statement
              | statement statements {result = val[0]+val[1]}
   iodeclaration: iotype type varname fdpsname EOL {result = [Iodeclaration.new(val)]}
                | iotype type varname EOL {result=[Iodeclaration.new(val[0..2]+[val[2]])]}
   iotype : IIN 
           | JIN
	   | IOUT
   nsdeclaration: nstype varname EOL {result = [Nsdeclaration.new(val)]}
   nstype : EPINAME
          | EPJNAME
          | FNAME
   type  : F64VEC
           | F64
	   | I64
	   | I32
    varname : TEXT
    fdpsname : TEXT
    statement : TEXT "=" expression EOL {result = [Statement.new([nil,val[0],val[2]])]}
              | type  TEXT "=" expression EOL{result = [Statement.new([val[0],val[1],val[3]])]}
     expression: expression '+' expression {result=Expression.new([:plus, val[0], val[2]])}
    | expression '-' expression           {result=Expression.new([:minus, val[0], val[2]])}
    | expression '*' expression           {result=Expression.new([:mult, val[0], val[2]])}
     | expression '/' expression        {result=Expression.new([:div, val[0], val[2]])}
     | '(' expression ')'               {result=val[1]}
     | TEXT '(' expression ')'          {result=Expression.new([:func, val[0], val[2]])}
     | '-' TEXT  =UMINUS                {result=Expression.new([:uminus, val[0],nil])}
     | TEXT                             {result=val[0]}
end
---- inner
require 'ripper'

def vartype(string)
  [[:F64VEC,:F64,:F32VEC,:F32,:I64, :I32][["f64vec","f64","f32vec","f32","i64","i32"].index(string)],string]
end
def iotype(string)
  [[:IIN,:JIN,:IOUT][["iin","jin","iout"].index(string)], string]
end
def nstype(string)
  [[:EPINAME, :EPJNAME, :FNAME][["epi","epj","fi"].index(string)], string]
end

def parse(program)
  @q=[]
  program.each_line{|str|
    #    p str
    a=Ripper.tokenize(str.chomp).select{|s| s=~/\S+/}
    #   p a
    if a[0] =~ /(iin|jin|iout)/
      #      print "vardecl\n"		    
      io=a[0]

      type=a[1]
      varname=a[2]
      fdpsname=nil
      fdpsname=a[3] if a[3]		    
      @q <<iotype(io)
      @q << vartype(type)
      @q<< [:TEXT, varname]
      @q<< [:TEXT, fdpsname] if fdpsname
      @q << [:EOL,:EOL]
    end
    if a[0] =~ /(epi|epj|fi)/
#      print "iodecl\n"		    
      ns=a[0]
      nsname=a[1]
      @q <<nstype(ns)
      @q<< [:TEXT, nsname]
      @q << [:EOL,:EOL]
    end
    if a[1]=="=" ||a[2]=="="
#      print "statement \n"		    
      if a[2]=="="
        type=a[0]
        @q << vartype(type)
        a.shift
      end
#      p a
      @q << [:TEXT,a.shift]      
      a.shift
#      p a
       @q << ['=', '=']
       symbols="()+-*/"
       a.each{|x|
       if symbols.index(x)
         @q << [x,x]
       else
           @q << [:TEXT,x]
       end
     }
      @q << [:EOL,:EOL]
    end
   }
  p @q
  do_parse
end
def next_token
  @q.shift
end
