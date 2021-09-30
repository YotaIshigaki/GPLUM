#!/usr/bin/python
# coding: UTF-8
#--------------------------------------------------------------------------------------------------
# This script use python
#   "./src/enum_model.hpp" file generator.
#   source: ***.mol2 files in ./model
#--------------------------------------------------------------------------------------------------

import subprocess
import os

target_locate = "./model/"
target_ext    = ".mol2"
result_file   = "./src/enum_model.hpp"

def load_atom_name(model, atom_list):
    model_file = model + target_ext
    f = open(model_file)
    line_file = f.readlines()
    mode = "header"
    for line in line_file:

        #--- skip comment line
        if len(line) <= 2:
            continue
        if line[0] == "!":
            continue
        if line[:2] == "//":
            continue

        #--- convert string to list
        line_list = line.split()

        #--- decode header
        if line_list[0] == "@<TRIPOS>MOLECULE":
            if mode != "header":
                print "ERROR: each ***.mol2 file must have 1 molecule."
                sys.exit()
            mode = "header"
            continue
        if line_list[0] == "@<TRIPOS>ATOM":
            mode = "atom"
            continue
        if line_list[0] == "@<TRIPOS>BOND":
            mode = "bond"
            continue

        #--- decode atom data
        if mode == "atom":
            atom = line_list[1]
            if not (atom in atom_list):
                atom_list.append(atom)

    f.close


#--- enum class definition
def make_enum_class(typename, value_list, result):
    result.append(" "*4 + "enum class " + typename + " : uint_fast32_t {")
    for value in value_list:
        result.append(" "*8 + value + ",")
    result.append(" "*4 + "};")
    result.append("")


#--- map between enum and string
def make_map(typename, value_list, result):
    width = 0
    for value in value_list:
        width = max(len(value), width)

    result.append("namespace {")
    result.append(" "*4 + "static const std::map<std::string, " + typename + "> table_str_" + typename + " {")
    for value in value_list:
        n_space = width - len(value)
        line =  " "*8 + "{"
        line += '"' + value + '"' + " "*n_space
        line += ", "
        line += typename + "::" + value + " "*n_space
        line += "},"
        result.append(line)
    result.append(" "*4 + "};")

    result.append("")

    result.append(" "*4 + "static const std::map<" + typename + ", std::string> table_" + typename + "_str {")
    for value in value_list:
        n_space = width - len(value)
        line =  " "*8 + "{"
        line += typename + "::" + value + " "*n_space
        line += ", "
        line += '"' + value + '"' + " "*n_space
        line += "},"
        result.append(line)
    result.append(" "*4 + "};")

    result.append("}")
    result.append("")


#--- "what(enum value)" function (return std::string)
def make_what(type_list, result):
    result.append("namespace ENUM {")

    for t in type_list:
        result.append(" "*4 + "inline std::string what(const " + t + " &e){")
        result.append(" "*8  + "if(table_" + t + "_str.find(e) != table_" + t + "_str.end()){")
        result.append(" "*12 +     "return table_" + t + "_str.at(e);")
        result.append(" "*8  + "} else {")
        result.append(" "*12 +     "using type_base = typename std::underlying_type<" + t + ">::type;")
        result.append(" "*12 +     'std::cerr << "  ' + t + ': input = " << static_cast<type_base>(e) << std::endl;')
        result.append(" "*12 +     'throw std::out_of_range("undefined enum value in ' + t + '.");')
        result.append(" "*8  + "}")
        result.append(" "*4  + "}")
        result.append("")

    result.append("}")
    result.append("")

    result.append("namespace std {")

    for t in type_list:
        result.append(" "*4 + "inline string to_string(const " + t + " &e){ return ENUM::what(e); }")
        result.append("")

    result.append("}")
    result.append("")


#--- "which_***(std::string)" function (return enum class value)
def make_which(type_list, result):
    result.append("namespace ENUM {")

    for t in type_list:
        result.append(" "*4  + "inline " + t + " which_" + t + "(const std::string &str){")
        result.append(" "*8  +     "if(table_str_" + t + ".find(str) != table_str_" + t + ".end()){")
        result.append(" "*12 +         "return table_str_" + t + ".at(str);")
        result.append(" "*8  +     "} else {")
        result.append(" "*12 +         'std::cerr << "  ' + t + ': input = " << str << std::endl;')
        result.append(" "*12 +         'throw std::out_of_range("undefined enum value in ' + t + '.");')
        result.append(" "*8  +     "}")
        result.append(" "*4  + "}")
        result.append("")

    result.append("}")
    result.append("")


#--- "size_***()" function (return the number of defined value)
def make_size(type_list, result):
    result.append("namespace ENUM {")

    for t in type_list:
        result.append(" "*4 + "inline size_t size_" + t + "(){")
        result.append(" "*8 + "return table_str_" + t + ".size();")
        result.append(" "*4 + "}")
        result.append("")

    result.append("}")
    result.append("")


#--- function of "std::cout << enum class::(value)"
def make_output_ostream(type_list, result):
    for t in type_list:
        result.append(        "inline std::ostream& operator << (std::ostream& s, const " + t + " &e){")
        result.append(" "*4 + "s << ENUM::what(e);")
        result.append(" "*4 + "return s;")
        result.append(        "}")
        result.append("")


#--- "what(std::tuple<...>)" function (return std::string)
def make_whatis_tuple(result):
    #--- string converter: int and double
    result.append("namespace ENUM {")
    result.append(" "*4 + "template<typename T>")
    result.append(" "*4 + "std::string whatis_string_Impl(T const &v){")
    result.append(" "*8 + "return std::to_string(v);")
    result.append(" "*4 + "}")
    result.append("")

    #--- recursive part
    result.append(" "*4  + "template<typename Tuple, size_t Index = std::tuple_size<Tuple>::value-1>")
    result.append(" "*4  + "struct whatis_Impl{")
    result.append(" "*8  + "static void apply(std::string &str, Tuple const &tuple){")
    result.append(" "*12 + "whatis_Impl<Tuple, Index-1>::apply(str, tuple);")
    result.append(" "*12 + 'str += ", " + whatis_string_Impl(std::get<Index>(tuple));')
    result.append(" "*8  + "}")
    result.append(" "*4  + "};")
    result.append("")

    #--- terminator part
    result.append(" "*4  + "template<typename Tuple>")
    result.append(" "*4  + "struct whatis_Impl<Tuple, 0>{")
    result.append(" "*8  + "static void apply(std::string &str, Tuple const &tuple){")
    result.append(" "*12 + "str = whatis_string_Impl(std::get<0>(tuple));")
    result.append(" "*8  + "}")
    result.append(" "*4  + "};")
    result.append("")

    #--- interface
    result.append(" "*4 + "template<typename Tuple>")
    result.append(" "*4 + "inline std::string what(Tuple const &tt){")
    result.append(" "*8 + "std::string str{""};")
    result.append(" "*8 + "whatis_Impl<Tuple>::apply(str, tt);")
    result.append(" "*8 + 'return "(" + str + ")";')
    result.append(" "*4 + "}")
    result.append("}")
    result.append("")


#--- specialize std::hash() for enum class (C++14)
#      support for GCC 4.x ~ 5.x. naturally supported by GCC 6.0 or later.
#      ref: http://qiita.com/taskie/items/479d649ea1b20bacbe03
def make_hash_func(type_list, result):
    result.append(         "namespace std {")

    for t in type_list:
        result.append(" "*4  + "template <>")
        result.append(" "*4  + "struct hash<" + t + "> {")
        result.append(" "*8  + "size_t operator() (" + t + " x) const noexcept {")
        result.append(" "*12 + "using type = typename underlying_type<" + t +">::type;")
        result.append(" "*12 + "return hash<type>{}(static_cast<type>(x));")
        result.append(" "*8  + "}")
        result.append(" "*4  + "};")
        result.append(         "")

    result.append(         "}")
    result.append(         "")


if __name__ == "__main__":

    #--- get result of ls
    pwd = os.getcwd()
    os.chdir(target_locate)
    cmd = "ls"
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    os.chdir(pwd)
    ls_result = p.communicate()[0]
    ls_files = ls_result.split()
    #print ls_files

    #--- get target file name
    model_list = []
    for f in ls_files:
        if f[(len(f)-len(target_ext)):] == target_ext:
            model_list.append(f[:(len(f)-len(target_ext))])

    #--- import all atom names
    atom_list = []
    print "  model list : " + str(model_list)
    for model in model_list:
        load_atom_name(target_locate + model, atom_list)

    print "  atom list  : " + str(atom_list)


    #--- output C++ header file (definition of enum class)
    result = []

    #------ file header
    result.append("//" + "-"*88)
    result.append("//  This file is enum class of atom types and molecular types.")
    result.append("//    generated by ./script/convert_model_indicator.py")
    result.append("//" + "-"*88)
    result.append("#pragma once")
    result.append("")
    result.append("#include <cstdint>")
    result.append("#include <string>")
    result.append("#include <map>")
    result.append("#include <tuple>")
    result.append("#include <stdexcept>")
    result.append("")
    result.append("")

    #------ define enum class
    make_enum_class("MolName", model_list, result)
    make_enum_class("AtomName", atom_list, result)

    #------ map data between string and enum class
    make_map("MolName", model_list, result)
    make_map("AtomName", atom_list, result)

    type_list = []
    type_list.append("MolName")
    type_list.append("AtomName")

    #------ interface
    result.append('//--- basic interface for enum class')
    make_what(  type_list, result)
    make_which( type_list, result)
    make_size(  type_list, result)

    #result.append('//--- convert std::tuple<...> to std::string')
    #make_whatis_tuple(result)

    #--- add ostream interface
    result.append('//--- output function as "std::cout << (enum class::value)"')
    make_output_ostream(type_list, result)
    result.append("")

    #--- add hash function
    result.append('//--- specialized hash function')
    result.append('//        support for gcc 4.x ~ 5.x. naturally supported by gcc 6.0 or later. (C++14 standard)')
    result.append('//        ref: http://qiita.com/taskie/items/479d649ea1b20bacbe03')
    make_hash_func(type_list, result)


    #--- output file
    print "  output file: " + result_file
    f = open(result_file, "w")
    for line in result:
        f.write(line + "\n")
    f.close
