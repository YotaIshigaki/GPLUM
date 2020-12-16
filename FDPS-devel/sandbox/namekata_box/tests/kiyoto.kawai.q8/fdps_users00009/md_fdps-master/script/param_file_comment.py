#!/usr/bin/python
# coding: UTF-8
#--------------------------------------------------------------------------------------------------
# This script use python
#   adding quick-reference comment in ./model/***.param files.
#--------------------------------------------------------------------------------------------------

import subprocess
import os

#--- target file property
target_dir = "./model/"
target_ext = ".param"

#--- added at the head of comment[].
comment_mark = "//"

#--- managed comment and keywords.
#       comment["keyword"] = """ comment """
#       comments are insert in ahead of keywords.
comment = {}
comment["@<PARAM>ATOM"] = """=======================================================================================
  definition of parameters on each atoms.
    atom_name, res_name, mass, vdw_d, vdw_r
      atom_name [-]              string, must be same to "atom_name" in ***.mol2 file.
      res_name  [-]              string, up to 3 characters for pdb files.
      mass      [atomic weight]
      vdw_d     [kcal/mol]       function: V(r) = vdw_d*((vdw_r/r)^12 - 2*(vdw_r/r)^6)
      vdw_r     [angstrom]
======================================================================================="""

comment["@<PARAM>BOND"] = """=======================================================================================
  definition of bond potential.
    i, j, form, r0, k, a
      i, j [-]          string, must be same to "atom_name" in ***.mol2 file.
      form [-]          must be "none", "harmonic", or "anharmonic".
      r0   [angstrom]   equivalent length
      k    [kcal/mol]
      a    [/angstrom]  used in "anharmonic" form.

    form == "none",       free stretching.
    form == "harmonic",   function: V(r) = 0.5*k*(r - r0)^2
    form == "anharmonic", function: V(r) = k*[ar^2 - ar^3 + 7/12*ar^4], ar = a*(r - r0)
======================================================================================="""

comment["@<PARAM>ANGLE"] = """=======================================================================================
  definition of angle potential.
    j, i, k, form, theta0, k
      j, i, k [-]              j-i-k shape, string, must be same to "atom_name" in ***.mol2 file.
      form    [-]              must be "none" or "harmonic". other form is not defined.
      theta0  [degree]         equivalent angle
      k       [kcal/mol·rad^2]

    form == "none",     free rotation.
    form == "harmonic", function: V(phi) = 0.5*k*[cos(phi) - cos(theta0)]^2/[sin(theta0)]^2
======================================================================================="""

comment["@<PARAM>TORSION"] = """=======================================================================================
  definition of torsion potential.
    shape, i, j, k, l, form, v1, v2, v3
      shape       [-]         must be "dihedral" or "improper"
      i, j, k, l  [-]         i-jk-l shape, string, must be same to "atom_name" in ***.mol2 file.
      form        [-]         must be "none", "cos", or "OPLS_3".

    form == "none" case, free rotation.
      v1, v2, v3  [-]         ignored.

    form == "cos" case, CHARMM style.
      v1 = theta0 [degree]    equivalent angle
      v2 = v      [kcal/mol]  function: V(phi) = 0.5*v*(1-cos(n*phi - theta0))
      v3 = n      [integer]   number of local minimum point

    form == "OPLS_3" case, OPLS_AA style.
      v1, v2, v3  [kcal/mol]  function: V(phi) = 0.5*v1*(1+cos(phi)) + 0.5*v2*(1-cos(2*phi)) + 0.5*v3*(1+cos(3*phi))
======================================================================================="""

comment["@<PARAM>SCALING"] = """=======================================================================================
  definition of scaling coefficient for intra-mask.
      scaling_LJ       1-2mask  1-3mask  1-4mask...
      scaling_coulomb  1-2mask  1-3mask  1-4mask...
          [-] numeric. accepts any order length (must be continous from 1-2 level).
======================================================================================="""

#------------------------------------------
#  all settings are written in above.
#------------------------------------------


def add_comment(file_name, mark, comment):
    f         = open(file_name, "r")
    line_file = f.readlines()
    f.close()

    processed_key = []
    result = ""
    buff   = ""
    for line in line_file:
        #--- void line
        if len(line) <= len(mark):
            if not(buff == ""):
                result += buff
                buff = ""
            result += line
            continue
        #--- comment line
        if line[:len(mark)] == mark:
            buff += line
            continue
        #--- valuable line
        else:
            #--- check keyword
            for key in comment.keys():
                #--- find keyword -> insert comment
                if line[:len(key)] == key:
                    result += comment[key]
                    buff = ""  #--- delete old comment
                    processed_key.append(key)
                    break
            #--- return buff comment
            if not(buff == ""):
                result += buff
                buff = ""
            #--- return line
            result += line

    #--- replace new model file
    f = open(file_name, "w")
    f.write(result)
    f.close()

    #--- report process error
    for key in comment.keys():
        if key in processed_key:
            pass
        else:
            print "    WORNING: keyword of '" + key + "' is not found."


def add_mark(str, mark):
    result = ""
    for line in str.split("\n"):
        result += mark + line + "\n"
    return result


if __name__ == "__main__":

    #--- get result of ls
    pwd = os.getcwd()
    os.chdir(target_dir)
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

    ##--- test data
    #print "  find_files:"
    #print model_list
    #model_list = ["AA_Ar"]

    #--- insert comment mark at head of comment lines
    for key in comment.keys():
        comment[key] = add_mark(comment[key], comment_mark)

    #--- update comment process
    for model in model_list:
        print "  update comment: " + model
        file_name = target_dir + model + target_ext
        add_comment(file_name, comment_mark, comment)
