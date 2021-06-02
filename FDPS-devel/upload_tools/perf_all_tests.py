#!/usr/bin/env python
# -*- coding: utf-8 -*-
#=================================
#   Module import
#=================================
try:
    import os
    import sys
    import re
    import struct
    import math
except ImportError:
    print("Module os,sys,struct,math are not found.")
    quit()

try:
    import argparse 
except ImportError:
    print("Module argparse is not found.")
    quit()

try:
    import itertools
except ImportError:
    print("Module itertools is not found.")
    quit()

try:
    import collections
except ImportError:
    print("Module collections is not found.")
    quit()

try:
    import textwrap
except ImportError:
    print("Module textwrap is not found.")
    quit()

try:
    import subprocess
except ImportError:
    print("Module subprocess is not found.")
    quit()


#================================
#   Class definition
#================================

#================================
#   Function definition
#================================
def analyze_CL_args():
    #-------------------------------------------------------------------------------------
    # More information about the module `argparser` can be found at the following sites:
    #  Python 3.5) http://docs.python.jp/3.5/library/argparse.html
    #  Python 2.x) http://docs.python.jp/2/library/argparse.html
    #-------------------------------------------------------------------------------------
    # Create an ArgumentParser object
    description = "Perform all the tests of the sample codes of FDPS." 
    parser = argparse.ArgumentParser(description=description)
    # Output option
    help_msg = "The PATH of the root directory of FDPS"
    parser.add_argument('fdps_root_dir', \
                        default="./",type=str, \
                        help=help_msg, metavar='DIRECTORY')
    # Return the result
    args = parser.parse_args()
    return args.fdps_root_dir 

def get_red_text(msg):
    return "\033[31;1m" + msg + "\033[m" # red

#================================
#   Main function
#================================
if __name__ == '__main__':
    # Analyze the command-line options and arguments
    fdps_root_dir = analyze_CL_args()
    cpp_samples_root_dir = fdps_root_dir + "/sample/c++"
    ftn_samples_root_dir = fdps_root_dir + "/sample/fortran"

    # Get the list of directories where the sample codes are placed
    target_dirs = []
    for root_dir in [cpp_samples_root_dir, ftn_samples_root_dir]:
        files = os.listdir(root_dir)
        for f in files:
            # [temporarily;start]
            if f == "nbody+sph" or f == "pmmm":
                continue
            # [temporarily;end]
            dir_path = os.path.join(root_dir, f)

            if os.path.isdir(dir_path):
                target_dirs.append(dir_path)
    # [debug]
    #print(target_dirs)

    # Perform the tests
    ret = 0
    current_path = os.getcwd()
    for d in target_dirs:
        # Change the working directory
        os.chdir(current_path)
        os.chdir(d)
        # Perform test.py
        print("Testing ... {0}".format(d))
        cmd = "./test.py; echo $?"
        outs,errs = subprocess.Popen(cmd, \
                                     stdout=subprocess.PIPE, \
                                     stderr=subprocess.PIPE, \
                                     shell=True).communicate()
        print("[STDOUT ({0})]\n {1}\n".format(d,outs))
        print("[STDERR ({0})]\n {1}\n".format(d,errs))
        # Extract the return value
        rettmp = re.findall(r"[+-]?\d+\n$",outs)
        rettmp = rettmp[-1]
        rettmp = rettmp.rstrip()
        rettmp = rettmp.strip()
        rettmp = rettmp.strip("\t")
        # Check the results
        if (rettmp != "0"): 
            msg = get_red_text("test.py in {} is failed.".format(d))
            print("{0}".format(msg))
            print("rettmp = {0}".format(rettmp))
            ret += 1
        # [DEBUG]
        #sys.exit(0)
 
    # Return 
    if (ret == 0):
        print("All the tests are completed!")
        sys.exit(ret)
    else:
        msg = get_red_text("Some test is failed! Please check STDOUT.")
        print("{0}".format(msg))
        sys.exit(-1)
