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
except ImportError:
    print("Module os,sys,struct are not found.")
    quit()

try:
    import argparse 
except ImportError:
    print("Module argparse is not found.")
    quit()

try:
    import subprocess
except ImportError:
    print("Module subprocess is not found.")
    quit()

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
    description = "Get the number of nodes described in a job script."
    parser = argparse.ArgumentParser(description=description)
    # Input options
    help_msg = "The PATH of job script file"
    parser.add_argument('job_script', nargs=1, \
                        type=str, \
                        help=help_msg, metavar='FILE')
    help_msg = "The PATH of original data file"
    parser.add_argument('data_file', nargs=1, \
                        type=str, \
                        help=help_msg, metavar='FILE')
    # Output option
    help_msg = "The PATH of output directory"
    parser.add_argument('output_dir', nargs=1, \
                        type=str, \
                        help=help_msg, metavar='FILE')
    # Return the result 
    args = parser.parse_args()
    return args.job_script[0], args.data_file[0], args.output_dir[0]

#================================
#   Main function
#================================
if __name__ == '__main__':
    job_script, data_file, output_dir = analyze_CL_args()

    # [1] Extract # of process from the job script file
    n_proc = -1
    fp = open(job_script,'r')
    pat = re.compile(r"^#PJM\s+--rsc-list\s+\"node=[0-9]+\"")
    for line in fp:
        if (pat.match(line)):
            val = re.match(r"^#PJM\s+--rsc-list\s+\"node=(?P<node>[0-9]+)\"",line).group('node')
            n_proc = int(val)
            break
    fp.close()
    if (n_proc <= 0):
        print("n_proc is invalid!")
        sys.exit(1)
    elif (n_proc == 1):
        # In this case, no split is needed. So, just copy.
        output_file = output_dir + "/ptcl_data_00000"
        cmd = "cp {0} {1}".format(data_file, output_file)
        subprocess.call(cmd,shell=True)
    else:
        # [2] Calculate byte number used to split file
        PTCL_SIZE = 66 # in bytes
        data_file_size = os.path.getsize(data_file)
        n_ptcl_tot = data_file_size / PTCL_SIZE
        n_ptcl_loc = n_ptcl_tot / n_proc + 1
        byte_num = n_ptcl_loc * PTCL_SIZE

        # [3] Split the original data and save them to the specified directory
        if (not os.path.exists(output_dir)):
            os.mkdir(output_dir)
        base_name = output_dir + "/ptcl_data_"
        cmd = "split -b {0} --numeric-suffixes --suffix-length=5 {1} {2}".format(byte_num,data_file,base_name)
        subprocess.call(cmd,shell=True)



