#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=================================
#   Module import
#=================================
try:
    import os
except ImportError:
    print("Modules os is not found.")
    quit()

try:
    import sys
except ImportError:
    print("Module sys is not found.")
    quit()

try:
    import re
except ImportError:
    print("Module re is not found.")
    quit()

try:
    import struct
except ImportError:
    print("Module struct is not found.")
    quit()

try:
    import math
except ImportError:
    print("Module math is not found.")
    quit()

try:
    import glob
except ImportError:
    print("Module glob is not found.")
    quit()

try:
    import string
except ImportError:
    print("Module string is not found.")
    quit()

try:
    import random
except ImportError:
    print("Module random is not found.")
    quit()

try:
    import shutil
except ImportError:
    print("Module shutil is not found.")
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

#try:
#    import numpy as np
#except ImportError:
#    print("Module numpy is not found.")
#    quit()

#================================
#   Class definition
#================================
class FileReader:
    def __init__(self):
        pass

    def reduc_data(self,model_name,snapshot_nums,output_dir):
        for snapshot_num in snapshot_nums:
            # Get the list of files that must be merged
            pat = "./{0}_{1:d}/{0}_{1:d}-*".format(model_name, int(snapshot_num))
            data_files = glob.glob(pat)
            data_files.sort()
            #print(data_files)

            # Set the name of a file for which all data output and open it
            output_file = "{0}/pos_{1:05d}.txt".format(output_dir, int(snapshot_num))
            #print(output_file)

            # Read files
            ofp = open(output_file,"w")
            nptcl_glb = 0
            for f in data_files:
                print("reading {}".format(f))
                ifp = open(f,"rb")
                # Read the header information
                mpi_nproc = 0
                mpi_rank = -1
                nptcl_loc = 0
                nptcl_tot = 0
                omegam = 0.0
                omegab = 0.0
                omegav = 0.0
                omeganu = 0.0
                hubble = 0.0
                anow = 0.0
                tnow = 0.0
                lunit = 0.0
                munit = 0.0
                tunit = 0.0
                try:
                    # See https://docs.python.org/ja/3/library/struct.html
                    mpi_nproc = struct.unpack(">i",ifp.read(4))[0];
                    mpi_rank  = struct.unpack(">i",ifp.read(4))[0];
                    nptcl_loc = struct.unpack(">q",ifp.read(8))[0];
                    nptcl_tot = struct.unpack(">q",ifp.read(8))[0];

                    omegam  = struct.unpack(">f",ifp.read(4))[0];
                    omegab  = struct.unpack(">f",ifp.read(4))[0];
                    omegav  = struct.unpack(">f",ifp.read(4))[0];
                    omeganu = struct.unpack(">f",ifp.read(4))[0];
                    #hubble  = strcut.unpack(">f",ifp.read(4))[0];
                    ifp.read(4) # skip to read `hubble` because a read error happens.
                    anow    = struct.unpack(">f",ifp.read(4))[0];
                    tnow    = struct.unpack(">f",ifp.read(4))[0];
                    
                    lunit = struct.unpack(">d",ifp.read(8))[0];
                    munit = struct.unpack(">d",ifp.read(8))[0];
                    tunit = struct.unpack(">d",ifp.read(8))[0];
                except:
                    print("fail to read the header information!")
                    ifp.close()
                    ofp.close()
                    sys.exit()
                # Reset the buffers for the particle data
                x = []
                vx = []
                y = []
                vy = []
                z = []
                vz = []
                mass = []
                id = []
                # Read the particle data
                for i in range(0,nptcl_loc):
                    try:
                        x.append(struct.unpack(">f",ifp.read(4))[0])
                        vx.append(struct.unpack(">f",ifp.read(4))[0])
                        y.append(struct.unpack(">f",ifp.read(4))[0])
                        vy.append(struct.unpack(">f",ifp.read(4))[0])
                        z.append(struct.unpack(">f",ifp.read(4))[0])
                        vz.append(struct.unpack(">f",ifp.read(4))[0])
                        mass.append(struct.unpack(">f",ifp.read(4))[0])
                        id.append(struct.unpack(">i",ifp.read(4))[0])
                    except:
                        print("fail to read the particle data!")
                        ifp.close()
                        ofp.close()
                        sys.exit()
                # Close the file
                ifp.close()
                # Output
                for i in range(0,nptcl_loc):
                    try:
                        msg = "{X}  {Y}  {Z}\n".format(X=x[i],Y=y[i],Z=z[i])
                        ofp.write(msg)
                    except:
                        print("cannot write into the target file!")
                        ofp.close()
                        sys.exit()
                # Update nptcl_glb
                nptcl_glb += nptcl_loc
            # Close the file 
            ofp.close()
            # Output the total number of particles
            print("nptcl_glb = {}".format(nptcl_glb))


#================================
#   Function definition
#================================
def analyze_CL_args():
    #-------------------------------------------------------------------------------------
    # More information about the module `argparser` can be found at the following sites:
    #  Python 3.5) http://docs.python.jp/3.5/library/argparse.html
    #-------------------------------------------------------------------------------------
    # Create an ArgumentParser object
    description = "Extract the data of positions of particles from output files."
    parser = argparse.ArgumentParser(description=description)
    # Input options
    help_msg = "Model name of simulation"
    parser.add_argument('--model_name', \
                        type=str,required=True, \
                        help=help_msg, metavar='MODEL NAME',
                        dest='model_name')
    help_msg = "A list of snapshot numbers"
    parser.add_argument('--snapshot_nums', nargs="+", \
                        type=int,required=True, \
                        help=help_msg, metavar='SNAPSHOT NUMBER',
                        dest='snapshot_nums')
    # Output option
    help_msg = "The PATH of output directory"
    parser.add_argument('-o','--output','--output_dir', \
                        default="./",type=str,required=False, \
                        help=help_msg, metavar='DIRECTORY', \
                        dest='output_dir')
    # Return the result 
    args = parser.parse_args()
    return args.model_name, args.snapshot_nums, args.output_dir 

    # [DEBUG] Check
    #print("model_name is {}".format(args.model_name))
    #print("snapshot_nums is {}".format(args.snapshot_nums))
    #print("output_dir is {}".format(args.output_dir))


#================================                                                                      
#   Main function                                                                                      
#================================                                                                      
if __name__ == '__main__':
    # Display the version of Python
    print(sys.version_info)
    # Analyze the command-line options and arguments
    model_name, snapshot_nums, output_dir = analyze_CL_args()
    print("model_name is {}".format(model_name))
    print("snapshot_nums is {}".format(snapshot_nums))
    print("output_dir is {}".format(output_dir))

    # Extract the data of positions of particles and output it to a file
    reader = FileReader()
    reader.reduc_data(model_name, snapshot_nums, output_dir)

