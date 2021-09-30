#!/usr/bin/env /home/x10225/gnu/python3/bin/python3
# -*- coding: utf-8 -*-
#=================================
#   Module import
#=================================
try:
    import os
    import sys
    import struct
    import math
    import re
    import glob
    import subprocess
except ImportError:
    print("Module os,sys,struct,math,re,subprocess are not found.")
    quit()

try:
    import numpy as np
except ImportError:
    print("Module numpy is not found.")
    quit()

#================================
#   Class definition
#================================
class FileConverter:
    def __init__(self):
        pass

    def convert(self,snapshot_num):
        # Get the list of files that must be merged
        pat = "snap{0:05d}-proc*.dat".format(snapshot_num)
        files = glob.glob(pat)
        files.sort()
        #print(files)

        # Set the name of a file for which all data output and open it
        output_file = "snap{0:05d}.txt".format(snapshot_num)
        ofp = open(output_file,"w")
        #print(output_file)
        
        # Read a file 
        nptcl_glb = 0
        for f in files:
            print("reading {}".format(f))
            ifp = open(f,"rb")
            # Read the header information
            try:
                snap_num_next  = struct.unpack("<i",ifp.read(4))[0]
                time_diag      = struct.unpack("<d",ifp.read(8))[0]
                time_snap_next = struct.unpack("<d",ifp.read(8))[0]
                time_sys       = struct.unpack("<d",ifp.read(8))[0]
                nptcl_loc      = struct.unpack("<i",ifp.read(4))[0]
            except:
                print("fail to read the header information!")
                sys.exit()
            # Reset the buffers for the particle data
            id = []
            mass = []
            x = []
            y = []
            z = []
            vx = []
            vy = []
            vz = []
            # Read the particle data
            for i in range(0,nptcl_loc):
                try:
                    id.append(struct.unpack("<q",ifp.read(8))[0])
                    mass.append(struct.unpack("<d",ifp.read(8))[0])
                    x.append(struct.unpack("<d",ifp.read(8))[0])
                    y.append(struct.unpack("<d",ifp.read(8))[0])
                    z.append(struct.unpack("<d",ifp.read(8))[0])
                    vx.append(struct.unpack("<d",ifp.read(8))[0])
                    vy.append(struct.unpack("<d",ifp.read(8))[0])
                    vz.append(struct.unpack("<d",ifp.read(8))[0])
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
        # Output the total number of particles
        print("nptcl_glb = {}".format(nptcl_glb))

#================================
#   Main function
#================================
conv = FileConverter()
conv.convert(0)


