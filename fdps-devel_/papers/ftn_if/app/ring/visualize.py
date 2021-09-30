#!/usr/bin/env /opt/local/bin/python3.6
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
    import subprocess
except ImportError:
    print("Module os,sys,struct,math,re,subprocess are not found.")
    quit()

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import matplotlib.cm as cm
except ImportError:
    print("Module matplotlib.pyplot,patches,cm are not found.")
    quit()

try:
    import numpy as np
except ImportError:
    print("Module numpy is not found.")
    quit()

#================================
#   Class definition
#================================
class Plotter:
    def __init__(self,xmin,xmax,ymin,ymax):
        km = 1.0e5
        self.xmin = xmin * km
        self.xmax = xmax * km
        self.ymin = ymin * km
        self.ymax = ymax * km
        self.x = []
        self.y = []

    def read_file(self,filename,skip_freq):
        # Read a file
        print("reading {0} (skip_freq = {1})".format(filename,skip_freq))
        fp = open(filename,"r")
        data_num = 0
        for line in fp:
            items = line.split()
            try:
                x = float(items[0])
                y = float(items[1])
                z = float(items[2])
                if ((self.xmin <= x) and (x <= self.xmax) and \
                    (self.ymin <= y) and (y <= self.ymax)):
                    if skip_freq == 0:
                        self.x.append(x)
                        self.y.append(y)
                    else:
                        data_num += 1
                        if (data_num == skip_freq):
                            self.x.append(x)
                            self.y.append(y)
                            data_num = 0
                            
            except:
                print("cannot convert to FP values!")
                sys.exit()
        fp.close()
        print("{} doubles are read.".format(2*len(self.x)))

        # Convert to numpy-format
        self.x = np.array(self.x)
        self.y = np.array(self.y)
        print("converted to numpy-format.")

        # Convert physical unit
        km = 1.0e5
        self.xmin = self.xmin/km
        self.xmax = self.xmax/km
        self.ymin = self.ymin/km
        self.ymax = self.ymax/km
        self.x = np.divide(self.x,km)
        self.y = np.divide(self.y,km)
        print("normalization completed.")

    def make_fig(self,basename,dpi,patch,save_fmt='ps'):
        # Set file names
        ps_file   = basename + '.ps'
        eps_file  = basename + '.eps'
        pdf_file  = basename + '.pdf'
        png_file  = basename + '.png'
        
        # Make a figure
        font = {'family' : 'Verdana',
                'weight' : 'normal',
                'size'   : '14'}
        plt.rc('font',**font)
        plt.rc('text',usetex=True)
        fig = plt.figure(1,figsize=(5,20)) # (b)
        #fig = plt.figure(1,figsize=(10,10)) # (a)(c)
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        
        # Plot column density of the reference models
        ax.set_xlabel(r'$x\, [\mathrm{km}]$',fontsize=20)
        ax.set_ylabel(r'$y\, [\mathrm{km}]$',fontsize=20)
        ax.tick_params(axis='x',labelsize=24)
        ax.tick_params(axis='y',labelsize=24)
        ax.set_xlim(self.xmin,self.xmax)
        ax.set_ylim(self.ymin,self.ymax)
        x_vals = np.resize(self.x,(np.size(self.x),1))
        y_vals = np.resize(self.y,(np.size(self.y),1))
        #ax.scatter(x_vals, y_vals, s=0.005, c='k',edgecolors='none') # (a)
        ax.scatter(x_vals, y_vals, s=0.25, c='k',edgecolors='none') # (b)
        #ax.scatter(x_vals, y_vals, s=1, c='k',edgecolors='none') # (c)
        #if (patch is not None):
        #    ax.add_patch(patch)
        ax.add_patch(patch)
        
        # Display figure 
        fig.tight_layout()

        # Save figure
        if (save_fmt == 'ps'):
            fig.savefig(ps_file)
            # Make eps,pdf,bb files
            cmd = 'ps2eps -B -l -g -a -f ' + ps_file
            subprocess.call(cmd,shell=True)
            cmd = 'eps2pdf -B -H -f ' + eps_file
            subprocess.call(cmd,shell=True)
            cmd = 'extractbb ' + pdf_file
            subprocess.call(cmd,shell=True)
        elif save_fmt == 'png':
            fig.savefig(png_file,dpi=dpi)
            cmd = 'extractbb ' + png_file
            subprocess.call(cmd,shell=True)
        

#================================
#   Main function
#================================
# Settings
figtype = 'b'
if figtype == 'a':
    # (a) entire picture
    skip_freq = 10
    basename = "./snapshots/snap00010a"
    dpi = 1600
    xmin = - 400.0;  xmax   = 400.0
    ymin = - 400.0;  ymax   = 400.0
    xmin_p = 384.0;  xmax_p = 396.0
    ymin_p = - 30.0; ymax_p = 30.0
    rect = patches.Rectangle(
               (xmin_p, ymin_p),   # (x,y)
               xmax_p-xmin_p,      # width
               ymax_p-ymin_p,      # height
               facecolor="none",   # background color
               edgecolor='r',      # border color
               linewidth=1         # border width
           )
elif figtype == 'b':
    # (b) ring segment
    skip_freq = 0
    basename = "./snapshots/snap00010b"
    dpi = 800
    xmin   = 384.0;  xmax   = 396.0
    ymin   = - 30.0; ymax   = 30.0
    xmin_p = 386.0;  xmax_p = 395.0
    ymin_p = - 4.5;  ymax_p = 4.5
    rect = patches.Rectangle(
               (xmin_p, ymin_p),   # (x,y)
               xmax_p-xmin_p,      # width
               ymax_p-ymin_p,      # height
               facecolor="none",   # background color
               edgecolor='r',      # border color
               linewidth=1         # border width
           )
elif figtype == 'c':
    # (c) most-enlarged view
    skip_freq = 0
    basename = "./snapshots/snap00010c"
    dpi = 800
    xmin = 386.0; xmax = 395.0
    ymin = - 4.5; ymax =   4.5
    rect = None

# 
P = Plotter(xmin,xmax,ymin,ymax)
filename = "./snapshots/snap00010.txt"
P.read_file(filename,skip_freq)
P.make_fig(basename, dpi, patch=rect, save_fmt='png')


