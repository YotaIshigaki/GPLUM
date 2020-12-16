#!/usr/bin/env /opt/local/bin/python3.5
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
    import matplotlib
    import matplotlib.pyplot as mpyplot
    import matplotlib.scale  as mscale
    import matplotlib.gridspec as gridspec
except ImportError:
    print("Module matplotlib.pyplot,scale is not found.")
    quit()

try:
    import numpy as np
except ImportError:
    print("Module numpy is not found.")
    quit()

try:
    import PhysicalConstants as constants
except ImportError:
    print("Module PhysicalConstants not found")
    quit()

#================================
#   Class definition
#================================
class MeasuredData:
    def __init__(self):
        self.nptcl = []
        self.wtcpp = []
        self.wtftn = []
        self.rdiff = []

    def read_file(self,file_name):
        # [1] Read a file
        print("reading {0}".format(file_name))
        fp = open(file_name,"r")
        for line in fp:
            items = line.split()
            if (len(items) != 0):
               ret = re.match(r'^#',items[0]) # To skip comments
               if ret != None:
                   continue
               try:
                   self.nptcl.append(float(items[0]))
                   self.wtcpp.append(float(items[1]))
                   self.wtftn.append(float(items[2]))
               except:
                   print("cannot convert to FP values!")
                   sys.exit()
        fp.close()

        # [2] Convert to numpy-format
        self.nptcl = np.array(self.nptcl)
        self.wtcpp = np.array(self.wtcpp)
        self.wtftn = np.array(self.wtftn)
        tmp = np.multiply([0.5],np.add(self.wtcpp,self.wtftn))
        tmp = np.divide(np.subtract(self.wtcpp,self.wtftn),tmp)
        self.rdiff = np.absolute(tmp)
        print("converted to numpy-format")

    def plot_panel(self,ax1,title='',panel_name='', \
                   with_xlabel=True,with_y1label=True,with_y2label=True):
        # Create ax2
        ax2 = ax1.twinx()
        
        # Set the title
        if (title != ''):
            ax1.set_title(title,fontsize=16)

        # Set the panel name
        if (panel_name != ''):
            x = 6.0e2
            y = 2.0e0
            ax1.text(x,y,panel_name,fontsize=10)

        # Set the labels, etc.
        if (with_xlabel == True):
            ax1.set_xlabel('Number of particles',fontsize=16)

        if (with_y1label == True):
            ax1.set_ylabel('Wall clock time per steps [s]',fontsize=16)
        else:
            ax1.tick_params(axis='y',labelleft=False)

        if (with_y2label == True):
            ax2.set_ylabel('Relative difference',fontsize=16)
        else:
            ax2.tick_params(axis='y',labelright=False)

        # Set the range of the axes
        ax1.set_xlim(5.12e2,5.24288e5)
        ax1.set_ylim(1.0e-3,1.0e1)
        ax2.set_ylim(1.0e-4,1.0e0)

        # TODO: investigate why the following setting does not work as expected.
        #major_ticks = [1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0e0]
        #ax2.set_yticks(major_ticks)
        #ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        #ax2.get_yaxis().get_major_formatter().labelOnlyBase = False
        
        # Set the aspect ratio using ax1
        #xmin,xmax = ax1.get_xlim()
        #ymin,ymax = ax1.get_ylim()
        #aspect_ratio = abs(math.log10(xmax)-math.log10(xmin)) \
        #             / abs(math.log10(ymax)-math.log10(ymin))
        #print("aspect_ratio = {0}".format(aspect_ratio))
        #ax1.set_aspect(aspect_ratio,adjustable='box-forced')
        #ax2.set_aspect(aspect_ratio,adjustable='box-forced')

        # Draw the grid lines
        ax1.grid(True,which='major',linestyle='-',color='#696969')
        ax1.grid(True,which='minor',linestyle=':',color='#8e8e8e')
        
        # Plot the graphs
        lns_ftn  = ax1.loglog(self.nptcl, self.wtftn,'r-', \
                              linewidth=2, marker='s', label="Fortran")
        lns_cpp  = ax1.loglog(self.nptcl, self.wtcpp,'b-', \
                              linewidth=2, marker='o', label="C++")
        lns_diff = ax2.loglog(self.nptcl, self.rdiff,'k-', \
                              linewidth=1, marker='*', label="Rel. diff.")
        
        # Plot the legends
        lns = lns_ftn + lns_cpp + lns_diff
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, fontsize=8, loc='lower right')
        
        # Adjust the size of the markers of the coordinate axes
        #SetMarkerAppearence(mpyplot.gca(),length=5,width=2)

#================================
#   Definitions of functions
#================================
def SetMarkerAppearence(ax,length,width):
    # ax := matplotlib.pyplot.gca()
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(length)      # major ticksの長さ
        line.set_markeredgewidth(width)  # major ticksの幅


#================================
#   Main function
#================================
# [1] Create instances
Kcomp  = MeasuredData()
ofp    = MeasuredData()
knl    = MeasuredData()
ajisai = MeasuredData()
xc30   = MeasuredData()

basename  = 'nbody_comp_serial'
ps_file   = basename + '.ps'
eps_file  = basename + '.eps'
pdf_file  = basename + '.pdf'

# [2] Read files
Kcomp.read_file("nbody_data_serial_K.txt")
ofp.read_file("nbody_data_serial_ofp.txt")
knl.read_file("nbody_data_serial_knl.txt")
ajisai.read_file("nbody_data_serial_ajisai.txt")
xc30.read_file("nbody_data_serial_xc30.txt")


# [3] Make a figure
# Basic setting
font = {'family' : 'Verdana',
        'weight' : 'normal',
        'size'   : '14'}
mpyplot.rc('font',**font)
mpyplot.rc('text',usetex=True)

# Create a drawing region
fig = mpyplot.figure(figsize=(8.267,2.923)) 
gs = gridspec.GridSpec(1,3)

# Plot K-computer data
ax = fig.add_subplot(gs[0,0])
panel_name = r'\textbf{(d) K Computer}' \
           + "\n" \
           + r'\hspace{1.5pc}\textbf{(normal case)}'
Kcomp.plot_panel(ax,panel_name=panel_name, \
                 with_xlabel=False,with_y2label=False)
# Plot Oakforest-PACS
ax = fig.add_subplot(gs[0,1])
panel_name = r'\textbf{(e) Oakforest-PACS}' \
           + "\n" \
           + r'\hspace{1.5pc}\textbf{(normal case)}'
ofp.plot_panel(ax,panel_name=panel_name, \
               with_y1label=False,with_y2label=False)
# Plot aterui
ax = fig.add_subplot(gs[0,2])
panel_name = r'\textbf{(f) Cray XC30}' \
           + "\n" \
           + r'\hspace{1.5pc}\textbf{(normal case)}'
xc30.plot_panel(ax,panel_name=panel_name, \
                with_xlabel=False,with_y1label=False)
# Plot ajisai
#ax = fig.add_subplot(gs[0,2])
#ajisai.plot_panel(ax,with_xlabel=False,with_y1label=False)

# Display figure 
#gs.tight_layout(fig) # TODO: why tight_layout does not working ?
fig.savefig(ps_file)

# Make eps,pdf,bb files
cmd = 'ps2eps -B -l -g -a -f ' + ps_file
subprocess.call(cmd,shell=True)
cmd = 'eps2pdf -B -H -f ' + eps_file
subprocess.call(cmd,shell=True)
cmd = 'extractbb ' + pdf_file
subprocess.call(cmd,shell=True)
