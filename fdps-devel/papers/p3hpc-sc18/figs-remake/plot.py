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
    import matplotlib.pyplot as mpyplot
    import matplotlib.scale  as mscale
    import matplotlib.cm     as cm
    from matplotlib.ticker import *
except ImportError:
    print("Module matplotlib.pyplot,scale,cm is not found.")
    quit()

try:
    import numpy as np
except ImportError:
    print("Module numpy is not found.")
    quit()

#================================
#   Class definition
#================================
class MeasuredData:
    def __init__(self):
        self.nproc = []
        self.total = [] # [s]
        self.force = [] # [s]
        self.comm  = [] # [s]
        self.misc  = [] # [s] 
        self.speed = [] # [Pflops]

    def set_data(self,nproc,total,force,comm,misc,speed):
        self.nproc = np.array(nproc)
        self.total = np.array(total)
        self.force = np.array(force)
        self.comm  = np.array(comm)
        self.misc  = np.array(misc)
        self.speed = np.multiply(np.array(speed),[1.0e-3])

#================================
#   Main function
#================================
# [1] Create instances
gyoukou_weak      = MeasuredData()
gyoukou_strong    = MeasuredData()
taihulight_weak   = MeasuredData()
taihulight_strong = MeasuredData()
shoubu_b_weak     = MeasuredData()

# [2] Set data
# Gyoukou (weak)
nproc = [1024, 2048, 4096, 8192]
total = [0.7270650734374997, \
         0.74774775, \
         0.7867435796874999, \
         0.8229698140625001]
force = [0.3322833125000001, \
         0.39179512499999997, \
         0.35494600000000004, \
         0.453257234375]
comm  = [0.11396353124999997, \
         0.12069094062499994, \
         0.14302774687500006, \
         0.14727031875]
misc  = [0.28081822968749964, \
         0.23526168437500014, \
         0.2887698328124998, \
         0.22244226093750016]
speed = [918.2438056658906, \
         2883.7799913139156, \
         3644.325894771016, \
         10634.98302178956]
gyoukou_weak.set_data(nproc,total,force,comm,misc,speed)
# Gyoukou (strong)
nproc = [1024, 2048, 4096, 8192]
total = [0.7270650734374997, \
         0.3837818745312501, \
         0.21932615546875006, \
         0.14768426921875]
force = [0.3322833125000001, \
         0.15120714062500004, \
         0.08535239843749999, \
         0.053817257812500005]
comm  = [0.11396353124999997, \
         0.06680079062499998, \
         0.041718070312499994, \
         0.03570155625]
misc  = [0.28081822968749964, \
         0.1657739432812501, \
         0.09225568671875009, \
         0.05816545515624999]
speed = [918.2438056658906, \
         1744.7280458928421, \
         3064.8373814028573, \
         4555.603677758406]
gyoukou_strong.set_data(nproc,total,force,comm,misc,speed)
# Taihulist (weak)
nproc = [10000, 20000, 40000, 80000, 160000]
total = [2.580628327941895, \
         2.2615068715836055, \
         2.67435802889122, \
         2.382595780539304, \
         2.8766929706248443]
force = [2.0738286050109305, \
         1.6314158783084165, \
         2.1326778060902143, \
         1.7008116195452203, \
         2.3113758577915178]
comm  = [0.04099083371738744, \
         0.04021286473107465, \
         0.06406917883668937, \
         0.05272958375053351, \
         0.08977260964184093]
misc  = [0.465808889213577, \
         0.5898781285441144, \
         0.4776110439643162, \
         0.6290545772435503, \
         0.4755445031914856]
speed  = [3204.3261547308553, \
          5165.548209091579, \
          12608.662351863843, \
          19952.045633082915, \
          47888.61627810115]
taihulight_weak.set_data(nproc,total,force,comm,misc,speed)
# Taihulight (strong)
nproc = [10000, 20000, 40000, 80000, 160000]
total = [2.580628327941895, \
         1.3178478957543114, \
         0.6988592633292683, \
         0.387629450194936, \
         0.24435274373445282]
force = [2.0738286050109305, \
         1.0498853363687886, \
         0.5565431241248007, \
         0.2991384625750016, \
         0.17652858685505632]
comm  = [0.04099083371738744, \
         0.025323317251604752, \
         0.029773664881531655, \
         0.023282531723907603, \
         0.03221608008516341]
misc = [0.465808889213577, \
        0.24263924213391808, \
        0.11254247432293596, \
        0.0652084558960268, \
        0.035608076794233084]
speed = [3204.3261547308553, \
         6290.965981159734, \
         11891.487147599888, \
         21346.740782566496, \
         34097.29450769898]
taihulight_strong.set_data(nproc,total,force,comm,misc,speed)
# Shoubu-B (weak)
nproc = [8, 16, 32, 64, 128, 256, 512]
total = [0.47708846874999994, \
         0.45259629531250006, \
         0.4968537953124998, \
         0.47271901875000005, \
         0.5067804046875001, \
         0.49608912499999996, \
         0.5251760640625]
force = [0.3270154375, \
         0.25166479687500004, \
         0.34432323437499995, \
         0.26950279687499995, \
         0.348261421875, \
         0.27892371874999994, \
         0.3597390468750001]
comm  = [0.018102199062500005, \
         0.020064457968750005, \
         0.018304459062499993, \
         0.023076479062500002, \
         0.026515634062500012, \
         0.030637049218750007, \
         0.030443944375000004]
misc  = [0.13197083218749994, \
         0.18086704046875002, \
         0.13422610187499984, \
         0.18013974281250011, \
         0.13200334875000008, \
         0.18652835703125, \
         0.13499307281249986]
speed = [16.06001926660484, \
         20.990317637135153, \
         63.66621388109619, \
         84.07108329402284, \
         256.8370813000585, \
         327.6870864685857, \
         1009.2482050684664]
shoubu_b_weak.set_data(nproc,total,force,comm,misc,speed)

# [3] Make figures
# Common setting
font = {'family' : 'Verdana',
        'weight' : 'normal',
        'size'   : '14'}
mpyplot.rc('font',**font)
mpyplot.rc('text',usetex=True)
# Printing rcParams is useful to check its content:
#print(mpyplot.rcParams)
mpyplot.rcParams['xtick.minor.size'] = 4
mpyplot.rcParams['xtick.minor.width'] = 1
mpyplot.rcParams['ytick.minor.size'] = 4
mpyplot.rcParams['ytick.minor.width'] = 1
##############################################################
# [Figure 3]
##############################################################
basename  = 'weak_scaling2'
ps_file   = basename + '.ps'
eps_file  = basename + '.eps'
pdf_file  = basename + '.pdf'
fig, ax = mpyplot.subplots(figsize=(5.5,5))
ax.set_xlabel("MPI processes",fontsize=18)
ax.set_ylabel("Time per one step (sec)",fontsize=18)
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
ax.grid(True,which='major',linestyle='-',color='#696969')
ax.grid(True,which='minor',linestyle=':',color='#696969')
ax.loglog(taihulight_weak.nproc, taihulight_weak.total, \
          'k-^',linewidth=1,markersize=6,label='Taihulight')
ax.loglog(gyoukou_weak.nproc, gyoukou_weak.total, \
          'k-+',linewidth=1,markersize=10,label='Gyoukou')
ax.loglog(shoubu_b_weak.nproc, shoubu_b_weak.total, \
          'k-D',linewidth=1,markersize=6,label='Shoubu-B')
ax.legend(loc='upper left',fontsize=12)
fig.tight_layout()
fig.savefig(ps_file)
cmd = 'ps2eps -B -l -g -a -f ' + ps_file
subprocess.call(cmd,shell=True)
cmd = 'eps2pdf -B -H -f ' + eps_file
subprocess.call(cmd,shell=True)
cmd = 'extractbb ' + pdf_file
subprocess.call(cmd,shell=True)
##############################################################
# [Figure 4]
##############################################################
basename  = 'weak_scaling_speed'
ps_file   = basename + '.ps'
eps_file  = basename + '.eps'
pdf_file  = basename + '.pdf'
fig, ax = mpyplot.subplots(figsize=(5.5,5))
ax.set_xlabel("MPI processes",fontsize=18)
ax.set_ylabel("Performance (PF)",fontsize=18)
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
ax.grid(True,which='major',linestyle='-',color='#696969')
ax.grid(True,which='minor',linestyle=':',color='#696969')
ax.loglog(taihulight_weak.nproc, taihulight_weak.speed, \
          'k-^',linewidth=1,markersize=6,label='Taihulight')
ax.loglog(gyoukou_weak.nproc, gyoukou_weak.speed, \
          'k-+',linewidth=1,markersize=10,label='Gyoukou')
ax.loglog(shoubu_b_weak.nproc, shoubu_b_weak.speed, \
          'k-D',linewidth=1,markersize=6,label='Shoubu-B')
ax.legend(loc='upper left',fontsize=12)
fig.tight_layout()
fig.savefig(ps_file)
cmd = 'ps2eps -B -l -g -a -f ' + ps_file
subprocess.call(cmd,shell=True)
cmd = 'eps2pdf -B -H -f ' + eps_file
subprocess.call(cmd,shell=True)
cmd = 'extractbb ' + pdf_file
subprocess.call(cmd,shell=True)
##############################################################
# [Figure 5]
##############################################################
basename  = 'strong_scaling'
ps_file   = basename + '.ps'
eps_file  = basename + '.eps'
pdf_file  = basename + '.pdf'
fig, ax = mpyplot.subplots(figsize=(5.5,5))
ax.set_xlabel("MPI processes",fontsize=18)
ax.set_ylabel("Time per one step (sec)",fontsize=18)
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
ax.grid(True,which='major',linestyle='-',color='#696969')
ax.grid(True,which='minor',linestyle=':',color='#696969')
ax.loglog(taihulight_strong.nproc, taihulight_strong.total, \
          'k-^',linewidth=1,markersize=6,label='Taihulight')
ax.loglog(gyoukou_strong.nproc, gyoukou_strong.total, \
          'k-+',linewidth=1,markersize=10,label='Gyoukou')
ax.legend(loc='upper right',fontsize=12)
fig.tight_layout()
fig.savefig(ps_file)
cmd = 'ps2eps -B -l -g -a -f ' + ps_file
subprocess.call(cmd,shell=True)
cmd = 'eps2pdf -B -H -f ' + eps_file
subprocess.call(cmd,shell=True)
cmd = 'extractbb ' + pdf_file
subprocess.call(cmd,shell=True)
