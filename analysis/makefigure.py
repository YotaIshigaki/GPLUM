from analysis import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys


def makeSnap(datfile, pngfile, \
             bHeader=True, \
             m_sun=1.0, \
             xrange=None, \
             yrange=None, \
             axis="a-e", \
             m_cut=0., \
             m_emp=0., \
             size=10., \
             color0="blue",\
             color1="midnightblue") :
    """
    Make a snapshot figure from a snapshot file.
    
    Arguments
    ----------
    datfile : character strings. snapshot filename
    pngfile : character string. snapshot figure filename
    bHeader : boolian. flag as to whether header exists in snapshot (default True)
    m_sun :   float. mass of a central star (default 1.0)
    xrange :  [float,float]. x direction range of a snapshot figure
    yrange :  [float,float]. y direction range of a snapshot figure
    axis :    character string. elements of x,y in a snapshot figure (default "a-e")
              "a-e"  -- semimajor axis - eccentricity
              "a-i"  -- semimajor axis - inclination
              "a-ei" -- semimajor axis - root mean squar of eccentricity and inclination
    m_cut :   float. minimum particle mass to draw in a snapshot figure (default 0.)
    m_emp :   float. minimum particle mass to emphasize in a snapshot figure (default 0.)
    size :    float. marker size of particles with mass m_size (default 10.)
    color0 :   character string. marker color (default "blue")
    color1 :   character string. marker color of emphasized particles (default "midnightblue")

    Returns
    ---------
    m_max :   float. largest particle mass
    m_ave :   float. mean particle mass
    ecc_rms : float. root mean square of eccentricities
    inc_rms : float. root mean square of inclinations
    """

    header, pp = readSnap(datfile, bHeader)

    time_str = " "
    if bHeader : time_str = "%0.f yr"%(header.time/(2*PI))

    m_max, m_ave, ecc_rms, inc_rms = clacOrbitalElements(pp, m_sun)

    m_size = 1.e23/M
    siz = []
    ax  = []
    ecc = []
    inc = []
    siz_emp = []
    ax_emp  = []
    ecc_emp = []
    inc_emp = []
    for ptcl in pp.values() :
        if ptcl.mass >= m_emp :
            siz_emp.append(((ptcl.mass/m_size)**(1./2.))*size)
            ax_emp.append(ptcl.ax)
            ecc_emp.append(ptcl.ecc)
            inc_emp.append(ptcl.inc)
        elif ptcl.mass >= m_cut :
            siz.append(((ptcl.mass/m_size)**(1./2.))*size)
            ax.append(ptcl.ax)
            ecc.append(ptcl.ecc)
            inc.append(ptcl.inc)

    # Make Figure
    if axis in ["a-e", "a-i", "a-ei"] and pngfile is not None: 
        plt.rcParams['font.family'] = "Times New Roman"
        plt.rcParams['font.size'] = 17
        
        plt.figure()
        if axis == "a-e" :
            plt.scatter(ax, ecc, s=siz,\
                        c=color0, alpha=0.08, linewidths=0.5)
            plt.scatter(ax_emp, ecc_emp, s=siz_emp,\
                        c=color1, alpha=0.16, linewidths=0.5)
            plt.xlabel("Semi-major Axis (AU)", fontsize=18)
            plt.ylabel("Eccentricity", fontsize=18)

            if xrange==None :
                xmin = min(ax+ax_emp)
                xmax = max(ax+ax_emp)
                xrange = [xmin-(xmax-xmin)*0.05, xmax+(xmax-xmin)*0.05]
            if yrange==None :
                ymax = max(ecc+ecc_emp)
                yrange = [0.-(ymax-0.)*0.05, ymax+(ymax-0.)*0.05]
            
        elif axis == "a-i" :
            plt.scatter(ax, inc, s=siz,\
                        c=color0, alpha=0.08, linewidths=0.5)
            plt.scatter(ax_emp, inc_emp, s=siz_emp,\
                        c=color1, alpha=0.16, linewidths=0.5)
            plt.xlabel("Semi-major Axis (AU)", fontsize=18)
            plt.ylabel("Inclination", fontsize=18)

            if xrange==None :
                xmin = min(ax+ax_emp)
                xmax = max(ax+ax_emp)
                xrange = [xmin-(xmax-xmin)*0.05, xmax+(xmax-xmin)*0.05]
            if yrange==None :
                ymax = max(inc+inc_emp)
                yrange = [0.-(ymax-0.)*0.05, ymax+(ymax-0.)*0.05]
            
        elif axis == "a-ei" :
            eccinc     = [sqrt(ecci**2 + inci**2) for (ecci, inci) in zip(ecc, inc)]
            eccinc_emp = [sqrt(ecci**2 + inci**2) for (ecci, inci) in zip(ecc_emp, inc_emp)]
            plt.scatter(ax, eccinc, s=siz, \
                        c=color0, alpha=0.08, linewidths=0.5)
            plt.scatter(ax_emp, eccinc_emp, s=siz_emp, \
                        c=color1, alpha=0.16, linewidths=0.5)
            plt.xlabel("$a$ (AU)", fontsize=18)
            plt.ylabel("$\sqrt{e^2+i^2}$", fontsize=18)

            if xrange==None :
                xmin = min(ax+ax_emp)
                xmax = max(ax+ax_emp)
                xrange = [xmin-(xmax-xmin)*0.05, xmax+(xmax-xmin)*0.05]
            if yrange==None :
                ymax = max(eccinc+eccinc_emp)
                yrange = [0.-(ymax-0.)*0.05, ymax+(ymax-0.)*0.05]

        plt.scatter([],[],alpha=0,label=' ')
        plt.legend(fontsize=10, loc='upper right', frameon=False, title=time_str)
        plt.subplots_adjust(bottom=0.13, left=0.15)

        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlim(xrange)
        plt.ylim(yrange)
    
        plt.savefig(pngfile, dpi=130)
        plt.close()

    return [header, m_max, m_ave, ecc_rms, inc_rms]


def makeCumulativeNumberDistribution(datfiles, pngfile, \
                                     bHeader=True, \
                                     var="mass" , \
                                     xrange=None, \
                                     yrange=None, \
                                     labels=None, \
                                     colors=["mediumorchid", "blueviolet", "blue", "midnightblue"], \
                                     styles=[":","-.","--","-"]) :
    """
    Make a figure of cumulative number distributuion
    
    Arguments
    ----------
    datfiles: list of character strings. snapshot filenames
    pngfile:  character string. snapshot figure filename
    bHeader : boolian. flag as to whether header exists in snapshot (default True)
    var :     charactor string. variable of number distribution (default mass)
              "mass" -- mass distribution
              "ecc"  -- eccentricity distribution
              "inc"  -- inclination distribution
    xrange:   [float,float]. x direction range of a figure
    yrange:   [float,float]. y direction range of a figure
    labels:   list of character strings. labels (default None)
    colors:   list of character strings.  marker colors 
              (default ["mediumorchid", "blueviolet", "blue", "midnightblue"])
    styles:   list of character strings. line styles (default ["-.","--","-"])

    Returns
    ---------
    """

    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = 17

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")

    i = 0
    xmin = sys.float_info.max
    xmax = 0.
    ymin = sys.float_info.max
    ymax = 0.
    for datfile in datfiles:
        header, pp = readSnap(datfile, bHeader)

        label = " "
        if labels == None :
            if bHeader :
                label = "%0.f yr"%(header.time/(2*PI))
        else :
            label = labels[i]

        cumulativeNumber = calcCumulativeNumberDistribution(pp, var=var)
        
        if xrange==None :
            mass = [i for i, j in cumulativeNumber]
            xmin = min(mass+[xmin])*M
            xmax = max(mass+[xmax])*M
        if yrange==None :
            Nc = [j for i, j in cumulativeNumber]
            ymin = min(Nc+[ymin])
            ymax = max(Nc+[ymax])

        mass = []
        Nc = []
        for m, nc in cumulativeNumber :
            mass.append(m*M)
            Nc.append(nc-1)
            mass.append(m*M)
            Nc.append(nc)
        mass.append(0.)
        Nc.append(cumulativeNumber[-1][1])       
            
        plt.plot(mass, Nc, linewidth=1, label=label, \
                 color=colors[i], ls=styles[i])

        i = i + 1

    if xrange==None : xrange = [xmin/1.5, xmax*1.5]
    if yrange==None : yrange = [ymin/1.5, ymax*1.5]
    if var == "ecc" :
        plt.xlabel("Eccentricity", fontsize=17)
    elif var == "inc" :
        plt.xlabel("Inclination", fontsize=17)
    else:
        plt.xlabel("Mass (g)", fontsize=17)
    plt.ylabel("Cumulative Number", fontsize=17)
    plt.legend(fontsize=10, loc='upper right', frameon=False)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(xrange)
    plt.ylim(yrange)

    plt.savefig(pngfile, dpi=100)
    plt.close()

    return


def makeEccentricityDistribution(datfiles, pngfile, \
                                 bHeader=True, \
                                 m_sun=1.0, \
                                 var="ecc",\
                                 xrange=None, \
                                 yrange=None, \
                                 w_b=2.,\
                                 labels=None, \
                                 colors=["mediumorchid", "blueviolet", "blue", "midnightblue"], \
                                 styles=[":","-.","--","-"],\
                                 markers=["^","s","o","D","h","*"]) :
    """
    Make a figure of distribution of eccentricity and inclination against mass
    
    Arguments
    ----------
    datfiles: list of character strings. snapshot filenames
    pngfile:  character string. snapshot figure filename
    bHeader : boolian. flag as to whether header exists in snapshot (default True)
    m_sun :   float. mass of a central star (default 1.0)
    var :     charactor string. variable of y axis (default ecc)
              "ecc"  -- rms of eccentricity
              "inc"  -- rms of inclination
              "num"  -- the number of particle
    xrange:   [float,float]. x direction range of a figure
    yrange:   [float,float]. y direction range of a figure
    w_b :     float. width of bin is set to logarithm this value in log scale (default 2.0)
    labels:   list of character strings. labels
    colors:   list of character strings.  marker colors 
              (default ["mediumorchid", "blueviolet", "blue", "midnightblue"])
    styles:   list of character strings. line styles 
              (default ["-.","--","-"])
    marks :   list of character strings. mark styles 
              (default ["^","s","o","D","h","*"])

    Returns
    ---------
    """

    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.unicode'] = True
    plt.rcParams['font.family'] = "Times"
    plt.rcParams['font.size'] = 17

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")
    if not var == "num" :
        plt.yticks([0.00005, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, \
                    0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1], \
                   ["0.00005", "0.0001", "0.0002", "0.0005", "0.001", "0.002", "0.005", \
                    "0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "1"])

    i = 0
    xmin = sys.float_info.max
    xmax = 0.
    ymin = sys.float_info.max
    ymax = 0.
    for datfile in datfiles:
        header, pp = readSnap(datfile, bHeader)

        label = " "
        if labels == None :
            if bHeader :
                label = "%0.f yr"%(header.time/(2*PI))
        else :
            label = labels[i]
            
        dist = calcEccentricityDistribution(pp, w_b, m_sun)

        mass = []
        ecc  = []
        mass_e = []
        ecc_e  = []
        decc0 = []
        decc1 = []       

        if var == "num" :
            for dmmy1, dmmy2, e, dmmy3, m, dmmy4, dmmy5, dmmy6, dmmy7, dmmy8, dmmy9 in dist :
                if e != None :
                    mass.append(m*M)
                    ecc.append(e)
        
        elif var == "inc" :
            for dmmy1, dmmy2, dmmy3, dmmy4, m, dmmy5, dmmy6, dmmy7, e, ebtm, etop in dist :
                if e != None :
                    mass.append(m*M)
                    ecc.append(e)
                    if not etop == None :
                        mass_e.append(m*M)
                        ecc_e.append(e)
                        decc0.append(e-ebtm)
                        decc1.append(etop-e)

                        #print(e-ebtm,etop-e)

        else :
            for dmmy1, dmmy2, dmmy3, dmmy4, m, e, ebtm, etop, dmmy5, dmmy6, dmmy7 in dist :
                if e != None :
                    mass.append(m*M)
                    ecc.append(e)
                    if not etop == None :
                        mass_e.append(m*M)
                        ecc_e.append(e)
                        decc0.append(e-ebtm)
                        decc1.append(etop-e)

                        #print(e-ebtm,etop-e)

        plt.plot(mass, ecc, linewidth=1, label=label, \
                 color=colors[i], marker=markers[i], ls=styles[i], fillstyle="none")
        if not var == "num" :
            plt.errorbar(mass_e, ecc_e, yerr=[decc0,decc1], color=colors[i], ls="None", capsize=3)

        if xrange==None :
            xmin = min(mass+[xmin])
            xmax = max(mass+[xmax])
        if yrange==None :
            ymin = min(ecc+[ymin])
            ymax = max(ecc+[ymax])

        i = i + 1

    if xrange==None : xrange = [xmin/2, xmax*1.5]
    if yrange==None : yrange = [ymin/2, ymax*1.5]
    plt.xlabel("Mass (g)", fontsize=17)
    if var == "num" :
        plt.ylabel("The Number of Particles", fontsize=17)
    elif var == "inc" :
        plt.ylabel(r"$\langle i^2\rangle^{1/2}$", fontsize=17)
    else :
        plt.ylabel(r"$\langle e^2\rangle^{1/2}$", fontsize=17)
    #plt.ylabel("Mean Eccentricity", fontsize=17)
    #plt.ylabel("Mean Inclination", fontsize=17)
    plt.legend(fontsize=8, loc='lower left', frameon=False)
    plt.subplots_adjust(bottom=0.15, left=0.18)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(xrange)
    plt.ylim(yrange)

    plt.savefig(pngfile, dpi=200)

    plt.rcParams['text.usetex'] = False
    plt.rcParams['text.latex.unicode'] = False
    plt.close()

    return

def makeEccentricityDistribution2(datfile, pngfile, \
                                  bHeader = True, \
                                  m_sun=1.0, \
                                  xrange=None, \
                                  yrange=None, \
                                  w_b=2.,\
                                  labels=None, \
                                  colors=["blue", "red"], \
                                  styles=["-","--"],\
                                  markers=["o","s"]) :
    """
    Make a figure of distribution of eccentricity and inclination against mass
    
    Arguments
    ----------
    datfile:  character strings. snapshot filename
    pngfile:  character string. snapshot figure filename
    bHeader : boolian. flag as to whether header exists in snapshot (default True)
    m_sun :   float. mass of a central star (default 1.0)
    xrange:   [float,float]. x direction range of a figure (default [3.e21, 1.26])
    yrange:   [float,float]. y direction range of a figure (default [1.e5, 0.5])
    w_b :     float. width of bin is set to logarithm this value in log scale (default 2.0)
    labels:   list of character strings. labels (default None)
    colors:   list of character strings.  marker colors 
              (default ["blue", "red"])
    styles:   list of character strings. line styles 
              (default ["-","--"])
    marks :   list of character strings. mark styles 
              (default ["o","s"])

    Returns
    ---------
    """

    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.unicode'] = True
    plt.rcParams['font.family'] = "Times"
    plt.rcParams['font.size'] = 17

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")
    plt.yticks([0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, \
                0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1], \
               [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, \
                0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1])
    
    header, pp = readSnap(datfile, bHeader)

    label1 = " "
    label2 = " "
    if labels == None :
        label1 = r"$\langle e^2\rangle^{1/2}$"
        label2 = r"$\langle i^2\rangle^{1/2}$"
    else :
        label1 = labels[0]
        label2 = labels[1]
        
    dist = calcEccentricityDistribution(pp)

    mass = []
    ecc = []
    inc = []
    decc = []
    dinc = []
    for dmmy1, dmmy2, dmmy3, dmmy4, m, e, de, i, di in dist :
        if e != None :
            mass.append(m*M)
            ecc.append(e)
            if de == None :
                decc.append(0.)
            else :
                decc.append(de)
            inc.append(i)
            if di == None :
                dinc.append(0.)
            else :
                dinc.append(di)
            
    plt.plot(mass, ecc, linewidth=1, label=label1, \
             color=colors[0], marker=markers[0], ls=styles[0], fillstyle="none")
    plt.errorbar(mass, ecc, yerr=decc, color=colors[0], ls=styles[0], capsize=3)
    plt.plot(mass, inc, linewidth=1, label=label2, \
             color=colors[1], marker=markers[1], ls=styles[1], fillstyle="none")
    plt.errorbar(mass, inc, yerr=dinc, color=colors[1], ls=styles[1], capsize=3)

    if xrange==None :
        xmin = min(mass)
        xmax = max(mass)
        xrange = [xmin/1.5, xmax*1.5]
    if yrange==None :
        ymin = min(ecc+inc)
        ymax = max(ecc+inc)
        yrange = [ymin/1.5, ymax*1.5]

    plt.xlabel("Mass (g)", fontsize=17)
    #plt.ylabel(r"$\langle e^2\rangle^{1/2}, \langle i^2\rangle^{1/2}$", fontsize=17)
    plt.ylabel("Eccentricity, Inclination", fontsize=17)
    plt.legend(fontsize=8, loc='lower left', frameon=False)
    plt.subplots_adjust(bottom=0.15, left=0.18)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(xrange)
    plt.ylim(yrange)

    plt.savefig(pngfile, dpi=100)

    plt.rcParams['text.usetex'] = False
    plt.rcParams['text.latex.unicode'] = False
    plt.close()

    return
