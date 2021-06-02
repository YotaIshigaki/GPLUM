import math as m
import numpy as np
import sys


def main() :
    makeinit(sys.argv[1], sys.argv[2], 1.)

    
def KeplerEq(u, ecc):
    return u - ecc * m.sin(u)

def solveKeplerEq(l, ecc):
    ecc2 = ecc *ecc
    ecc3 = ecc2*ecc
    ecc4 = ecc2*ecc2
    ecc5 = ecc3*ecc2
    ecc6 = ecc3*ecc3
    u = l \
        + (ecc     - ecc3/8. + ecc5/192.)*m.sin(l) \
        + (ecc2/2. - ecc4/6. + ecc6/48. )*m.sin(2.*l) \
        + (3.*ecc3/8. - 27.*ecc5/128.)*m.sin(3.*l) \
        + (   ecc4/3. -  4.*ecc6/15. )*m.sin(4.*l) \
        + 125.*ecc5/384.*m.sin(5.*l) \
        +  27.*ecc6/ 80.*m.sin(6.*l)

    if m.fabs(KeplerEq(u,ecc)-l) > 1.e-15 :
        u0 = u
        loop = 0
        while True : 
            u0 = u
            sinu0 = m.sin(u0);
            cosu0 = m.cos(u0);
            u = u0 - ((u0 - ecc*sinu0 - l)/(1. - ecc*cosu0))
            loop = loop + 1

            if m.fabs(u - u0) < 1.e-15 or loop > 10 : break

    return u



class Particle:
    def __init__(self, mass, r_p, f, ax, ecc, inc, Omega, omega, l):
        self.mass  = mass
        self.r_p   = r_p
        self.f     = f
        self.ax    = ax
        self.ecc   = ecc
        self.inc   = inc
        self.Omega = Omega
        self.omega = omega
        self.l     = l

    def getPosVel(self, mu=1.):
        u = solveKeplerEq(self.l, self.ecc)
        Px =  m.cos(self.omega)*m.cos(self.Omega) - m.sin(self.omega)*m.sin(self.Omega)*m.cos(self.inc)
        Py =  m.cos(self.omega)*m.sin(self.Omega) + m.sin(self.omega)*m.cos(self.Omega)*m.cos(self.inc)
        Pz =  m.sin(self.omega)*m.sin(self.inc)
        Qx = -m.sin(self.omega)*m.cos(self.Omega) - m.cos(self.omega)*m.sin(self.Omega)*m.cos(self.inc)
        Qy = -m.sin(self.omega)*m.sin(self.Omega) + m.cos(self.omega)*m.cos(self.Omega)*m.cos(self.inc)
        Qz =  m.cos(self.omega)*m.sin(self.inc)

        
        cosu = m.cos(u)
        sinu = m.sin(u)
        ecc_sq = m.sqrt(1.-self.ecc*self.ecc)
        self.x = self.ax * ((cosu-self.ecc) * Px + ecc_sq*sinu * Qx)
        self.y = self.ax * ((cosu-self.ecc) * Py + ecc_sq*sinu * Qy)
        self.z = self.ax * ((cosu-self.ecc) * Pz + ecc_sq*sinu * Qz)
        
        r2 = self.x*self.x + self.y*self.y + self.z*self.z
        n  = m.sqrt(mu / (self.ax*self.ax*self.ax));
        rinv = m.sqrt(1./r2)
        self.vx = self.ax*self.ax*n*rinv* (-sinu * Px + ecc_sq*cosu * Qx)
        self.vy = self.ax*self.ax*n*rinv* (-sinu * Py + ecc_sq*cosu * Qy)
        self.vz = self.ax*self.ax*n*rinv* (-sinu * Pz + ecc_sq*cosu * Qz)


def makedist(Mass, n, f, dens, p, a0, a1, ecc_hill, inc_hill, mu=1.):
    
    mass = Mass / n
    r_p = (3.*mass/(4.*m.pi*dens))**(1./3.)
    
    system = []
    
    for i in range(n) :
        R = np.random.rand(1)[0]

        ax = 0.
        if p != 2 :
            ax = ( (a1**(2.-p)-a0**(2.-p)) * R + a0**(2.-p) )**(1./(2.-p))
        else :
            ax = m.exp( (m.log(a1) - m.log(a0)) * R + m.log(a0) )

        h = ( mass/(3.*mu) )**(1./3.)
        ecc = m.fabs(np.random.normal(0.,ecc_hill*h,1)[0])
        inc = m.fabs(np.random.normal(0.,inc_hill*h,1)[0])

        Omega, omega, l = np.random.rand(3) * m.pi * 2.

        ptcl = Particle(mass, r_p, f, ax, ecc, inc, Omega, omega, l)
        ptcl.getPosVel(mu)

        system.append(ptcl)

    return system


def makeinit(filename0, filename1, mu=1.) :

    n0 = 0
    n1 = 0
    idx = 0

    f1 = open(filename1, mode='w')

    i = 0
    with open(filename0) as f0:
        for line in f0:
            part = [par.strip() for par in line.split()]
            part = [par for par in part if par != '']

            if i == 0 :
                n0 = int(part[0])
                n1 = int(part[1])
                np.random.seed(int(part[2]))

            elif i <= n0 :
                mass, r_p, f, ax, ecc, inc = float(part[0]), float(part[1]), float(part[2]), float(part[3]), float(part[4]), float(part[5])
                Omega, omega, l = np.random.rand(3) * m.pi * 2.
                
                p = Particle(mass, r_p, f, ax, ecc, inc, Omega, omega, l)
                p.getPosVel(mu)

                line1 = "{:d}".format(idx) + "\t" \
                        + "{:.15e}".format(p.mass) + "\t" \
                        + "{:.15e}".format(p.r_p) + "\t" \
                        + "{:.15e}".format(p.f) + "\t" \
                        + "{:.15e}".format(p.x) + "\t" \
                        + "{:.15e}".format(p.y) + "\t" \
                        + "{:.15e}".format(p.z) + "\t" \
                        + "{:.15e}".format(p.vx) + "\t" \
                        + "{:.15e}".format(p.vy) + "\t" \
                        + "{:.15e}".format(p.vz) + "\t" \
                        + "{:d}".format(0) + "\t" \
                        + "{:d}".format(0) + "\n"
                f1.write(line1)

                idx = idx + 1

            elif i <= n0+n1 :

                Mass, n, f, dens, p, a0, a1, ecc_hill, inc_hill = float(part[0]), int(part[1]), float(part[2]), float(part[3]), float(part[4]), float(part[5]), float(part[6]), float(part[7]), float(part[8])

                pp = makedist(Mass, n, f, dens, p, a0, a1, ecc_hill, inc_hill, mu)

                for p in pp :
                    line1 = "{:d}".format(idx) + "\t" \
                            + "{:.15e}".format(p.mass) + "\t" \
                            + "{:.15e}".format(p.r_p) + "\t" \
                            + "{:.15e}".format(p.f) + "\t" \
                            + "{:.15e}".format(p.x) + "\t" \
                            + "{:.15e}".format(p.y) + "\t" \
                            + "{:.15e}".format(p.z) + "\t" \
                            + "{:.15e}".format(p.vx) + "\t" \
                            + "{:.15e}".format(p.vy) + "\t" \
                            + "{:.15e}".format(p.vz) + "\t" \
                            + "{:d}".format(0) + "\t" \
                            + "{:d}".format(0) + "\n"
                    f1.write(line1)

                    idx = idx + 1

            i = i + 1

        
if __name__ == "__main__":
    main()

