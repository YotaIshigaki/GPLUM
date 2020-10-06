import math as m
#import numpy as np

class Particle:
    """
    Class to store data of particle
    """
    #def __init__(self, mass, x, y, z, vx, vy, vz):
    #    self.mass = mass
    #    self.x  = x
    #    self.y  = y
    #    self.z  = z
    #    self.vx = vx
    #    self.vy = vy
    #    self.vz = vz

    def __init__(self, mass, r_p, f, x, y, z, vx, vy, vz, ngb, flag):
        self.mass = mass
        self.r_p  = r_p
        self.f  = f
        self.x  = x
        self.y  = y
        self.z  = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.ngb = ngb
        self.flag = flag
        
    def pos(self, target=None):
        if target is None:
            return [self.x, self.y, self.z]
        else:
            return [self.x-target.x, self.y-target.y, self.z-target.z]

    def pos2(self, target=None):
        x, y, z = self.pos(target)
        return x**2 + y**2 + z**2
        
    def vel(self, target=None):
        if target is None:
            return [self.vx, self.vy, self.vz]
        else:
            return [self.vx-target.vx, self.vy-target.vy, self.vz-target.vz]

    def vel2(self, target=None):
        vx, vy, vz = self.vel(target)
        return vx**2 + vy**2 + vz**2

    def getOrbitalElement(self, mu=1.):
        r = m.sqrt(self.pos2())
        v2 = self.vel2()
        rv = self.x*self.vx + self.y*self.vy + self.z*self.vz
        rxv = [self.y*self.vz - self.z*self.vy,\
               self.z*self.vx - self.x*self.vz,\
               self.x*self.vy - self.y*self.vx\
        ]
        
        self.ax  = 1./(2./r - v2/mu)
        self.ecc = m.sqrt( (1.-r/self.ax)**2 + (rv)**2/(mu*self.ax) )
        self.inc = m.atan2(m.sqrt(rxv[0]**2+rxv[1]**2), rxv[2])

        return [self.ax, self.ecc, self.inc]

class Header:
    """
    Class to store data of header in snapshot
    """
    def __init__(self, time, n_body, id_next,\
                 etot0, ekin0, ephi_sun0, ephi_planet0, edisp0, \
                 etot1, ekin1, ephi_sun1, ephi_planet1, edisp1):
        self.time    = time
        self.n_body  = n_body
        self.id_next = id_next
        self.etot0        = etot0
        self.ekin0        = ekin0
        self.ephi_sun0    = ephi_sun0
        self.ephi_planet0 = ephi_planet0
        self.edisp0       = edisp0
        self.etot1        = etot1
        self.ekin1        = ekin1
        self.ephi_sun1    = ephi_sun1
        self.ephi_planet1 = ephi_planet1
        self.edisp1       = edisp1
