import numpy as np
import matplotlib.pylab as plt

def Func0(a, x):
    if a>x:
        ret = 0.5*(x*np.sqrt(a*a-x*x) + a*a*np.arctan(x/np.sqrt(a*a-x*x)))
    else:
        ret = 0.5*a*a*0.5*np.pi 
    return ret

def CalcS0(a0, a1, x0, x1):
    S0 = (Func0(a1, x1) - Func0(a1, x0)) - (Func0(a0, x1) - Func0(a0, x0))
    return S0

def CalcX(dS, a0, a1, x0):
    xmin = x0
    xmax = a1
    #while True:
    for i in range(30):
        x1 = (xmin+xmax)*0.5
        dSnow = CalcS0(a0, a1, x0, x1)
        eps = dSnow - dS
        if( (np.fabs(eps))/dS < 1e-10):
            break
        elif ( eps > 0.0 ):
            xmax = x1
        else:
            xmin = x1
    return x1

def CalcDx(y_a0x0, y_a1x1, a0, a1, x0, x1, y):
    if y > y_a0x0:
        xmin = x0
    else:
        xmin = np.sqrt(a0*a0-y*y)
    if y > y_a1x1:
        xmax = np.sqrt(a1*a1-y*y)
    else:
        xmax = x1
    return xmax - xmin

def CalcS1(y_a0x0, y_a1x1, a0, a1, x0, x1, y0, y1):
    n_div = 20
    y_array = np.linspace(y0, y1, n_div)
    dy = y_array[1] - y_array[0]
    S1 = 0.0
    for y_coord in y_array:
        dx = CalcDx(y_a0x0, y_a1x1, a0, a1, x0, x1, y_coord)
        S1 += dx*dy
    S1 -= 0.5*dy*(CalcDx(y_a0x0, y_a1x1, a0, a1, x0, x1, y_array[0])+CalcDx(y_a0x0, y_a1x1, a0, a1, x0, x1, y_array[n_div-1]))
    return S1

def CalcY(dS, a0, a1, x0, x1, y1):
    # intersection
    y_a0x0 = 0.0
    if x0 < a0:
        y_a0x0 = np.sqrt(a0*a0-x0*x0)
    y_a1x1 = np.sqrt(a1*a1-x1*x1)
    y_a0x1 = 0.0
    if a0 > x1:
        y_a0x1 = np.sqrt(a0*a0 - x1*x1)
    ymax = y1
    ymin = y_a0x1
    #while True:
    for i in range(30):
        y0 = (ymin+ymax)*0.5
        dSnow = CalcS1(y_a0x0, y_a1x1, a0, a1, x0, x1, y0, y1)
        eps = dSnow - dS
        if( (np.fabs(eps))/dS < 1e-10):
            break
        elif ( eps > 0.0 ):
            ymin = y0
        else:
            ymax = y0
    #print("dS, CalcS1", dS, CalcS1(y_a0x0, y_a1x1, a0, a1, x0, x1, y0, y1))
    return y0

nx = 8
ny = 8
nproc = nx*ny
ax = 1.0
delta_ax = 0.1

ax_out = ax + 0.5*delta_ax
ax_in  = ax - 0.5*delta_ax
Stot = (np.pi * (ax_out*ax_out - ax_in*ax_in)) / 4.0
x_coord = np.zeros(nx+1)
x_prev = 0.0
dS0 = Stot / nx
dS1 = Stot / nproc
for ix in range(1,len(x_coord)):
    x_coord[ix] = CalcX(dS0, ax_in, ax_out, x_prev)
    x_now = x_coord[ix]
    y_coord = np.zeros(ny+1)
    y_prev = np.sqrt(ax_out*ax_out - x_prev*x_prev) # y_high y_a1x0
    y_coord[-1] = y_prev
    y_low = 0.0
    if ax_in > x_now:
        y_low = np.sqrt(ax_in*ax_in - x_now*x_now)
        
    """       
    plt.plot([x_prev, x_prev], [y_low, y_prev], color="k")
    plt.plot([x_now, x_now], [y_low, y_prev], color="k")
    plt.plot([x_prev, x_now], [y_prev, y_prev], color="k")
    """
    
    plt.plot([x_prev, x_prev], [y_low, y_prev], color="k")
    plt.plot([-x_prev, -x_prev], [y_low, y_prev], color="k")
    plt.plot([x_prev, x_prev], [-y_low, -y_prev], color="k")
    plt.plot([-x_prev, -x_prev], [-y_low, -y_prev], color="k")

    plt.plot([x_now, x_now],   [y_low, y_prev], color="k")
    plt.plot([-x_now, -x_now], [y_low, y_prev], color="k")
    plt.plot([x_now, x_now],   [-y_low, -y_prev], color="k")
    plt.plot([-x_now, -x_now], [-y_low, -y_prev], color="k")
    
    plt.plot([x_prev, x_now], [y_prev, y_prev], color="k")
    plt.plot([-x_prev, -x_now], [y_prev, y_prev], color="k")
    plt.plot([x_prev, x_now], [-y_prev, -y_prev], color="k")
    plt.plot([-x_prev, -x_now], [-y_prev, -y_prev], color="k")

    
    for iy in reversed(range(0, len(y_coord)-1)):
        y_coord[iy] = CalcY(dS1, ax_in, ax_out, x_prev, x_now, y_prev)
        plt.plot([x_prev, x_now], [y_coord[iy], y_coord[iy]], color="k")
        plt.plot([-x_prev, -x_now], [y_coord[iy], y_coord[iy]], color="k")
        plt.plot([-x_prev, -x_now], [-y_coord[iy], -y_coord[iy]], color="k")
        plt.plot([x_prev, x_now],  [-y_coord[iy], -y_coord[iy]], color="k")
        #print("iy, y_coord[iy]", iy, y_coord[iy])
        y_prev = y_coord[iy]
    #print("y_coord", y_coord)
    x_prev = x_coord[ix]

#print("x_coord", x_coord)
    
circle_out = plt.Circle((0,0), ax_out, fill=False)
circle_in  = plt.Circle((0,0), ax_in,  fill=False)

#plt.xlim([0.0, ax_out])
#plt.ylim([0.0, ax_out])
plt.xlim([-ax_out, ax_out])
plt.ylim([-ax_out, ax_out])
plt.xlabel("x")
plt.ylabel("y")
plt.axes().set_aspect(1.0)
plt.gca().add_patch(circle_in)
plt.gca().add_patch(circle_out)

#plt.show()
#plt.savefig("domain_cart.eps", bbox_inches='tight', pad_inches=0.1)
plt.savefig("domain_cart.eps")

plt.clf()

n_theta = 128
n_r = 2

r_array = np.linspace(ax_in, ax_out, n_r+1)
for r_coord in r_array:
    circle_tmp = plt.Circle((0,0), r_coord, fill=False)
    plt.gca().add_patch(circle_tmp)

#theta_array = np.linspace(0.0, np.pi*0.5, n_theta+1)
theta_array = np.linspace(0.0, 2.0*np.pi, n_theta+1)
for theta_coord in theta_array:
    x0 = ax_in*np.cos(theta_coord)
    x1 = ax_out*np.cos(theta_coord)
    y0 = ax_in*np.sin(theta_coord)
    y1 = ax_out*np.sin(theta_coord)
    #y0 = np.sqrt(ax_in*ax_in-x0*x0)
    #y1 = np.sqrt(ax_out*ax_out-x1*x1)
    plt.plot([x0, x1], [y0, y1], color="k")
    
circle_out = plt.Circle((0,0), ax_out, fill=False)
circle_in  = plt.Circle((0,0), ax_in,  fill=False)
#plt.xlim([0.0,ax_out])
#plt.ylim([0.0,ax_out])
plt.xlim([-ax_out, ax_out])
plt.ylim([-ax_out, ax_out])
plt.xlabel("x")
plt.ylabel("y")
plt.axes().set_aspect(1.0)
plt.gca().add_patch(circle_in)
plt.gca().add_patch(circle_out)
    
#plt.show()
#plt.savefig("domain_cyl.org.eps")
plt.savefig("domain_cyl.eps")

# after make eps files, do the below steps
# ps2pdf  -dEPSCrop xxx.org.eps
# pdf2ps xxx.org.pdf xxx.eps
