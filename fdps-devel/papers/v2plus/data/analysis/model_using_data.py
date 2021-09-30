import sys
import math
import numpy as np
import matplotlib.pylab as plt

BW = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 8.0])

# N=4M, n_list=3649, n_grp=234, n_reuse=2
wtime_BW_reuse2 = np.array([0.970536315355,
                            0.533930584623,
                            0.315627719257,
                            0.206476286574,
                            0.151900570233,
                            0.124612712062])

# N=4M, n_list=3649, n_grp=234, n_reuse=4
wtime_BW_reuse4 = np.array([0.817605525207,
                            0.446640652075,
                            0.261158215509,
                            0.168416997226,
                            0.122046388085,
                            0.0988610835139])

# N=4M, n_list=3649, n_grp=234, n_reuse=8
wtime_BW_reuse8 = np.array([0.741140130133,
                            0.402995685801,
                            0.233923463635,
                            0.149387352552,
                            0.107119297011,
                            0.0859852692398])

# N=4M, n_list=3649, n_grp=234, n_reuse=16
wtime_BW_reuse16 = np.array([0.702907432596,
                            0.381173202664,
                            0.220306087698,
                            0.139872530215,
                            0.0996557514734,
                            0.0795473621027])

# N=4M, n_list=3649, n_grp=234, n_reuse=32
wtime_BW_reuse32 = np.array([0.683791083827,
                            0.370261961095,
                            0.213497399729,
                            0.135115119046,
                            0.0959239787049,
                            0.0763284085342])



SPEED = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 8.0])

# N=4M, n_list=3649, n_grp=234, n_reuse=2
wtime_SPEED_reuse2 =  np.array([0.607602280932,
                                0.412952573149,
                                0.315627719257,
                                0.266965292312,
                                0.242634078839,
                                0.230468472102])

# N=4M, n_list=3649, n_grp=234, n_reuse=8
wtime_SPEED_reuse8 =  np.array([0.428477188042,
                                0.298774705104,
                                0.233923463635,
                                0.201497842901,
                                0.185285032533,
                                0.17717862735])

# N=4M, n_list=3649, n_grp=234, n_reuse=32
wtime_SPEED_reuse32 =  np.array([0.38369591482,
                                 0.270230238093,
                                 0.213497399729,
                                 0.185130980548,
                                 0.170947770957,
                                 0.163856166161])

n_proc = ([32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384])
# Ntot=1G, n_smp=500, n_dc=8, n_grp=256, n_reuse=8, bw=1, speed=1
wtime_para_bw1_speed1 = ([2.06480148514,
                          1.03312063153,
                          0.51757192055,
                          0.260293024102,
                          0.132505557895,
                          0.070095710481,
                          0.0415094602088,
                          0.0318995829267,
                          0.0355814343849,
                          0.052999174976])

# Ntot=1G, n_smp=500, n_dc=8, n_grp=256, n_reuse=8, bw=8, speed=1
wtime_para_bw8_speed1 = ([0.844674461643,
                          0.423030245622,
                          0.212499767781,
                          0.107729920998,
                          0.0561969273457,
                          0.0319142752827,
                          0.0223915905766,
                          0.022313470833,
                          0.0307611811774,
                          0.050561835528])

# Ntot=1G, n_smp=500, n_dc=8, n_grp=256, n_reuse=8, bw=1, speed=8
wtime_para_bw1_speed8 = ([1.4785599289,
                          0.739779327269,
                          0.370680039433,
                          0.186625305534,
                          0.095449491615,
                          0.0513451345124,
                          0.031911365909,
                          0.026877522307,
                          0.0328472274478,
                          0.0514087661823])

plt.plot(n_proc, wtime_para_bw1_speed1,  linewidth=3.0, linestyle='solid',  label="N=1G, n_grp=256, theta=0.5")
plt.plot(n_proc, wtime_para_bw8_speed1,  linewidth=3.0, linestyle='dashed', label="N=1G, n_grp=256, theta=0.5, bw*8")
plt.plot(n_proc, wtime_para_bw1_speed8,  linewidth=3.0, linestyle='dashed', label="N=1G, n_grp=256, theta=0.5, speed*8")
plt.grid()
plt.rcParams["font.size"] = 16
plt.legend(loc="top right")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# of processes")
plt.ylabel("wall clock time [sec]")
plt.xlim([30, 2e4])
plt.savefig("parallel.eps", dpi=150)
sys.exit()

"""
plt.plot(SPEED, wtime_SPEED_reuse2,  linewidth=3.0, linestyle='solid', label="n_reuse=2")
plt.plot(SPEED, wtime_SPEED_reuse8,  linewidth=3.0, linestyle='dashed', label="n_reuse=8")
plt.plot(SPEED, wtime_SPEED_reuse32, linewidth=3.0, linestyle='dotted', label="n_reuse=32")

plt.grid()
plt.rcParams["font.size"] = 16
plt.legend(loc="top right")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("speed up factor of floating operation")
plt.ylabel("wall clock time [sec]")
plt.xlim([0.1, 10])
plt.ylim([0.1,1])
#plt.show()
plt.savefig("perf_model_speed.eps", dpi=150)
"""

plt.plot(BW, wtime_BW_reuse2,  linewidth=3.0, linestyle='solid', label="n_reuse=2")
plt.plot(BW, wtime_BW_reuse8,  linewidth=3.0, linestyle='dashed', label="n_reuse=8")
plt.plot(BW, wtime_BW_reuse32, linewidth=3.0, linestyle='dotted', label="n_reuse=32")

plt.grid()
plt.rcParams["font.size"] = 16
plt.legend(loc="lower left")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("speed up factor of bandwidth")
plt.ylabel("wall clock time [sec]")
plt.xlim([0.1, 10])
plt.ylim([5e-2,1])
#plt.show()
plt.savefig("perf_model_bw.eps", dpi=150)

