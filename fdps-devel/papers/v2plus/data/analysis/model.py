import sys
import math
import numpy as np
import matplotlib.pylab as plt
import matplotlib.colors as clr
from matplotlib import ticker

#FLAG_PLOT = "COMPARE_WITH_MEASUREMENT"
#FLAG_PLOT = "COMPARE_SINGLE"
#FLAG_PLOT = "COMPARE_N_REUSE"
#FLAG_PLOT = "2D_MAP_TO_1D"
FLAG_PLOT = "2D_MAP"
#FLAG_PLOT = "BREAK_DOWN_PARALLEL"

#FLAG_PLOT = "BREAK_DOWN_NORMAL"
#FLAG_PLOT = "BREAK_DOWN_FORCE"
#FLAG_PLOT = "BREAK_DOWN_CHANGING_NGR"


FLAG_OUT_PUT_FILE = True
#FLAG_OUT_PUT_FILE = False

#FLAG_BOOST_CPU_MM = True
FLAG_BOOST_CPU_MM = False

#FLAG_BOOST_GPU_MM = True
FLAG_BOOST_GPU_MM = False

#FLAG_BOOST_GPU_SPEED = True
FLAG_BOOST_GPU_SPEED = False

#FLAG_BOOST_TRANS = True
FLAG_BOOST_TRANS = False

N_LEAF = 4.0

if(FLAG_PLOT == "COMPARE_WITH_MEASUREMENT"):
    N_LOC = 2**22
    BW_CPU_MM  = 4.0e10  # 40 GB/s
    BW_TRANS   = 1.2e10  # 12 GB/s
    S_GPU      = 1.38e13 # 13.8Tflops
    BW_GPU_MM  = 5.5e11  # 550 GB/s

if(FLAG_PLOT == "BREAK_DOWN_NORMAL" or FLAG_PLOT == "BREAK_DOWN_CHANGING_NGR" or FLAG_PLOT == "BREAK_DOWN_FORCE" or FLAG_PLOT == "BREAK_DOWN_PARALLEL" \
   or FLAG_PLOT == "COMPARE_SINGLE" or FLAG_PLOT == "2D_MAP_TO_1D" \
   or FLAG_PLOT == "2D_MAP" or FLAG_PLOT == "COMPARE_N_REUSE"):
    N_LOC = 2**22
    BW_CPU_MM  = 1.0e11  # 100 GB/s
    BW_TRANS   = 1.0e10  # 10 GB/s
    S_GPU      = 1.0e13  # 10 Tflops
    BW_GPU_MM  = 5.0e11  # 500 GB/s
    
BW_INJ = 8.0e9  # 8[GB/s] ## injection bandwidth of k


if(FLAG_BOOST_CPU_MM == True):
    BW_CPU_MM = 4.0*BW_CPU_MM
if(FLAG_BOOST_GPU_MM == True):
    BW_GPU_MM = 4.0*BW_GPU_MM
if(FLAG_BOOST_GPU_SPEED == True):
    S_GPU     = 4.0*S_GPU
if(FLAG_BOOST_TRANS == True):
    BW_TRANS = 4.0*BW_TRANS
    
N_OP = 23.0
# 3sub 3mad 3mul 3mad 1sub 1rsqrt
# = 19 [+ 4(rsqrt)]
# = 23/26

TAU_CPU_MM     = 1.0 / BW_CPU_MM
TAU_TRANS      = 1.0 / BW_TRANS
TAU_GPU_KERNEL = 1.0 / S_GPU
TAU_GPU_MM     = 1.0 / BW_GPU_MM
TAU_INJ        = 1.0 / BW_INJ

def SetTau(f_bw_cpu_mm, f_bw_tans, f_s_gpu, f_bw_gpu_mm):
    global TAU_CPU_MM, TAU_TRANS, TAU_GPU_KERNEL, TAU_GPU_MM
    TAU_CPU_MM     = 1.0 / (BW_CPU_MM * f_bw_cpu_mm)
    TAU_TRANS      = 1.0 / (BW_TRANS * f_bw_tans)
    TAU_GPU_KERNEL = 1.0 / (S_GPU * f_s_gpu)
    TAU_GPU_MM     = 1.0 / (BW_GPU_MM * f_bw_gpu_mm)

b_FP     = 88
b_EPI    = 24
b_EPJ    = 32
b_FORCE  = 32
b_key    = 16
b_EPI_buf= 16
b_EPJ_buf= 16
b_FORCE_buf = 16
b_EPJ_GPU= 16
b_ID     = 4
b_ID_buf = 4
b_pos = 24
#def GetSGpu(ni):
#    return N_OP / (2.4 + b_EPJ_GPU/(BW_GPU_MM*1e-12)/ni) * 1e12

# for T_{root}
beta_root     = b_FP
alpha_root = 0.7

# for T_{const_lt}
beta_key_lt   = b_FP + b_key
alpha_key  = 0.85
beta_sort     = 24.0*b_key
alpha_sort = 1.1
beta_reorder_lt     = 2.0*b_FP + b_EPI + b_EPJ + b_key
alpha_reorder_lt = 1.1
beta_link     = 16
alpha_link = 3.4
#beta_link     = math.log(8*N_LEAF, 2)/N_LEAF*b_key
#alpha_link = 6.5

# for T_{mom_lt}
beta_mom     = b_EPJ
alpha_mom = 1.8

# for T_{const_gt}
beta_key_gt  = b_EPJ + b_key
beta_reorder_gt_const = 2.0*b_EPJ + 4
alpha_reorder_gt_const = 5.0
beta_reorder_gt_reuse = 2.0*b_EPJ + 4
alpha_reorder_gt_reuse = 1.0


# for T_{mom_gt}

# for T_{force}
beta_cp_all      = b_EPJ + b_EPJ_buf
alpha_cp_all  = 0.87
beta_H2D_all     = b_EPJ_buf
alpha_H2D_all = 1.0

beta_cp_EPI      = b_EPI + b_EPI_buf
alpha_cp_EPI  = 1.0
beta_H2D_EPI     = b_EPI_buf
alpha_H2D_EPI = 1.0

beta_cp_list      = b_ID
alpha_cp_list  = 1.0
beta_H2D_list     = b_ID_buf
alpha_H2D_list = 1.0

beta_cp_FORCE   = b_FORCE + b_FORCE_buf
alpha_cp_FORCE  = 1.2
beta_D2H_FORCE  = b_FORCE_buf
alpha_D2H_FORCE = 1.2

n_op = 23
alpha_gpu_kernel   = 1.7

beta_gpu_mm    = b_EPJ_GPU + b_ID
alpha_gpu_mm   = 2.7

beta_const_list  = b_ID
alpha_const_list = 19.0

beta_wb_int_cp   = b_FORCE + b_FP + b_EPI + b_EPJ
alpha_wb_int_cp  = 1.4


# for T_{dc}
#tau_dc_sort          = 24.0 / BW_MM
beta_dc_sort          = b_pos
alpha_dc_sort         = 7.7

tau_p2p_startup       = 1.0e-5
tau_allgather_startup = 1.2e-5
#tau_inj_ptcl          = 32.0 / BW_INJ

beta_allgather_ptcl  = b_EPJ
alpha_allgather_ptcl = 1.0

beta_p2p_ptcl        = b_EPJ
alpha_p2p_ptcl        = 1.7

def GetTRoot(n_loc, reuse):
    time = 0.0
    if(reuse==False):
        time = n_loc*alpha_root*beta_root*TAU_CPU_MM
    return time

def GetTConstLt(n_loc, reuse):
    time = 0.0
    if(reuse==False):
        #time = n_loc*(alpha_key*beta_key_lt + alpha_sort*beta_sort + alpha_reorder_lt*beta_reorder_lt + alpha_link*beta_link)*TAU_CPU_MM
        time = (n_loc*(alpha_key*beta_key_lt + alpha_sort*beta_sort + alpha_reorder_lt*beta_reorder_lt) + alpha_link*beta_link*n_loc/4.0*math.log2(n_loc/4.0))*TAU_CPU_MM
    return time

def GetTMom(n, reuse):
    return n*(alpha_mom*beta_mom)*TAU_CPU_MM

def GetTConstGt(n_loc, n_let, reuse):
    n_glb = n_loc + n_let
    if(reuse==True):
        time = n_glb*alpha_reorder_gt_reuse*beta_reorder_gt_reuse*TAU_CPU_MM
    if(reuse==False):
        #time = (n_let*alpha_key*beta_key_gt + n_glb*(alpha_link*beta_link + alpha_sort*beta_sort + alpha_reorder_gt_const*beta_reorder_gt_const))*TAU_CPU_MM
        time = (n_let*alpha_key*beta_key_gt + n_glb*(alpha_sort*beta_sort + alpha_reorder_gt_const*beta_reorder_gt_const) + alpha_link*beta_link*n_glb/4.0*math.log2(n_glb/4.0) )*TAU_CPU_MM
    return time

def GetTForce(n_loc, n_let, n_list, ni, reuse):
        n_glb = n_loc + n_let
        t_cp_all   = n_glb * alpha_cp_all  * beta_cp_all * TAU_CPU_MM
        t_H2D_all  = n_glb * alpha_H2D_all * beta_H2D_all * TAU_TRANS
        t_cp_EPI   = n_loc * alpha_cp_EPI  * beta_cp_EPI * TAU_CPU_MM
        t_H2D_EPI  = n_loc * alpha_H2D_EPI * beta_H2D_EPI * TAU_TRANS
        t_cp_list  = n_loc*n_list/ni * alpha_cp_list  * beta_cp_list * TAU_CPU_MM
        t_H2D_list = n_loc*n_list/ni * alpha_H2D_list * beta_H2D_list * TAU_TRANS
        t_cp_FORCE  = n_loc * alpha_cp_FORCE  * beta_cp_FORCE  * TAU_CPU_MM
        t_D2H_FORCE = n_loc * alpha_D2H_FORCE * beta_D2H_FORCE * TAU_TRANS
        t_kernel    = n_loc*n_list*(alpha_gpu_kernel*n_op*TAU_GPU_KERNEL + alpha_gpu_mm*beta_gpu_mm*TAU_GPU_MM/ni)
        t_kernel_trans = n_loc*n_list/ni*alpha_gpu_mm*beta_gpu_mm*TAU_GPU_MM
        t_kernel_calc  = n_loc*n_list*alpha_gpu_kernel*n_op*TAU_GPU_KERNEL
        t_const_list = n_loc*n_list/ni * alpha_const_list * beta_const_list * TAU_CPU_MM
        t_wb_int_cp = n_loc * alpha_wb_int_cp * beta_wb_int_cp * TAU_CPU_MM
        if(reuse==True):
            time = max(t_cp_all, t_H2D_all) + max(t_kernel, t_H2D_EPI, t_D2H_FORCE, t_cp_EPI+t_wb_int_cp) + t_cp_FORCE
        if(reuse==False):
            time = max(t_cp_all, t_H2D_all) + max(t_kernel, t_H2D_EPI+t_H2D_list, t_D2H_FORCE, t_const_list+t_cp_EPI+t_cp_list+t_wb_int_cp) + t_cp_FORCE
        return time, t_cp_all, t_H2D_all, t_cp_EPI, t_H2D_EPI, t_cp_list, t_H2D_list, t_cp_FORCE, t_D2H_FORCE, t_const_list, t_kernel, t_wb_int_cp, t_kernel_trans, t_kernel_calc
    
def GetTstep(n_loc, n_list, ni, n_let, reuse):
    time_root        = GetTRoot(n_loc, reuse)
    time_const_lt    = GetTConstLt(n_loc, reuse)
    time_mom_lt      = GetTMom(n_loc, reuse)
    time_const_gt    = GetTConstGt(n_loc, n_let, reuse)
    time_mom_gt      = GetTMom(n_loc+n_let, reuse)
    time_force_ar    = GetTForce(n_loc, n_let, n_list, ni, reuse)
    time_force = time_force_ar[0]
    time_total = time_root + time_const_lt + time_mom_lt + time_const_gt + time_mom_gt + time_force
    #if(reuse==False):
    #    print("time_root= ", time_root, "time_const_lt= ", time_const_lt, "time_mom_lt= ", time_mom_lt, "time_const_gt= ", time_const_gt, "time_mom_gt= ", time_mom_gt, "time_force_ar[0]= ", time_force_ar[0])
    return time_total, time_root, time_const_lt, time_mom_lt, time_const_gt, time_mom_gt, time_force

def GetNlist(n_tot, theta, ni):
    inv_theta = 1.0 / theta
    inv_theta_2 = inv_theta * inv_theta
    inv_theta_3 = inv_theta * inv_theta_2
    n_list = ni + 14.0*inv_theta*pow(ni, (2.0/3.0)) \
             + 21.0*inv_theta_2*math.pi*pow(ni, (1.0/3.0)) \
             + 28.0/3.0*inv_theta_3*math.pi*math.log(theta/2.8*(pow(n_tot,(1.0/3.0)) - pow(ni, 1.0/3.0)), 2)
    return n_list

#### FOR PARALLEL MODEL
def GetParallelTime(n_loc, n_p, n_p_close, n_smp, n_let, n_dc):
    time_dc_gather = pow(n_p, (1.0/6.0))*tau_allgather_startup + n_smp*pow(n_p, (2.0/3.0))*alpha_allgather_ptcl*b_pos*TAU_INJ
    time_dc_sort = n_smp*pow(n_p, (2.0/3.0))*math.log(pow(n_smp, 3.0)*pow(n_p, (5.0/3.0)), 2)*alpha_dc_sort*beta_dc_sort*TAU_CPU_MM
    #time_dc = (time_dc_gather + time_dc_sort) / n_dc
    time_dc = time_dc_gather + time_dc_sort
    time_const_let = (n_let-n_p+n_p_close)*alpha_const_list*beta_const_list*TAU_CPU_MM
    time_allgather_let = pow(n_p, (1.0/4.0))*tau_allgather_startup + n_p*alpha_allgather_ptcl*beta_allgather_ptcl*TAU_INJ
    time_p2p_let = n_p_close*tau_p2p_startup + (n_let-n_p+n_p_close) * alpha_p2p_ptcl* beta_p2p_ptcl * TAU_INJ
    time_let = time_const_let + time_allgather_let + time_p2p_let
    time = time_dc + time_let
    return time, time_dc_gather, time_dc, time_const_let, time_allgather_let, time_p2p_let
    #return time, time_dc_gather, time_dc_sort, time_const_let, time_allgather_let, time_p2p_let
    
def GetNProcClose(theta):
    inv_theta = 1.0 / theta
    return 6*(inv_theta+1) + 3.0*math.pi*pow((inv_theta+1), 2.0) + 4.0*math.pi/3.0*pow((inv_theta+1), 3.0)
    
def GetNlet(n_loc, theta, n_p):
    inv_theta = 1.0 / theta
    inv_theta_2 = inv_theta * inv_theta
    inv_theta_3 = inv_theta * inv_theta_2
    n_p_close = GetNProcClose(theta)
    n_list = 14.0*inv_theta*pow(n_loc, (2.0/3.0)) \
             + 21.0*inv_theta_2*math.pi*pow(n_loc, (1.0/3.0)) \
             + 28.0/3.0*inv_theta_3*math.pi*math.log(theta/2.8*(pow(n_loc*n_p,(1.0/3.0)) - pow(n_loc, 1.0/3.0)), 2)
    n_list = n_list + n_p - n_p_close # new line
    return n_list

#### FOR PARALLEL MODEL

#####################################
# compare model with estimaed value #
#####################################
if(FLAG_PLOT == "COMPARE_WITH_MEASUREMENT"):
    ni_ar = np.logspace(0, 20, 100, base=2)
    
    ni_measure_ar = np.array([30, 96, 168, 225, 438, 892, 1477, 3498, 6864, 10280])
    wtime_total_const_measure_ar = np.array([1.02057448215783, 0.58665977511555, 0.517638379707932, 0.425352575723082, 0.385704821906984, 0.363689037039876, 0.368042776826769, 0.504200789146125, 0.62384789576754, 0.652912518940866])
    wtime_total_reuse_measure_ar = np.array([0.118371980730444, 0.10267551895231, 0.096680029295384, 0.094520645681768, 0.134411729872227, 0.151146794203669, 0.157689225394279, 0.295281784143299, 0.402100322302431, 0.436437553726137])

    #wtime_force_const_measure_ar = np.array([0.809013470076025,  0.26522418949753,   0.30648113694042,   0.21132396813482,   0.171324928756803,  ])
    #wtime_force_reuse_measure_ar = np.array([0.0962629169225693, 0.0807113572955132, 0.0746883279643953, 0.0724211800843477, 0.111294223926961,  ])

    
    wtime_total_const_ar    = np.zeros(len(ni_ar))
    wtime_root_const_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_const_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_const_ar   = np.zeros(len(ni_ar))
    wtime_force_const_ar     = np.zeros(len(ni_ar))

    wtime_total_reuse_ar    = np.zeros(len(ni_ar))
    wtime_root_reuse_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_force_reuse_ar    = np.zeros(len(ni_ar))
    
    for id in range(len(ni_ar)):
        n_loc = float(N_LOC)
        theta = 0.5
        ni = ni_ar[id]
        n_let = 0
        n_list = GetNlist(n_loc, theta, ni)
        #print(n_loc, theta, ni, n_list, n_loc*n_list/ni)
        wtime_total_const_ar[id], wtime_root_const_ar[id], wtime_const_lt_const_ar[id], wtime_mom_lt_const_ar[id], wtime_const_gt_const_ar[id], \
            wtime_mom_gt_const_ar[id], wtime_force_const_ar[id] = GetTstep(n_loc, n_list, ni, n_let, False)
        wtime_total_reuse_ar[id], wtime_root_reuse_ar[id], wtime_const_lt_reuse_ar[id], wtime_mom_lt_reuse_ar[id], wtime_const_gt_reuse_ar[id], \
            wtime_mom_gt_reuse_ar[id], wtime_force_reuse_ar[id] = GetTstep(n_loc, n_list, ni, n_let, True)

    plt.plot(ni_ar, wtime_total_const_ar, linewidth=5.0, linestyle='solid',   color="k", label=r"$T_{\rm step, const}$")
    plt.scatter(ni_measure_ar, wtime_total_const_measure_ar, color="k", marker="o", s=200, label=r"$T_{\rm step, const}$(measurment)")
    plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=5.0, linestyle='dashed',   color="r", label=r"$T_{\rm step, reuse}$")
    plt.scatter(ni_measure_ar, wtime_total_reuse_measure_ar, color="r", marker="^", s=200, label=r"$T_{\rm step, reuse}$(measurment)")
    
    plt.grid()
    plt.rcParams["font.size"] = 16
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\langle n_{\rm i} \rangle$", fontsize=16)
    plt.ylabel("wall clock time [sec]", fontsize=16)
    plt.xlim([1e1,2e4])
    plt.ylim([0.06,4])
    plt.tick_params(which='major', width = 2, length = 14, labelsize=16)
    plt.tick_params(which='minor', width = 1, length = 7)
    plt.legend(bbox_to_anchor=(1, 1), ncol=1, labelspacing=0.0, columnspacing=0.1, fontsize=16)

    output_file_name = "compare.eps"
    if(FLAG_OUT_PUT_FILE == True):
        plt.savefig(output_file_name, dpi=150)
    else:
        plt.show()
            
    sys.exit()
#####################################
# compare model with estimaed value #
#####################################


###################
# normal perf fig #
###################
if(FLAG_PLOT == "BREAK_DOWN_NORMAL"):
    ni_ar = np.logspace(0, 20, 100, base=2)

    wtime_total_const_ar    = np.zeros(len(ni_ar))
    wtime_root_const_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_const_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_const_ar   = np.zeros(len(ni_ar))
    wtime_force_const_ar     = np.zeros(len(ni_ar))

    wtime_total_reuse_ar    = np.zeros(len(ni_ar))
    wtime_root_reuse_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_force_reuse_ar    = np.zeros(len(ni_ar))
    
    for id in range(len(ni_ar)):
        n_loc = float(N_LOC)
        theta = 0.5
        ni = ni_ar[id]
        n_let = 0
        n_list = GetNlist(n_loc, theta, ni)
        wtime_total_const_ar[id], wtime_root_const_ar[id], wtime_const_lt_const_ar[id], wtime_mom_lt_const_ar[id], wtime_const_gt_const_ar[id], \
            wtime_mom_gt_const_ar[id], wtime_force_const_ar[id] = GetTstep(n_loc, n_list, ni, n_let, False)
        wtime_total_reuse_ar[id], wtime_root_reuse_ar[id], wtime_const_lt_reuse_ar[id], wtime_mom_lt_reuse_ar[id], wtime_const_gt_reuse_ar[id], \
            wtime_mom_gt_reuse_ar[id], wtime_force_reuse_ar[id] = GetTstep(n_loc, n_list, ni, n_let, True)

    ni_left_tmp_ar     = np.array([10, 20])
    ni_right_tmp_ar     = np.array([5e3, 1e4])
    
    wtime_const_lt_const_tmp_ar  = np.array([wtime_const_lt_const_ar[0], wtime_const_lt_const_ar[0]])
    plt.plot(ni_left_tmp_ar, wtime_const_lt_const_tmp_ar,     linewidth=10.0, linestyle='solid', color="m")
    plt.text(12, wtime_const_lt_const_tmp_ar[0]*0.65, r"$T_{\rm const\, lt}$, construction step", color="m", fontsize=16)

    wtime_const_gt_const_tmp_ar  = np.array([wtime_const_gt_const_ar[0], wtime_const_gt_const_ar[0]])
    plt.plot(ni_right_tmp_ar, wtime_const_gt_const_tmp_ar,     linewidth=10.0, linestyle='solid', color="g")
    plt.text(4e2, wtime_const_gt_const_tmp_ar[0]*0.7, r"$T_{\rm const\, gt}$, construction step", color="g", fontsize=16)

    wtime_other_const_tmp_ar = np.array([wtime_total_const_ar[0]-wtime_force_const_ar[0]-wtime_const_gt_const_ar[0]-wtime_const_lt_const_ar[0],
    wtime_total_const_ar[0]-wtime_force_const_ar[0]-wtime_const_gt_const_ar[0]-wtime_const_lt_const_ar[0]])
    plt.plot(ni_left_tmp_ar, wtime_other_const_tmp_ar, linewidth=10.0, linestyle='solid', color="b")
    plt.text(12, wtime_other_const_tmp_ar[0]*1.2, r"$T_{\rm root}+T_{\rm mom\, lt}+T_{\rm mom\, gt}$, construction step", color="b", fontsize=16)

    #print(wtime_const_gt_reuse_ar[0])
    wtime_const_gt_reuse_tmp_ar  = np.array([wtime_const_gt_reuse_ar[0], wtime_const_gt_reuse_ar[0]])
    plt.plot(ni_right_tmp_ar, wtime_const_gt_reuse_tmp_ar,     linewidth=3.0, linestyle='solid', color="g")
    plt.text(5e2, wtime_const_gt_reuse_tmp_ar[0]*1.1, r"$T_{\rm const\, gt}$, reusing step", color="g", fontsize=16)

    wtime_other_reuse_tmp_ar = np.array([wtime_total_reuse_ar[0]-wtime_force_reuse_ar[0]-wtime_const_gt_reuse_ar[0]-wtime_const_lt_reuse_ar[0],
    wtime_total_reuse_ar[0]-wtime_force_reuse_ar[0]-wtime_const_gt_reuse_ar[0]-wtime_const_lt_reuse_ar[0]])
    plt.plot(ni_right_tmp_ar, wtime_other_reuse_tmp_ar, linewidth=3.0, linestyle='solid', color="b")
    plt.text(3e2, wtime_other_reuse_tmp_ar[0]*1.1, r"$T_{\rm mom\, lt}+T_{\rm mom\, gt}$, reusing step", color="b", fontsize=16)
    
    plt.plot(ni_ar, wtime_total_const_ar,          linewidth=10.0, linestyle='solid',   color="k", label=r"$T_{\rm step, single}$, construction step")
    plt.plot(ni_ar, wtime_force_const_ar,          linewidth=10.0, linestyle='dashed',   color="r", label=r"$T_{\rm force}$, construction step")

    plt.plot(ni_ar, wtime_total_reuse_ar,  linewidth=3.0, linestyle='solid',  color="k", label=r"$T_{\rm step, single}$, reusing step")
    plt.plot(ni_ar, wtime_force_reuse_ar,  linewidth=3.0, linestyle='dashed', color="r", label=r"$T_{\rm force}$, reusing step")

    plt.grid()
    plt.rcParams["font.size"] = 16
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\langle n_{\rm i} \rangle$", fontsize=16)
    plt.ylabel("wall clock time [sec]", fontsize=16)
    plt.xlim([1e1,1e4])
    plt.ylim([2e-3,4])
    plt.tick_params(which='major', width = 2, length = 12)
    plt.tick_params(which='minor', width = 1, length = 7)    
    plt.legend(bbox_to_anchor=(0.9, 1), ncol=1, labelspacing=0.0, columnspacing=0.1, fontsize=16)

    str0 = "total_break_down"

    output_file_name = str0 + ".eps"        
    
    if(FLAG_OUT_PUT_FILE == True):
        print(output_file_name)
        plt.savefig(output_file_name, dpi=150)
    else:
        plt.show()
            
    sys.exit()
###################
# normal perf fig #
###################


####################
# FORCE BREAK DOWN #
####################
if(FLAG_PLOT == "BREAK_DOWN_FORCE"):
    ni_ar = np.logspace(0, 20, 100, base=2)
    
    wtime_force_total_const_ar     = np.zeros(len(ni_ar))
    wtime_force_cp_all_const_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_all_const_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_EPI_const_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_EPI_const_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_list_const_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_list_const_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_FORCE_const_ar  = np.zeros(len(ni_ar))
    wtime_force_D2H_FORCE_const_ar  = np.zeros(len(ni_ar))
    wtime_force_const_list_const_ar  = np.zeros(len(ni_ar))
    wtime_force_kernel_const_ar      = np.zeros(len(ni_ar))
    wtime_force_wb_int_cp_const_ar   = np.zeros(len(ni_ar))
    wtime_force_kernel_trans_const_ar = np.zeros(len(ni_ar))
    wtime_force_kernel_calc_const_ar = np.zeros(len(ni_ar))
    
    wtime_force_total_reuse_ar     = np.zeros(len(ni_ar))
    wtime_force_cp_all_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_all_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_EPI_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_EPI_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_list_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_list_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_FORCE_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_D2H_FORCE_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_const_list_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_kernel_reuse_ar      = np.zeros(len(ni_ar))
    wtime_force_wb_int_cp_reuse_ar   = np.zeros(len(ni_ar))
    wtime_force_kernel_trans_reuse_ar = np.zeros(len(ni_ar))
    wtime_force_kernel_calc_reuse_ar = np.zeros(len(ni_ar))
    
    for id in range(len(ni_ar)):
        n_loc = float(N_LOC)
        theta = 0.5
        ni = ni_ar[id]
        n_let = 0
        n_list = GetNlist(n_loc, theta, ni)
        wtime_force_total_const_ar[id], wtime_force_cp_all_const_ar[id], wtime_force_H2D_all_const_ar[id], \
            wtime_force_cp_EPI_const_ar[id], wtime_force_H2D_EPI_const_ar[id], \
            wtime_force_cp_list_const_ar[id], wtime_force_H2D_list_const_ar[id], \
            wtime_force_cp_FORCE_const_ar[id], wtime_force_D2H_FORCE_const_ar[id], \
            wtime_force_const_list_const_ar[id], wtime_force_kernel_const_ar[id], \
            wtime_force_wb_int_cp_const_ar[id], \
            wtime_force_kernel_trans_const_ar[id], wtime_force_kernel_calc_const_ar[id]\
            = GetTForce(n_loc, n_let, n_list, ni, False)

        wtime_force_total_reuse_ar[id], wtime_force_cp_all_reuse_ar[id], wtime_force_H2D_all_reuse_ar[id], \
            wtime_force_cp_EPI_reuse_ar[id], wtime_force_H2D_EPI_reuse_ar[id], \
            wtime_force_cp_list_reuse_ar[id], wtime_force_H2D_list_reuse_ar[id], \
            wtime_force_cp_FORCE_reuse_ar[id], wtime_force_D2H_FORCE_reuse_ar[id], \
            wtime_force_const_list_reuse_ar[id], wtime_force_kernel_reuse_ar[id], \
            wtime_force_wb_int_cp_reuse_ar[id], \
            wtime_force_kernel_trans_reuse_ar[id], wtime_force_kernel_calc_reuse_ar[id]\
            = GetTForce(n_loc, n_let, n_list, ni, True)
        
    plt.plot(ni_ar, wtime_force_total_const_ar, linewidth=10.0, linestyle='solid',   color="k", label=r"$T_{\rm force}$, construction step")
    plt.plot(ni_ar, wtime_force_total_reuse_ar, linewidth=3.0, linestyle='solid',   color="k", label=r"$T_{\rm force}$, reusing step")
    plt.plot(ni_ar, wtime_force_kernel_const_ar,      linewidth=10.0, linestyle='dashed',  color="r", label=r"$T_{\rm kernel}$")
    plt.plot(ni_ar, wtime_force_const_list_const_ar,  linewidth=3.0, linestyle='dashed',  color="r", label=r"$T_{\rm const\, list}$")
    plt.plot(ni_ar, wtime_force_cp_list_const_ar,  linewidth=3.0, linestyle='dotted',  color="g", label=r"$T_{\rm cp\, list}$")
    plt.plot(ni_ar, wtime_force_H2D_list_reuse_ar, linewidth=3.0, linestyle='dashdot', color="b", label=r"$T_{\rm H2D\, list}$")
    
    #plt.plot(ni_ar, wtime_force_H2D_EPI_const_ar+wtime_force_H2D_list_reuse_ar, linewidth=10.0, linestyle='dotted', color="g", label=r"$T_{\rm H2D\, EPI}+T_{\rm H2D\, list}$")

    #plt.plot(ni_ar, wtime_force_cp_EPI_reuse_ar+wtime_force_wb_int_cp_reuse_ar,  linewidth=3.0, linestyle='dashdot',  color="b", label=r"$T_{\rm cp\, EPI}+T_{\rm wb+int+cp}$, reusing step")
    
    #plt.plot(ni_ar, wtime_force_const_list_const_ar,  linewidth=10.0, linestyle='dashed',  color="r", label=r"$T_{\rm const\, list}$, construction step")
    #plt.plot(ni_ar, wtime_force_H2D_list_const_ar,  linewidth=10.0, linestyle='dotted',  color="g", label=r"$T_{\rm H2D\, list}$, construction step")
    #plt.plot(ni_ar, wtime_force_kernel_trans_const_ar, linewidth=3.0, linestyle='dotted',  color="r", label=r"$T_{\rm kernel\, trans}$")
    #plt.plot(ni_ar, wtime_force_kernel_calc_const_ar,  linewidth=3.0, linestyle='dashdot',  color="g", label=r"$T_{\rm kernel\, calc}$")
    
    #plt.plot(ni_ar, wtime_force_kernel_calc_const_ar+wtime_force_kernel_trans_const_ar,  linewidth=1.0, linestyle='dashdot',  color="m", label="xxx")

    

    ni_left_tmp_ar     = np.array([10, 20])
    ni_right_tmp_ar     = np.array([5e3, 1e4])

    wtime_force_cp_EPI_reuse_tmp_ar  = np.array([wtime_force_cp_EPI_reuse_ar[0], wtime_force_cp_EPI_reuse_ar[0]])
    wtime_force_wb_int_cp_reuse_tmp_ar  = np.array([wtime_force_wb_int_cp_reuse_ar[0], wtime_force_wb_int_cp_reuse_ar[0]])
    plt.plot(ni_left_tmp_ar, wtime_force_cp_EPI_reuse_tmp_ar+wtime_force_wb_int_cp_reuse_tmp_ar, linewidth=3.0, linestyle='solid', color="k")
    plt.text(12, (wtime_force_cp_EPI_reuse_tmp_ar[0]+wtime_force_wb_int_cp_reuse_tmp_ar[0])*1.2, r"$T_{\rm cp\, EPI}+T_{\rm wb+int+cp}$", color="k", fontsize=16)
    
    wtime_force_D2H_FORCE_const_tmp_ar  = np.array([wtime_force_D2H_FORCE_const_ar[0], wtime_force_D2H_FORCE_const_ar[0]])
    plt.plot(ni_left_tmp_ar, wtime_force_D2H_FORCE_const_tmp_ar, linewidth=3.0, linestyle='solid', color="r")
    plt.text(20, wtime_force_D2H_FORCE_const_tmp_ar[0]*0.9, r"$T_{\rm D2H\, FORCE}$", color="r", fontsize=16)

    wtime_force_H2D_all_reuse_tmp_ar  = np.array([wtime_force_H2D_all_reuse_ar[0], wtime_force_H2D_all_reuse_ar[0]])
    plt.plot(ni_left_tmp_ar, wtime_force_H2D_all_reuse_tmp_ar, linewidth=3.0, linestyle='solid', color="g")
    plt.text(12, wtime_force_H2D_all_reuse_tmp_ar[0]*0.7, r"$T_{\rm H2D\, all}$", color="g", fontsize=16)

    wtime_force_cp_all_reuse_tmp_ar  = np.array([wtime_force_cp_all_reuse_ar[0], wtime_force_cp_all_reuse_ar[0]])
    plt.plot(ni_left_tmp_ar, wtime_force_cp_all_reuse_tmp_ar, linewidth=3.0, linestyle='solid', color="b")
    plt.text(12, wtime_force_cp_all_reuse_tmp_ar[0]*1.2, r"$T_{\rm cp\, all}$", color="b", fontsize=16)
    
    wtime_force_H2D_EPI_reuse_tmp_ar  = np.array([wtime_force_H2D_EPI_reuse_ar[0], wtime_force_H2D_EPI_reuse_ar[0]])
    plt.plot(ni_right_tmp_ar, wtime_force_H2D_EPI_reuse_tmp_ar, linewidth=3.0, linestyle='solid', color="k")
    plt.text(2e3, wtime_force_H2D_EPI_reuse_tmp_ar[0]*0.9, r"$T_{\rm H2D\, EPI}$", color="k", fontsize=16)
    
    wtime_force_cp_FORCE_reuse_tmp_ar  = np.array([wtime_force_cp_FORCE_reuse_ar[0], wtime_force_cp_FORCE_reuse_ar[0]])
    plt.plot(ni_right_tmp_ar, wtime_force_cp_FORCE_reuse_tmp_ar, linewidth=3.0, linestyle='solid', color="r")
    plt.text(1.8e3, wtime_force_cp_FORCE_reuse_tmp_ar[0]*0.9, r"$T_{\rm cp\, FORCE}$", color="r", fontsize=16)




    
    #plt.plot(ni_ar, wtime_force_cp_EPI_reuse_ar+wtime_force_wb_int_cp_reuse_ar,  linewidth=3.0, linestyle='dashdot',  color="b", label=r"$T_{\rm cp\, EPI}+T_{\rm wb+int+cp}$, reusing step")
    #plt.plot(ni_ar, wtime_force_H2D_EPI_reuse_ar, linewidth=3.0, linestyle='dotted', color="g", label=r"$T_{\rm H2D\, EPI}$, reusing. step")
    #plt.plot(ni_ar, wtime_force_D2H_FORCE_const_ar, linewidth=1.0, linestyle='solid', color="k", label=r"$T_{\rm D2H\, FORCE}$") # const value
    #plt.plot(ni_ar, wtime_force_cp_all_const_ar, linewidth=1.0, linestyle='solid', color="r", label=r"$T_{\rm cp\, all}$") # const value
    #plt.plot(ni_ar, wtime_force_H2D_all_const_ar, linewidth=1.0, linestyle='solid', color="g", label=r"$T_{\rm H2D\, all}$") # const value
    #plt.plot(ni_ar, wtime_force_cp_FORCE_const_ar, linewidth=1.0, linestyle='solid', color="b", label=r"$T_{\rm cp\, FORCE}$") # const value
    
    
    plt.grid()
    plt.rcParams["font.size"] = 16
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\langle n_{\rm i} \rangle$", fontsize=16)
    plt.ylabel("wall clock time [sec]", fontsize=16)
    plt.xlim([1e1,1e4])
    plt.ylim([1e-3,4.0])
    plt.tick_params(which='major', width = 2, length = 14, labelsize=16)
    plt.tick_params(which='minor', width = 1, length = 7)
    plt.legend(bbox_to_anchor=(0.8, 1), ncol=1, labelspacing=0.0, columnspacing=0.5, fontsize=16)

    str0 = "force_break_down"

    output_file_name = str0 + ".eps"
    
    if(FLAG_OUT_PUT_FILE == True):
        print(output_file_name)
        plt.savefig(output_file_name, dpi=150)
    else:
        plt.show()

    sys.exit()
####################
# FORCE BREAK DOWN #
####################

################################
# PEFORMANCE FIG CHANGING NI#
################################
if(FLAG_PLOT == "BREAK_DOWN_CHANGING_NGR"):
    PARAM = "CPU_MM"
    #PARAM = "GPU_SPEED"
    #PARAM = "GPU_MM"
    #PARAM = "TRANS"

    ni_ar = np.logspace(0, 20, 100, base=2)
    size_of_ar = len(ni_ar)

    wtime_total_const_ar    = np.zeros(len(ni_ar))
    wtime_root_const_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_const_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_const_ar   = np.zeros(len(ni_ar))
    wtime_force_const_ar     = np.zeros(len(ni_ar))

    wtime_total_reuse_ar    = np.zeros(len(ni_ar))
    wtime_root_reuse_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_force_reuse_ar    = np.zeros(len(ni_ar))
    
    f_cpu_mm = f_trans = f_s_gpu = f_gpu_mm = 1.0

    for fact in 1.0, 2.0, 4.0, 8.0, 16.0, 1000000000.0:
        F_CPU_MM     = f_cpu_mm
        F_TRANS      = f_trans
        F_GPU_SPEED = f_s_gpu
        F_GPU_MM     = f_gpu_mm
        if(PARAM == "CPU_MM"):
            F_CPU_MM    = fact
        if(PARAM == "TRANS"):
            F_TRANS     = fact
        if(PARAM == "GPU_MM"):
            F_GPU_MM    = fact
        if(PARAM == "GPU_SPEED"):
            F_GPU_SPEED = fact
        for id in range(size_of_ar):
            SetTau(F_CPU_MM, F_TRANS, F_GPU_SPEED, F_GPU_MM)
            n_loc = float(N_LOC)
            ni = ni_ar[id]
            n_let = 0
            theta = 0.5
            n_list = GetNlist(n_loc, theta, ni)
            wtime_total_const_ar[id], wtime_root_const_ar[id], wtime_const_lt_const_ar[id], wtime_mom_lt_const_ar[id], wtime_const_gt_const_ar[id], \
                wtime_mom_gt_const_ar[id], wtime_force_const_ar[id] = GetTstep(n_loc, n_list, ni, n_let, False)
            wtime_total_reuse_ar[id], wtime_root_reuse_ar[id], wtime_const_lt_reuse_ar[id], wtime_mom_lt_reuse_ar[id], wtime_const_gt_reuse_ar[id], \
                wtime_mom_gt_reuse_ar[id], wtime_force_reuse_ar[id] = GetTstep(n_loc, n_list, ni, n_let, True)
    
        # FOR CPU MM
        if(PARAM == "CPU_MM"):
            if(fact==1.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=3.0, linestyle='solid',   color="k", \
                    label=r"$\tau_{\rm host\, mm}=10^{-11}{\rm s/B}\left( {\rm BW_{host\, mm}=100GB/s} \right)$")
                #plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=1.0, linestyle='solid',   color="k", label=r"${\rm BW_{host\, mm}=100GB/s}$, reusing step")
            if(fact==2.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=3.0, linestyle='dashed',  color="r", \
                    label=r"${\tau_{\rm host\, mm}=5 \times 10^{-12}{\rm s/B}}\left( {\rm BW_{host\, mm}=200GB/s} \right)$")
                #plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=1.0, linestyle='dashed',  color="r", label=r"${\rm BW_{host\, mm}=200GB/s}$, reusing step")
                
            if(fact==4.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=3.0, linestyle='dotted',  color="g", \
                    label=r"${\tau_{\rm host\, mm}=2.5 \times 10^{-12}{\rm s/B}}\left( {\rm BW_{host\, mm}=400GB/s} \right)$")
                #plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=1.0, linestyle='dotted',  color="g", label=r"${\rm BW_{host\, mm}=400GB/s}$, reusing step")

            if(fact>100.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=3.0, linestyle='dashdot',  color="b", \
                    label=r"${\tau_{\rm host\, mm}=0{\rm s/B}}\left( {\rm BW_{host\, mm}=\infty GB/s} \right)$")
                #plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=1.0, linestyle='dashdot',  color="b", label=r"${\rm BW_{host\, mm}=\infty GB/s}$, reusing step")
                
            #if(fact==8.0):
            #    plt.plot(ni_ar, wtime_total_const_ar, linewidth=10.0, linestyle='dashdot',  color="b", label=r"${\rm BW_{host\, mm}=800GB/s}$, construction step")

        # FOR SPEED
        if(PARAM == "GPU_SPEED"):
            if(fact==1.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=10.0, linestyle='solid',   color="k", \
                    label=r"${\tau_{\rm GPU\, speed}=10^{-13}{\rm s/B}} \left( S_{\rm GPU} =10{\rm Tflops} \right) $")
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='solid',   color="k")      
            if(fact==2.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=10.0, linestyle='dashed',  color="r", \
                    label=r"${\tau_{\rm GPU\, speed}=5 \times 10^{-14}{\rm s/B}} \left( S_{\rm GPU} =20{\rm Tflops} \right) $")
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='dashed',  color="r")
                
            if(fact==4.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=10.0, linestyle='dashed',  color="g", \
                    label=r"${\tau_{\rm GPU\, speed}=2.5 \times 10^{-14}{\rm s/B}} \left( S_{\rm GPU} =40{\rm Tflops} \right) $")
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='dashed',  color="g")
                
            #if(fact==8.0):
            #    plt.plot(ni_ar, wtime_total_const_ar, linewidth=10.0, linestyle='dotted',  color="g", label=r"${\tau_{GPU\, speed}=1.25 \times 10^{-14}}$, construction step")
            #    plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='dotted',  color="g", label=r"${\tau_{GPU\, speed}=1.25 \times 10^{-14}}$, reusing step")

            #if(fact==8.0):
            #    plt.plot(ni_ar, wtime_total_const_ar, linewidth=10.0, linestyle='dotted',  color="g", label=r"${\tau_{GPU\, speed}=6.25 \times 10^{-15}}$, construction step")
            #    plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='dotted',  color="g", label=r"${\tau_{GPU\, speed}=6.25 \times 10^{-15}}$, reusing step") 
            
            if(fact>100.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=10.0, linestyle='dashdot',  color="b", \
                         label=r"$\tau_{\rm GPU\, speed}=0.0{\rm s/B} \left( S_{\rm GPU} = \infty {\rm Tflops} \right) $")
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0,  linestyle='dashdot',  color="b")
        
        # FOR GPU MM
        if(PARAM == "GPU_MM"):
            if(fact==1.0):
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='solid',   color="k", \
                         label=r"$\tau_{\rm GPU\, mm}=2 \times 10^{-12}{\rm s/B} \left( {\rm BW_{GPU\, mm}=500GB/s} \right) $")
            if(fact==2.0):
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='dashed',  color="r", \
                         label=r"$\tau_{\rm GPU\, mm}=10^{-12}{\rm s/B} \left( {\rm BW_{GPU\, mm}=1TB/s} \right) $")
            if(fact==4.0):
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='dotted',  color="g", \
                         label=r"$\tau_{\rm GPU\, mm}=5 \times 10^{-13}{\rm s/B} \left( {\rm BW_{GPU\, mm}=2TB/s} \right)$")
            if(fact>100.0):
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='dashdot',  color="b", \
                         label=r"$\tau_{\rm GPU\, mm}=0.0{\rm s/B} \left( {\rm BW_{GPU\, mm}=\infty TB/s} \right)$")


        # FOR TRANS
        if(PARAM == "TRANS"):
            if(fact==1.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=10.0, linestyle='solid', color="k", label=r"${\rm BW_{trans}=10GB/s}$, construction step")
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0,  linestyle='solid', color="k", label=r"${\rm BW_{trans}=10GB/s}$, reusing step")
            if(fact==2.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=10.0, linestyle='solid', color="r", label=r"${\rm BW_{trans}=20GB/s}$, construction step")
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='dashed', color="r", label=r"${\rm BW_{trans}=20GB/s}$, reusing step")
            if(fact==4.0):
                plt.plot(ni_ar, wtime_total_const_ar, linewidth=10.0, linestyle='solid', color="g", label=r"${\rm BW_{trans}=40GB/s}$, construction step")
                plt.plot(ni_ar, wtime_total_reuse_ar, linewidth=3.0, linestyle='dotted', color="g", label=r"${\rm BW_{trans}=40GB/s}$, reusing step")
                

                
    plt.grid()
    plt.rcParams["font.size"] = 16
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\langle {\rm n_{i}} \rangle$", fontsize=16)
    plt.ylabel("wall clock time [sec]", fontsize=16)
    plt.tick_params(which='major', width = 2, length = 14, labelsize=16)
    plt.tick_params(which='minor', width = 1, length = 7)
    
    #FOE CPU MM
    if(PARAM == "CPU_MM"):
        plt.xlim([10,1e4])
        plt.ylim([6e-2,1.2])
        #plt.ylim([0.1,2.0])
        plt.legend(bbox_to_anchor=(0.95, 1), ncol=1, labelspacing=0.0, columnspacing=0.2, fontsize=16)
        output_file_name = "total_cpu_mm.eps"

    #FOE SPEED
    if(PARAM == "GPU_SPEED"):
        plt.xlim([10,1e4])
        plt.ylim([0.02, 2.0])
        plt.legend(bbox_to_anchor=(0.95, 1), ncol=1, labelspacing=0.0, columnspacing=0.2, fontsize=16)
        output_file_name = "total_gpu_speed.eps"
        
    #FOE GPU MM
    if(PARAM == "GPU_MM"):
        plt.xlim([1,1e4])
        plt.ylim([4e-2,1.0])
        plt.legend(bbox_to_anchor=(0.95, 1), ncol=1, labelspacing=0.0, columnspacing=0.2, fontsize=16)
        output_file_name = "total_gpu_mm.eps"
        
    #FOE PCI
    if(PARAM == "TRANS"):
        plt.xlim([10,1e4])
        plt.ylim([0.1,2.0])
        plt.legend(bbox_to_anchor=(0.8, 1), ncol=1, labelspacing=0.0, columnspacing=0.2, fontsize=16)
        output_file_name = "total_trans.eps"
        


    if(FLAG_OUT_PUT_FILE == True):
        print(output_file_name)
        plt.savefig(output_file_name, dpi=150)
    else:
        plt.show()
        
    sys.exit()
################################
# PEFORMANCE FIG CHANGING NI#
################################


#######################
# COMPARE SINGLE NODE #
#######################
if(FLAG_PLOT == "COMPARE_SINGLE" or FLAG_PLOT == "COMPARE_N_REUSE"):
    print("COMPARE_SINGLE or COMPARE_N_REUSE")
    
    ni_ar = np.logspace(0, 20, 100, base=2)
    size_of_ar = len(ni_ar)
    n_list_ar = np.zeros(len(ni_ar))
    
    wtime_total_const_ar    = np.zeros(len(ni_ar))
    wtime_root_const_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_const_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_const_ar   = np.zeros(len(ni_ar))
    wtime_force_const_ar     = np.zeros(len(ni_ar))

    wtime_total_reuse_ar    = np.zeros(len(ni_ar))
    wtime_root_reuse_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_force_reuse_ar    = np.zeros(len(ni_ar))

    wtime_force_total_const_ar     = np.zeros(len(ni_ar))
    wtime_force_cp_all_const_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_all_const_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_EPI_const_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_EPI_const_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_list_const_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_list_const_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_FORCE_const_ar  = np.zeros(len(ni_ar))
    wtime_force_D2H_FORCE_const_ar  = np.zeros(len(ni_ar))
    wtime_force_const_list_const_ar  = np.zeros(len(ni_ar))
    wtime_force_kernel_const_ar      = np.zeros(len(ni_ar))
    wtime_force_wb_int_cp_const_ar   = np.zeros(len(ni_ar))
    wtime_force_kernel_trans_const_ar = np.zeros(len(ni_ar))
    wtime_force_kernel_calc_const_ar = np.zeros(len(ni_ar))
    
    wtime_force_total_reuse_ar     = np.zeros(len(ni_ar))
    wtime_force_cp_all_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_all_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_EPI_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_EPI_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_list_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_list_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_FORCE_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_D2H_FORCE_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_const_list_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_kernel_reuse_ar      = np.zeros(len(ni_ar))
    wtime_force_wb_int_cp_reuse_ar   = np.zeros(len(ni_ar))
    wtime_force_kernel_trans_reuse_ar = np.zeros(len(ni_ar))
    wtime_force_kernel_calc_reuse_ar = np.zeros(len(ni_ar))
    
    wtime_total_ave_ar    = np.zeros(len(ni_ar))
    
    f_cpu_mm = f_trans = f_s_gpu = f_gpu_mm = 1.0

    t_min_gpu_speed_ar = np.zeros(10)
    arg_t_min_gpu_speed_ar = np.zeros(10)
    id_t_min_gpu_speed = 0

    #target_ele = "GPU_SPEED"
    #target_ele = "CPU_MM"
    #target_ele = "GPU_MM"
    #target_ele = "TRANS"
    target_ele = "GPU_SPEED-TRANS"
    
    #n_reuse = 16.0
    for n_reuse in 1.0, 4,0, 16.0:
        #for fact in 1.0, 10.0, 1000000000.0:
        #for fact in 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 10.0, 1000000000.0:
        for boost_ele in "GPU_SPEED", "CPU_MM", "GPU_MM", "CPU_MM-GPU_SPEED", "GPU_SPEED-GPU_MM", "TRANS", "GPU_SPEED-TRANS", "CPU_MM-GPU_SPEED-GPU_MM", "CPU_MM-GPU_SPEED-TRANS":
            for fact in 1.0, 2.0, 4.0, 0.5, 8.0, 16.0:
                if(boost_ele != target_ele and (fact != 1.0 and fact != 2.0 and fact != 4.0 and fact != 0.5) ):
                    continue
                
                F_CPU_MM     = f_cpu_mm
                F_TRANS      = f_trans
                F_GPU_SPEED = f_s_gpu
                F_GPU_MM     = f_gpu_mm
                if(boost_ele == "CPU_MM"):
                    F_CPU_MM    = fact
                elif(boost_ele == "TRANS"):
                    F_TRANS     = fact
                elif(boost_ele == "GPU_MM"):
                    F_GPU_MM    = fact
                elif(boost_ele == "GPU_SPEED"):
                    F_GPU_SPEED = fact
                elif(boost_ele == "CPU_MM-GPU_SPEED"):
                    F_GPU_SPEED = F_CPU_MM = fact
                elif(boost_ele == "GPU_SPEED-GPU_MM"):
                    F_GPU_SPEED = F_GPU_MM = fact                    
                elif(boost_ele == "CPU_MM-GPU_SPEED-GPU_MM"):
                    F_GPU_SPEED = F_CPU_MM = F_GPU_MM = fact
                elif(boost_ele == "CPU_MM-GPU_SPEED-TRANS"):
                    F_GPU_SPEED = F_CPU_MM = F_TRANS = fact
                elif(boost_ele == "GPU_SPEED-TRANS"):
                    F_GPU_SPEED = F_TRANS = fact
                    
                for id in range(size_of_ar):
                    SetTau(F_CPU_MM, F_TRANS, F_GPU_SPEED, F_GPU_MM)
                    n_loc = float(N_LOC)
                    ni = ni_ar[id]
                    n_let = 0
                    theta = 0.5
                    n_list = GetNlist(n_loc, theta, ni)
                    n_list_ar[id] = n_list
                    wtime_total_const_ar[id], wtime_root_const_ar[id], wtime_const_lt_const_ar[id], wtime_mom_lt_const_ar[id], wtime_const_gt_const_ar[id], \
                        wtime_mom_gt_const_ar[id], wtime_force_const_ar[id] = GetTstep(n_loc, n_list, ni, n_let, False)
                    wtime_total_reuse_ar[id], wtime_root_reuse_ar[id], wtime_const_lt_reuse_ar[id], wtime_mom_lt_reuse_ar[id], wtime_const_gt_reuse_ar[id], \
                        wtime_mom_gt_reuse_ar[id], wtime_force_reuse_ar[id] = GetTstep(n_loc, n_list, ni, n_let, True)
                    wtime_total_ave_ar[id] = (wtime_total_const_ar[id] + wtime_total_reuse_ar[id]*(n_reuse-1)) / n_reuse

                    wtime_force_total_const_ar[id], wtime_force_cp_all_const_ar[id], wtime_force_H2D_all_const_ar[id], \
                        wtime_force_cp_EPI_const_ar[id], wtime_force_H2D_EPI_const_ar[id], \
                        wtime_force_cp_list_const_ar[id], wtime_force_H2D_list_const_ar[id], \
                        wtime_force_cp_FORCE_const_ar[id], wtime_force_D2H_FORCE_const_ar[id], \
                        wtime_force_const_list_const_ar[id], wtime_force_kernel_const_ar[id], \
                        wtime_force_wb_int_cp_const_ar[id], \
                        wtime_force_kernel_trans_const_ar[id], wtime_force_kernel_calc_const_ar[id]\
                        = GetTForce(n_loc, n_let, n_list, ni, False)
                    
                    wtime_force_total_reuse_ar[id], wtime_force_cp_all_reuse_ar[id], wtime_force_H2D_all_reuse_ar[id], \
                        wtime_force_cp_EPI_reuse_ar[id], wtime_force_H2D_EPI_reuse_ar[id], \
                        wtime_force_cp_list_reuse_ar[id], wtime_force_H2D_list_reuse_ar[id], \
                        wtime_force_cp_FORCE_reuse_ar[id], wtime_force_D2H_FORCE_reuse_ar[id], \
                        wtime_force_const_list_reuse_ar[id], wtime_force_kernel_reuse_ar[id], \
                        wtime_force_wb_int_cp_reuse_ar[id], \
                        wtime_force_kernel_trans_reuse_ar[id], wtime_force_kernel_calc_reuse_ar[id]\
                        = GetTForce(n_loc, n_let, n_list, ni, True)
                    

                #if(boost_ele == "CPU_MM" and fact==1.0 and n_reuse==4.0):
                #    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='solid',   color="k", \
                #             label=r"standard model, $n_{\rm reuse}= 4$")
                #    
                #if(boost_ele == "CPU_MM" and fact==10.0 and n_reuse==4.0):
                #    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dashed',   color="r", \
                #             label=r"${\rm BW_{host\, mm}=1TB/s}$, $n_{\rm reuse}= 4$")
                #
                #if(boost_ele == "GPU_MM" and fact==10.0 and n_reuse==4.0):
                #    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dotted',   color="g", \
                #             label=r"${\rm BW_{GPU\, mm}=5TB/s}$, $n_{\rm reuse}= 4$")
                #
                #if(boost_ele == "GPU_SPEED" and fact==10.0 and n_reuse==4.0):
                #    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dashdot',   color="b", \
                #             label=r"${\rm S_{GPU}=100Tflops}$, $n_{\rm reuse}= 4$")


                """
                if(boost_ele == "CPU_MM" and fact==1.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='solid',   color="k", \
                             label=r"ST")
                    
                if(boost_ele == "GPU_SPEED" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dashed',   color="r", \
                             label=r"SG")                    
                    
                if(boost_ele == "CPU_MM" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dotted',   color="g", \
                             label=r"BH")
            
                if(boost_ele == "GPU_MM" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dashdot',   color="b", \
                             label=r"BG")
            
                if(boost_ele == "CPU_MM-GPU_SPEED" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=10.0, linestyle='solid',   color="k", \
                             label=r"SGBH")

                if(boost_ele == "CPU_MM-GPU_SPEED-GPU_MM" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=10.0, linestyle='dashed',   color="r", \
                             label=r"SGBHBG")

                if(boost_ele == "GPU_SPEED-TRANS" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=10.0, linestyle='dotted',   color="g", \
                             label=r"SGTR")
                
                    
                if(boost_ele == "CPU_MM-GPU_SPEED-TRANS" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=10.0, linestyle='dashdot',   color="b", \
                             label=r"SGBHTR")
                """

                """
                if(boost_ele == "CPU_MM" and fact==1.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=10.0, linestyle='solid',   color="k", \
                             label=r"ST")
                    
                if(boost_ele == "GPU_SPEED" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='solid',   color="k", \
                             label=r"FG")                    
            
                if(boost_ele == "CPU_MM" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dashed',   color="r", \
                             label=r"BH")
            
                if(boost_ele == "GPU_MM" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dotted',   color="g", \
                             label=r"BG")

                if(boost_ele == "TRANS" and fact==2.0 and n_reuse==16.0):
                    plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dashdot',   color="b", \
                             label=r"BT")

                
                if(boost_ele == target_ele and n_reuse == 16):
                    arg_t_min_gpu_speed_ar[id_t_min_gpu_speed] = np.argmin(wtime_total_ave_ar)
                    t_min_gpu_speed_ar[id_t_min_gpu_speed] = np.min(wtime_total_ave_ar)
                    wtime_force_const_list_const_ar
                    id_t_min_gpu_speed += 1
                """

                if(FLAG_PLOT == "COMPARE_SINGLE"):
                    if(boost_ele == "CPU_MM" and fact==1.0 and n_reuse==16.0):
                        plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='solid',   color="k", \
                                 label=r"reference")
                    
                    if(boost_ele == "GPU_SPEED-GPU_MM" and fact==2.0 and n_reuse==16.0):
                        plt.plot(ni_ar, wtime_total_ave_ar, linewidth=5.0, linestyle='solid',   color="k", \
                                 label=r"GPU2X")
            
                    if(boost_ele == "GPU_SPEED-GPU_MM" and fact==4.0 and n_reuse==16.0):
                        plt.plot(ni_ar, wtime_total_ave_ar, linewidth=5.0, linestyle='dashed',   color="r", \
                                 label=r"GPU4X")
            
                    if(boost_ele == "TRANS" and fact==4.0 and n_reuse==16.0):
                        plt.plot(ni_ar, wtime_total_ave_ar, linewidth=5.0, linestyle='dotted',   color="g", \
                                 label=r"LINK4X")

                    if(boost_ele == "TRANS" and fact==0.5 and n_reuse==16.0):
                        plt.plot(ni_ar, wtime_total_ave_ar, linewidth=5.0, linestyle='dashdot',   color="b", \
                                 label=r"LINK0.5X")

                        
                if(FLAG_PLOT == "COMPARE_N_REUSE"):
                    if(boost_ele == "CPU_MM" and fact==1.0 and n_reuse==1.0):
                        plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='solid',   color="k", \
                                 label=r"$n_{\rm reuse}=1$")
                    if(boost_ele == "CPU_MM" and fact==1.0 and n_reuse==4.0):
                        plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dashed',   color="r", \
                                label=r"$n_{\rm reuse}=4$")
                    if(boost_ele == "CPU_MM" and fact==1.0 and n_reuse==16.0):
                        plt.plot(ni_ar, wtime_total_ave_ar, linewidth=3.0, linestyle='dotted',   color="g", \
                                 label=r"$n_{\rm reuse}=16$")

    plt.grid()
    plt.rcParams["font.size"] = 16
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\langle {\rm n_{i}} \rangle$", fontsize=16)
    plt.ylabel("wall clock time [sec]", fontsize=16)
    plt.tick_params(which='major', width = 2, length = 14, labelsize=16)
    plt.tick_params(which='minor', width = 1, length = 7)
    
    plt.xlim([10, 1e4])
    #plt.ylim([3e-2, 0.3])
    #plt.ylim([1e-2, 0.6])
    #plt.ylim([1e-2, 0.2])
    if(FLAG_PLOT == "COMPARE_SINGLE"):
        plt.ylim([3e-2, 0.6])
    else:
        plt.ylim([6e-2, 0.6])
        
    #plt.legend(bbox_to_anchor=(0.5, 1.0), ncol=2, labelspacing=0.0, columnspacing=0.2, fontsize=16)
    plt.legend(bbox_to_anchor=(0.5, 1.0), ncol=1, labelspacing=0.0, columnspacing=0.2, fontsize=16)

    print(t_min_gpu_speed_ar)
    print(arg_t_min_gpu_speed_ar)
    for i in arg_t_min_gpu_speed_ar:
        #print(ni_ar[int(i)])
        #print(n_list_ar[int(i)])
        print("i, ni, n_list: ", i, ni_ar[int(i)], n_list_ar[int(i)])

    output_file_name = "compare_total_hyp.eps"
    if(FLAG_PLOT == "COMPARE_N_REUSE"):
        output_file_name = "compare_total_n_reuse.eps"
        
    if(FLAG_OUT_PUT_FILE == True):
        print(output_file_name)
        plt.savefig(output_file_name, dpi=150)
    else:
        plt.show()

    sys.exit()
#######################
# COMPARE SINGLE NODE #
#######################


###############
# 2D MAP TO 1D#
###############
if(FLAG_PLOT == "2D_MAP_TO_1D"):

    ni_ar = np.logspace(0, 20, 100, base=2)
    size_of_ar = len(ni_ar)
    n_list_ar = np.zeros(len(ni_ar))

    target_bwgpu = 0.05
    #target_bwgpu = 0.0125

    type_fig = "tstep-bh"
    #type_fig = "tstep-bt"

    """
    bwhost_sgpu_ar = np.zeros(8)
    bwhost_sgpu_ar[0] = 0.0003125
    bwhost_sgpu_ar[1] = 0.000625
    bwhost_sgpu_ar[2] = 0.00125
    bwhost_sgpu_ar[3] = 0.0025
    bwhost_sgpu_ar[4] = 0.005
    bwhost_sgpu_ar[5] = 0.01
    bwhost_sgpu_ar[6] = 0.02
    bwhost_sgpu_ar[7] = 0.04

    bwtrans_sgpu_ar = np.zeros(8)
    bwtrans_sgpu_ar[0] = 0.00003125
    bwtrans_sgpu_ar[1] = 0.0000625
    bwtrans_sgpu_ar[2] = 0.000125
    bwtrans_sgpu_ar[3] = 0.00025
    bwtrans_sgpu_ar[4] = 0.0005
    bwtrans_sgpu_ar[5] = 0.001
    bwtrans_sgpu_ar[6] = 0.002
    bwtrans_sgpu_ar[7] = 0.004
    """

    bwhost_sgpu_ar = np.zeros(7)
    bwhost_sgpu_ar[0] = 0.000625
    bwhost_sgpu_ar[1] = 0.00125
    bwhost_sgpu_ar[2] = 0.0025
    bwhost_sgpu_ar[3] = 0.005
    bwhost_sgpu_ar[4] = 0.01
    bwhost_sgpu_ar[5] = 0.02
    bwhost_sgpu_ar[6] = 0.04

    bwtrans_sgpu_ar = np.zeros(7)
    bwtrans_sgpu_ar[0] = 0.0000625
    bwtrans_sgpu_ar[1] = 0.000125
    bwtrans_sgpu_ar[2] = 0.00025
    bwtrans_sgpu_ar[3] = 0.0005
    bwtrans_sgpu_ar[4] = 0.001
    bwtrans_sgpu_ar[5] = 0.002
    bwtrans_sgpu_ar[6] = 0.004
    
    wtime_total_const_ar    = np.zeros(len(ni_ar))
    wtime_root_const_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_const_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_const_ar   = np.zeros(len(ni_ar))
    wtime_force_const_ar     = np.zeros(len(ni_ar))

    wtime_total_reuse_ar    = np.zeros(len(ni_ar))
    wtime_root_reuse_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_force_reuse_ar    = np.zeros(len(ni_ar))

    wtime_force_total_const_ar     = np.zeros(len(ni_ar))
    wtime_force_cp_all_const_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_all_const_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_EPI_const_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_EPI_const_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_list_const_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_list_const_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_FORCE_const_ar  = np.zeros(len(ni_ar))
    wtime_force_D2H_FORCE_const_ar  = np.zeros(len(ni_ar))
    wtime_force_const_list_const_ar  = np.zeros(len(ni_ar))
    wtime_force_kernel_const_ar      = np.zeros(len(ni_ar))
    wtime_force_wb_int_cp_const_ar   = np.zeros(len(ni_ar))
    wtime_force_kernel_trans_const_ar = np.zeros(len(ni_ar))
    wtime_force_kernel_calc_const_ar = np.zeros(len(ni_ar))
    
    wtime_force_total_reuse_ar     = np.zeros(len(ni_ar))
    wtime_force_cp_all_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_all_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_EPI_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_EPI_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_list_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_H2D_list_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_cp_FORCE_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_D2H_FORCE_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_const_list_reuse_ar  = np.zeros(len(ni_ar))
    wtime_force_kernel_reuse_ar      = np.zeros(len(ni_ar))
    wtime_force_wb_int_cp_reuse_ar   = np.zeros(len(ni_ar))
    wtime_force_kernel_trans_reuse_ar = np.zeros(len(ni_ar))
    wtime_force_kernel_calc_reuse_ar = np.zeros(len(ni_ar))
    
    wtime_total_ave_ar    = np.zeros(len(ni_ar))
    
    f_cpu_mm = f_trans = f_s_gpu = f_gpu_mm = 1.0

    
    t_min_bg0p1_bt0p00025_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p1_bt0p00025_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p1_bt0p00025 = 0
    
    t_min_bg0p1_bt0p0005_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p1_bt0p0005_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p1_bt0p0005 = 0
    
    t_min_bg0p1_bt0p001_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p1_bt0p001_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p1_bt0p001 = 0
    
    t_min_bg0p1_bt0p002_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p1_bt0p002_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p1_bt0p002 = 0


    t_min_bg0p05_bt0p00003125_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p05_bt0p00003125_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p05_bt0p00003125 = 0
    
    t_min_bg0p05_bt0p0000625_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p05_bt0p0000625_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p05_bt0p0000625 = 0
    
    t_min_bg0p05_bt0p000125_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p05_bt0p000125_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p05_bt0p000125 = 0
    
    t_min_bg0p05_bt0p00025_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p05_bt0p00025_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p05_bt0p00025 = 0
    
    t_min_bg0p05_bt0p0005_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p05_bt0p0005_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p05_bt0p0005 = 0
    
    t_min_bg0p05_bt0p001_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p05_bt0p001_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p05_bt0p001 = 0
    
    t_min_bg0p05_bt0p002_ar = np.zeros(len(bwhost_sgpu_ar))
    arg_t_min_bg0p05_bt0p002_ar = np.zeros(len(bwhost_sgpu_ar))
    id_t_min_bg0p05_bt0p002 = 0


    # BG=0.05
    t_min_bg0p05_bh0p0003125_ar = np.zeros(len(bwtrans_sgpu_ar))
    arg_t_min_bg0p05_bh0p0003125_ar = np.zeros(len(bwtrans_sgpu_ar))
    id_t_min_bg0p05_bh0p0003125 = 0
    
    t_min_bg0p05_bh0p000625_ar = np.zeros(len(bwtrans_sgpu_ar))
    arg_t_min_bg0p05_bh0p000625_ar = np.zeros(len(bwtrans_sgpu_ar))
    id_t_min_bg0p05_bh0p000625 = 0
    
    t_min_bg0p05_bh0p00125_ar = np.zeros(len(bwtrans_sgpu_ar))
    arg_t_min_bg0p05_bh0p00125_ar = np.zeros(len(bwtrans_sgpu_ar))
    id_t_min_bg0p05_bh0p00125 = 0
    
    t_min_bg0p05_bh0p0025_ar = np.zeros(len(bwtrans_sgpu_ar))
    arg_t_min_bg0p05_bh0p0025_ar = np.zeros(len(bwtrans_sgpu_ar))
    id_t_min_bg0p05_bh0p0025 = 0
    
    t_min_bg0p05_bh0p005_ar = np.zeros(len(bwtrans_sgpu_ar))
    arg_t_min_bg0p05_bh0p005_ar = np.zeros(len(bwtrans_sgpu_ar))
    id_t_min_bg0p05_bh0p005 = 0
    
    t_min_bg0p05_bh0p01_ar = np.zeros(len(bwtrans_sgpu_ar))
    arg_t_min_bg0p05_bh0p01_ar = np.zeros(len(bwtrans_sgpu_ar))
    id_t_min_bg0p05_bh0p01 = 0
    
    t_min_bg0p05_bh0p02_ar = np.zeros(len(bwtrans_sgpu_ar))
    arg_t_min_bg0p05_bh0p02_ar = np.zeros(len(bwtrans_sgpu_ar))
    id_t_min_bg0p05_bh0p02 = 0
    

    for n_reuse in 4,0, 16.0:
        for bwgpu_sgpu in 0.0125, 0.05, 0.1:
            for bwtrans_sgpu in bwtrans_sgpu_ar:
                for bwhost_sgpu in bwhost_sgpu_ar:
                    TAU_CPU_MM  = 1.0/(bwhost_sgpu  * 1.0e13)
                    TAU_GPU_MM  = 1.0/(bwgpu_sgpu   * 1.0e13)
                    TAU_TRANS   = 1.0/(bwtrans_sgpu * 1.0e13)
                    TAU_GPU_KERNEL = 1.0/1.0e13
                    for id in range(size_of_ar):
                        n_loc = float(N_LOC)
                        ni = ni_ar[id]
                        n_let = 0
                        theta = 0.5
                        n_list = GetNlist(n_loc, theta, ni)
                        n_list_ar[id] = n_list
                        wtime_total_const_ar[id], wtime_root_const_ar[id], wtime_const_lt_const_ar[id], wtime_mom_lt_const_ar[id], wtime_const_gt_const_ar[id], \
                        wtime_mom_gt_const_ar[id], wtime_force_const_ar[id] = GetTstep(n_loc, n_list, ni, n_let, False)
                        wtime_total_reuse_ar[id], wtime_root_reuse_ar[id], wtime_const_lt_reuse_ar[id], wtime_mom_lt_reuse_ar[id], wtime_const_gt_reuse_ar[id], \
                        wtime_mom_gt_reuse_ar[id], wtime_force_reuse_ar[id] = GetTstep(n_loc, n_list, ni, n_let, True)
                        wtime_total_ave_ar[id] = (wtime_total_const_ar[id] + wtime_total_reuse_ar[id]*(n_reuse-1)) / n_reuse
                        #wtime_total_ave_ar[id] = wtime_total_reuse_ar[id]

                        wtime_force_total_const_ar[id], wtime_force_cp_all_const_ar[id], wtime_force_H2D_all_const_ar[id], \
                            wtime_force_cp_EPI_const_ar[id], wtime_force_H2D_EPI_const_ar[id], \
                            wtime_force_cp_list_const_ar[id], wtime_force_H2D_list_const_ar[id], \
                            wtime_force_cp_FORCE_const_ar[id], wtime_force_D2H_FORCE_const_ar[id], \
                            wtime_force_const_list_const_ar[id], wtime_force_kernel_const_ar[id], \
                            wtime_force_wb_int_cp_const_ar[id], \
                            wtime_force_kernel_trans_const_ar[id], wtime_force_kernel_calc_const_ar[id]\
                            = GetTForce(n_loc, n_let, n_list, ni, False)
                    
                        wtime_force_total_reuse_ar[id], wtime_force_cp_all_reuse_ar[id], wtime_force_H2D_all_reuse_ar[id], \
                            wtime_force_cp_EPI_reuse_ar[id], wtime_force_H2D_EPI_reuse_ar[id], \
                            wtime_force_cp_list_reuse_ar[id], wtime_force_H2D_list_reuse_ar[id], \
                            wtime_force_cp_FORCE_reuse_ar[id], wtime_force_D2H_FORCE_reuse_ar[id], \
                            wtime_force_const_list_reuse_ar[id], wtime_force_kernel_reuse_ar[id], \
                            wtime_force_wb_int_cp_reuse_ar[id], \
                            wtime_force_kernel_trans_reuse_ar[id], wtime_force_kernel_calc_reuse_ar[id]\
                            = GetTForce(n_loc, n_let, n_list, ni, True)

                    """    
                    if(bwgpu_sgpu==0.1 and bwtrans_sgpu == 0.002 and n_reuse == 16):
                        arg_t_min_bg0p1_bt0p002_ar[id_t_min_bg0p1_bt0p002] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p1_bt0p002_ar[id_t_min_bg0p1_bt0p002] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p1_bt0p002 += 1

                    if(bwgpu_sgpu==0.1 and bwtrans_sgpu == 0.001 and n_reuse == 16):
                        arg_t_min_bg0p1_bt0p001_ar[id_t_min_bg0p1_bt0p001] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p1_bt0p001_ar[id_t_min_bg0p1_bt0p001] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p1_bt0p001 += 1
                        
                    if(bwgpu_sgpu==0.1 and bwtrans_sgpu == 0.0005 and n_reuse == 16):
                        arg_t_min_bg0p1_bt0p0005_ar[id_t_min_bg0p1_bt0p0005] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p1_bt0p0005_ar[id_t_min_bg0p1_bt0p0005] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p1_bt0p0005 += 1

                    if(bwgpu_sgpu==0.1 and bwtrans_sgpu == 0.00025 and n_reuse == 16):
                        arg_t_min_bg0p1_bt0p00025_ar[id_t_min_bg0p1_bt0p00025] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p1_bt0p00025_ar[id_t_min_bg0p1_bt0p00025] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p1_bt0p00025 += 1
                    """
                    
                    # BG=0.05
                    if(bwgpu_sgpu==target_bwgpu and bwtrans_sgpu == 0.002 and n_reuse == 16):
                        arg_t_min_bg0p05_bt0p002_ar[id_t_min_bg0p05_bt0p002] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bt0p002_ar[id_t_min_bg0p05_bt0p002] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bt0p002 += 1

                    if(bwgpu_sgpu==target_bwgpu and bwtrans_sgpu == 0.001 and n_reuse == 16):
                        arg_t_min_bg0p05_bt0p001_ar[id_t_min_bg0p05_bt0p001] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bt0p001_ar[id_t_min_bg0p05_bt0p001] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bt0p001 += 1
                        
                    if(bwgpu_sgpu==target_bwgpu and bwtrans_sgpu == 0.0005 and n_reuse == 16):
                        arg_t_min_bg0p05_bt0p0005_ar[id_t_min_bg0p05_bt0p0005] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bt0p0005_ar[id_t_min_bg0p05_bt0p0005] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bt0p0005 += 1

                    if(bwgpu_sgpu==target_bwgpu and bwtrans_sgpu == 0.00025 and n_reuse == 16):
                        arg_t_min_bg0p05_bt0p00025_ar[id_t_min_bg0p05_bt0p00025] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bt0p00025_ar[id_t_min_bg0p05_bt0p00025] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bt0p00025 += 1
                        
                    if(bwgpu_sgpu==target_bwgpu and bwtrans_sgpu == 0.000125 and n_reuse == 16):
                        arg_t_min_bg0p05_bt0p000125_ar[id_t_min_bg0p05_bt0p000125] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bt0p000125_ar[id_t_min_bg0p05_bt0p000125] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bt0p000125 += 1

                    if(bwgpu_sgpu==target_bwgpu and bwtrans_sgpu == 0.0000625 and n_reuse == 16):
                        arg_t_min_bg0p05_bt0p0000625_ar[id_t_min_bg0p05_bt0p0000625] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bt0p0000625_ar[id_t_min_bg0p05_bt0p0000625] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bt0p0000625 += 1

                    if(bwgpu_sgpu==target_bwgpu and bwtrans_sgpu == 0.00003125 and n_reuse == 16):
                        arg_t_min_bg0p05_bt0p00003125_ar[id_t_min_bg0p05_bt0p00003125] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bt0p00003125_ar[id_t_min_bg0p05_bt0p00003125] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bt0p00003125 += 1
                        

                    # BG=0.05
                    if(bwgpu_sgpu==target_bwgpu and bwhost_sgpu == 0.02 and n_reuse == 16):
                        arg_t_min_bg0p05_bh0p02_ar[id_t_min_bg0p05_bh0p02] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bh0p02_ar[id_t_min_bg0p05_bh0p02] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bh0p02 += 1

                    if(bwgpu_sgpu==target_bwgpu and bwhost_sgpu == 0.01 and n_reuse == 16):
                        arg_t_min_bg0p05_bh0p01_ar[id_t_min_bg0p05_bh0p01] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bh0p01_ar[id_t_min_bg0p05_bh0p01] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bh0p01 += 1
                        
                    if(bwgpu_sgpu==target_bwgpu and bwhost_sgpu == 0.005 and n_reuse == 16):
                        arg_t_min_bg0p05_bh0p005_ar[id_t_min_bg0p05_bh0p005] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bh0p005_ar[id_t_min_bg0p05_bh0p005] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bh0p005 += 1

                    if(bwgpu_sgpu==target_bwgpu and bwhost_sgpu == 0.0025 and n_reuse == 16):
                        arg_t_min_bg0p05_bh0p0025_ar[id_t_min_bg0p05_bh0p0025] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bh0p0025_ar[id_t_min_bg0p05_bh0p0025] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bh0p0025 += 1
                        
                    if(bwgpu_sgpu==target_bwgpu and bwhost_sgpu == 0.00125 and n_reuse == 16):
                        arg_t_min_bg0p05_bh0p00125_ar[id_t_min_bg0p05_bh0p00125] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bh0p00125_ar[id_t_min_bg0p05_bh0p00125] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bh0p00125 += 1

                    if(bwgpu_sgpu==target_bwgpu and bwhost_sgpu == 0.000625 and n_reuse == 16):
                        arg_t_min_bg0p05_bh0p000625_ar[id_t_min_bg0p05_bh0p000625] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bh0p000625_ar[id_t_min_bg0p05_bh0p000625] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bh0p000625 += 1

                    if(bwgpu_sgpu==target_bwgpu and bwhost_sgpu == 0.0003125 and n_reuse == 16):
                        arg_t_min_bg0p05_bh0p0003125_ar[id_t_min_bg0p05_bh0p0003125] = np.argmin(wtime_total_ave_ar)
                        t_min_bg0p05_bh0p0003125_ar[id_t_min_bg0p05_bh0p0003125] = np.min(wtime_total_ave_ar)
                        id_t_min_bg0p05_bh0p0003125 += 1

    plt.grid()
    plt.rcParams["font.size"] = 16
    plt.xscale("log")
    plt.yscale("log")

    plt.tick_params(which='major', width = 2, length = 14, labelsize=16)
    plt.tick_params(which='minor', width = 1, length = 7)

    if(target_bwgpu == 0.05):
        if(type_fig == "tstep-bh"):
            output_file_name = "wtime-bh_bg0.05.eps"
            xmin = 6e-4
            xmax = 1e-3
            y_at_xmin = 2e12
            data_inv_x_ar = np.array([xmin, xmax])
            plt.xlim([5e-4, 5e-2])
        elif(type_fig == "tstep-bt"):
            output_file_name = "wtime-bt_bg0.05.eps"
            xmin = 3e-5
            xmax = 8e-5
            y_at_xmin = 3e12
            data_inv_x_ar = np.array([xmin, xmax])
            plt.xlim([5e-5, 4e-3])
        data_inv_y_ar = np.array([y_at_xmin, (xmin*y_at_xmin)/xmax])


    elif(target_bwgpu == 0.0125):
        if(type_fig == "tstep-bh"):
            output_file_name = "wtime-bh_bg0.0125.eps"
            xmin = 3e-4
            xmax = 8e-4
            y_at_xmin = 3e12
            data_inv_x_ar = np.array([xmin, xmax])
            plt.xlim([5e-4, 5e-2])
        elif(type_fig == "tstep-bt"):
            output_file_name = "wtime-bt_bg0.0125.eps"
            xmin = 2e-5
            xmax = 7e-5
            y_at_xmin = 3e12
            data_inv_x_ar = np.array([xmin, xmax])
            plt.xlim([5e-5, 4e-3])
        data_inv_y_ar = np.array([y_at_xmin, (xmin*y_at_xmin)/xmax])

    plt.ylim([6e11, 1e13])
    
    if(type_fig == "tstep-bt"):
        plt.xlabel(r"$BW_{\rm transfer}/F_{\rm GPU}$", fontsize=16)
        plt.ylabel(r"$T_{\rm step, single} \cdot F_{\rm GPU}$", fontsize=16)
        plt.plot(bwtrans_sgpu_ar, t_min_bg0p05_bh0p01_ar*S_GPU, linewidth=3.0, linestyle='solid',   color="k", \
                 label=r"$BW_{\rm host}/F_{\rm GPU}=10^{-2}$",  marker='o', markersize=12)
        plt.plot(bwtrans_sgpu_ar, t_min_bg0p05_bh0p005_ar*S_GPU, linewidth=3.0, linestyle='dashed', color="r", \
                 label=r"$BW_{\rm host}/F_{\rm GPU}=5 \cdot 10^{-3}$", marker='^', markersize=12)
        plt.plot(bwtrans_sgpu_ar, t_min_bg0p05_bh0p0025_ar*S_GPU, linewidth=3.0, linestyle='dotted',   color="g", \
                 label=r"$BW_{\rm host}/F_{\rm GPU}=2.5 \cdot 10^{-3}$", marker='s', markersize=12)
        plt.plot(bwtrans_sgpu_ar, t_min_bg0p05_bh0p00125_ar*S_GPU, linewidth=3.0, linestyle='dashdot',   color="b", \
                 label=r"$BW_{\rm host}/F_{\rm GPU}=1.25 \cdot 10^{-3}$", marker='*', markersize=12)
        plt.plot(data_inv_x_ar, data_inv_y_ar, linewidth=3.0, linestyle='solid',   color="k")
        plt.text(data_inv_x_ar[1]*0.8, data_inv_y_ar[1]*1.2, r"$\propto \frac{F_{\rm GPU}}{BW_{\rm transfer}}$", color="k", fontsize=16)
        #plt.plot(data_inv_x_ar, data_inv_y_ar, linewidth=3.0, linestyle='dashdot',   color="b", label=r"$\propto (BW_{\rm host\, mm}/F_{\rm GPU})^{-1}$")

    if(type_fig == "tstep-bh"):
        plt.xlabel(r"$BW_{\rm host}/F_{\rm GPU}$", fontsize=16)
        plt.ylabel(r"$T_{\rm step, single} \cdot F_{\rm GPU}$", fontsize=16)
        #plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p001_ar*S_GPU, linewidth=3.0, linestyle='solid',   color="k", \
        #         label=r"$BW_{\rm transfer}/F_{\rm GPU}=10^{-3}$",  marker='o', markersize=12)
        #plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p0005_ar*S_GPU, linewidth=3.0, linestyle='dashed',   color="r", \
        #         label=r"$BW_{\rm transfer}/F_{\rm GPU}=5\cdot 10^{-4}$",  marker='^', markersize=12)
        #plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p00025_ar*S_GPU, linewidth=3.0, linestyle='dotted',   color="g", \
        #         label=r"$BW_{\rm transfer}/F_{\rm GPU}=2.5\cdot 10^{-4}$", marker='s', markersize=12)
        #plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p000125_ar*S_GPU, linewidth=3.0, linestyle='dashdot',   color="b", \
        #         label=r"$BW_{\rm transfer}/F_{\rm GPU}=1.25\cdot 10^{-4}$", marker='*', markersize=12)

        #plt.scatter(bwhost_sgpu_ar, t_min_bg0p05_bt0p002_ar*S_GPU, color="k", \
        #            label=r"$BW_{\rm transfer}/F_{\rm GPU}=2 \cdot 10^{-3}$",  marker='o', s=200)
        #plt.scatter(bwhost_sgpu_ar, t_min_bg0p05_bt0p001_ar*S_GPU, color="r", \
        #            label=r"$BW_{\rm transfer}/F_{\rm GPU}=10^{-3}$",  marker='^', s=200)
        #plt.scatter(bwhost_sgpu_ar, t_min_bg0p05_bt0p0005_ar*S_GPU, color="g", \
        #            label=r"$BW_{\rm transfer}/F_{\rm GPU}=5\cdot 10^{-4}$",  marker='s', s=200)
        #plt.scatter(bwhost_sgpu_ar, t_min_bg0p05_bt0p00025_ar*S_GPU, color="b", \
        #            label=r"$BW_{\rm transfer}/F_{\rm GPU}=2.5\cdot 10^{-4}$", marker='p', s=200)
        #plt.scatter(bwhost_sgpu_ar, t_min_bg0p05_bt0p000125_ar*S_GPU, color="m", \
        #            label=r"$BW_{\rm transfer}/F_{\rm GPU}=1.25\cdot 10^{-4}$", marker='h', s=200)
        #plt.scatter(bwhost_sgpu_ar, t_min_bg0p05_bt0p0000625_ar*S_GPU, color="c", \
        #            label=r"$BW_{\rm transfer}/F_{\rm GPU}=6.25\cdot 10^{-5}$", marker='v', s=200)
        #plt.scatter(bwhost_sgpu_ar, t_min_bg0p05_bt0p00003125_ar*S_GPU, color="y", \
        #            label=r"$BW_{\rm transfer}/F_{\rm GPU}=3.125\cdot 10^{-5}$", marker='*', s=200)

        plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p002_ar*S_GPU, color="k", \
                 label=r"$BW_{\rm transfer}/F_{\rm GPU}=2 \cdot 10^{-3}$",  marker='o', markersize=12)
        plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p001_ar*S_GPU, color="r", \
                 label=r"$BW_{\rm transfer}/F_{\rm GPU}=10^{-3}$",  marker='^', markersize=12)
        plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p0005_ar*S_GPU, color="g", \
                 label=r"$BW_{\rm transfer}/F_{\rm GPU}=5\cdot 10^{-4}$",  marker='s', markersize=12)
        plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p00025_ar*S_GPU, color="b", \
                 label=r"$BW_{\rm transfer}/F_{\rm GPU}=2.5\cdot 10^{-4}$", marker='p', markersize=12)
        plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p000125_ar*S_GPU, color="m", \
                 label=r"$BW_{\rm transfer}/F_{\rm GPU}=1.25\cdot 10^{-4}$", marker='h', markersize=12)
        plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p0000625_ar*S_GPU, color="c", \
                 label=r"$BW_{\rm transfer}/F_{\rm GPU}=6.25\cdot 10^{-5}$", marker='v', markersize=12)
        #plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p00003125_ar*S_GPU, color="y", \
        #         label=r"$BW_{\rm transfer}/F_{\rm GPU}=3.125\cdot 10^{-5}$", marker='*', markersize=12)
        
        plt.plot(data_inv_x_ar, data_inv_y_ar, linewidth=3.0, linestyle='solid',   color="k")
        plt.text(data_inv_x_ar[1]*0.8, data_inv_y_ar[1]*1.2, r"$\propto \frac{F_{\rm GPU}}{BW_{\rm host}}$", color="k", fontsize=16)

        
        #plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p001_ar*S_GPU, linewidth=3.0, linestyle='solid',   color="k", \
        #         label=r"$BW_{\rm transfer}/F_{\rm GPU}=0.001$",  marker='s', markersize=12)
        #plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p00025_ar*S_GPU, linewidth=3.0, linestyle='dashed',   color="r", \
        #         label=r"$BW_{\rm transfer}/F_{\rm GPU}=0.00025$",  marker='*', markersize=12)
        #plt.plot(bwhost_sgpu_ar, t_min_bg0p05_bt0p0000625_ar*S_GPU, linewidth=3.0, linestyle='dotted',   color="g", \
        #         label=r"$BW_{\rm transfer}/F_{\rm GPU}=0.0000625$", marker='*', markersize=12)
        #plt.plot(data_inv_x_ar, data_inv_y_ar, linewidth=3.0, linestyle='dashdot',   color="b", label=r"$\propto (BW_{\rm transfer}/F_{\rm GPU})^{-1}$")


    plt.legend(loc='upper right', ncol=1, labelspacing=0.0, columnspacing=0.2, fontsize=16)        
    if(FLAG_OUT_PUT_FILE == True):
        print(output_file_name)
        plt.savefig(output_file_name, dpi=150)
    else:
        plt.show()
    
    sys.exit()
###############
# 2D MAP TO 1D#
###############

##########
# 2D MAP #
##########
if(FLAG_PLOT == "2D_MAP"):
    print("2D_MAP")
    
    #FLAG_NORMALIZED = False
    FLAG_NORMALIZED = True
    
    NORMALIZED_FACTOR_X = S_GPU / BW_CPU_MM
    NORMALIZED_FACTOR_Y = S_GPU / BW_TRANS
    
    ni_ar = np.logspace(0, 20, 100, base=2)
    size_of_ar = len(ni_ar)
    n_list_ar = np.zeros(len(ni_ar))

    n_reuse = 16
    #bwgpu_sgpu = 0.2
    #bwgpu_sgpu = 0.1
    bwgpu_sgpu = 0.05
    #bwgpu_sgpu = 0.025
    #bwgpu_sgpu = 0.0125
    #bwgpu_sgpu = 0.00625
    #bwgpu_sgpu = 0.01
    #bwgpu_sgpu = 0.005

    bwhost_sgpu_ar = np.zeros(40)
    bwtrans_sgpu_ar = np.zeros(40)
    for i in range(len(bwhost_sgpu_ar)):
        bwhost_sgpu_ar[i] = 1.2**i * (3e-4)
        bwtrans_sgpu_ar[i] = 1.2**i * (3e-5)

    bwhost_sgpu_x_ar, bwtrans_sgpu_y_ar = np.meshgrid(bwhost_sgpu_ar, bwtrans_sgpu_ar)
    print(bwhost_sgpu_x_ar)
    print(bwtrans_sgpu_y_ar)
    tmin_z_ar = np.zeros(len(bwhost_sgpu_ar)*len(bwtrans_sgpu_ar)).reshape((len(bwhost_sgpu_ar),len(bwtrans_sgpu_ar)))
    print(tmin_z_ar)
    
    wtime_total_ave_ar    = np.zeros(len(ni_ar))
    f_cpu_mm = f_trans = f_s_gpu = f_gpu_mm = 1.0

    wtime_total_const_ar    = np.zeros(len(ni_ar))
    wtime_root_const_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_const_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_const_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_const_ar   = np.zeros(len(ni_ar))
    wtime_force_const_ar     = np.zeros(len(ni_ar))

    wtime_total_reuse_ar    = np.zeros(len(ni_ar))
    wtime_root_reuse_ar     = np.zeros(len(ni_ar))
    wtime_const_lt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_lt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_const_gt_reuse_ar = np.zeros(len(ni_ar))
    wtime_mom_gt_reuse_ar   = np.zeros(len(ni_ar))
    wtime_force_reuse_ar    = np.zeros(len(ni_ar))

    for id in range(size_of_ar):
        n_loc = float(N_LOC)
        ni = ni_ar[id]
        n_let = 0
        theta = 0.5
        n_list = GetNlist(n_loc, theta, ni)
        n_list_ar[id] = n_list
        wtime_total_const_ar[id], wtime_root_const_ar[id], wtime_const_lt_const_ar[id], wtime_mom_lt_const_ar[id], wtime_const_gt_const_ar[id], \
            wtime_mom_gt_const_ar[id], wtime_force_const_ar[id] = GetTstep(n_loc, n_list, ni, n_let, False)
        wtime_total_reuse_ar[id], wtime_root_reuse_ar[id], wtime_const_lt_reuse_ar[id], wtime_mom_lt_reuse_ar[id], wtime_const_gt_reuse_ar[id], \
            wtime_mom_gt_reuse_ar[id], wtime_force_reuse_ar[id] = GetTstep(n_loc, n_list, ni, n_let, True)
        wtime_total_ave_ar[id] = (wtime_total_const_ar[id] + wtime_total_reuse_ar[id]*(n_reuse-1)) / n_reuse
        if(id==0):
            tmin_tmp = wtime_total_ave_ar[id]
            id_tmin_tmp = id
        if(tmin_tmp > wtime_total_ave_ar[id]):
            tmin_tmp = wtime_total_ave_ar[id]
            id_tmin_tmp = id
    NORMALIZED_FACTOR_Z = 1.0 / (tmin_tmp*S_GPU)
    print("tmin_tmp*S_GPU=%f", tmin_tmp*S_GPU)
    
    for id_bwhost_sgpu in range(len(bwhost_sgpu_ar)):
        bwhost_sgpu = bwhost_sgpu_ar[id_bwhost_sgpu]    
        for id_bwtrans_sgpu in range(len(bwtrans_sgpu_ar)):
            bwtrans_sgpu = bwtrans_sgpu_ar[id_bwtrans_sgpu]
            tmin_tmp = 99999999.9
            id_tmin_tmp = -1
            TAU_CPU_MM  = 1.0/(bwhost_sgpu  * 1.0e13)
            TAU_GPU_MM  = 1.0/(bwgpu_sgpu   * 1.0e13)
            TAU_TRANS   = 1.0/(bwtrans_sgpu * 1.0e13)
            TAU_GPU_KERNEL = 1.0/1.0e13
            for id in range(size_of_ar):
                n_loc = float(N_LOC)
                ni = ni_ar[id]
                n_let = 0
                theta = 0.5
                n_list = GetNlist(n_loc, theta, ni)
                n_list_ar[id] = n_list
                wtime_total_const_ar[id], wtime_root_const_ar[id], wtime_const_lt_const_ar[id], wtime_mom_lt_const_ar[id], wtime_const_gt_const_ar[id], \
                    wtime_mom_gt_const_ar[id], wtime_force_const_ar[id] = GetTstep(n_loc, n_list, ni, n_let, False)
                wtime_total_reuse_ar[id], wtime_root_reuse_ar[id], wtime_const_lt_reuse_ar[id], wtime_mom_lt_reuse_ar[id], wtime_const_gt_reuse_ar[id], \
                    wtime_mom_gt_reuse_ar[id], wtime_force_reuse_ar[id] = GetTstep(n_loc, n_list, ni, n_let, True)
                wtime_total_ave_ar[id] = (wtime_total_const_ar[id] + wtime_total_reuse_ar[id]*(n_reuse-1)) / n_reuse
                if(id==0):
                    tmin_tmp = wtime_total_ave_ar[id]
                    id_tmin_tmp = id
                if(tmin_tmp > wtime_total_ave_ar[id]):
                    tmin_tmp = wtime_total_ave_ar[id]
                    id_tmin_tmp = id

            tmin_z_ar[id_bwhost_sgpu][id_bwtrans_sgpu] = tmin_tmp*S_GPU

    print(tmin_z_ar)
    plt.rcParams["font.size"] = 16
    plt.xscale("log")
    plt.yscale("log")
    if(FLAG_NORMALIZED == True):
        plt.xlabel(r"Relative  $B_{\rm host}$", fontsize=16)
        plt.ylabel(r"Relative  $B_{\rm transfer}$", fontsize=16)
    else:
        plt.xlabel(r"$B_{\rm host}/F_{\rm GPU}$", fontsize=16)
        plt.ylabel(r"$B_{\rm transfer}/F_{\rm GPU}$", fontsize=16)
    plt.tick_params(which='major', width = 2, length = 14, labelsize=16)
    plt.tick_params(which='minor', width = 1, length = 7)

    #plt.clim(tmin_z_ar.min(), tmin_z_ar.max())
    #plt.pcolor(bwhost_sgpu_x_ar, bwtrans_sgpu_y_ar, tmin_z_ar, norm=clr.LogNorm(vmin=tmin_z_ar.min(), vmax=tmin_z_ar.max()))
    if(FLAG_NORMALIZED == True):
        print("FLAG_NORMALIZED: ", FLAG_NORMALIZED)
        #plt.pcolor(bwhost_sgpu_x_ar*NORMALIZED_FACTOR_X, bwtrans_sgpu_y_ar*NORMALIZED_FACTOR_Y, tmin_z_ar*NORMALIZED_FACTOR_Z, norm=clr.LogNorm(vmin=4e11*NORMALIZED_FACTOR_Z, vmax=1.6e13*NORMALIZED_FACTOR_Z))
        #cbar = plt.colorbar(ticks=[0.01, 0.1, 1.0, 10.0, 100.0])
        #cbar.ax.set_yticklabels([r"$10^{-2}$", r"$2\cdot 10^{12}$", r"$4\cdot 10^{12}$", r"$8\cdot 10^{12}$", r"$1.6\cdot 10^{13}$", r"$3.2\cdot 10^{13}$"])
    else:
        plt.pcolor(bwhost_sgpu_x_ar, bwtrans_sgpu_y_ar, tmin_z_ar, norm=clr.LogNorm(vmin=4e11, vmax=1.6e13))
    
    #levels = [1e12, 2e12, 4e12, 8e12, 16e12, 32e12]
    #cbar = plt.colorbar(ticks=[1e12, 2e12, 4e12, 8e12, 16e12, 32e12])
    #cbar.ax.set_yticklabels([r"$10^{12}$", r"$2\cdot 10^{12}$", r"$4\cdot 10^{12}$", r"$8\cdot 10^{12}$", r"$1.6\cdot 10^{13}$", r"$3.2\cdot 10^{13}$"])

    if(bwgpu_sgpu == 0.2):
        plt.title(r"$B_{\rm GPU}/F_{\rm GPU}=0.2$", fontsize=16)
        levels = [5e11, 1.0/math.sqrt(2.0)*1e12, 1e12, math.sqrt(2.0)*1e12, 2e12, 2*math.sqrt(2.0)*1e12, 4e12, 4*math.sqrt(2.0)*1e12, 8e12, 8*math.sqrt(2.0)*1e12, 16e12, 32e12]
    elif(bwgpu_sgpu == 0.1):
        plt.title(r"$B_{\rm GPU}/F_{\rm GPU}=0.1$", fontsize=16)
        levels = [1.0/math.sqrt(2.0)*1e12, 1e12, math.sqrt(2.0)*1e12, 2e12, 2*math.sqrt(2.0)*1e12, 4e12, 4*math.sqrt(2.0)*1e12, 8e12, 8*math.sqrt(2.0)*1e12, 16e12, 32e12]        
    elif(bwgpu_sgpu == 0.05):
        if(FLAG_NORMALIZED == True):
            levels = [0.73, 1.0, 2.0, 4.0, 8.0]
        else:
            plt.title(r"$B_{\rm GPU}/F_{\rm GPU}=0.05$", fontsize=16)
            levels = [1.0/math.sqrt(2.0)*1e12, 1e12, math.sqrt(2.0)*1e12, 2e12, 2*math.sqrt(2.0)*1e12, 4e12, 4*math.sqrt(2.0)*1e12, 8e12, 8*math.sqrt(2.0)*1e12, 16e12, 32e12]
    elif(bwgpu_sgpu == 0.025):
        plt.title(r"$B_{\rm GPU}/F_{\rm GPU}=0.025$", fontsize=16)
        levels = [1.0/math.sqrt(2.0)*1e12, 1e12, math.sqrt(2.0)*1e12, 2e12, 2*math.sqrt(2.0)*1e12, 4e12, 4*math.sqrt(2.0)*1e12, 8e12, 8*math.sqrt(2.0)*1e12, 16e12, 32e12]
    elif(bwgpu_sgpu == 0.0125):
        plt.title(r"$B_{\rm GPU}/F_{\rm GPU}=0.0125$", fontsize=16)
        levels = [1.0/math.sqrt(2.0)*1e12, 1e12, math.sqrt(2.0)*1e12, 2e12, 2*math.sqrt(2.0)*1e12, 4e12, 4*math.sqrt(2.0)*1e12, 8e12, 8*math.sqrt(2.0)*1e12, 16e12, 32e12]
    elif(bwgpu_sgpu == 0.00625):
        plt.title(r"$B_{\rm GPU}/F_{\rm GPU}=0.00625$", fontsize=16)
        levels = [1.0/math.sqrt(2.0)*1e12, 1e12, math.sqrt(2.0)*1e12, 2e12, 2*math.sqrt(2.0)*1e12, 4e12, 4*math.sqrt(2.0)*1e12, 8e12, 8*math.sqrt(2.0)*1e12, 16e12, 32e12]
    elif(bwgpu_sgpu == 0.01):
        levels = [1e12, math.sqrt(2.0)*1e12, 2e12, 2*math.sqrt(2.0)*1e12, 4e12, 4*math.sqrt(2.0)*1e12, 8e12, 8*math.sqrt(2.0)*1e12, 16e12, 32e12]
    elif(bwgpu_sgpu == 0.005):
        levels = [math.sqrt(2.0)*1e12, 2e12, 2*math.sqrt(2.0)*1e12, 4e12, 4*math.sqrt(2.0)*1e12, 8e12, 8*math.sqrt(2.0)*1e12, 16e12, 32e12]




    if(FLAG_NORMALIZED == True):
        #cbar = plt.colorbar(ticks=[0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
        #cbar.ax.set_yticklabels([r"$0.25$", r"$0.5$", r"$1.0$", r"$2.0$", r"$4.0$", r"$8.0$", r"$16.0$"])
        #cbar.set_label(r"$T_{\rm step, single}\cdot F_{\rm GPU}$")
        cont = plt.contour(bwhost_sgpu_x_ar*NORMALIZED_FACTOR_X, bwtrans_sgpu_y_ar*NORMALIZED_FACTOR_Y, tmin_z_ar*NORMALIZED_FACTOR_Z, levels, locator=ticker.LogLocator(), colors='black', linewidth=5.0)
        cont.clabel(fontsize='16', norm=clr.LogNorm(vmin=tmin_z_ar.min(), vmax=tmin_z_ar.max()), fmt="%2.1e")
    else:
        cbar = plt.colorbar(ticks=[5e11, 1e12, 2e12, 4e12, 8e12, 16e12, 32e12])
        cbar.ax.set_yticklabels([r"$5\cdot 10^{11}$", r"$10^{12}$", r"$2\cdot 10^{12}$", r"$4\cdot 10^{12}$", r"$8\cdot 10^{12}$", r"$1.6\cdot 10^{13}$", r"$3.2\cdot 10^{13}$"])
        cbar.set_label(r"$T_{\rm step, single}\cdot F_{\rm GPU}$")
        cont = plt.contour(bwhost_sgpu_x_ar, bwtrans_sgpu_y_ar, tmin_z_ar, levels, locator=ticker.LogLocator(), colors='black', linewidth=5.0)
        cont.clabel(fontsize='12', norm=clr.LogNorm(vmin=tmin_z_ar.min(), vmax=tmin_z_ar.max()), fmt="%2.1e")
        #cont.clabel(fontsize='16', norm=clr.LogNorm(), fmt=ticker.LogFormatterMathtext())

    if(bwgpu_sgpu == 0.05):
        if(FLAG_NORMALIZED == True):
            x_tmp = np.array([1.0])
            y_tmp = np.array([1.0])
            #plt.scatter(x_tmp, y_tmp, color='k', marker='*', s=300)
        else:
            x_tmp = np.array([1e-2])
            y_tmp = np.array([1e-3])
            plt.scatter(x_tmp, y_tmp, color='k', marker='*', s=300)

        
    #x_tmp2 = [8e-4, 8e-3]
    #y_tmp2 = [2e-5, 1.5e-4]
    #c0 = (y_tmp2[1]-y_tmp2[0]) / (x_tmp2[1]-x_tmp2[0])
    #c1 = y_tmp2[1] - c0*x_tmp2[1]
    #x_tmp3 = [1e-3, 1e-2]
    #y_tmp3 = [c0*x_tmp3[0]+c1, c0*x_tmp3[1]+c1]
    #plt.plot(x_tmp3, y_tmp3)
    #plt.plot(x_tmp2, y_tmp2)

    if(FLAG_NORMALIZED == True):
        plt.xlim([0.03,30])
        plt.ylim([0.03,30])
        plt.axes().set_aspect('equal')
    else:
        plt.xlim([3e-4,0.1])
        plt.ylim([3e-5,0.01])
        
    if(bwgpu_sgpu == 0.2):
        output_file_name = "2d_map_bg0.2.eps"
    if(bwgpu_sgpu == 0.1):
        output_file_name = "2d_map_bg0.1.eps"
    if(bwgpu_sgpu == 0.05):
        output_file_name = "2d_map_bg0.05.eps"
    if(bwgpu_sgpu == 0.025):
        output_file_name = "2d_map_bg0.025.eps"
    if(bwgpu_sgpu == 0.0125):
        output_file_name = "2d_map_bg0.0125.eps"
    if(bwgpu_sgpu == 0.00625):
        output_file_name = "2d_map_bg0.00625.eps"
    elif(bwgpu_sgpu == 0.01):
        output_file_name = "2d_map_bg0.01.eps"
    elif(bwgpu_sgpu == 0.005):
        output_file_name = "2d_map_bg0.005.eps"

    if(FLAG_OUT_PUT_FILE == True):
        print(output_file_name)
        plt.savefig(output_file_name, dpi=150)
    else:
        plt.show()    
    
    sys.exit()
##########
# 2D MAP #
##########

###########################
# PARALLEL PEFORMANCE FIG #
###########################
if(FLAG_PLOT == "BREAK_DOWN_PARALLEL"):
    #PARAM = "WEAK"
    PARAM = "STRONG"
    if(PARAM == "STRONG"):
        #n_tot = float(2**30)
        n_tot = float(2**40)
    elif(PARAM == "WEAK"):
        n_loc = float(2**20)
        
    n_p_ar = np.logspace(5, 20, 100, base=2)
    time_total_ar = np.zeros(len(n_p_ar))
    time_single_ar = np.zeros(len(n_p_ar))
    time_dc_sort_ar = np.zeros(len(n_p_ar))
    time_const_let_ar = np.zeros(len(n_p_ar))
    time_allgather_let_ar = np.zeros(len(n_p_ar))
    time_p2p_let_ar = np.zeros(len(n_p_ar))
    for id in range(len(n_p_ar)):
        n_p = n_p_ar[id]
        n_smp = 30.0
        n_dc = 16.0
        n_reuse = 16.0
        theta = 0.5
        ni = 230.0
        if(PARAM == "STRONG"):
            n_loc = n_tot / n_p
        elif(PARAM == "WEAK"):
            n_tot = n_loc * n_p
            
        n_list = GetNlist(n_tot, theta, ni)
        n_let = GetNlet(n_loc, theta, n_p)
        n_p_close = GetNProcClose(theta)
        time_single_const_tmp_ar  = GetTstep(n_loc, n_list, ni, n_let, False)
        time_single_reuse_tmp_ar  = GetTstep(n_loc, n_list, ni, n_let, True)
        time_single_ar[id]        = (time_single_const_tmp_ar[0] + (n_reuse-1.0)*time_single_reuse_tmp_ar[0]) / n_reuse
        time_para_tmp_ar          = GetParallelTime(n_loc, n_p, n_p_close, n_smp, n_let, n_dc)
        time_total_ar[id]         = time_single_ar[id] + time_para_tmp_ar[0]
        time_dc_sort_ar[id]       = time_para_tmp_ar[2] / n_dc
        time_const_let_ar[id]     = time_para_tmp_ar[3] / n_reuse
        time_allgather_let_ar[id] = time_para_tmp_ar[4]
        time_p2p_let_ar[id]       = time_para_tmp_ar[5]

    plt.plot(n_p_ar, time_total_ar,         linewidth=5.0, linestyle='solid',   color="k")
    plt.plot(n_p_ar, time_single_ar,        linewidth=5.0, linestyle='dashed',  color="r")
    plt.plot(n_p_ar, time_dc_sort_ar,       linewidth=5.0, linestyle='dotted',  color="g")
    plt.plot(n_p_ar, time_const_let_ar,     linewidth=5.0, linestyle='solid',   color="b")
    plt.plot(n_p_ar, time_allgather_let_ar, linewidth=5.0, linestyle='dashed',  color="m")
    plt.plot(n_p_ar, time_p2p_let_ar,       linewidth=5.0, linestyle='dotted',  color="c")

    if(PARAM == "WEAK"):
        plt.annotate(s=r"$T_{\rm step, para}$", xy=(1000, 0.04), xytext=(1000, 0.15), color='k', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='k', edgecolor='k'))
        plt.annotate(s=r"$T_{\rm step, single}$", xy=(300, 0.03), xytext=(300, 0.007), color='r', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='r', edgecolor='r'))
        plt.annotate(s=r"$T_{\rm dc}/n_{\rm dc}$", xy=(8e4, 3e-4), xytext=(1.5e5, 1e-4), color='g', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='g', edgecolor='g'))
        plt.annotate(s=r"$T_{\rm LET, const}/n_{\rm reuse}$", xy=(7e4, 1.5e-5), xytext=(7e4, 5e-5), color='b', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='b', edgecolor='b'))

        plt.annotate(s=r"$T_{\rm LET, allgather}$", xy=(1e4, 2.0e-4), xytext=(1e4, 1e-3), color='m', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='m', edgecolor='m'))

        plt.annotate(s=r"$T_{\rm LET, p2p}$", xy=(1e3, 4.0e-3), xytext=(1e3, 1e-3), color='c', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='c', edgecolor='c'))
    elif(PARAM == "STRONG"):
        if(n_tot == float(2**40)):
            plt.annotate(s=r"$T_{\rm step, para}$", xy=(8e5, 0.1), xytext=(2e5, 1), color='k', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='k', edgecolor='k'))
            plt.annotate(s=r"$T_{\rm step, single}$", xy=(8e5, 4e-2), xytext=(3e4, 0.1), color='r', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='r', edgecolor='r'))
            #plt.annotate(s=r"$T_{\rm step, single}$", xy=(8e5, 4e-2), xytext=(8e4, 0.02), color='r', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='r', edgecolor='r'))
            plt.annotate(s=r"$T_{\rm dc}/n_{\rm dc}$", xy=(2e5, 4e-4), xytext=(3e5, 7e-5), color='g', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='g', edgecolor='g'))
            plt.annotate(s=r"$T_{\rm LET, const}/n_{\rm reuse}$", xy=(2e2, 3e-3), xytext=(2e2, 5e-2), color='b', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='b', edgecolor='b'))
            plt.annotate(s=r"$T_{\rm LET, allgather}$", xy=(2e2, 5.0e-5), xytext=(2e2, 3e-4), color='m', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='m', edgecolor='m'))
            plt.annotate(s=r"$T_{\rm LET, p2p}$", xy=(1e3, 0.2), xytext=(1e3, 2.0), color='c', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='c', edgecolor='c'))
        elif(n_tot == float(2**30)):
            plt.annotate(s=r"$T_{\rm step, para}$", xy=(3000, 2e-2), xytext=(10000, 0.1), color='k', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='k', edgecolor='k'))
            plt.annotate(s=r"$T_{\rm step, single}$", xy=(30000, 1e-3), xytext=(7e4, 3e-2), color='r', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='r', edgecolor='r'))
            plt.annotate(s=r"$T_{\rm dc}/n_{\rm dc}$", xy=(8e4, 3e-4), xytext=(1.5e5, 1e-4), color='g', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='g', edgecolor='g'))
            plt.annotate(s=r"$T_{\rm LET, const}/n_{\rm reuse}$", xy=(1.1e2, 8e-5), xytext=(2e2, 2e-4), color='b', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='b', edgecolor='b'))
            plt.annotate(s=r"$T_{\rm LET, allgather}$", xy=(1e4, 2.0e-4), xytext=(5e3, 1e-3), color='m', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='m', edgecolor='m'))
            plt.annotate(s=r"$T_{\rm LET, p2p}$", xy=(1e3, 4.0e-3), xytext=(1e3, 1e-3), color='c', xycoords='data', fontsize=16, arrowprops=dict(width=1, facecolor='c', edgecolor='c'))
            
            
    plt.grid()
    plt.rcParams["font.size"] = 16
    plt.xscale("log")
    plt.yscale("log")
    plt.tick_params(which='major', width = 2, length = 14, labelsize=16)
    plt.tick_params(which='minor', width = 1, length = 7)
    plt.xlabel(r"${\rm n_{p}}$", fontsize=16)
    plt.ylabel("wall clock time [sec]", fontsize=16)    
    if(PARAM == "STRONG"):
        plt.legend(loc="upper left")
        plt.xlim([100,1e6])
        if(n_tot == float(2**40)):
            plt.ylim([1e-5,1000])
            output_file_name = "total_parallel_break_down_strong_n1T.eps"
        elif(n_tot == float(2**30)):
            plt.ylim([1e-5,0.5])
            #plt.ylim([1e-5,10])
            output_file_name = "total_parallel_break_down_strong_n1G.eps"

        plt.legend(loc='upper right', ncol=2, labelspacing=0.0, columnspacing=0.0, fontsize=16)

    if(PARAM == "WEAK"):
        plt.xlim([100,1e6])
        #plt.ylim([1e-6,0.05])
        #plt.ylim([1e-5,0.05])
        plt.ylim([1e-5,0.3])
        output_file_name = "total_parallel_break_down_weak_n1M.eps"
        #plt.legend(bbox_to_anchor=(1.0, 0.25), ncol=2, labelspacing=0.0, columnspacing=0.0, fontsize=16)
        plt.legend(bbox_to_anchor=(0.7, 0.7), ncol=2, labelspacing=0.0, columnspacing=0.0, fontsize=16)
            
    if(FLAG_OUT_PUT_FILE == True):
        print(output_file_name)
        plt.savefig(output_file_name, dpi=150)
    else:
        plt.show()
    sys.exit()
###########################
# PARALLEL PEFORMANCE FIG #
###########################
    
    


    
    

