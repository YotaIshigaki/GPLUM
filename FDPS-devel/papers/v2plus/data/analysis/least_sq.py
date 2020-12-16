import sys
from scipy.optimize import curve_fit
import numpy as np

def linear_fit(x, a):
    return a*x
    
#from sklearn.linear_model import LinearRegression

#file_name = sys.argv[1]
#print(file_name)
#file = open(file_name)

key_ward_ar = [ "set_particle_local_tree",
                "make_local_tree",
                "calc_moment_local_tree",
                "make_global_tree",
                "calc_moment_global_tree",
                "calc_force",
                "write_back",
                "WTIME_KICK_DRIFT",
                "morton_key_local_tree",
                "morton_sort_local_tree",
                "link_cell_local_tree",
                "morton_key_global_tree",
                "morton_sort_global_tree",
                "link_cell_global_tree",
                "WTIME_KERNEL",
                "calc_force__core__walk_tree",
                "WTIME_COPY_JP",
                "WTIME_SEND_JP",
                "WTIME_COPY_IP",
                "WTIME_SEND_IP",
                "WTIME_COPY_ID",
                "WTIME_SEND_ID",
                "WTIME_COPY_FORCE",
                "WTIME_RECV_FORCE"
]

                #"WTIME_COMM_IP",
                #"WTIME_COPY_IP",
                #"WTIME_COMM_JP",
                #"WTIME_COPY_JP",
                #"WTIME_COMM_ID",
                #"WTIME_COPY_ID"

for id in range(len(key_ward_ar)):
    key = key_ward_ar[id]
    print(id, key)

key_id = int(input("input number: "))
key_word = key_ward_ar[key_id]
print(key_word)

#dir_name = "../data_overlap_stream_01"
dir_name = "../data_for_tau_00"

file_name_ar = [
    dir_name+"/gpu_multi_index_reuse_N512k_l16_n1024.dat",
    dir_name+"/gpu_multi_index_reuse_N1m_l16_n1024.dat",
    dir_name+"/gpu_multi_index_reuse_N2m_l16_n1024.dat",
    dir_name+"/gpu_multi_index_reuse_N4m_l16_n1024.dat",
    dir_name+"/gpu_multi_index_reuse_N8m_l16_n1024.dat"
]


print(file_name_ar)

"""
file_name_ar = [
    "../data_overlap_01/gpu_multi_index_reuse_N512k_l16_n1024.dat",
    "../data_overlap_01/gpu_multi_index_reuse_N1m_l16_n1024.dat",
    "../data_overlap_01/gpu_multi_index_reuse_N2m_l16_n1024.dat",
    "../data_overlap_01/gpu_multi_index_reuse_N4m_l16_n1024.dat",
    "../data_overlap_01/gpu_multi_index_reuse_N8m_l16_n1024.dat"
    ]
"""

"""
file_name_ar = [
    "../data_no_overlap/gpu_multi_index_reuse_N512k_l16_n1024.dat",
    "../data_no_overlap/gpu_multi_index_reuse_N1m_l16_n1024.dat",
    "../data_no_overlap/gpu_multi_index_reuse_N2m_l16_n1024.dat",
    "../data_no_overlap/gpu_multi_index_reuse_N4m_l16_n1024.dat",
    "../data_no_overlap/gpu_multi_index_reuse_N8m_l16_n1024.dat"
    ]
"""
ni_ave = np.array([221, 419, 605, 234, 450])
nj_ave = np.array([3227, 4704, 6552, 3649, 5204])
n = np.array([2**19, 2**20, 2**21, 2**22, 2**23])

y8 = np.zeros(5)
y9 = np.zeros(5)
n_loop = 0

for id_file in range(len(file_name_ar)):
    file_name = file_name_ar[id_file]
    print(file_name)
    file = open(file_name)
    for str_line in file:
        str_line = str_line.replace("\n", "")
        str_list = str_line.replace("=", "").split()
        if(len(str_list) <= 0):
            continue
        if(str_list[0] == "n_loop"):
            n_loop = int(str_list[1])
        for id in range(len(str_list)):
            if( str_list[id] == key_word and (n_loop==8 or n_loop==9) ):
                print("n_loop= %d, %s = %f \n", n_loop, str_list[id], str_list[id+1])
                if(n_loop==8):
                    y8[id_file] = float(str_list[id+1])
                if(n_loop==9):
                    y9[id_file] = float(str_list[id+1])

x = n
if(key_word == "WTIME_KERNEL"):
    x = n*nj_ave
elif(key_word == "calc_force__core__walk_tree" or key_word == "WTIME_SEND_ID" or key_word == "WTIME_COPY_ID"):
    x = n*nj_ave/ni_ave
else:
    print("x=", x)

param8, cov8 = curve_fit(linear_fit, x, y8)
print(param8, cov8)
param9, cov9 = curve_fit(linear_fit, x, y9)
print(param9, cov9)

#print(param)
#print(optimize.leastsq(func1, a, args=(x,y8)))
#print(optimize.leastsq(func1, a, args=(x,y9)))

#lr = LinearRegression(fit_intercept=False)
#lr.fit(x.reshape(-1, 1), y)

#print(lr.get_params())
#print(lr.coef_)
#print(lr.intercept_)
#print(lr.score(x.reshape(-1, 1), y))

#print(y8)
#fit8 = np.polyfit(x, y8, 1)
#print(fit8)
#print(y9)
#fit9 = np.polyfit(x, y9, 1)
#print(fit9)
