import sys
import numpy as np
import matplotlib.pylab as plt
#import matplotlib.pyplot as plt

plt.rcParams["font.size"] = 12

file_name = sys.argv[1]
print(file_name)
file_strong = open(file_name)

file_name = sys.argv[2]
print(file_name)
file_weak = open(file_name)

plt.subplots_adjust(hspace=0.4)

def read_file(file, data):
    n_cnt = 0
    for str_line in file:
        str_list = str_line.split()
        #print(str_line)
        if(len(str_list) == 0):
            continue
        elif(str_list[0] == "#"):
            continue
        else:
            #for i in range(n_col):
            for i in range(len(str_list)):
                print(i, n_cnt, len(str_list))
                data[n_cnt][i] = str_list[i]
            print(data)
            n_cnt += 1

line_styles = [["k", "o"], ["r", "^"], ["g", "s"], ["b", "*"]]
line_labels = ["total", "interaction", "communication", "others"]

##################
###  weak scaling
n_col  = 6
n_line = 8
data = [[0 for i in range(n_col)] for j in range(n_line)]

read_file(file_weak, data)

xs = []
wtime_tot = []
wtime_int = []
wtime_comm = []
wtime_other = []

for i in range(n_line):
    xs.append( int(data[i][0]) )
    wtime_tot.append( float(data[i][1]) )
    wtime_int.append( float(data[i][2]) )
    wtime_comm.append( float(data[i][3]) )
    wtime_other.append( float(data[i][4]) )
    """
    plt.subplot(2,1,1)
    if(i==0):
        plt.plot(xs, wtime_tot,   ms = 9, color=line_styles[0][0], marker=line_styles[0][1], label=line_labels[0])
        plt.plot(xs, wtime_int,   ms = 9, color=line_styles[1][0], marker=line_styles[1][1], label=line_labels[1])
        plt.plot(xs, wtime_comm,  ms = 9, color=line_styles[2][0], marker=line_styles[2][1], label=line_labels[2])
        plt.plot(xs, wtime_other, ms = 9, color=line_styles[3][0], marker=line_styles[3][1], label=line_labels[3])
    else:
        plt.plot(xs, wtime_tot,   ms = 9, color=line_styles[0][0], marker=line_styles[0][1])
        plt.plot(xs, wtime_int,   ms = 9, color=line_styles[1][0], marker=line_styles[1][1])
        plt.plot(xs, wtime_comm,  ms = 9, color=line_styles[2][0], marker=line_styles[2][1])
        plt.plot(xs, wtime_other, ms = 9, color=line_styles[3][0], marker=line_styles[3][1])
    """
plt.subplot(2,1,1)
plt.plot(xs, wtime_tot,   ms = 9, color=line_styles[0][0], marker=line_styles[0][1], label=line_labels[0])
plt.plot(xs, wtime_int,   ms = 9, color=line_styles[1][0], marker=line_styles[1][1], label=line_labels[1])
plt.plot(xs, wtime_comm,  ms = 9, color=line_styles[2][0], marker=line_styles[2][1], label=line_labels[2])
plt.plot(xs, wtime_other, ms = 9, color=line_styles[3][0], marker=line_styles[3][1], label=line_labels[3])
plt.grid()
#plt.rcParams["font.size"] = 18
plt.legend(loc="upper left", ncol=2)
plt.xlim([100,20000])
plt.ylim([3e-3,20])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# of processes")
plt.ylabel("wall clock time [sec]")
#plt.savefig("comm_np-wtime.eps", dpi=150)

##################
### strong scaling
n_col  = 6
n_line = 4
data = [[0 for i in range(n_col)] for j in range(n_line)]

read_file(file_strong, data)

xs = []
wtime_tot = []
wtime_int = []
wtime_comm = []
wtime_other = []

for i in range(n_line):
    xs.append( int(data[i][0]) )
    wtime_tot.append( float(data[i][1]) )
    wtime_int.append( float(data[i][2]) )
    wtime_comm.append( float(data[i][3]) )
    wtime_other.append( float(data[i][4]) )
    plt.subplot(2,1,2)
    """
    if(i==0):
        plt.plot(xs, wtime_tot,   ms = 9, color=line_styles[0][0], marker=line_styles[0][1], label=line_labels[0])
        plt.plot(xs, wtime_int,   ms = 9, color=line_styles[1][0], marker=line_styles[1][1], label=line_labels[1])
        plt.plot(xs, wtime_comm,  ms = 9, color=line_styles[2][0], marker=line_styles[2][1], label=line_labels[2])
        plt.plot(xs, wtime_other, ms = 9, color=line_styles[3][0], marker=line_styles[3][1], label=line_labels[3])
    else:
        plt.plot(xs, wtime_tot,   ms = 9, color=line_styles[0][0], marker=line_styles[0][1])
        plt.plot(xs, wtime_int,   ms = 9, color=line_styles[1][0], marker=line_styles[1][1])
        plt.plot(xs, wtime_comm,  ms = 9, color=line_styles[2][0], marker=line_styles[2][1])
        plt.plot(xs, wtime_other, ms = 9, color=line_styles[3][0], marker=line_styles[3][1])
    """
plt.subplot(2,1,2)
plt.plot(xs, wtime_tot,   ms = 9, color=line_styles[0][0], marker=line_styles[0][1], label=line_labels[0])
plt.plot(xs, wtime_int,   ms = 9, color=line_styles[1][0], marker=line_styles[1][1], label=line_labels[1])
plt.plot(xs, wtime_comm,  ms = 9, color=line_styles[2][0], marker=line_styles[2][1], label=line_labels[2])
plt.plot(xs, wtime_other, ms = 9, color=line_styles[3][0], marker=line_styles[3][1], label=line_labels[3])
plt.grid()
#plt.rcParams["font.size"] = 18
plt.legend(loc="upper right", ncol=2)
plt.xlim([1000,9000])
plt.ylim([5e-3,10])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# of processes")
plt.ylabel("wall clock time [sec]")
plt.savefig("scaling.eps", dpi=150)
#plt.savefig("comm_np-wtime.eps", dpi=150)
#plt.show()
#plt.show()
