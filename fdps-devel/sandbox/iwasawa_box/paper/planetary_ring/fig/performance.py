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
line_labels = ["speed", "interaction", "communication", "others"]

##################
###  weak scaling
n_col  = 6
n_line = 8
data = [[0 for i in range(n_col)] for j in range(n_line)]

read_file(file_weak, data)

xs = []
wtime_tot = []
speed = []

for i in range(n_line):
    xs.append( int(data[i][0]) )
    speed.append( float(data[i][5]) )
    """
    if(i==0):
        plt.plot(xs, speed, ms = 9, color=line_styles[0][0], marker=line_styles[0][1], label="measurement")
        plt.plot(xs, 0.754*xs*0.35, label="35% of theoretical peack")
    else:
        plt.plot(xs, speed, ms = 9, color=line_styles[0][0], marker=line_styles[0][1])
        plt.plot(xs, 0.754*xs*0.35)
    """
    
plt.subplot(2,1,1)
plt.plot(xs, speed, ms = 9, color=line_styles[0][0], marker=line_styles[0][1], label="measurement")
t = np.arange(0, 100000, 100)
plt.plot(t, 0.754*t*0.35, color=line_styles[0][0], linestyle="dashed", label="35% of theoretical peack")
plt.grid()
#plt.rcParams["font.size"] = 18
plt.legend(loc="upper left")
plt.xlim([100,20000])
plt.ylim([20,5000])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# of processes")
plt.ylabel("performance [Tflops]")
#plt.savefig("comm_np-wtime.eps", dpi=150)


##################
### strong scaling
n_col  = 6
n_line = 4
data = [[0 for i in range(n_col)] for j in range(n_line)]

read_file(file_strong, data)

xs = []
speed = []


for i in range(n_line):
    xs.append( int(data[i][0]) )
    speed.append( float(data[i][5]) )

    """
    if(i==0):
        plt.plot(xs, speed, ms = 9, color=line_styles[0][0], marker=line_styles[0][1])
        plt.plot(xs, 0.754*float(xs[i])*0.35, label="35% of theoretical peack")
    else:
        plt.plot(xs, speed, ms = 9, color=line_styles[0][0])
        plt.plot(xs, 0.754*float(xs[i])*0.35)
    """
plt.subplot(2,1,2)
plt.plot(xs, speed, ms = 9, color=line_styles[0][0], marker=line_styles[0][1], label="measurement")
t = np.arange(0, 100000, 100)
plt.plot(t, 0.754*t*0.35, color=line_styles[0][0], linestyle="dashed", label="35% of theoretical peack")
plt.grid()
#plt.rcParams["font.size"] = 18
plt.legend(loc="upper left")
plt.xlim([1000,9000])
plt.ylim([200,3000])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# of processes")
plt.ylabel("performance [Tflops]")
plt.savefig("performance.eps", dpi=150)
#plt.savefig("comm_np-wtime.eps", dpi=150)
#plt.show()
    

#plt.show()
