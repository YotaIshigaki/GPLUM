import sys
import numpy as np
import matplotlib.pylab as plt
#import matplotlib.pyplot as plt

file_name = sys.argv[1]
print(file_name)
file = open(file_name)

plt.rcParams["font.size"] = 12

ni = 8*11
nj = 8

data = [[0 for i in range(nj)] for j in range(ni)]
cnt = 0

for str_line in file:
    str_list = str_line.split()
    if(len(str_list) == 0):
        continue
    elif(str_list[0] == "#"):
        n_proc = int(str_list[3])
        if(n_proc < 128):
            break
    else:
        #print("cnt", cnt)
        #print(data)
        data[cnt][0] = n_proc
        data[cnt][1] = str_list[1]
        for i in range (2, 8):
            data[cnt][i] = str_list[i*2]
        cnt += 1

#print(data)

line_styles = [["k", "o"], ["r", "^"], ["g", "s"], ["b", "t"]]

plt.subplots_adjust(wspace=0.4)

sizes = [8, 64, 512]
xs = []
ys = []

cnt = 0
for size in sizes:
    for i in range(ni):
        if( int(data[i][1]) == size):
            xs.append( int(data[i][0]) )
            ys.append( float(data[i][2]) )
    print(xs)
    print(ys)
    plt.subplot(1,2,1)
    plt.plot(xs, ys, ms = 9, color=line_styles[cnt][0], marker=line_styles[cnt][1], label=(str(size)+"B"))
    cnt += 1
    xs = []
    ys = []

plt.subplot(1,2,1)
plt.grid()

plt.legend(loc="upper left")
plt.xlim([100,20000])
plt.ylim([1e-4,1])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# of processes")
plt.ylabel("wall clock time [sec]")
#plt.savefig("comm_np-wtime.eps", dpi=150)
#plt.show()

sizes = [8, 64, 512]
xs = []
ys = []

cnt = 0
for size in sizes:
    for i in range(ni):
        if( int(data[i][1]) == size):
            xs.append( int(data[i][0]) )
            ys.append( float(data[i][4]) )
    print(xs)
    print(ys)
    plt.subplot(1,2,2)
    plt.plot(xs, ys, ms = 9, color=line_styles[cnt][0], marker=line_styles[cnt][1], label=(str(size)+"B"))
    cnt += 1
    xs = []
    ys = []
    
plt.subplot(1,2,2)
plt.grid()
#plt.rcParams["font.size"] = 18
plt.legend(loc="upper left")
plt.xlim([100,20000])
plt.ylim([1e-3,300])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# of processes")
plt.ylabel("wall clock time [sec]")


plt.savefig("comm_np-wtime.eps", dpi=150)
plt.clf()

#plt.show()

###################
# wtime against size 

sizes = [256, 2048, 16384]
xs = []
ys = []

cnt = 0
for size in sizes:
    for i in range(ni):
        if( int(data[i][0]) == size):
            xs.append( int(data[i][1]) )
            ys.append( float(data[i][2]) )
    print(xs)
    print(ys)
    plt.subplot(1,2,1)
    plt.plot(xs, ys, ms = 9, color=line_styles[cnt][0], marker=line_styles[cnt][1], label=("np="+str(size)))
    cnt += 1
    xs = []
    ys = []

plt.grid()
#plt.rcParams["font.size"] = 18
plt.legend(loc="upper left")
plt.xlim([1,2000])
plt.ylim([3e-4,3])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("message size [B]")
plt.ylabel("wall clock time [sec]")
#plt.savefig("comm_msize-wtime.eps", dpi=150)


sizes = [256, 2048, 16384]
xs = []
ys = []

cnt = 0
for size in sizes:
    for i in range(ni):
        if( int(data[i][0]) == size):
            xs.append( int(data[i][1]) )
            ys.append( float(data[i][4]) )
    print(xs)
    print(ys)
    plt.subplot(1,2,2)
    plt.plot(xs, ys, ms = 9, color=line_styles[cnt][0], marker=line_styles[cnt][1], label=("np="+str(size)))
    cnt += 1
    xs = []
    ys = []

plt.grid()
plt.legend(loc="upper left")
plt.xlim([1,2000])
plt.ylim([1e-3,5000])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("message size [B]")
plt.ylabel("wall clock time [sec]")
plt.savefig("comm_msize-wtime.eps", dpi=150)

plt.show()
