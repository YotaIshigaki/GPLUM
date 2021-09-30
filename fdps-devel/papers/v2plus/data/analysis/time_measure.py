import sys
import difflib


file_name = sys.argv[1]
print(file_name, file=sys.stderr)

file = open(file_name)

n_merge = 0
n_block = 0
n_block_stop = 2
key_sentence_0 = 'BIG STEP START'
key_sentence_1 = 'n_loop'

#d = difflib.Differ()

class ExeTime:
    def __init__(self, name, level=1 ):
        self.name = name
        # break down level 0: total, the sum of times with level 1 is total
        self.level = level 
        
    def setTime(self, time):
        self.time = time

    def accumulateTime(self, time):
        self.time += time

    def clearTime(self):
        self.time = 0.0
        
exe_time_ar = []
#exe_time_ar.append(ExeTime(""))
exe_time_ar.append(ExeTime("wtime_loop", 0))
exe_time_ar.append(ExeTime("collect_sample_particle", 2))
exe_time_ar.append(ExeTime("decompose_domain", 2))
exe_time_ar.append(ExeTime("set_particle_local_tree"))
exe_time_ar.append(ExeTime("set_particle_global_tree", 2))
exe_time_ar.append(ExeTime("make_local_tree"))
exe_time_ar.append(ExeTime("make_global_tree"))
exe_time_ar.append(ExeTime("set_root_cell", 2))
exe_time_ar.append(ExeTime("calc_force"))
exe_time_ar.append(ExeTime("calc_moment_local_tree"))
exe_time_ar.append(ExeTime("calc_moment_global_tree"))
exe_time_ar.append(ExeTime("write_back"))
exe_time_ar.append(ExeTime("add_moment_as_sp_glboal", 2))
exe_time_ar.append(ExeTime("WTIME_KICK_DRIFT"))
exe_time_ar.append(ExeTime("morton_key_local_tree", 2))
exe_time_ar.append(ExeTime("morton_key_global_tree", 2))
exe_time_ar.append(ExeTime("morton_sort_local_tree", 2))
exe_time_ar.append(ExeTime("morton_sort_global_tree", 2))
exe_time_ar.append(ExeTime("link_cell_local_tree", 2))
exe_time_ar.append(ExeTime("link_cell_global_tree", 2))
exe_time_ar.append(ExeTime("calc_force__core__walk_tree", 2))

exe_time_ar.append(ExeTime("WTIME_COPY_IP", 2))
exe_time_ar.append(ExeTime("WTIME_COPY_JP", 2))
exe_time_ar.append(ExeTime("WTIME_COPY_ID", 2))
exe_time_ar.append(ExeTime("WTIME_COPY_IP_ID", 2))
exe_time_ar.append(ExeTime("WTIME_COPY_IP_JP", 2))
exe_time_ar.append(ExeTime("WTIME_COPY_FORCE", 2))

class Param:
    def __init__(self, name, level=1 ):
        self.name = name
        # break down level 0: total, the sum of times with level 1 is total
        self.level = level 
        
    def setVal(self, val):
        self.val = val

    def clearVal(self):
        self.val = 0.0
        
    def dump(self):
        print(self.name, self.val, file=sys.stderr)

val_ar = []
val_ar.append(Param("ni_ave", 0))

for str_line in file:
    str_line = str_line.replace("\n", "")
    if(str_line == key_sentence_0 and n_block == n_block_stop):
        break
    if(str_line == key_sentence_0):
        n_block += 1
        for inst in exe_time_ar:
            inst.clearTime()
        n_merge = 0
        continue
    str_list = str_line.replace("=", "").split()
    if(len(str_list) <= 0):
        continue
    if(str_list[0] == key_sentence_1):
        n_merge += 1
    for id in range(len(str_list)):
        for inst in exe_time_ar:
            if(str_list[id] == inst.name):
                inst.accumulateTime( float(str_list[id+1]) )
        for inst in val_ar:
            if(str_list[id] == inst.name):
                inst.setVal( float(str_list[id+1]) )

time_sum = 0.0
time_loop = 0.0
time_array = []
time_const_lt_tot = 0.0
time_const_gt_tot = 0.0

print("n_merge= ", n_merge, file=sys.stderr)
for inst in val_ar:
    inst.dump()
        
for inst in exe_time_ar:
    #print(inst.name)
    #print(inst.time / n_merge)
    if(inst.level == 1):
        #print(inst.name)
        #print(inst.time / n_merge)
        print(inst.name, (inst.time / n_merge), file=sys.stderr)
        time_sum += inst.time
    if(inst.level == 2):
        print("     ", inst.name, (inst.time / n_merge), file=sys.stderr)
    if(inst.name == "wtime_loop"):
        time_loop += inst.time
    if(inst.name == "make_local_tree"  or inst.name == "calc_moment_local_tree" ):
        time_const_lt_tot += inst.time
    if(inst.name == "make_global_tree" or inst.name == "calc_moment_global_tree" ):
        time_const_gt_tot += inst.time

print("time_const_lt_tot/n_merge= ", time_const_lt_tot/n_merge, file=sys.stderr)
print("time_const_gt_tot/n_merge= ", time_const_gt_tot/n_merge, file=sys.stderr)
print("time_sum/n_merge= ", time_sum/n_merge, file=sys.stderr)
print("time_loop/n_merge= ", time_loop/n_merge, file=sys.stderr)
print("(time_loop - time_sum)/n_merge= ", (time_loop - time_sum)/n_merge, file=sys.stderr)



key_ward_list = []
key_ward_list.append("set_particle_local_tree")
key_ward_list.append("make_local_tree")
key_ward_list.append("calc_moment_local_tree")
key_ward_list.append("make_global_tree")
key_ward_list.append("calc_moment_global_tree")
key_ward_list.append("calc_force")
key_ward_list.append("write_back")
key_ward_list.append("WTIME_KICK_DRIFT")

print(time_loop/n_merge, end='  ', file=sys.stderr)
time_sum = 0.0
for key in key_ward_list:
    for inst in exe_time_ar:
        if(inst.name == key):
            time_sum += inst.time
            print(inst.time/n_merge, end='  ', file=sys.stderr)
print("time_misc=", (time_loop-time_sum)/n_merge, file=sys.stderr)
print("\n", file=sys.stderr)

n_grp = val_ar[0].val

print(n_merge, n_grp, time_loop/n_merge)
