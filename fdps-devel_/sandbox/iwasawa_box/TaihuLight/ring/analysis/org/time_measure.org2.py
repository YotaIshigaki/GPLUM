import sys

file_name = sys.argv[1]
print(file_name)

file = open(file_name)

n_merge = 64.0
n_merge_cum = 64.0

"""
def sum_wtime(i, str_list, dict):
    for key in dict:
        if(str_list[i] == key):
            #dict[key] += float(str_list[i+1])
            return float(str_list[i+1])
        else:
            return 0.0
"""

wtime_tot = 0.0
wtime_tot_cum = 0.0
dict_tot = {"wtime_calc_shift_z": 0.0,
            "exchange_particle": 0.0,
            "calc_force_all": 0.0,
            "wtime_copy_epi_to_psys": 0.0,
            "wtime_rotate": 0.0}

wtime_copy = 0.0
wtime_copy_cum = 0.0
dict_copy = {"set_particle_local_tree": 0.0,
             "add_moment_as_sp_local_tree": 0.0,
             "set_particle_global_tree": 0.0,
             "add_moment_as_sp_global_tree": 0.0,
             "wtime_copy_epi_to_psys": 0.0,
             "wtime_rotate": 0.0}

wtime_tree = 0.0
wtime_tree_cum = 0.0
dict_tree = {"morton_sort_local_tree": 0.0,
             "link_cell_local_tree": 0.0,
             "calc_moment_local_tree": 0.0,
             "morton_sort_global_tree": 0.0,
             "link_cell_global_tree": 0.0,
             "calc_moment_global_tree": 0.0}

wtime_list = 0.0
wtime_list_cum = 0.0
dict_list = {"make_all_interaction_list_id": 0.0}


wtime_force = 0.0
wtime_force_cum = 0.0
dict_force = {"calc_force": 0.0}

wtime_comm = 0.0
wtime_comm_cum = 0.0
dict_comm = {"exchange_particle": 0.0, "exchange_LET_tot": 0.0}

n_op_tot = 0.0
n_op_tot_cum = 0.0

n_loop = 0

for str_line in file:
    str_list = str_line.replace("=", "").split()
    for i in range(len(str_list)):
        for key in dict_tot:
            if(str_list[i] == key):
                dict_tot[key] += float(str_list[i+1])
                wtime_tot += float(str_list[i+1])

        for key in dict_copy:
            if(str_list[i] == key):
                dict_copy[key] += float(str_list[i+1])
                wtime_copy += float(str_list[i+1])

        for key in dict_tree:
            if(str_list[i] == key):
                dict_tree[key] += float(str_list[i+1])
                wtime_tree += float(str_list[i+1])

        for key in dict_list:
            if(str_list[i] == key):
                dict_list[key] += float(str_list[i+1])
                wtime_list += float(str_list[i+1])

        for key in dict_force:
            if(str_list[i] == key):
                dict_force[key] += float(str_list[i+1])
                wtime_force += float(str_list[i+1])

        for key in dict_comm:
            if(str_list[i] == key):
                dict_comm[key] += float(str_list[i+1])
                wtime_comm += float(str_list[i+1])

        if(str_list[i] == "n_op_tot"):
            n_op_tot += float(str_list[i+1])
                
        """
        wtime_copy += sum_wtime(dict_copy, str_list, i)
        wtime_tree += sum_wtime(dict_tree, str_list, i)
        wtime_list += sum_wtime(dict_list, str_list, i)
        wtime_force += sum_wtime(dict_force, str_list, i)
        wtime_comm += sum_wtime(dict_comm, str_list, i)
        """

        if("n_loop"==str_list[i]):
            if(n_loop == 0):
                n_loop = 1
                continue
            
            print("n_loop", n_loop)
            
            print("speed", n_op_tot/wtime_tot*1e-12, "[Tflops]")
            print("wtime_tot", wtime_tot/n_merge)
            print("wtime_copy", wtime_copy/n_merge)
            print("wtime_tree", wtime_tree/n_merge)
            print("wtime_list", wtime_list/n_merge)
            print("wtime_force", wtime_force/n_merge)
            print("wtime_comm", wtime_comm/n_merge)
            print("wtime_tot-wtime_foce-wtime_comm", (wtime_tot-wtime_force-wtime_comm)/n_merge)

            n_merge_cum     += n_merge
            n_op_tot_cum    += n_op_tot
            wtime_tot_cum   += wtime_tot
            wtime_copy_cum  += wtime_copy
            wtime_tree_cum  += wtime_tree
            wtime_list_cum  += wtime_list
            wtime_force_cum += wtime_force
            wtime_comm_cum  += wtime_comm
            
            print("speed(cum)",       n_op_tot_cum/wtime_tot_cum*1e-12, "[Tflops]")
            print("wtime_tot(cum)",   wtime_tot_cum/n_merge_cum)
            print("wtime_copy(cum)",  wtime_copy_cum/n_merge_cum)
            print("wtime_tree(cum)",  wtime_tree_cum/n_merge_cum)
            print("wtime_list(cum)",  wtime_list_cum/n_merge_cum)
            print("wtime_force(cum)", wtime_force_cum/n_merge_cum)
            print("wtime_comm(cum)",  wtime_comm_cum/n_merge_cum)
            print("wtime_tot-wtime_foce-wtime_comm(cum)", (wtime_tot_cum-wtime_force_cum-wtime_comm_cum)/n_merge_cum)

            n_loop += 1
            
            #print("speed(cum, 1st)",       n_op_tot_cum_1st/wtime_tot_cum_1st*1e-12, "[Tflops]")
            #print("wtime_tot(cum, 1st)",   wtime_tot_cum_1st/n_loop)
            #print("wtime_copy(cum, 1st)",  wtime_copy_cum_1st/n_loop)
            #print("wtime_tree(cum, 1st)",  wtime_tree_cum_1st/n_loop)
            #print("wtime_list(cum, 1st)",  wtime_list_cum_1st/n_loop)
            #print("wtime_force(cum, 1st)", wtime_force_cum_1st/n_loop)
            #print("wtime_comm(cum, 1st)",  wtime_comm_cum_1st/n_loop)
            #print("wtime_tot-wtime_foce-wtime_comm(cum, 1st)", (wtime_tot_cum_1st-wtime_force_cum_1st-wtime_comm_cum_1st)/n_loop)



            
            n_op_tot = 0.0
            wtime_tot = 0.0
            wtime_copy = 0.0
            wtime_tree = 0.0
            wtime_list = 0.0
            wtime_force = 0.0
            wtime_comm = 0.0
            
            """            
            for key in dict_tot:
                #print(key, dict_tot[key]/n_merge)
                dict_tot[key] = 0.0
            """
            print()
