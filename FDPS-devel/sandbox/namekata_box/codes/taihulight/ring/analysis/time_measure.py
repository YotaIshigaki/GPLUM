import sys

file_name = sys.argv[1]
print(file_name)

file = open(file_name)

n_merge = 64.0
n_merge_cum = 0.0

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
wtime_tot_1st = 0.0
wtime_tot_cum = 0.0
wtime_tot_cum_1st = 0.0
dict_tot = {"wtime_calc_shift_z": 0.0,
            "exchange_particle": 0.0,
            "calc_force_all": 0.0,
            "wtime_copy_epi_to_psys": 0.0,
            "wtime_rotate": 0.0}
dict_tot_1st = {"wtime_calc_shift_z": 0.0,
                "exchange_particle": 0.0,
                "calc_force_all": 0.0,
                "wtime_copy_epi_to_psys": 0.0,
                "wtime_rotate": 0.0}


wtime_copy = 0.0
wtime_copy_1st = 0.0
wtime_copy_cum = 0.0
wtime_copy_cum_1st = 0.0
dict_copy = {"set_particle_local_tree": 0.0,
             "add_moment_as_sp_local_tree": 0.0,
             "set_particle_global_tree": 0.0,
             "add_moment_as_sp_global_tree": 0.0,
             "wtime_copy_epi_to_psys": 0.0,
             "wtime_rotate": 0.0}
dict_copy_1st = {"set_particle_local_tree": 0.0,
                 "add_moment_as_sp_local_tree": 0.0,
                 "set_particle_global_tree": 0.0,
                 "add_moment_as_sp_global_tree": 0.0,
                 "wtime_copy_epi_to_psys": 0.0,
                 "wtime_rotate": 0.0}


wtime_tree = 0.0
wtime_tree_1st = 0.0
wtime_tree_cum = 0.0
wtime_tree_cum_1st = 0.0
dict_tree = {"morton_sort_local_tree": 0.0,
             "link_cell_local_tree": 0.0,
             "calc_moment_local_tree": 0.0,
             "morton_sort_global_tree": 0.0,
             "link_cell_global_tree": 0.0,
             "calc_moment_global_tree": 0.0}
dict_tree_1st = {"morton_sort_local_tree": 0.0,
                     "link_cell_local_tree": 0.0,
                     "calc_moment_local_tree": 0.0,
                     "morton_sort_global_tree": 0.0,
                     "link_cell_global_tree": 0.0,
                     "calc_moment_global_tree": 0.0}


wtime_let = 0.0
wtime_let_1st = 0.0
wtime_let_cum = 0.0
wtime_let_cum_1st = 0.0
wtime_let_comm = 0.0
wtime_let_const = 0.0
wtime_let_comm_cum = 0.0
wtime_let_const_cum = 0.0
wtime_let_comm_1st = 0.0
wtime_let_const_1st = 0.0
dict_let = {"wtime_ex_let_sd_2_allgather_bd": 0.0,
            "wtime_ex_let_sd_2_ex_n_d_d": 0.0,
            "wtime_ex_let_sd_2_ex_n_d_sd_1d_ring": 0.0,
            "wtime_ex_let_sd_2_ex_n_d_sd_sub": 0.0,
            "wtime_ex_let_sd_2_allgather_sd": 0.0,
            "wtime_ex_let_sd_2_ep_sp_among_sd": 0.0,
            "wtime_ex_let_sd_2_ep_sp_in_sd": 0.0,
            "wtime_ex_let_sd_2_make_list": 0.0}
dict_let_1st = {"wtime_ex_let_sd_2_allgather_bd": 0.0,
                "wtime_ex_let_sd_2_ex_n_d_d": 0.0,
                "wtime_ex_let_sd_2_ex_n_d_sd_1d_ring": 0.0,
                "wtime_ex_let_sd_2_ex_n_d_sd_sub": 0.0,
                "wtime_ex_let_sd_2_allgather_sd": 0.0,
                "wtime_ex_let_sd_2_ep_sp_among_sd": 0.0,
                "wtime_ex_let_sd_2_ep_sp_in_sd": 0.0,
                "wtime_ex_let_sd_2_make_list": 0.0}

wtime_list = 0.0
wtime_list_1st = 0.0
wtime_list_cum = 0.0
wtime_list_cum_1st = 0.0
dict_list = {"make_all_interaction_list_id": 0.0}
dict_list_1st = {"make_all_interaction_list_id": 0.0}

wtime_force = 0.0
wtime_force_1st = 0.0
wtime_force_cum = 0.0
wtime_force_cum_1st = 0.0
dict_force = {"calc_force": 0.0}
dict_force_1st = {"calc_force": 0.0}

wtime_comm = 0.0
wtime_comm_1st = 0.0
wtime_comm_cum = 0.0
wtime_comm_cum_1st = 0.0
dict_comm = {"exchange_particle": 0.0, "exchange_LET_tot": 0.0}
dict_comm_1st = {"exchange_particle": 0.0, "exchange_LET_tot": 0.0}

n_op_tot = 0.0
n_op_tot_cum = 0.0

n_loop = 0

def ClearWtimeTot():
    #global wtime_tot = 0.0
    global wtime_tot, wtime_copy, wtime_tree, wtime_list, wtime_force, wtime_comm, \
    wtime_let_comm, wtime_let_const, \
    wtime_tot_1st, wtime_copy_1st, wtime_tree_1st, wtime_list_1st, wtime_force_1st, \
    wtime_comm_1st, wtime_let_comm_1st, wtime_let_const_1st
    wtime_tot = 0.0
    wtime_copy = 0.0
    wtime_tree = 0.0
    wtime_list = 0.0
    wtime_force = 0.0
    wtime_comm = 0.0
    wtime_let_comm = 0.0
    wtime_let_const = 0.0
    wtime_tot_1st = 0.0
    wtime_copy_1st = 0.0
    wtime_tree_1st = 0.0
    wtime_list_1st = 0.0
    wtime_force_1st = 0.0
    wtime_comm_1st = 0.0
    wtime_let_comm_1st = 0.0
    wtime_let_const_1st = 0.0    

def ClearDict(dict):
    for key in dict:
        dict[key] = 0.0

wtime_ex_ptcl_cum = wtime_const_local_tree_cum = wtime_calc_moment_local_tree_cum \
                  = wtime_const_let_cum = wtime_comm_let_cum = wtime_const_global_tree_cum \
                  = wtime_calc_moment_global_tree_cum = wtime_const_list_cum \
                  = wtime_calc_force_cum = 0.0

wtime_ex_ptcl_cum_1st = wtime_const_local_tree_cum_1st = wtime_calc_moment_local_tree_cum_1st \
                  = wtime_const_let_cum_1st = wtime_comm_let_cum_1st = wtime_const_global_tree_cum_1st \
                  = wtime_calc_moment_global_tree_cum_1st = wtime_const_list_cum_1st \
                  = wtime_calc_force_cum_1st = 0.0
        
for str_line in file:
    str_list = str_line.replace("=", "").split()
    for i in range(len(str_list)):
        for key in dict_tot:
            if(str_list[i] == key):
                if(dict_tot[key] == 0.0):
                    dict_tot_1st[key] += float(str_list[i+1])
                    wtime_tot_1st += float(str_list[i+1])
                dict_tot[key] += float(str_list[i+1])
                wtime_tot += float(str_list[i+1])

        for key in dict_copy:
            if(str_list[i] == key):
                if(dict_copy[key] == 0.0):
                    dict_copy_1st[key] += float(str_list[i+1])
                    wtime_copy_1st += float(str_list[i+1])
                dict_copy[key] += float(str_list[i+1])
                wtime_copy += float(str_list[i+1])

        for key in dict_tree:
            if(str_list[i] == key):
                if(dict_tree[key] == 0.0):
                    dict_tree_1st[key] += float(str_list[i+1])
                    wtime_tree_1st += float(str_list[i+1])
                dict_tree[key] += float(str_list[i+1])
                wtime_tree += float(str_list[i+1])

        for key in dict_list:
            if(str_list[i] == key):
                if(dict_list[key] == 0.0):
                    dict_list_1st[key] += float(str_list[i+1])
                    wtime_list_1st += float(str_list[i+1])
                dict_list[key] += float(str_list[i+1])
                wtime_list += float(str_list[i+1])

        for key in dict_force:
            if(str_list[i] == key):
                if(dict_force[key] == 0.0):
                    dict_force_1st[key] += float(str_list[i+1])
                    wtime_force_1st += float(str_list[i+1])
                dict_force[key] += float(str_list[i+1])
                wtime_force += float(str_list[i+1])

        for key in dict_comm:
            if(str_list[i] == key):
                if(dict_comm[key] == 0.0):
                    dict_comm_1st[key] += float(str_list[i+1])
                    wtime_comm_1st += float(str_list[i+1])
                dict_comm[key] += float(str_list[i+1])
                wtime_comm += float(str_list[i+1])

        for key in dict_let:
            if(str_list[i] == key):
                if(dict_let[key] == 0.0):
                    dict_let_1st[key] += float(str_list[i+1])
                    wtime_let_1st += float(str_list[i+1])
                dict_let[key] += float(str_list[i+1])
                wtime_let += float(str_list[i+1])
                

        if(str_list[i] == "n_op_tot"):
            n_op_tot += float(str_list[i+1])
                
        """
        wtime_copy += sum_wtime(dict_copy, str_list, i)
        wtime_tree += sum_wtime(dict_tree, str_list, i)
        wtime_list += sum_wtime(dict_list, str_list, i)
        wtime_force += sum_wtime(dict_force, str_list, i)
        wtime_comm += sum_wtime(dict_comm, str_list, i)
        """

        if("n_loop"==str_list[i] and n_loop == 0):
            n_op_tot = 0.0
            ClearWtimeTot()
            n_loop = 1
        elif( ("n_loop"==str_list[i] and n_loop != 0)):
            print("n_loop", n_loop)
            print("speed", n_op_tot/wtime_tot*1e-12, "[Tflops]")
            print("wtime_tot", wtime_tot/n_merge)
            print("wtime_copy", wtime_copy/n_merge)
            print("wtime_tree", wtime_tree/n_merge)
            print("wtime_list", wtime_list/n_merge)
            print("wtime_force", wtime_force/n_merge)
            print("wtime_comm", wtime_comm/n_merge)
            print("wtime_tot-wtime_foce-wtime_comm", (wtime_tot-wtime_force-wtime_comm)/n_merge)

            # averaget (per 64 steps)
            wtime_ex_ptcl = dict_comm["exchange_particle"]
            wtime_const_local_tree = dict_tree["morton_sort_local_tree"] \
                                     + dict_tree["link_cell_local_tree"]
            wtime_calc_moment_local_tree = dict_tree["calc_moment_local_tree"]
            wtime_const_let = dict_let["wtime_ex_let_sd_2_make_list"]
            wtime_comm_let = dict_let["wtime_ex_let_sd_2_allgather_bd"] \
                             + dict_let["wtime_ex_let_sd_2_ex_n_d_d"] \
                             + dict_let["wtime_ex_let_sd_2_ex_n_d_sd_1d_ring"] \
                             + dict_let["wtime_ex_let_sd_2_ex_n_d_sd_sub"] \
                             + dict_let["wtime_ex_let_sd_2_allgather_sd"] \
                             + dict_let["wtime_ex_let_sd_2_ep_sp_among_sd"] \
                             + dict_let["wtime_ex_let_sd_2_ep_sp_in_sd"]
            wtime_const_global_tree = dict_tree["morton_sort_global_tree"] \
                                     + dict_tree["link_cell_global_tree"]
            wtime_calc_moment_global_tree = dict_tree["calc_moment_global_tree"]
            wtime_const_list = dict_list["make_all_interaction_list_id"]
            wtime_calc_force = dict_force["calc_force"]
            print()
            print("wtime_ex_ptcl", wtime_ex_ptcl/n_merge)
            print("wtime_const_local_tree", wtime_const_local_tree/n_merge)
            print("wtime_calc_moment_local_tree", wtime_calc_moment_local_tree/n_merge)
            print("wtime_const_let", wtime_const_let/n_merge)
            print("wtime_comm_let",  wtime_comm_let/n_merge)
            print("wtime_const_global_tree", wtime_const_global_tree/n_merge)
            print("wtime_calc_moment_global_tree", wtime_calc_moment_global_tree/n_merge)
            print("wtime_const_list", wtime_const_list/n_merge)
            print("wtime_calc_force", wtime_calc_force/n_merge)
            print("sum", (wtime_ex_ptcl+wtime_const_local_tree+wtime_calc_moment_local_tree \
                          +wtime_const_let+wtime_comm_let+wtime_const_global_tree \
                          +wtime_calc_moment_global_tree+wtime_const_list+wtime_calc_force)/n_merge)
            print()

            #1st step only 
            wtime_ex_ptcl_1st = dict_comm_1st["exchange_particle"]
            wtime_const_local_tree_1st = dict_tree_1st["morton_sort_local_tree"] \
                                         + dict_tree_1st["link_cell_local_tree"]
            wtime_calc_moment_local_tree_1st = dict_tree_1st["calc_moment_local_tree"]
            wtime_const_let_1st = dict_let_1st["wtime_ex_let_sd_2_make_list"]
            wtime_comm_let_1st = dict_let_1st["wtime_ex_let_sd_2_allgather_bd"] \
                                 + dict_let_1st["wtime_ex_let_sd_2_ex_n_d_d"] \
                                 + dict_let_1st["wtime_ex_let_sd_2_ex_n_d_sd_1d_ring"] \
                                 + dict_let_1st["wtime_ex_let_sd_2_ex_n_d_sd_sub"] \
                                 + dict_let_1st["wtime_ex_let_sd_2_allgather_sd"] \
                                 + dict_let_1st["wtime_ex_let_sd_2_ep_sp_among_sd"] \
                                 + dict_let_1st["wtime_ex_let_sd_2_ep_sp_in_sd"]
            wtime_const_global_tree_1st = dict_tree_1st["morton_sort_global_tree"] \
                                          + dict_tree_1st["link_cell_global_tree"]
            wtime_calc_moment_global_tree_1st = dict_tree_1st["calc_moment_global_tree"]
            wtime_const_list_1st = dict_list_1st["make_all_interaction_list_id"]
            wtime_calc_force_1st = dict_force_1st["calc_force"]
            print()
            print("wtime_ex_ptcl_1st", wtime_ex_ptcl_1st)
            print("wtime_const_local_tree_1st", wtime_const_local_tree_1st)
            print("wtime_calc_moment_local_tree_1st", wtime_calc_moment_local_tree_1st)
            print("wtime_const_let_1st", wtime_const_let_1st)
            print("wtime_comm_let_1st",  wtime_comm_let_1st)
            print("wtime_const_global_tree_1st", wtime_const_global_tree_1st)
            print("wtime_calc_moment_global_tree_1st", wtime_calc_moment_global_tree_1st)
            print("wtime_const_list_1st", wtime_const_list_1st)
            print("wtime_calc_force_1st", wtime_calc_force_1st)
            print("sum", (wtime_ex_ptcl_1st+wtime_const_local_tree_1st+wtime_calc_moment_local_tree_1st \
                          +wtime_const_let_1st+wtime_comm_let_1st+wtime_const_global_tree_1st \
                          +wtime_calc_moment_global_tree_1st+wtime_const_list_1st+wtime_calc_force_1st))
            print()
            
            n_merge_cum     += n_merge
            n_op_tot_cum    += n_op_tot
            wtime_tot_cum   += wtime_tot
            wtime_copy_cum  += wtime_copy
            wtime_tree_cum  += wtime_tree
            wtime_list_cum  += wtime_list
            wtime_force_cum += wtime_force
            wtime_comm_cum  += wtime_comm
            wtime_let_comm_cum += wtime_let_comm
            wtime_let_const_cum += wtime_let_const

            # average (over all steps)
            wtime_ex_ptcl_cum += wtime_ex_ptcl
            wtime_const_local_tree_cum += wtime_const_local_tree
            wtime_calc_moment_local_tree_cum += wtime_calc_moment_local_tree
            wtime_const_let_cum += wtime_const_let
            wtime_comm_let_cum += wtime_comm_let
            wtime_const_global_tree_cum += wtime_const_global_tree
            wtime_calc_moment_global_tree_cum += wtime_calc_moment_global_tree
            wtime_const_list_cum += wtime_const_list
            wtime_calc_force_cum += wtime_calc_force
            print("wtime_ex_ptcl_cum", wtime_ex_ptcl_cum/n_merge_cum)
            print("wtime_const_local_tree_cum", wtime_const_local_tree_cum/n_merge_cum)
            print("wtime_calc_moment_local_tree_cum", wtime_calc_moment_local_tree_cum/n_merge_cum)
            print("wtime_const_let_cum", wtime_const_let_cum/n_merge_cum)
            print("wtime_comm_let_cum",  wtime_comm_let_cum/n_merge_cum)
            print("wtime_const_global_tree_cum", wtime_const_global_tree_cum/n_merge_cum)
            print("wtime_calc_moment_global_tree_cum", wtime_calc_moment_global_tree_cum/n_merge_cum)
            print("wtime_const_list_cum", wtime_const_list_cum/n_merge_cum)
            print("wtime_calc_force_cum", wtime_calc_force_cum/n_merge_cum)
            print("sum", (wtime_ex_ptcl_cum+wtime_const_local_tree_cum+wtime_calc_moment_local_tree_cum \
                          +wtime_const_let_cum+wtime_comm_let_cum+wtime_const_global_tree_cum \
                          +wtime_calc_moment_global_tree_cum+wtime_const_list_cum+wtime_calc_force_cum)/n_merge_cum)
            print()

            # 1st step only (over all steps)
            wtime_ex_ptcl_cum_1st += wtime_ex_ptcl_1st
            wtime_const_local_tree_cum_1st += wtime_const_local_tree_1st
            wtime_calc_moment_local_tree_cum_1st += wtime_calc_moment_local_tree_1st
            wtime_const_let_cum_1st += wtime_const_let_1st
            wtime_comm_let_cum_1st += wtime_comm_let_1st
            wtime_const_global_tree_cum_1st += wtime_const_global_tree_1st
            wtime_calc_moment_global_tree_cum_1st += wtime_calc_moment_global_tree_1st
            wtime_const_list_cum_1st += wtime_const_list_1st
            wtime_calc_force_cum_1st += wtime_calc_force_1st
            wtime_tot_cum_1st += wtime_tot_1st
            print("wtime_ex_ptcl_cum_1st", wtime_ex_ptcl_cum_1st/n_loop)
            print("wtime_const_local_tree_cum_1st", wtime_const_local_tree_cum_1st/n_loop)
            print("wtime_calc_moment_local_tree_cum_1st", wtime_calc_moment_local_tree_cum_1st/n_loop)
            print("wtime_const_let_cum_1st", wtime_const_let_cum_1st/n_loop)
            print("wtime_comm_let_cum_1st",  wtime_comm_let_cum_1st/n_loop)
            print("wtime_const_global_tree_cum_1st", wtime_const_global_tree_cum_1st/n_loop)
            print("wtime_calc_moment_global_tree_cum_1st", wtime_calc_moment_global_tree_cum_1st/n_loop)
            print("wtime_const_list_cum_1st", wtime_const_list_cum_1st/n_loop)
            print("wtime_calc_force_cum_1st", wtime_calc_force_cum_1st/n_loop)
            print("sum", (wtime_ex_ptcl_cum_1st+wtime_const_local_tree_cum_1st+wtime_calc_moment_local_tree_cum_1st \
                          +wtime_const_let_cum_1st+wtime_comm_let_cum_1st+wtime_const_global_tree_cum_1st \
                          +wtime_calc_moment_global_tree_cum_1st+wtime_const_list_cum_1st+wtime_calc_force_cum_1st)/n_loop)
            print("wtime_tot_cum_1st", wtime_tot_cum_1st/n_loop)
            print()
            




            """
            n_merge_cum     += n_merge
            n_op_tot_cum    += n_op_tot
            wtime_tot_cum   += wtime_tot
            wtime_copy_cum  += wtime_copy
            wtime_tree_cum  += wtime_tree
            wtime_list_cum  += wtime_list
            wtime_force_cum += wtime_force
            wtime_comm_cum  += wtime_comm
            wtime_let_comm_cum += wtime_let_comm
            wtime_let_const_cum += wtime_let_const
            
            wtime_ex_ptcl_cum += wtime_ex_ptcl
            wtime_const_local_tree_cum += wtime_const_local_tree
            wtime_calc_moment_local_tree_cum += wtime_calc_moment_local_tree
            wtime_const_let_cum += wtime_const_let
            wtime_comm_let_cum += wtime_comm_let
            wtime_const_global_tree_cum += wtime_const_global_tree
            wtime_calc_moment_global_tree_cum += wtime_calc_moment_global_tree
            wtime_const_list_cum += wtime_const_list
            wtime_calc_force_cum += wtime_calc_force
            print("wtime_ex_ptcl_cum", wtime_ex_ptcl_cum/n_merge_cum)
            print("wtime_const_local_tree_cum", wtime_const_local_tree_cum/n_merge_cum)
            print("wtime_calc_moment_local_tree_cum", wtime_calc_moment_local_tree_cum/n_merge_cum)
            print("wtime_const_let_cum", wtime_const_let_cum/n_merge_cum)
            print("wtime_comm_let_cum",  wtime_comm_let_cum/n_merge_cum)
            print("wtime_const_global_tree_cum", wtime_const_global_tree_cum/n_merge_cum)
            print("wtime_calc_moment_global_tree_cum", wtime_calc_moment_global_tree_cum/n_merge_cum)
            print("wtime_const_list_cum", wtime_const_list_cum/n_merge_cum)
            print("wtime_calc_force_cum", wtime_calc_force_cum/n_merge_cum)
            print("sum", (wtime_ex_ptcl_cum+wtime_const_local_tree_cum+wtime_calc_moment_local_tree_cum \
                          +wtime_const_let_cum+wtime_comm_let_cum+wtime_const_global_tree_cum \
                          +wtime_calc_moment_global_tree_cum+wtime_const_list_cum+wtime_calc_force_cum)/n_merge_cum)
            print()            
            """            




            
            print("speed(cum)",       n_op_tot_cum/wtime_tot_cum*1e-12, "[Tflops]")
            print("wtime_tot(cum)",   wtime_tot_cum/n_merge_cum)
            print("wtime_copy(cum)",  wtime_copy_cum/n_merge_cum)
            print("wtime_tree(cum)",  wtime_tree_cum/n_merge_cum)
            print("wtime_list(cum)",  wtime_list_cum/n_merge_cum)
            print("wtime_force(cum)", wtime_force_cum/n_merge_cum)
            print("wtime_comm(cum)",  wtime_comm_cum/n_merge_cum)
            print("wtime_let_comm(cum)",  wtime_let_comm_cum/n_merge_cum)
            print("wtime_let_const(cum)",  wtime_let_const_cum/n_merge_cum)
            print("wtime_tot-wtime_foce-wtime_comm(cum)", (wtime_tot_cum-wtime_force_cum-wtime_comm_cum)/n_merge_cum)
            print(wtime_tot_cum/n_merge_cum, wtime_force_cum/n_merge_cum, wtime_comm_cum/n_merge_cum, (wtime_tot_cum-wtime_force_cum-wtime_comm_cum)/n_merge_cum, n_op_tot_cum/wtime_tot_cum*1e-12)
            
            n_loop += 1

            ClearDict(dict_tot)
            ClearDict(dict_copy)
            ClearDict(dict_tree)
            ClearDict(dict_list)
            ClearDict(dict_force)
            ClearDict(dict_comm)
            ClearDict(dict_let)

            ClearDict(dict_tot_1st)
            ClearDict(dict_copy_1st)
            ClearDict(dict_tree_1st)
            ClearDict(dict_list_1st)
            ClearDict(dict_force_1st)
            ClearDict(dict_comm_1st)
            ClearDict(dict_let_1st)
                        
            """
            print("speed(cum, 1st)",       n_op_tot_cum_1st/wtime_tot_cum_1st*1e-12, "[Tflops]")
            print("wtime_tot(cum, 1st)",   wtime_tot_cum_1st/n_loop)
            print("wtime_copy(cum, 1st)",  wtime_copy_cum_1st/n_loop)
            print("wtime_tree(cum, 1st)",  wtime_tree_cum_1st/n_loop)
            print("wtime_list(cum, 1st)",  wtime_list_cum_1st/n_loop)
            print("wtime_force(cum, 1st)", wtime_force_cum_1st/n_loop)
            print("wtime_comm(cum, 1st)",  wtime_comm_cum_1st/n_loop)
            print("wtime_tot-wtime_foce-wtime_comm(cum, 1st)", (wtime_tot_cum_1st-wtime_force_cum_1st-wtime_comm_cum_1st)/n_loop)
            """
            """            
            for key in dict_tot:
                #print(key, dict_tot[key]/n_merge)
                dict_tot[key] = 0.0
            """
            print()
            n_op_tot = 0.0
            ClearWtimeTot()


print("n_loop", n_loop)
print("speed", n_op_tot/wtime_tot*1e-12, "[Tflops]")
print("wtime_tot", wtime_tot/n_merge)
print("wtime_copy", wtime_copy/n_merge)
print("wtime_tree", wtime_tree/n_merge)
print("wtime_list", wtime_list/n_merge)
print("wtime_force", wtime_force/n_merge)
print("wtime_comm", wtime_comm/n_merge)
print("wtime_tot-wtime_foce-wtime_comm", (wtime_tot-wtime_force-wtime_comm)/n_merge)

print("wtime_tot_1st", wtime_tot_1st)
print("wtime_copy_1st", wtime_copy_1st)
print("wtime_tree_1st", wtime_tree_1st)
print("wtime_list_1st", wtime_list_1st)
print("wtime_force_1st", wtime_force_1st)
print("wtime_comm_1st", wtime_comm_1st)
print("wtime_tot_1st-wtime_foce_1st-wtime_comm_1st", \
      (wtime_tot_1st-wtime_force_1st-wtime_comm_1st))

# average 
wtime_ex_ptcl = dict_comm["exchange_particle"]
wtime_const_local_tree = dict_tree["morton_sort_local_tree"] \
                         + dict_tree["link_cell_local_tree"]
wtime_calc_moment_local_tree = dict_tree["calc_moment_local_tree"]
wtime_const_let = dict_let["wtime_ex_let_sd_2_make_list"]
wtime_comm_let = dict_let["wtime_ex_let_sd_2_allgather_bd"] \
                 + dict_let["wtime_ex_let_sd_2_ex_n_d_d"] \
                 + dict_let["wtime_ex_let_sd_2_ex_n_d_sd_1d_ring"] \
                 + dict_let["wtime_ex_let_sd_2_ex_n_d_sd_sub"] \
                 + dict_let["wtime_ex_let_sd_2_allgather_sd"] \
                 + dict_let["wtime_ex_let_sd_2_ep_sp_among_sd"] \
                 + dict_let["wtime_ex_let_sd_2_ep_sp_in_sd"]
wtime_const_global_tree = dict_tree["morton_sort_global_tree"] \
                          + dict_tree["link_cell_global_tree"]
wtime_calc_moment_global_tree = dict_tree["calc_moment_global_tree"]
wtime_const_list = dict_list["make_all_interaction_list_id"]
wtime_calc_force = dict_force["calc_force"]
print()
print("wtime_ex_ptcl", wtime_ex_ptcl/n_merge)
print("wtime_const_local_tree", wtime_const_local_tree/n_merge)
print("wtime_calc_moment_local_tree", wtime_calc_moment_local_tree/n_merge)
print("wtime_const_let", wtime_const_let/n_merge)
print("wtime_comm_let",  wtime_comm_let/n_merge)
print("wtime_const_global_tree", wtime_const_global_tree/n_merge)
print("wtime_calc_moment_global_tree", wtime_calc_moment_global_tree/n_merge)
print("wtime_const_list", wtime_const_list/n_merge)
print("wtime_calc_force", wtime_calc_force/n_merge)
print("sum", (wtime_ex_ptcl+wtime_const_local_tree+wtime_calc_moment_local_tree \
              +wtime_const_let+wtime_comm_let+wtime_const_global_tree \
              +wtime_calc_moment_global_tree+wtime_const_list+wtime_calc_force)/n_merge)
print()

#1st step only 
wtime_ex_ptcl_1st = dict_comm_1st["exchange_particle"]
wtime_const_local_tree_1st = dict_tree_1st["morton_sort_local_tree"] \
                             + dict_tree_1st["link_cell_local_tree"]
wtime_calc_moment_local_tree_1st = dict_tree_1st["calc_moment_local_tree"]
wtime_const_let_1st = dict_let_1st["wtime_ex_let_sd_2_make_list"]
wtime_comm_let_1st = dict_let_1st["wtime_ex_let_sd_2_allgather_bd"] \
                     + dict_let_1st["wtime_ex_let_sd_2_ex_n_d_d"] \
                     + dict_let_1st["wtime_ex_let_sd_2_ex_n_d_sd_1d_ring"] \
                     + dict_let_1st["wtime_ex_let_sd_2_ex_n_d_sd_sub"] \
                     + dict_let_1st["wtime_ex_let_sd_2_allgather_sd"] \
                     + dict_let_1st["wtime_ex_let_sd_2_ep_sp_among_sd"] \
                     + dict_let_1st["wtime_ex_let_sd_2_ep_sp_in_sd"]
wtime_const_global_tree_1st = dict_tree_1st["morton_sort_global_tree"] \
                              + dict_tree_1st["link_cell_global_tree"]
wtime_calc_moment_global_tree_1st = dict_tree_1st["calc_moment_global_tree"]
wtime_const_list_1st = dict_list_1st["make_all_interaction_list_id"]
wtime_calc_force_1st = dict_force_1st["calc_force"]
print()
print("wtime_ex_ptcl_1st", wtime_ex_ptcl_1st)
print("wtime_const_local_tree_1st", wtime_const_local_tree_1st)
print("wtime_calc_moment_local_tree_1st", wtime_calc_moment_local_tree_1st)
print("wtime_const_let_1st", wtime_const_let_1st)
print("wtime_comm_let_1st",  wtime_comm_let_1st)
print("wtime_const_global_tree_1st", wtime_const_global_tree_1st)
print("wtime_calc_moment_global_tree_1st", wtime_calc_moment_global_tree_1st)
print("wtime_const_list_1st", wtime_const_list_1st)
print("wtime_calc_force_1st", wtime_calc_force_1st)
print("sum", (wtime_ex_ptcl_1st+wtime_const_local_tree_1st+wtime_calc_moment_local_tree_1st \
              +wtime_const_let_1st+wtime_comm_let_1st+wtime_const_global_tree_1st \
              +wtime_calc_moment_global_tree_1st+wtime_const_list_1st+wtime_calc_force_1st))
print()

n_merge_cum     += n_merge
n_op_tot_cum    += n_op_tot
wtime_tot_cum   += wtime_tot
wtime_copy_cum  += wtime_copy
wtime_tree_cum  += wtime_tree
wtime_list_cum  += wtime_list
wtime_force_cum += wtime_force
wtime_comm_cum  += wtime_comm
wtime_let_comm_cum += wtime_let_comm
wtime_let_const_cum += wtime_let_const


wtime_ex_ptcl_cum += wtime_ex_ptcl
wtime_const_local_tree_cum += wtime_const_local_tree
wtime_calc_moment_local_tree_cum += wtime_calc_moment_local_tree
wtime_const_let_cum += wtime_const_let
wtime_comm_let_cum += wtime_comm_let
wtime_const_global_tree_cum += wtime_const_global_tree
wtime_calc_moment_global_tree_cum += wtime_calc_moment_global_tree
wtime_const_list_cum += wtime_const_list
wtime_calc_force_cum += wtime_calc_force
print("wtime_ex_ptcl_cum", wtime_ex_ptcl_cum/n_merge_cum)
print("wtime_const_local_tree_cum", wtime_const_local_tree_cum/n_merge_cum)
print("wtime_calc_moment_local_tree_cum", wtime_calc_moment_local_tree_cum/n_merge_cum)
print("wtime_const_let_cum", wtime_const_let_cum/n_merge_cum)
print("wtime_comm_let_cum",  wtime_comm_let_cum/n_merge_cum)
print("wtime_const_global_tree_cum", wtime_const_global_tree_cum/n_merge_cum)
print("wtime_calc_moment_global_tree_cum", wtime_calc_moment_global_tree_cum/n_merge_cum)
print("wtime_const_list_cum", wtime_const_list_cum/n_merge_cum)
print("wtime_calc_force_cum", wtime_calc_force_cum/n_merge_cum)
print("sum", (wtime_ex_ptcl_cum+wtime_const_local_tree_cum+wtime_calc_moment_local_tree_cum \
              +wtime_const_let_cum+wtime_comm_let_cum+wtime_const_global_tree_cum \
              +wtime_calc_moment_global_tree_cum+wtime_const_list_cum+wtime_calc_force_cum)/n_merge_cum)
print()


# 1st step only (over all steps)
wtime_ex_ptcl_cum_1st += wtime_ex_ptcl_1st
wtime_const_local_tree_cum_1st += wtime_const_local_tree_1st
wtime_calc_moment_local_tree_cum_1st += wtime_calc_moment_local_tree_1st
wtime_const_let_cum_1st += wtime_const_let_1st
wtime_comm_let_cum_1st += wtime_comm_let_1st
wtime_const_global_tree_cum_1st += wtime_const_global_tree_1st
wtime_calc_moment_global_tree_cum_1st += wtime_calc_moment_global_tree_1st
wtime_const_list_cum_1st += wtime_const_list_1st
wtime_calc_force_cum_1st += wtime_calc_force_1st
wtime_tot_cum_1st += wtime_tot_1st
print("wtime_ex_ptcl_cum_1st", wtime_ex_ptcl_cum_1st/n_loop)
print("wtime_const_local_tree_cum_1st", wtime_const_local_tree_cum_1st/n_loop)
print("wtime_calc_moment_local_tree_cum_1st", wtime_calc_moment_local_tree_cum_1st/n_loop)
print("wtime_const_let_cum_1st", wtime_const_let_cum_1st/n_loop)
print("wtime_comm_let_cum_1st",  wtime_comm_let_cum_1st/n_loop)
print("wtime_const_global_tree_cum_1st", wtime_const_global_tree_cum_1st/n_loop)
print("wtime_calc_moment_global_tree_cum_1st", wtime_calc_moment_global_tree_cum_1st/n_loop)
print("wtime_const_list_cum_1st", wtime_const_list_cum_1st/n_loop)
print("wtime_calc_force_cum_1st", wtime_calc_force_cum_1st/n_loop)
print("sum", (wtime_ex_ptcl_cum_1st+wtime_const_local_tree_cum_1st+wtime_calc_moment_local_tree_cum_1st \
              +wtime_const_let_cum_1st+wtime_comm_let_cum_1st+wtime_const_global_tree_cum_1st \
              +wtime_calc_moment_global_tree_cum_1st+wtime_const_list_cum_1st+wtime_calc_force_cum_1st)/n_loop)
print("wtime_tot_cum_1st", wtime_tot_cum_1st/n_loop)
print()

print("speed(cum)",       n_op_tot_cum/wtime_tot_cum*1e-12, "[Tflops]")
print("wtime_tot(cum)",   wtime_tot_cum/n_merge_cum)
print("wtime_copy(cum)",  wtime_copy_cum/n_merge_cum)
print("wtime_tree(cum)",  wtime_tree_cum/n_merge_cum)
print("wtime_list(cum)",  wtime_list_cum/n_merge_cum)
print("wtime_force(cum)", wtime_force_cum/n_merge_cum)
print("wtime_comm(cum)",  wtime_comm_cum/n_merge_cum)
print("wtime_let_comm(cum)",  wtime_let_comm_cum/n_merge_cum)
print("wtime_let_const(cum)",  wtime_let_const_cum/n_merge_cum)
print("wtime_tot-wtime_foce-wtime_comm(cum)", (wtime_tot_cum-wtime_force_cum-wtime_comm_cum)/n_merge_cum)
print(wtime_tot_cum/n_merge_cum, wtime_force_cum/n_merge_cum, wtime_comm_cum/n_merge_cum, (wtime_tot_cum-wtime_force_cum-wtime_comm_cum)/n_merge_cum, n_op_tot_cum/wtime_tot_cum*1e-12)
n_loop += 1
n_op_tot = 0.0
ClearWtimeTot()
