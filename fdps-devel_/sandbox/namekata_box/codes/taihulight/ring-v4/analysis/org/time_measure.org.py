file = open("../code-5.2.0_alpha_ver13/debug_00000.dat", "r")
n_merge = 8.0

def sum_wtime(i, str_list, dict):
    for key in dict:
        if(str_list[i] == key):
            #dict[key] += float(str_list[i+1])
            return float(str_list[i+1])
        else:
            return 0.0

    
wtime_tot = 0.0
dict_tot = {"wtime_calc_shift_z": 0.0,
            "exchange_particle": 0.0,
            "calc_force_all": 0.0,
            "wtime_copy_epi_to_psys": 0.0,
            "wtime_rotate": 0.0}

wtime_copy = 0.0
dict_copy = {"set_particle_local_tree": 0.0,
             "add_moment_as_sp_local_tree": 0.0,
             "set_particle_global_tree": 0.0,
             "add_moment_as_sp_global_tree": 0.0,
             "wtime_copy_epi_to_psys": 0.0,
             "wtime_rotate": 0.0}

wtime_tree = 0.0
dict_tree = {"morton_sort_local_tree": 0.0,
             "link_cell_local_tree": 0.0,
             "calc_moment_local_tree": 0.0,
             "morton_sort_global_tree": 0.0,
             "link_cell_global_tree": 0.0,
             "calc_moment_global_tree": 0.0}

wtime_list = 0.0
dict_list = {"make_all_interaction_list_id": 0.0}


wtime_force = 0.0
dict_force = {"calc_force": 0.0}

wtime_comm = 0.0
dict_comm = {"exchange_particle": 0.0, "exchange_LET_tot": 0.0}

for str_line in file:
    str_list = str_line.replace("=", "").split()
    for i in range(len(str_list)):
        wtime_tot += sum_wtime(i, str_list, dict_tot)
        """
        wtime_copy += sum_wtime(dict_copy, str_list, i)
        wtime_tree += sum_wtime(dict_tree, str_list, i)
        wtime_list += sum_wtime(dict_list, str_list, i)
        wtime_force += sum_wtime(dict_force, str_list, i)
        wtime_comm += sum_wtime(dict_comm, str_list, i)
        """

        if("n_loop"==str_list[i]):
            print("wtime_tot", wtime_tot/n_merge)
            #wtime_tot = 0.0
            
            """
            print("wtime_copy", wtime_copy/n_merge)
            wtime_copy = 0.0

            print("wtime_tree", wtime_tree/n_merge)
            wtime_tree = 0.0

            print("wtime_list", wtime_list/n_merge)
            wtime_list = 0.0

            print("wtime_force", wtime_force/n_merge)
            wtime_force = 0.0

            print("wtime_comm", wtime_comm/n_merge)
            wtime_comm = 0.0

            
            for key in dict_tot:
                #print(key, dict_tot[key]/n_merge)
                dict_tot[key] = 0.0
            print()
            """
