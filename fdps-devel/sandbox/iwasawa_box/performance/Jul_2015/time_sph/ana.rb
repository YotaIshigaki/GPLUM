
#input = File.open("timing_plummer2/t-tcal_76544_00000.dat")
input = File.open("time_sph_nx644.dat")

n_cnt_wtime = 0
flag = 0
n_ep_ep = 0
n_ep_sp = 0
time_DDex = 0
time_tot = 0
time_grav = 0
time_hydr = 0
n_op_dens = 0
n_op_hydr = 0
n_op_grav_ep = 0
n_op_grav_sp = 0

while line = input.gets("\n")
  line_array = line.split
  if flag == 2
    #p line_array
    time_tot += line_array[3].to_f()
    #p time_tot
  end
  if flag == 6 || flag == 7
    #p line_array
    time_DDex += line_array[7].to_f()
    #p time_DDex
  end
  if flag == 9 || flag == 12
    #p line_array
    time_hydr += line_array[7].to_f()
    #p time_hydr
  end
  if flag == 13
    #p line_array
    time_grav += line_array[7].to_f()
    #p time_grav
  end
  if flag == 22
    #p line_array
    line_array2 = line_array[0].split("=")
    n_op_dens += line_array2[1].to_i()
    #p n_op_dens
    #p line_array2
  end
  if flag == 37
    #p line_array
    line_array2 = line_array[0].split("=")
    n_op_hydr += line_array2[1].to_i()
    #p n_op_hydr
    #p line_array2
  end
  if flag == 52
    #p line_array
    line_array2 = line_array[0].split("=")
    line_array3 = line_array[1].split("=")
    #p line_array2[1].chomp(",")
    #p line_array2[3]
    n_op_grav_ep += line_array2[1].chomp(",").to_i()
    n_op_grav_sp += line_array3[1].to_i()
    #p line_array2
    #p n_op_grav_ep
    #p n_op_grav_sp
  end
  if flag != 0
    flag += 1
  end
  if line_array[0] =~ /wtime_tot/
    n_cnt_wtime += 1
    flag = 0
    if n_cnt_wtime >= 97
      flag = 1
    end
    #p n_cnt_wtime
  end
end

flops = (n_op_grav_ep*29 + n_op_grav_sp*29 + n_op_dens*(42.0*3.0+74.0) + n_op_hydr*119.0) / time_tot * 1e-9
printf("%f  %f  %f  %f %f \n", flops, time_tot, time_DDex, time_grav, time_hydr)
