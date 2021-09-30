
#input = File.open("timing_plummer2/t-tcal_76544_00000.dat")
input = File.open("timing_plummer2/t-tcal_00128_00000.dat")

flag = 0
n_ep_ep = 0
n_ep_sp = 0
time_tot = 0
time_force = 0
time_DDex = 0

while line = input.gets("\n")
  #  line_array = line.split(/\s+/)
  line_array = line.split
  if flag == 2
    n_ep_ep += line_array[9].to_i()
    n_ep_sp += line_array[11].to_i()
  end
  if flag == 11 || flag == 12
    time_DDex += line_array[10].to_f()
    p time_DDex
  end
  if flag == 24
    time_force += line_array[10].to_f()
    #p time_force
  end
  if flag == 26
    time_tot += line_array[10].to_f()
    p time_tot
    flag = 0
  end
  if flag != 0
    flag += 1
  end
  if /n_loop=29/ =~ line_array[2] || /n_loop=3[0-2]/ =~ line_array[2]
    flag = 1
  end
end

flops = (n_ep_ep * 29 + n_ep_sp * 64) / time_tot * 1e-9
printf("%f  %f  %f  %f \n", flops, time_tot, time_DDex, time_force)
