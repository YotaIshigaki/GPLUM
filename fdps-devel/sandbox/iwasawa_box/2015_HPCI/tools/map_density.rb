
xmin = ymin = 0.0
xmax = ymax = 1.0
#nx = ny = 1024
nx = ny = 512
#nx = ny = 128

dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny

dens = Array.new(nx)
dens.length.times{ |i|
  dens[i] = Array.new(ny, 0)
}

def getID(x, y, dx, dy)
  idx = (x / dx).to_i
  idy = (y / dy).to_i
  return idx, idy
end

dir_name = "../data/sb256_np512_j2991936/"
Dir.foreach(dir_name) { |file|
  #if( file  =~ /snap_00000_00512*/ ) # z = 60
  if( file  =~ /snap_00006_00512*/ )  # z = 9
  #if( file  =~ /snap_00008_00512*/ )  # z = 2
  #if( file  =~ /snap_00009_00512*/ )  # z = 1
  #if( file  =~ /snap_00017_00512*/ )  # z = 0
    STDERR.printf("file: %s, dir_name: %s, dens[0][0]=%d\n", file, dir_name, dens[0][0])
    input = File.open(dir_name+file)
    line = input.gets("\n")
    while line = input.gets("\n")
      line_array = line.split
      x = line_array[2].to_f < 1.0 ? line_array[2].to_f : 0.0
      y = line_array[3].to_f < 1.0 ? line_array[3].to_f : 0.0
      idx, idy = getID(x, y, dx, dy)
      dens[idx][idy] += 1
    end
  end
}

ncnt = 0

for i in 0..nx-1 do
  for j in 0..ny-1 do
    ncnt += dens[i][j]
  end
end
STDERR.printf("ncnt=%d\n", ncnt)

factor_x = nx.to_f/(nx-1).to_f
factor_y = ny.to_f/(ny-1).to_f

STDERR.printf("nx=%d, ny=%d", nx, ny)

for i in 0..nx-1 do
  for j in 0..ny-1 do
    printf("%e %e %e \n", i*dx*factor_x, j*dx*factor_y, (dens[i][j].to_f)/(dx*dy*ncnt.to_f) )
  end
  printf("\n")
end

#dir_name = "../test_128_7/"
#dir = Dir.open(dir_name)

#Dir::foreach(dir_name) { |f|
#  input = File.open(dir_name+f)
#  while line = input.gets("\n")
#  end
#}

#input = File.open("timing_plummer2/t-tcal_76544_00000.dat")
#input = File.open("timing_plummer2/t-tcal_00128_00000.dat")

#n_snap = 64



#flag = 0
#n_ep_ep = 0
#n_ep_sp = 0
#time_tot = 0
#time_force = 0
#time_DDex = 0

#while line = input.gets("\n")
#  #  line_array = line.split(/\s+/)
#  line_array = line.split
#  if flag == 2
#    n_ep_ep += line_array[9].to_i()
#    n_ep_sp += line_array[11].to_i()
#  end
#  if flag == 11 || flag == 12
#    time_DDex += line_array[10].to_f()
#    p time_DDex
#  end
#  if flag == 24
#    time_force += line_array[10].to_f()
#    #p time_force
#  end
#  if flag == 26
#    time_tot += line_array[10].to_f()
#    p time_tot
#    flag = 0
#  end
#  if flag != 0
#    flag += 1
#  end
#  if /n_loop=29/ =~ line_array[2] || /n_loop=3[0-2]/ =~ line_array[2]
#    flag = 1
#  end
#end

#flops = (n_ep_ep * 29 + n_ep_sp * 64) / time_tot * 1e-9
#printf("%f  %f  %f  %f \n", flops, time_tot, time_DDex, time_force)
