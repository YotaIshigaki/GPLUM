rm log.txt
mpirun -q --bind-to none -n  1 -host t1n121 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n  2 -host t1n121 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n  3 -host t1n121 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n  4 -host t1n121 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n  5 -host t1n121,t1n123 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n  6 -host t1n121,t1n123 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n  7 -host t1n121,t1n123 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n  8 -host t1n121,t1n123 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n  9 -host t1n121,t1n123,t1n124 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n 10 -host t1n121,t1n123,t1n124 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n 11 -host t1n121,t1n123,t1n124 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n 12 -host t1n121,t1n123,t1n124 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n 13 -host t1n121,t1n111,t1n123,t1n124 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n 14 -host t1n121,t1n111,t1n123,t1n124 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n 15 -host t1n121,t1n111,t1n123,t1n124 ./SPzyH.out >> log.txt
mpirun -q --bind-to none -n 16 -host t1n121,t1n111,t1n123,t1n124 ./SPzyH.out >> log.txt
