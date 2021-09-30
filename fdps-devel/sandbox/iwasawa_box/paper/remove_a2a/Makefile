SRC_FDPS_DIR  = /home/masaki/project/fdps.git/src
SRC_PHANTOM_DIR = /home/masaki/project/git/phantom-grape
SRC_PM_DIR = /home/masaki/project/fdps.git/sandbox/iwasawa_box/paper/remove_a2a/particle_mesh

#DEST = hoge
#DEST = xc
#DEST = hokusai
DEST = fugaku_test

ifeq ($(findstring xc, $(DEST)), xc)
USER = iwaswams
DEST_SYS = $(DEST)
DEST_DIR = /home/iwaswams/work
DEST_FDPS_DIR = /home/iwaswams/work/fdps/src
DEST_PHANTOM_DIR = /home/iwaswams/work/remove_a2a
DEST_PM_DIR = /home/iwaswams/work/remove_a2a
COMPILE = "module switch PrgEnv-cray PrgEnv-gnu; module switch gcc/7.3.0 gcc/8.3.0; cd $(DEST_DIR)/remove_a2a/nbody/; make -B"
else ifeq ($(findstring hokusai, $(DEST)), hokusai)
USER = iwasawa
DEST_SYS = $(DEST)
DEST_DIR = /home/iwasawa/work/
DEST_FDPS_DIR = $(DEST_DIR)/fdps/src
DEST_PHANTOM_DIR = $(DEST_DIR)/remove_a2a
DEST_PM_DIR = $(DEST_DIR)/remove_a2a
COMPILE = "module load intel; cd $(DEST_DIR)/remove_a2a/nbody/; make -B"
else ifeq ($(findstring fugaku_test, $(DEST)), fugaku_test)
USER = go0013
DEST_SYS = $(DEST)
DEST_DIR = /home/users/go19/go0013/work/
DEST_FDPS_DIR = $(DEST_DIR)/fdps/src
DEST_PHANTOM_DIR = $(DEST_DIR)/remove_a2a
DEST_PM_DIR = $(DEST_DIR)/remove_a2a
COMPILE = "module load lang; cd $(DEST_DIR)/remove_a2a/nbody/; make -B"
else ifeq ($(findstring jenever0, $(DEST)), jenever0)
USER = masaki
DEST_SYS = jenever0
DEST_DIR = /nfshome/masaki/work
DEST_FDPS_DIR = /nfshome/masaki/work/fdps/src
#COMPILE = "cd $(DEST_DIR)/remove_a2a/treepm/; make -B"
else ifeq ($(findstring Tsukumo00, $(DEST)), Tsukumo00)
USER = masaki
DEST_SYS = 172.19.8.68
DEST_DIR = /home/masaki/work
DEST_FDPS_DIR = /home/masaki/work/fdps/src
else
$(error "specify DEST_SYS")
endif

send_src:
	ssh $(USER)@$(DEST_SYS) "cd ./work; mkdir -pv fdps; cd fdps; mkdir -pv src;"
	rsync -av --include='*.hpp' --include='*.cpp' --exclude=$(SRC_FDPS_DIR)/* $(SRC_FDPS_DIR)/ $(USER)@$(DEST_SYS):$(DEST_FDPS_DIR)
	rsync -av $(SRC_PHANTOM_DIR) $(USER)@$(DEST_SYS):$(DEST_PHANTOM_DIR)
	rsync -av $(SRC_PM_DIR) $(USER)@$(DEST_SYS):$(DEST_PM_DIR)
	rsync -av --include='*.hpp' --include='*.cpp' --include='Makefile' --include='input.para' --include='IC' \
	--exclude='.git' --exclude='result/*' --exclude='utils/*' --exclude='*~' --exclude='*.eps' ../remove_a2a $(USER)@$(DEST_SYS):$(DEST_DIR)

compile_src:
	ssh $(USER)@$(DEST_SYS) $(COMPILE)

recv_result:
	rsync -av $(USER)@$(DEST_SYS):$(DEST_DIR)/remove_a2a/nbody/result_n1M_np64_m1_M2_l16_n256 ./nbody/

#ssh $(USER)@$(DEST_SYS) "$(COMPILE); cd $(DEST_DIR)/remove_a2a/treepm/; make -B"
#compile:
#	ssh $(XC) "module switch PrgEnv-cray PrgEnv-gnu; modlue load fftw3; cd $(DEST_DIR)/remove_a2a/treepm/; make -B"

#send_src_jenever0:
#	ssh $(JENEVER0) "cd ./work; mkdir -pv fdps; cd fdps; mkdir -pv src"
#	rsync -av --include='*.hpp' --include='*.cpp' --exclude=$(SRC_FDPS_DIR)/* $(SRC_FDPS_DIR)/ $(JENEVER0):/nfshome/masaki/work/fdps/src
#	rsync -av --include='*.hpp' --include='*.cpp' --include='Makefile' --include='IC' \
#	--exclude='.git' --exclude='result/*' --exclude='utils/*' ../remove_a2a $(JENEVER0):/nfshome/masaki/work
#
#compile_jenever0:
#	ssh $(JENEVER0) "cd /nfshome/masaki/work/remove_a2a/treepm/; make -B"
#
#send_src_tsukumo00:
#	ssh $(JENEVER0) "cd ./work; mkdir -pv fdps; cd fdps; mkdir -pv src"
#	rsync -av --include='*.hpp' --include='*.cpp' --exclude=$(SRC_FDPS_DIR)/* $(SRC_FDPS_DIR)/ $(JENEVER0):/nfshome/masaki/work/fdps/src
#	rsync -av --include='*.hpp' --include='*.cpp' --include='Makefile' --include='IC' \
#	--exclude='.git' --exclude='result/*' --exclude='utils/*' ../remove_a2a $(JENEVER0):/nfshome/masaki/work
#
#compile_tsukumo00:
#	ssh $(JENEVER0) "cd /nfshome/masaki/work/remove_a2a/treepm/; make -B"

#send_src_xc:
#	rsync -av --include='*.hpp' --include='*.cpp' --exclude=$(SRC_FDPS_DIR)/* $(SRC_FDPS_DIR)/ $(DEST_SYS):$(DEST_FDPS_DIR)
#	rsync -av --include='*.hpp' --include='*.cpp' --include='Makefile' --include='IC' \
#	--exclude='.git' --exclude='result/*' --exclude='utils/*' ../remove_a2a $(DEST_SYS):$(DEST_DIR)

#	ssh $(XC) "cd ./work; mkdir -pv fdps; cd fdps; mkdir -pv src" \

#compile_xc:
#	ssh $(XC) "module switch PrgEnv-cray PrgEnv-gnu; modlue load fftw3; cd $(DEST_DIR)/remove_a2a/treepm/; make -B"

#	ssh $(DEST_SYS) "cd $(DEST_DIR);"
