TRT=ms
LTX=latex
#LTX=pdflatex
DVI=dvipdfmx
VIEW=evince
FIGS0 = src_reuse.eps src_dispatch.eps src_dispatch_index.eps src_retrieve.eps
FIGS1 = wtime_ngrp_th0.5.eps theta_wtime.eps compare.eps compare_total_n_reuse.eps compare_total_hyp.eps 2d_map_bg0.05.eps total_parallel_break_down_weak_n1M.eps total_parallel_break_down_strong_n1G.eps total_parallel_break_down_strong_n1T.eps wtime_sim_force.eps wtime_allgather_startup.eps wtime_allgather.eps

#VIEW=acroread
#VIEW=/cygdrive/c/Program\ Files\ \(x86\)/Adobe/Acrobat\ Reader\ DC/Reader/AcroRd32.exe

all:
	make tex
	make dvi
	#make copy
	#make view

tex:
	$(LTX) $(TRT).tex
	$(LTX) $(TRT).tex

dvi:
	$(DVI) $(TRT).dvi

copy:
	cp $(TRT).pdf ~/windows_dir/

view:
	$(VIEW) $(TRT).pdf

clean:
	rm -f *~

comp:
	mkdir -p src
	mkdir -p src/fig/
	mkdir -p src/fig/snipet/
	cp $(TRT).tex src/$(TRT).tex
	cd fig/snipet/; cp $(FIGS0) ../../src/fig/snipet/
	cd fig/; cp $(FIGS1) ../src/fig/
	tar czvf src.tar.gz src/

distclean:
	make clean
	rm -f $(TRT).aux $(TRT).dvi $(TRT).log $(TRT).toc

