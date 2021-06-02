all:
	# Update PDFs
	cd specs/c++/ja && $(MAKE)
	cd specs/c++/en && $(MAKE)
	cd specs/fortran/ja && $(MAKE)
	cd specs/fortran/en && $(MAKE)
	cd tutorial/c++/ja && $(MAKE)
	cd tutorial/c++/en && $(MAKE)
	cd tutorial/fortran/ja && $(MAKE)
	cd tutorial/fortran/en && $(MAKE)
	cd tutorial/c/ja && $(MAKE)
	cd tutorial/c/en && $(MAKE)
	# Delete symbolic links
	rm -f doc_specs.pdf
	rm -f doc_specs_en.pdf
	rm -f doc_tutorial.pdf
	rm -f doc_tutorial_e.pdf
	rm -f doc_specs_c_ja.pdf
	rm -f doc_specs_c_en.pdf
	# Remake symbolic links
	ln -s doc_specs_cpp_ja.pdf doc_specs.pdf
	ln -s doc_specs_cpp_en.pdf doc_specs_en.pdf
	ln -s doc_tutorial_cpp_ja.pdf doc_tutorial.pdf
	ln -s doc_tutorial_cpp_en.pdf doc_tutorial_e.pdf
	ln -s doc_specs_ftn_ja.pdf doc_specs_c_ja.pdf
	ln -s doc_specs_ftn_en.pdf doc_specs_c_en.pdf
