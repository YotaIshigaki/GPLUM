This is the working directory of FDPS on pure-C project.

The current code is still more like a proof of concept but this sample
at least can be compiled and run without relying on fortran sample code.

If you want to compile the current code, you first need to download
FDPS. Better to test if the fortran sample works or not.  Then you
place the source files in some directory, edit Makefile so that
FDPS_LOC points to reasonable place, and then try

   make fdpsc
   


