#!/bin/sh

# Sample make process for linking against surfpack_c_interface
/usr/bin/gcc -I/apps/surfpack/trunk/include -o eval_model.c.o -c eval_model.c

/usr/bin/gcc eval_model.c.o -o eval_model -rdynamic -L/apps/boost/1.49/lib -L/apps/surfpack/trunk/lib -lsurfpack_c_interface -lsurfpack -lsurfpack_fortran -lncsuopt -lconmin -lboost_serialization -llapack -lblas -lstdc++ -Wl,-rpath,/apps/boost/1.49/lib:/apps/surfpack/trunk/lib 
