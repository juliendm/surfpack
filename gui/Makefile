TOP=/home/mricha/Surfpack
sp: spgui.cpp
	g++ -g -c -o spgui.o spgui.cpp -I/home/mricha/software/fltk-1.1.7 -I$(TOP) -I$(TOP)/src -I$(TOP)/src/surfaces -I$(TOP)/src/interpreter
	libtool --mode=link g++ -g -o sp spgui.o -L/usr/X11R6/lib -L/home/mricha/software/fltk-1.1.7/lib -L$(TOP)/src $(TOP)/src/libsurfpack.la -llapack -lblas -L/usr/lib/gcc-lib/i386-redhat-linux/3.3.2 -L/usr/lib/gcc-lib/i386-redhat-linux/3.3.2/../../.. -lfrtbegin -lg2c -lgcc_s -lfl -lfltk_gl -lfltk -lGLU -lGL -lm -lXext -lX11 
