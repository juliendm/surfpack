%module firstswig

%{
#include "firstswig.h"
%}

%include "firstswig.h"
%include "carrays.i"
%array_class(double, doubleArray);
