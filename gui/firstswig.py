# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _firstswig
reload(_firstswig)

def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class FirstClass(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FirstClass, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FirstClass, name)
    def __repr__(self):
        return "<C FirstClass instance at %s>" % (self.this,)
    def __init__(self, *args):
        _swig_setattr(self, FirstClass, 'this', _firstswig.new_FirstClass(*args))
        _swig_setattr(self, FirstClass, 'thisown', 1)
    def __del__(self, destroy=_firstswig.delete_FirstClass):
        try:
            if self.thisown: destroy(self)
        except: pass
    def printVal(*args): return _firstswig.FirstClass_printVal(*args)
    def shiftArray(*args): return _firstswig.FirstClass_shiftArray(*args)
    def myevaluate(*args): return _firstswig.FirstClass_myevaluate(*args)

class FirstClassPtr(FirstClass):
    def __init__(self, this):
        _swig_setattr(self, FirstClass, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FirstClass, 'thisown', 0)
        _swig_setattr(self, FirstClass,self.__class__,FirstClass)
_firstswig.FirstClass_swigregister(FirstClassPtr)

class doubleArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, doubleArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, doubleArray, name)
    def __repr__(self):
        return "<C doubleArray instance at %s>" % (self.this,)
    def __init__(self, *args):
        _swig_setattr(self, doubleArray, 'this', _firstswig.new_doubleArray(*args))
        _swig_setattr(self, doubleArray, 'thisown', 1)
    def __del__(self, destroy=_firstswig.delete_doubleArray):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __getitem__(*args): return _firstswig.doubleArray___getitem__(*args)
    def __setitem__(*args): return _firstswig.doubleArray___setitem__(*args)
    def cast(*args): return _firstswig.doubleArray_cast(*args)
    __swig_getmethods__["frompointer"] = lambda x: _firstswig.doubleArray_frompointer
    if _newclass:frompointer = staticmethod(_firstswig.doubleArray_frompointer)

class doubleArrayPtr(doubleArray):
    def __init__(self, this):
        _swig_setattr(self, doubleArray, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, doubleArray, 'thisown', 0)
        _swig_setattr(self, doubleArray,self.__class__,doubleArray)
_firstswig.doubleArray_swigregister(doubleArrayPtr)

doubleArray_frompointer = _firstswig.doubleArray_frompointer


