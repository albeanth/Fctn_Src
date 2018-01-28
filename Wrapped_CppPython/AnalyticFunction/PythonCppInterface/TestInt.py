# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.10
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_TestInt', [dirname(__file__)])
        except ImportError:
            import _TestInt
            return _TestInt
        if fp is not None:
            try:
                _mod = imp.load_module('_TestInt', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _TestInt = swig_import_helper()
    del swig_import_helper
else:
    import _TestInt
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _TestInt.delete_SwigPyIterator
    __del__ = lambda self : None;
    def value(self): return _TestInt.SwigPyIterator_value(self)
    def incr(self, n=1): return _TestInt.SwigPyIterator_incr(self, n)
    def decr(self, n=1): return _TestInt.SwigPyIterator_decr(self, n)
    def distance(self, *args): return _TestInt.SwigPyIterator_distance(self, *args)
    def equal(self, *args): return _TestInt.SwigPyIterator_equal(self, *args)
    def copy(self): return _TestInt.SwigPyIterator_copy(self)
    def next(self): return _TestInt.SwigPyIterator_next(self)
    def __next__(self): return _TestInt.SwigPyIterator___next__(self)
    def previous(self): return _TestInt.SwigPyIterator_previous(self)
    def advance(self, *args): return _TestInt.SwigPyIterator_advance(self, *args)
    def __eq__(self, *args): return _TestInt.SwigPyIterator___eq__(self, *args)
    def __ne__(self, *args): return _TestInt.SwigPyIterator___ne__(self, *args)
    def __iadd__(self, *args): return _TestInt.SwigPyIterator___iadd__(self, *args)
    def __isub__(self, *args): return _TestInt.SwigPyIterator___isub__(self, *args)
    def __add__(self, *args): return _TestInt.SwigPyIterator___add__(self, *args)
    def __sub__(self, *args): return _TestInt.SwigPyIterator___sub__(self, *args)
    def __iter__(self): return self
SwigPyIterator_swigregister = _TestInt.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class LineInt(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, LineInt, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, LineInt, name)
    __repr__ = _swig_repr
    def iterator(self): return _TestInt.LineInt_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _TestInt.LineInt___nonzero__(self)
    def __bool__(self): return _TestInt.LineInt___bool__(self)
    def __len__(self): return _TestInt.LineInt___len__(self)
    def pop(self): return _TestInt.LineInt_pop(self)
    def __getslice__(self, *args): return _TestInt.LineInt___getslice__(self, *args)
    def __setslice__(self, *args): return _TestInt.LineInt___setslice__(self, *args)
    def __delslice__(self, *args): return _TestInt.LineInt___delslice__(self, *args)
    def __delitem__(self, *args): return _TestInt.LineInt___delitem__(self, *args)
    def __getitem__(self, *args): return _TestInt.LineInt___getitem__(self, *args)
    def __setitem__(self, *args): return _TestInt.LineInt___setitem__(self, *args)
    def append(self, *args): return _TestInt.LineInt_append(self, *args)
    def empty(self): return _TestInt.LineInt_empty(self)
    def size(self): return _TestInt.LineInt_size(self)
    def clear(self): return _TestInt.LineInt_clear(self)
    def swap(self, *args): return _TestInt.LineInt_swap(self, *args)
    def get_allocator(self): return _TestInt.LineInt_get_allocator(self)
    def begin(self): return _TestInt.LineInt_begin(self)
    def end(self): return _TestInt.LineInt_end(self)
    def rbegin(self): return _TestInt.LineInt_rbegin(self)
    def rend(self): return _TestInt.LineInt_rend(self)
    def pop_back(self): return _TestInt.LineInt_pop_back(self)
    def erase(self, *args): return _TestInt.LineInt_erase(self, *args)
    def __init__(self, *args): 
        this = _TestInt.new_LineInt(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _TestInt.LineInt_push_back(self, *args)
    def front(self): return _TestInt.LineInt_front(self)
    def back(self): return _TestInt.LineInt_back(self)
    def assign(self, *args): return _TestInt.LineInt_assign(self, *args)
    def resize(self, *args): return _TestInt.LineInt_resize(self, *args)
    def insert(self, *args): return _TestInt.LineInt_insert(self, *args)
    def reserve(self, *args): return _TestInt.LineInt_reserve(self, *args)
    def capacity(self): return _TestInt.LineInt_capacity(self)
    __swig_destroy__ = _TestInt.delete_LineInt
    __del__ = lambda self : None;
LineInt_swigregister = _TestInt.LineInt_swigregister
LineInt_swigregister(LineInt)

class LineDouble(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, LineDouble, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, LineDouble, name)
    __repr__ = _swig_repr
    def iterator(self): return _TestInt.LineDouble_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _TestInt.LineDouble___nonzero__(self)
    def __bool__(self): return _TestInt.LineDouble___bool__(self)
    def __len__(self): return _TestInt.LineDouble___len__(self)
    def pop(self): return _TestInt.LineDouble_pop(self)
    def __getslice__(self, *args): return _TestInt.LineDouble___getslice__(self, *args)
    def __setslice__(self, *args): return _TestInt.LineDouble___setslice__(self, *args)
    def __delslice__(self, *args): return _TestInt.LineDouble___delslice__(self, *args)
    def __delitem__(self, *args): return _TestInt.LineDouble___delitem__(self, *args)
    def __getitem__(self, *args): return _TestInt.LineDouble___getitem__(self, *args)
    def __setitem__(self, *args): return _TestInt.LineDouble___setitem__(self, *args)
    def append(self, *args): return _TestInt.LineDouble_append(self, *args)
    def empty(self): return _TestInt.LineDouble_empty(self)
    def size(self): return _TestInt.LineDouble_size(self)
    def clear(self): return _TestInt.LineDouble_clear(self)
    def swap(self, *args): return _TestInt.LineDouble_swap(self, *args)
    def get_allocator(self): return _TestInt.LineDouble_get_allocator(self)
    def begin(self): return _TestInt.LineDouble_begin(self)
    def end(self): return _TestInt.LineDouble_end(self)
    def rbegin(self): return _TestInt.LineDouble_rbegin(self)
    def rend(self): return _TestInt.LineDouble_rend(self)
    def pop_back(self): return _TestInt.LineDouble_pop_back(self)
    def erase(self, *args): return _TestInt.LineDouble_erase(self, *args)
    def __init__(self, *args): 
        this = _TestInt.new_LineDouble(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _TestInt.LineDouble_push_back(self, *args)
    def front(self): return _TestInt.LineDouble_front(self)
    def back(self): return _TestInt.LineDouble_back(self)
    def assign(self, *args): return _TestInt.LineDouble_assign(self, *args)
    def resize(self, *args): return _TestInt.LineDouble_resize(self, *args)
    def insert(self, *args): return _TestInt.LineDouble_insert(self, *args)
    def reserve(self, *args): return _TestInt.LineDouble_reserve(self, *args)
    def capacity(self): return _TestInt.LineDouble_capacity(self)
    __swig_destroy__ = _TestInt.delete_LineDouble
    __del__ = lambda self : None;
LineDouble_swigregister = _TestInt.LineDouble_swigregister
LineDouble_swigregister(LineDouble)

class ArrayInt(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ArrayInt, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ArrayInt, name)
    __repr__ = _swig_repr
    def iterator(self): return _TestInt.ArrayInt_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _TestInt.ArrayInt___nonzero__(self)
    def __bool__(self): return _TestInt.ArrayInt___bool__(self)
    def __len__(self): return _TestInt.ArrayInt___len__(self)
    def pop(self): return _TestInt.ArrayInt_pop(self)
    def __getslice__(self, *args): return _TestInt.ArrayInt___getslice__(self, *args)
    def __setslice__(self, *args): return _TestInt.ArrayInt___setslice__(self, *args)
    def __delslice__(self, *args): return _TestInt.ArrayInt___delslice__(self, *args)
    def __delitem__(self, *args): return _TestInt.ArrayInt___delitem__(self, *args)
    def __getitem__(self, *args): return _TestInt.ArrayInt___getitem__(self, *args)
    def __setitem__(self, *args): return _TestInt.ArrayInt___setitem__(self, *args)
    def append(self, *args): return _TestInt.ArrayInt_append(self, *args)
    def empty(self): return _TestInt.ArrayInt_empty(self)
    def size(self): return _TestInt.ArrayInt_size(self)
    def clear(self): return _TestInt.ArrayInt_clear(self)
    def swap(self, *args): return _TestInt.ArrayInt_swap(self, *args)
    def get_allocator(self): return _TestInt.ArrayInt_get_allocator(self)
    def begin(self): return _TestInt.ArrayInt_begin(self)
    def end(self): return _TestInt.ArrayInt_end(self)
    def rbegin(self): return _TestInt.ArrayInt_rbegin(self)
    def rend(self): return _TestInt.ArrayInt_rend(self)
    def pop_back(self): return _TestInt.ArrayInt_pop_back(self)
    def erase(self, *args): return _TestInt.ArrayInt_erase(self, *args)
    def __init__(self, *args): 
        this = _TestInt.new_ArrayInt(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _TestInt.ArrayInt_push_back(self, *args)
    def front(self): return _TestInt.ArrayInt_front(self)
    def back(self): return _TestInt.ArrayInt_back(self)
    def assign(self, *args): return _TestInt.ArrayInt_assign(self, *args)
    def resize(self, *args): return _TestInt.ArrayInt_resize(self, *args)
    def insert(self, *args): return _TestInt.ArrayInt_insert(self, *args)
    def reserve(self, *args): return _TestInt.ArrayInt_reserve(self, *args)
    def capacity(self): return _TestInt.ArrayInt_capacity(self)
    __swig_destroy__ = _TestInt.delete_ArrayInt
    __del__ = lambda self : None;
ArrayInt_swigregister = _TestInt.ArrayInt_swigregister
ArrayInt_swigregister(ArrayInt)

class ArrayDouble(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ArrayDouble, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ArrayDouble, name)
    __repr__ = _swig_repr
    def iterator(self): return _TestInt.ArrayDouble_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _TestInt.ArrayDouble___nonzero__(self)
    def __bool__(self): return _TestInt.ArrayDouble___bool__(self)
    def __len__(self): return _TestInt.ArrayDouble___len__(self)
    def pop(self): return _TestInt.ArrayDouble_pop(self)
    def __getslice__(self, *args): return _TestInt.ArrayDouble___getslice__(self, *args)
    def __setslice__(self, *args): return _TestInt.ArrayDouble___setslice__(self, *args)
    def __delslice__(self, *args): return _TestInt.ArrayDouble___delslice__(self, *args)
    def __delitem__(self, *args): return _TestInt.ArrayDouble___delitem__(self, *args)
    def __getitem__(self, *args): return _TestInt.ArrayDouble___getitem__(self, *args)
    def __setitem__(self, *args): return _TestInt.ArrayDouble___setitem__(self, *args)
    def append(self, *args): return _TestInt.ArrayDouble_append(self, *args)
    def empty(self): return _TestInt.ArrayDouble_empty(self)
    def size(self): return _TestInt.ArrayDouble_size(self)
    def clear(self): return _TestInt.ArrayDouble_clear(self)
    def swap(self, *args): return _TestInt.ArrayDouble_swap(self, *args)
    def get_allocator(self): return _TestInt.ArrayDouble_get_allocator(self)
    def begin(self): return _TestInt.ArrayDouble_begin(self)
    def end(self): return _TestInt.ArrayDouble_end(self)
    def rbegin(self): return _TestInt.ArrayDouble_rbegin(self)
    def rend(self): return _TestInt.ArrayDouble_rend(self)
    def pop_back(self): return _TestInt.ArrayDouble_pop_back(self)
    def erase(self, *args): return _TestInt.ArrayDouble_erase(self, *args)
    def __init__(self, *args): 
        this = _TestInt.new_ArrayDouble(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _TestInt.ArrayDouble_push_back(self, *args)
    def front(self): return _TestInt.ArrayDouble_front(self)
    def back(self): return _TestInt.ArrayDouble_back(self)
    def assign(self, *args): return _TestInt.ArrayDouble_assign(self, *args)
    def resize(self, *args): return _TestInt.ArrayDouble_resize(self, *args)
    def insert(self, *args): return _TestInt.ArrayDouble_insert(self, *args)
    def reserve(self, *args): return _TestInt.ArrayDouble_reserve(self, *args)
    def capacity(self): return _TestInt.ArrayDouble_capacity(self)
    __swig_destroy__ = _TestInt.delete_ArrayDouble
    __del__ = lambda self : None;
ArrayDouble_swigregister = _TestInt.ArrayDouble_swigregister
ArrayDouble_swigregister(ArrayDouble)

GaussianIntegration_h = _TestInt.GaussianIntegration_h
class QuadParams(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, QuadParams, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, QuadParams, name)
    __repr__ = _swig_repr
    __swig_setmethods__["nw"] = _TestInt.QuadParams_nw_set
    __swig_getmethods__["nw"] = _TestInt.QuadParams_nw_get
    if _newclass:nw = _swig_property(_TestInt.QuadParams_nw_get, _TestInt.QuadParams_nw_set)
    __swig_setmethods__["xw"] = _TestInt.QuadParams_xw_set
    __swig_getmethods__["xw"] = _TestInt.QuadParams_xw_get
    if _newclass:xw = _swig_property(_TestInt.QuadParams_xw_get, _TestInt.QuadParams_xw_set)
    __swig_setmethods__["w"] = _TestInt.QuadParams_w_set
    __swig_getmethods__["w"] = _TestInt.QuadParams_w_get
    if _newclass:w = _swig_property(_TestInt.QuadParams_w_get, _TestInt.QuadParams_w_set)
    def __init__(self): 
        this = _TestInt.new_QuadParams()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _TestInt.delete_QuadParams
    __del__ = lambda self : None;
QuadParams_swigregister = _TestInt.QuadParams_swigregister
QuadParams_swigregister(QuadParams)

class GaussianIntegration(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, GaussianIntegration, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, GaussianIntegration, name)
    __repr__ = _swig_repr
    def GaussInt_1D(self, *args): return _TestInt.GaussianIntegration_GaussInt_1D(self, *args)
    def GaussInt_2D_Serial(self, *args): return _TestInt.GaussianIntegration_GaussInt_2D_Serial(self, *args)
    def GaussInt_2D_Parallel(self, *args): return _TestInt.GaussianIntegration_GaussInt_2D_Parallel(self, *args)
    def GaussInt_3D_Serial(self, *args): return _TestInt.GaussianIntegration_GaussInt_3D_Serial(self, *args)
    def GaussInt_3D_Parallel(self, *args): return _TestInt.GaussianIntegration_GaussInt_3D_Parallel(self, *args)
    def __init__(self): 
        this = _TestInt.new_GaussianIntegration()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _TestInt.delete_GaussianIntegration
    __del__ = lambda self : None;
GaussianIntegration_swigregister = _TestInt.GaussianIntegration_swigregister
GaussianIntegration_swigregister(GaussianIntegration)

# This file is compatible with both classic and new-style classes.


