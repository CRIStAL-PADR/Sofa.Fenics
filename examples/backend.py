"""Compilation backend for Sofa.Fenics"""

import ffcx.main

def compile_to_c_code():
    import ffcx
    import ffcx.main

    ffcx.main.main(["fenics-hyperelasticity.py"])

def compile_to_binary():
    from cffi import FFI
    ffibuilder = FFI()

    ffibuilder.cdef("float pi_approx(int n);")
    ffibuilder.cdef("void* get_fenics_function();")

    ffibuilder.set_source("_fenics2",  # name of the output C extension
        """
        #include "fenics.h"
        #include "ufcx.h"
        """,
        sources=['fenics.c'],   # includes pi.c as additional sources
        libraries=['m'])        # on Unix, link with the math library

    ffibuilder.compile(verbose=True, debug=False)

    import _fenics2
    return _fenics2


compile_to_c_code()
compile_to_binary()
