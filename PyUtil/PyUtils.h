#ifndef __LARCV_PYUTILS_H__
#define __LARCV_PYUTILS_H__

struct _object;
typedef _object PyObject;

//#ifndef __CLING__
//#ifndef __CINT__
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#include <numpy/ndarrayobject.h>
#include "numpy/arrayobject.h"
//#endif
//#endif

#include <vector>
namespace sproc {
  namespace pyutil {
    /// Utility function: call one-time-only numpy module initialization (you don't
    /// have to call)
    void SetPyUtil();
    /// A function to convert 1D float pyarray to std::vector
    std::vector<float> as_float32_vector(PyObject* pyarray);
    /// A function to convert 2D float pyarray to std::vector
    std::vector<std::vector<float> > as_float32_vector_2d(PyObject* pyarray);
    
    class load_pyutil {
      load_pyutil(){}
    };
  }
}

#endif
