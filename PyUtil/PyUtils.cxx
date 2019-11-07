#ifndef __LARCV_PYUTILS_CXX__
#define __LARCV_PYUTILS_CXX__

#include "PyUtils.h"
#include <iostream>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#include <numpy/ndarrayobject.h>
#include "numpy/arrayobject.h"
#include <cassert>

namespace sproc {
  namespace pyutil {
    
    void SetPyUtil() {
      static bool once = false;
      if (!once) {
	_import_array();
	once = true;
      }
    }
    
    std::vector<float> as_float32_vector(PyObject* pyarray) {
      SetPyUtil();
      float *carray;
      const int dtype = NPY_FLOAT;
      PyArray_Descr *descr = PyArray_DescrFromType(dtype);
      npy_intp dims[1];
      if (PyArray_AsCArray(&pyarray, (void *)&carray, dims, 1, descr) < 0) {
	std::cerr<<"ERROR: cannot convert pyarray to 1D C-array"<<std::endl;
	throw std::exception();
      }
      
      size_t npts = dims[0];
      std::vector<float> data(npts);
      
      for(size_t i=0; i<npts; ++i)
	data[i] = carray[i];
      PyArray_Free(pyarray,  (void *)carray);
      
      return data;
    }

    std::vector<std::vector<float> > as_float32_vector_2d(PyObject* pyarray) {
      SetPyUtil();
      float **carray;
      const int dtype = NPY_FLOAT;
      PyArray_Descr *descr = PyArray_DescrFromType(dtype);
      npy_intp dims[2];
      if (PyArray_AsCArray(&pyarray, (void *)&carray, dims, 2, descr) < 0) {
	std::cerr<<"ERROR: cannot convert pyarray to 2D C-array"<<std::endl;
	throw std::exception();
      }
      
      std::vector<std::vector<float > > data(dims[0],std::vector<float>(dims[1]));

      for(int i=0; i<dims[0]; ++i) {
	for(int j=0; j<dims[1]; ++j) {
	  data[i][j] = carray[i][j];
	}
      }
      
      PyArray_Free(pyarray,  (void *)carray);
      
      return data;
    }

    
  }
}
#endif
