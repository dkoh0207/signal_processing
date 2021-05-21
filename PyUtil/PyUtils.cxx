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

    //
    // Float std=>py
    //
    
    /*
    PyObject* as_ndarray(const std::vector<float>& data) {
      SetPyUtil();
      npy_intp dim_data[1];
      dim_data[0] = data.size();
      return PyArray_SimpleNewFromData(1, dim_data, NPY_FLOAT, (char *)&(data[0]));
    }
    */

    PyObject* as_ndarray(const std::vector<float>& data) {
      SetPyUtil();
      npy_intp dims_[1];
      dims_[0] = data.size();
      PyObject* res = PyArray_ZEROS(1,dims_,NPY_FLOAT32,0);
      PyArrayObject *ptr = (PyArrayObject*)(res);
      npy_intp loc_[1];
      loc_[0] = 0;
      auto fptr = (float*)(PyArray_GetPtr(ptr, loc_));
      for (size_t i = 0; i < data.size(); ++i)
	fptr[i] = data[i];
      
      //PyArray_INCREF(ptr);
      //PyArray_Free(ptr, (void*)fptr);
      return res;
    }


    PyObject* as_ndarray(const std::vector<int>& data) {
      SetPyUtil();
      npy_intp dims_[1];
      dims_[0] = data.size();
      PyObject* res = PyArray_ZEROS(1,dims_,NPY_INT,0);
      PyArrayObject *ptr = (PyArrayObject*)(res);
      npy_intp loc_[1];
      loc_[0] = 0;
      auto fptr = (int*)(PyArray_GetPtr(ptr, loc_));
      for (size_t i = 0; i < data.size(); ++i)
	fptr[i] = data[i];
      
      //PyArray_INCREF(ptr);
      //PyArray_Free(ptr, (void*)fptr);
      return res;
    }


    PyObject* as_ndarray(const std::vector<std::vector<float> >& data) {
      SetPyUtil();
      npy_intp dim_data[2];
      dim_data[0] = data.size();
      dim_data[1] = data[0].size();
      auto result = (PyArrayObject*)(PyArray_SimpleNew(2,dim_data,NPY_FLOAT32));
      auto data_ptr = (float*)(PyArray_DATA(result));
      for (size_t i=0; i<data.size(); ++i) {
	if( (int)(data[i].size()) != dim_data[1]) {
	  Py_DECREF(result); // delete
	  std::cerr<<"data size mismatch among rows: " << dim_data[1]
		   << " vs " << data[i].size() << std::endl;
	  throw std::exception();
	}
	std::copy(data[i].begin(),data[i].end(),data_ptr+i*dim_data[1]);
      }
      return (PyObject*)result;
    }

    /*
    // This implementation assumes underlying data stays as is (hence static)
    PyObject* as_ndarray(const std::vector<std::vector<float> >& data) {
      SetPyUtil();
      npy_intp dim_data[2];
      dim_data[0] = data.size();
      dim_data[1] = data[0].size();
      for(auto const& d : data) { assert((int)(d.size()) == dim_data[1]); }
      static std::vector<float> source;
      source.resize(data.size() * data[0].size(),0.);
      size_t ctr=0;
      for(int i=0; i<dim_data[0]; ++i) {
	for(int j=0; j<dim_data[1]; ++j) {
	  source[ctr] = data[i][j];
	  ctr++;
	}
      }
      return PyArray_SimpleNewFromData(2, dim_data, NPY_FLOAT, (char *)&(source[0]));
    }
    */

    /*
    PyObject* as_ndarray(const std::vector<std::vector<float> >& data) {
      SetPyUtil();
      npy_intp dims_[2];
      dims_[0] = data.size();
      dims_[1] = data[0].size();
      // make sure it's NxM
      for(auto const& d : data) { assert((int)(d.size()) == dims_[1]); }
      PyObject* res = PyArray_ZEROS(2,dims_,NPY_FLOAT32,0);
      //PyArrayObject *ptr = (PyArrayObject*)(res);
      float **carray;
      //const int dtype = NPY_FLOAT;
      PyArray_Descr *descr = PyArray_DescrFromType(NPY_FLOAT32);
      if (PyArray_AsCArray(&res, (void *)&carray, dims_, 2, descr) < 0) {
	std::cerr<<"ERROR: cannot convert pyarray to 2D C-array"<<std::endl;
	throw std::exception();
      }
      for (int i = 0; i < dims_[0]; ++i) {
	for (int j = 0; j < dims_[1]; ++j)
	  carray[i][j] = data[i][j];
      }
      //PyArray_INCREF(ptr);
      PyArray_Free(res, (void*)carray);
      return res;
    }
    */

    //
    // Bool std=>py
    //
    PyObject* as_ndarray(const std::vector<bool>& data) {
      SetPyUtil();
      npy_intp dims_[1];
      dims_[0] = data.size();
      PyObject* res = PyArray_ZEROS(1,dims_,NPY_BOOL,0);
      PyArrayObject *ptr = (PyArrayObject*)(res);
      npy_intp loc_[1];
      loc_[0] = 0;
      auto fptr = (bool*)(PyArray_GetPtr(ptr, loc_));
      for (size_t i = 0; i < data.size(); ++i)
	fptr[i] = data[i];
      
      //PyArray_INCREF(ptr);
      //PyArray_Free(ptr, (void*)fptr);
      return res;
    }


    PyObject* as_ndarray(const std::vector<std::vector<bool> >& data) {
      SetPyUtil();
      npy_intp dim_data[2];
      dim_data[0] = data.size();
      dim_data[1] = data[0].size();
      auto result = (PyArrayObject*)(PyArray_SimpleNew(2,dim_data,NPY_BOOL));
      auto data_ptr = (bool*)(PyArray_DATA(result));
      for (size_t i=0; i<data.size(); ++i) {
	if( (int)(data[i].size()) != dim_data[1]) {
	  Py_DECREF(result); // delete
	  std::cerr<<"data size mismatch among rows: " << dim_data[1]
		   << " vs " << data[i].size() << std::endl;
	  throw std::exception();
	}
	std::copy(data[i].begin(),data[i].end(),data_ptr+i*dim_data[1]);
      }
      return (PyObject*)result;
    }

    //
    // Float py=>std
    //
    
    std::vector<float> as_float32_vector(PyObject* pyarray) {
      ::sproc::pyutil::SetPyUtil();
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
      ::sproc::pyutil::SetPyUtil();
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

    std::vector<std::vector<int> > as_int_vector_2d(PyObject* pyarray) {
      ::sproc::pyutil::SetPyUtil();
      int **carray;
      const int dtype = NPY_INT;
      PyArray_Descr *descr = PyArray_DescrFromType(dtype);
      npy_intp dims[2];
      if (PyArray_AsCArray(&pyarray, (void *)&carray, dims, 2, descr) < 0) {
	std::cerr<<"ERROR: cannot convert pyarray to 2D C-array"<<std::endl;
	throw std::exception();
      }
      
      std::vector<std::vector<int > > data(dims[0],std::vector<int>(dims[1]));

      for(int i=0; i<dims[0]; ++i) {
	for(int j=0; j<dims[1]; ++j) {
	  data[i][j] = carray[i][j];
	}
      }
      
      PyArray_Free(pyarray,  (void *)carray);
      
      return data;
    }

    std::vector<bool> as_bool_vector(PyObject* pyarray) {
      ::sproc::pyutil::SetPyUtil();
      bool *carray;
      const int dtype = NPY_BOOL;
      PyArray_Descr *descr = PyArray_DescrFromType(dtype);
      npy_intp dims[1];
      if (PyArray_AsCArray(&pyarray, (void *)&carray, dims, 1, descr) < 0) {
	std::cerr<<"ERROR: cannot convert pyarray to 1D C-array"<<std::endl;
	throw std::exception();
      }
      
      size_t npts = dims[0];
      std::vector<bool> data(npts);
      
      for(size_t i=0; i<npts; ++i)
	data[i] = carray[i];
      PyArray_Free(pyarray,  (void *)carray);
      
      return data;
    }

    std::vector<std::vector<bool> > as_bool_vector_2d(PyObject* pyarray) {
      ::sproc::pyutil::SetPyUtil();
      bool **carray;
      const int dtype = NPY_BOOL;
      PyArray_Descr *descr = PyArray_DescrFromType(dtype);
      npy_intp dims[2];
      if (PyArray_AsCArray(&pyarray, (void *)&carray, dims, 2, descr) < 0) {
	std::cerr<<"ERROR: cannot convert pyarray to 2D C-array"<<std::endl;
	throw std::exception();
      }
      
      std::vector<std::vector<bool > > data(dims[0],std::vector<bool>(dims[1]));

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
