/**
 * \file SigprocParallelDefs.h
 *
 * \ingroup sigproc_multithreading
 * 
 * \brief Class def header for a class SigprocParallelDefs
 *
 * @author koh0207
 */

/** \addtogroup sigproc_multithreading

    @{*/
#ifndef __SIGPROC_MULTITHREADING_SIGPROCPARALLELDEFS_H__
#define __SIGPROC_MULTITHREADING_SIGPROCPARALLELDEFS_H__

#include <vector>
#include "tbb/concurrent_vector.h"
#include <map>

namespace sigproc_multithreading {

      template <class T> 
      using Vector        = std::vector<T>;

      using VectorShort   = Vector<short>;
      using VectorInt     = Vector<int>;
      using VectorFloat   = Vector<float>;
      using VectorDouble  = Vector<double>;
      using VectorBool    = Vector<bool>;

      template <class T> 
      using Array2D        = std::vector<Vector<T>>;

      using ArrayShort    = Vector<VectorShort>;
      using ArrayInt      = Vector<VectorInt>;
      using ArrayFloat    = Vector<VectorFloat>;
      using ArrayDouble   = Vector<VectorDouble>;
      using ArrayBool     = Vector<VectorBool>;

      template <class T>
      using ConcurrentVector = tbb::concurrent_vector<T>;

      template <class T>
      using ConcurrentArray2D = ConcurrentVector<ConcurrentVector<T>>;
}

#endif
/** @} */ // end of doxygen group 

