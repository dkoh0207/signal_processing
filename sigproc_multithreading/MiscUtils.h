/**
 * \file MiscUtils.h
 *
 * \ingroup sigproc_multithreading
 * 
 * \brief Class def header for a class MiscUtils
 *
 * @author koh0207
 */

/** \addtogroup sigproc_multithreading

    @{*/
#ifndef __SIGPROC_MULTITHREADING_MISCUTILS_H__
#define __SIGPROC_MULTITHREADING_MISCUTILS_H__

#include "SigprocParallelDefs.h"

namespace sigproc_multithreading {

  /**
     \class MiscUtils
     User defined class MiscUtils ... these comments are used to generate
     doxygen documentation!
  */
  class MiscUtils{
    
    public:
      
      /// Default constructor
      MiscUtils(){}

      ConcurrentArray2D<float> converttoConcurrent(
        const Array2D<float>& inputSTL) const;

      ConcurrentArray2D<bool> converttoConcurrent(
        const Array2D<bool>& inputSTL) const;

      Array2D<float> converttoSTL(
        const ConcurrentArray2D<float>& inputTBB) const;

      Array2D<bool> converttoSTL(
        const ConcurrentArray2D<bool>& inputTBB) const;
      
      /// Default destructor
      ~MiscUtils(){}

    private:
      
      template <typename T>
      ConcurrentArray2D<T> converttoConcurrent(
        const Array2D<T>& inputSTL) const;

      template <typename T>
      Array2D<T> converttoSTL(const ConcurrentArray2D<T>& inputTBB) const;
    
  };
}

#endif
/** @} */ // end of doxygen group 

