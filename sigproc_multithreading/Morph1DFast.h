/**
 * \file Morph1DFast.h
 *
 * \ingroup sigproc_multithreading
 * 
 * \brief Class def header for a class Morph1DFast
 *
 * @author koh0207
 */

/** \addtogroup sigproc_multithreading

    @{*/
#ifndef __SIGPROC_MULTITHREADING_MORPH1DFAST_H__
#define __SIGPROC_MULTITHREADING_MORPH1DFAST_H__

#include <iostream>
#include <iterator>
#include <numeric>
#include <cmath>
#include <limits>
#include <assert.h>

#include "SigprocParallelDefs.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"

namespace sigproc_multithreading {

  /**
     \class Morph1DFast
     User defined class Morph1DFast ... these comments are used to generate
     doxygen documentation!
  */
  class Morph1DFast{
    
    public:
      
      /// Default constructor
      Morph1DFast(){}

      void getDilation(const VectorBool&,
                       const unsigned int,
                       VectorBool&) const;

      void getDilation(const VectorShort&,
                       const unsigned int,
                       VectorShort&) const;

      void getDilation(const VectorFloat&,
                       const unsigned int,
                       VectorFloat&) const;

      void getDilation(const VectorDouble&,
                       const unsigned int,
                       VectorDouble&) const;

      // Column Major Operations

      void getDilation(const Array2D<bool>&,
                       const unsigned int,
                       Array2D<bool>&,
                       const unsigned int) const;

      void getDilation(const Array2D<short>&,
                       const unsigned int,
                       Array2D<short>&,
                       const unsigned int) const;

      void getDilation(const Array2D<float>&,
                       const unsigned int,
                       Array2D<float>&,
                       const unsigned int) const;

      void getDilation(const Array2D<double>&,
                       const unsigned int,
                       Array2D<double>&,
                       const unsigned int) const;


      void getErosion(const VectorBool&,
                      const unsigned int,
                      VectorBool&) const;

      void getErosion(const VectorShort&,
                      const unsigned int,
                      VectorShort&) const;

      void getErosion(const VectorFloat&,
                      const unsigned int,
                      VectorFloat&) const;

      void getErosion(const VectorDouble&,
                      const unsigned int,
                      VectorDouble&) const;

      // Column Major Operations

      void getErosion(const Array2D<bool>&,
                      const unsigned int,
                      Array2D<bool>&,
                      const unsigned int) const;

      void getErosion(const Array2D<short>&,
                      const unsigned int,
                      Array2D<short>&,
                      const unsigned int) const;

      void getErosion(const Array2D<float>&,
                      const unsigned int,
                      Array2D<float>&,
                      const unsigned int) const;

      void getErosion(const Array2D<double>&,
                      const unsigned int,
                      Array2D<double>&,
                      const unsigned int) const;


      void getGradient(const VectorShort&,
                      const unsigned int,
                      VectorShort&) const;

      void getGradient(const VectorFloat&,
                      const unsigned int,
                      VectorFloat&) const;

      void getGradient(const VectorDouble&,
                      const unsigned int,
                      VectorDouble&) const;


      void getOpening(const VectorBool&,
                      const unsigned int,
                      VectorBool&) const;

      void getOpening(const VectorShort&,
                      const unsigned int,
                      VectorShort&) const;

      void getOpening(const VectorFloat&,
                      const unsigned int,
                      VectorFloat&) const;

      void getOpening(const VectorDouble&,
                      const unsigned int,
                      VectorDouble&) const;


      void getClosing(const VectorBool&,
                      const unsigned int,
                      VectorBool&) const;            

      void getClosing(const VectorShort&,
                      const unsigned int,
                      VectorShort&) const;

      void getClosing(const VectorFloat&,
                      const unsigned int,
                      VectorFloat&) const;

      void getClosing(const VectorDouble&,
                      const unsigned int,
                      VectorDouble&) const;

    
    /// Default destructor
      ~Morph1DFast(){}

    private:

      // Row Major Operations

      template <typename T> 
      void getDilation(
        const Vector<T>& inputVector,
        const unsigned int fStructuringElement,
        Vector<T>& dilationVec) const;

      template <typename T> 
      void getErosion(
        const Vector<T>& inputVector,
        const unsigned int fStructuringElement,
        Vector<T>& erosionVec) const;

      template <typename T>
      void getOpening(
        const Vector<T>& inputVector,
        const unsigned int fStructuringElement,
        Vector<T>& openingVec) const;

      template <typename T>
      void getClosing(
        const Vector<T>& inputVector,
        const unsigned int fStructuringElement,
        Vector<T>& closingVec) const;

      template <typename T> 
      void getGradient(
        const Vector<T>& inputVector,
        const unsigned int fStructuringElement,
        Vector<T>& gradientVec) const;

      template <typename T>
      void getDilation(
        const Array2D<T>& inputArray2D,
        const unsigned int fStructuringElement,
        Array2D<T>& dilation2D,
        const int columnNum) const;

      template <typename T>
      void getErosion(
        const Array2D<T>& inputArray2D,
        const unsigned int fStructuringElement,
        Array2D<T>& erosion2D,
        const int columnNum) const;

      
    
  };
}

#endif
/** @} */ // end of doxygen group 

