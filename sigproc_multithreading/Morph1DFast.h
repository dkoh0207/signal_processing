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

      void getDilation(const ConcurrentVector<bool>&,
                       const unsigned int,
                       ConcurrentVector<bool>&) const;

      void getDilation(const ConcurrentVector<short>&,
                       const unsigned int,
                       ConcurrentVector<short>&) const;

      void getDilation(const ConcurrentVector<float>&,
                       const unsigned int,
                       ConcurrentVector<float>&) const;

      void getDilation(const ConcurrentVector<double>&,
                       const unsigned int,
                       ConcurrentVector<double>&) const;

      // Column Major Operations

      void getDilation(const ConcurrentArray2D<bool>&,
                       const unsigned int,
                       ConcurrentArray2D<bool>&,
                       const unsigned int) const;

      void getDilation(const ConcurrentArray2D<short>&,
                       const unsigned int,
                       ConcurrentArray2D<short>&,
                       const unsigned int) const;

      void getDilation(const ConcurrentArray2D<float>&,
                       const unsigned int,
                       ConcurrentArray2D<float>&,
                       const unsigned int) const;

      void getDilation(const ConcurrentArray2D<double>&,
                       const unsigned int,
                       ConcurrentArray2D<double>&,
                       const unsigned int) const;


      void getErosion(const ConcurrentVector<bool>&,
                      const unsigned int,
                      ConcurrentVector<bool>&) const;

      void getErosion(const ConcurrentVector<short>&,
                      const unsigned int,
                      ConcurrentVector<short>&) const;

      void getErosion(const ConcurrentVector<float>&,
                      const unsigned int,
                      ConcurrentVector<float>&) const;

      void getErosion(const ConcurrentVector<double>&,
                      const unsigned int,
                      ConcurrentVector<double>&) const;

      // Column Major Operations

      void getErosion(const ConcurrentArray2D<bool>&,
                      const unsigned int,
                      ConcurrentArray2D<bool>&,
                      const unsigned int) const;

      void getErosion(const ConcurrentArray2D<short>&,
                      const unsigned int,
                      ConcurrentArray2D<short>&,
                      const unsigned int) const;

      void getErosion(const ConcurrentArray2D<float>&,
                      const unsigned int,
                      ConcurrentArray2D<float>&,
                      const unsigned int) const;

      void getErosion(const ConcurrentArray2D<double>&,
                      const unsigned int,
                      ConcurrentArray2D<double>&,
                      const unsigned int) const;


      void getGradient(const ConcurrentVector<short>&,
                      const unsigned int,
                      ConcurrentVector<short>&) const;

      void getGradient(const ConcurrentVector<float>&,
                      const unsigned int,
                      ConcurrentVector<float>&) const;

      void getGradient(const ConcurrentVector<double>&,
                      const unsigned int,
                      ConcurrentVector<double>&) const;


      void getOpening(const ConcurrentVector<bool>&,
                      const unsigned int,
                      ConcurrentVector<bool>&) const;

      void getOpening(const ConcurrentVector<short>&,
                      const unsigned int,
                      ConcurrentVector<short>&) const;

      void getOpening(const ConcurrentVector<float>&,
                      const unsigned int,
                      ConcurrentVector<float>&) const;

      void getOpening(const ConcurrentVector<double>&,
                      const unsigned int,
                      ConcurrentVector<double>&) const;


      void getClosing(const ConcurrentVector<bool>&,
                      const unsigned int,
                      ConcurrentVector<bool>&) const;            

      void getClosing(const ConcurrentVector<short>&,
                      const unsigned int,
                      ConcurrentVector<short>&) const;

      void getClosing(const ConcurrentVector<float>&,
                      const unsigned int,
                      ConcurrentVector<float>&) const;

      void getClosing(const ConcurrentVector<double>&,
                      const unsigned int,
                      ConcurrentVector<double>&) const;

    
    /// Default destructor
      ~Morph1DFast(){}

    private:

      // Row Major Operations

      template <typename T> 
      void getDilation(
        const ConcurrentVector<T>& inputVector,
        const unsigned int structuringElement,
        ConcurrentVector<T>& dilationVec) const;

      template <typename T> 
      void getErosion(
        const ConcurrentVector<T>& inputVector,
        const unsigned int structuringElement,
        ConcurrentVector<T>& erosionVec) const;

      template <typename T>
      void getOpening(
        const ConcurrentVector<T>& inputVector,
        const unsigned int structuringElement,
        ConcurrentVector<T>& openingVec) const;

      template <typename T>
      void getClosing(
        const ConcurrentVector<T>& inputVector,
        const unsigned int structuringElement,
        ConcurrentVector<T>& closingVec) const;

      template <typename T> 
      void getGradient(
        const ConcurrentVector<T>& inputVector,
        const unsigned int structuringElement,
        ConcurrentVector<T>& gradientVec) const;


      // Column Major Operations 

      template <typename T> 
      void getDilation(
        const ConcurrentArray2D<T>& inputArray2D,
        const unsigned int structuringElementy,
        ConcurrentArray2D<T>& dilation2D,
        const unsigned int columnNum) const;

      template <typename T> 
      void getErosion(
        const ConcurrentArray2D<T>& inputArray2D,
        const unsigned int structuringElementy,
        ConcurrentArray2D<T>& erosion2D,
        const unsigned int columnNum) const;
    
  };
}

#endif
/** @} */ // end of doxygen group 

