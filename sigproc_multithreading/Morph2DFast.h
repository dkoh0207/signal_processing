/**
 * \file Morph2DFast.h
 *
 * \ingroup sigproc_multithreading
 * 
 * \brief Class def header for a class Morph2DFast
 *
 * @author koh0207
 */

/** \addtogroup sigproc_multithreading

    @{*/
#ifndef __SIGPROC_MULTITHREADING_MORPH2DFAST_H__
#define __SIGPROC_MULTITHREADING_MORPH2DFAST_H__

#include "SigprocParallelDefs.h"
#include "Morph1DFast.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"

#include <chrono>

namespace sigproc_multithreading {

  /**
     \class Morph2DFast
     User defined class Morph2DFast ... these comments are used to generate
     doxygen documentation!
  */
  class Morph2DFast{
    
    public:

      /// Default constructor
      Morph2DFast(){}

        void getDilation(const ConcurrentArray2D<bool>&,
                        const unsigned int,
                        const unsigned int,
                        ConcurrentArray2D<bool>&) const;

        void getDilation(const ConcurrentArray2D<float>&,
                        const unsigned int,
                        const unsigned int,
                        ConcurrentArray2D<float>&) const;


        void getErosion(const ConcurrentArray2D<bool>&,
                        const unsigned int,
                        const unsigned int,
                        ConcurrentArray2D<bool>&) const;

        void getErosion(const ConcurrentArray2D<float>&,
                        const unsigned int,
                        const unsigned int,
                        ConcurrentArray2D<float>&) const;


        void getGradient(const ConcurrentArray2D<float>&,
                        const unsigned int,
                        const unsigned int,
                        ConcurrentArray2D<float>&) const;


        void getOpening(const ConcurrentArray2D<bool>&,
                        const unsigned int,
                        const unsigned int,
                        ConcurrentArray2D<bool>&) const;

        void getOpening(const ConcurrentArray2D<float>&,
                        const unsigned int,
                        const unsigned int,
                        ConcurrentArray2D<float>&) const;


        void getClosing(const ConcurrentArray2D<bool>&,
                        const unsigned int,
                        const unsigned int,
                        ConcurrentArray2D<bool>&) const;

        void getClosing(const ConcurrentArray2D<float>&,
                        const unsigned int,
                        const unsigned int,
                        ConcurrentArray2D<float>&) const;

      ConcurrentArray2D<float> converttoConcurrent(
        const Array2D<float>& inputSTL) const;

      ConcurrentArray2D<bool> converttoConcurrent(
        const Array2D<bool>& inputSTL) const;

      Array2D<float> converttoSTL(
        const ConcurrentArray2D<float>& inputTBB) const;

      Array2D<bool> converttoSTL(
        const ConcurrentArray2D<bool>& inputTBB) const;

        // void contrastStretching(const ConcurrentArray2D<float>&,
        //                         const unsigned int,
        //                         const unsigned int,
        //                         ConcurrentArray2D<float>&,
        //                         const float) const;
      /// Default destructor
      ~Morph2DFast(){}

    private:

      template <typename T>
      void getDilation(
        const ConcurrentArray2D<T>& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        ConcurrentArray2D<T>& dilation2D) const;

      template <typename T>
      void getErosion(
        const ConcurrentArray2D<T>& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        ConcurrentArray2D<T>& erosion2D) const;

      template <typename T>
      void getClosing(
        const ConcurrentArray2D<T>& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        ConcurrentArray2D<T>& closing2D) const;

      template <typename T>
      void getOpening(
        const ConcurrentArray2D<T>& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        ConcurrentArray2D<T>& opening2D) const;

      template <typename T>
      ConcurrentArray2D<T> converttoConcurrent(
        const Array2D<T>& inputSTL) const;

      template <typename T>
      Array2D<T> converttoSTL(const ConcurrentArray2D<T>& inputTBB) const;
    
  };
}

#endif
/** @} */ // end of doxygen group 

