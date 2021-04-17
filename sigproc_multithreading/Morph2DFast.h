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

        void getDilation(const Array2D<bool>&,
                        const unsigned int,
                        const unsigned int,
                        Array2D<bool>&) const;

        void getDilation(const Array2D<float>&,
                        const unsigned int,
                        const unsigned int,
                        Array2D<float>&) const;

        void getErosion(const Array2D<bool>&,
                        const unsigned int,
                        const unsigned int,
                        Array2D<bool>&) const;

        void getErosion(const Array2D<float>&,
                        const unsigned int,
                        const unsigned int,
                        Array2D<float>&) const;

        void getGradient(const Array2D<short>&,
                         const unsigned int,
                         const unsigned int,
                         Array2D<short>&) const;
                      
        void getGradient(const Array2D<float>&,
                         const unsigned int,
                         const unsigned int,
                         Array2D<float>&) const;

        void getGradient(const Array2D<double>&,
                         const unsigned int,
                         const unsigned int,
                         Array2D<double>&) const;

        void getOpening(const Array2D<bool>&,
                        const unsigned int,
                        const unsigned int,
                        Array2D<bool>&) const;

        void getOpening(const Array2D<float>&,
                        const unsigned int,
                        const unsigned int,
                        Array2D<float>&) const;


        void getClosing(const Array2D<bool>&,
                        const unsigned int,
                        const unsigned int,
                        Array2D<bool>&) const;

        void getClosing(const Array2D<float>&,
                        const unsigned int,
                        const unsigned int,
                        Array2D<float>&) const;

        template <typename T>
        void getDilation(
          const Array2D<T>& waveform2D,
          const unsigned int structuringElementx,
          const unsigned int structuringElementy,
          Array2D<T>& dilation2D) const;

        template <typename T>
        void getErosion(
          const Array2D<T>& waveform2D,
          const unsigned int structuringElementx,
          const unsigned int structuringElementy,
          Array2D<T>& erosion2D) const;

        template <typename T>
        void getGradient(const Array2D<T>&,
                         const unsigned int,
                         const unsigned int,
                         Array2D<T>&) const;

        template <typename T>
        void getClosing(
          const Array2D<T>& waveform2D,
          const unsigned int structuringElementx,
          const unsigned int structuringElementy,
          Array2D<T>& closing2D) const;

        template <typename T>
        void getOpening(
          const Array2D<T>& waveform2D,
          const unsigned int structuringElementx,
          const unsigned int structuringElementy,
          Array2D<T>& opening2D) const;

      /// Default destructor
      ~Morph2DFast(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

