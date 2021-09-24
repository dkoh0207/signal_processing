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

#include "Morph1DFast.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"
#include "tbb/concurrent_vector.h"


namespace sigproc_multithreading {

  /**
     \class Morph2DFast
     User defined class Morph2DFast ... these comments are used to generate
     doxygen documentation!
  */
    template <class T>
    using ConcurrentVector = tbb::concurrent_vector<T>;

    template <class T>
    using ConcurrentArray2D = ConcurrentVector<ConcurrentVector<T>>;
  class Morph2DFast{
    
    public:

      /// Default constructor
        Morph2DFast(){};

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

        void getGradient(const ConcurrentArray2D<short>&,
                         const unsigned int,
                         const unsigned int,
                         ConcurrentArray2D<short>&) const;
                      
        void getGradient(const ConcurrentArray2D<float>&,
                         const unsigned int,
                         const unsigned int,
                         ConcurrentArray2D<float>&) const;

        void getGradient(const ConcurrentArray2D<double>&,
                         const unsigned int,
                         const unsigned int,
                         ConcurrentArray2D<double>&) const;

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
        void getGradient(const ConcurrentArray2D<T>&,
                         const unsigned int,
                         const unsigned int,
                         ConcurrentArray2D<T>&) const;

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

      /// Default destructor
      ~Morph2DFast(){};
    
  };
}

#endif
/** @} */ // end of doxygen group 

