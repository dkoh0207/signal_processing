/**
 * \file CoherentNoiseCorrection.h
 *
 * \ingroup sigproc_multithreading
 * 
 * \brief Class def header for a class CoherentNoiseCorrection
 *
 * @author koh0207
 */

/** \addtogroup sigproc_multithreading

    @{*/
#ifndef __SIGPROC_MULTITHREADING_COHERENTNOISECORRECTION_H__
#define __SIGPROC_MULTITHREADING_COHERENTNOISECORRECTION_H__

#include "SigprocParallelDefs.h"
#include "Morph2DFast.h"

namespace sigproc_multithreading {

  /**
     \class CoherentNoiseCorrection
     User defined class CoherentNoiseCorrection ... these comments are used to generate
     doxygen documentation!
  */
  class CoherentNoiseCorrection{
    
    public:

      short computeMedian(const ConcurrentVector<short> &cVector) const;

      float computeMedian(const ConcurrentVector<float> &cVector) const;

      double computeMedian(const ConcurrentVector<double> &cVector) const;


      void getSelectVals(const Array2D<short> &morphedWaveforms,
                          Array2D<bool> &selectVals,
                          const short threshold) const;

      void getSelectVals(const Array2D<float> &morphedWaveforms,
                          Array2D<bool> &selectVals,
                          const float threshold) const;

      void getSelectVals(const Array2D<double> &morphedWaveforms,
                          Array2D<bool> &selectVals,
                          const double threshold) const;

      void ParallelDenoiseMorph2D(
          const Array2D<short> &inputArray2D,
          Array2D<bool> &selectVals,
          const char filterName,
          Array2D<short> &waveLessCoherent,
          const size_t fStructuringElementx,
          const size_t fStructuringElementy,
          const size_t grouping,
          const size_t groupingOffset,
          const short threshold) const;

      void ParallelDenoiseMorph2D(
          const Array2D<float> &inputArray2D,
          Array2D<bool> &selectVals,
          const char filterName,
          Array2D<float> &waveLessCoherent,
          const size_t fStructuringElementx,
          const size_t fStructuringElementy,
          const size_t grouping,
          const size_t groupingOffset,
          const float threshold) const;

      void ParallelDenoiseMorph2D(
          const Array2D<double> &inputArray2D,
          Array2D<bool> &selectVals,
          const char filterName,
          Array2D<double> &waveLessCoherent,
          const size_t fStructuringElementx,
          const size_t fStructuringElementy,
          const size_t grouping,
          const size_t groupingOffset,
          const double threshold) const;
      
      /// Default constructor
      CoherentNoiseCorrection(){}
      
      /// Default destructor
      ~CoherentNoiseCorrection(){}

    private:

      template <typename T>
      T computeMedian(const ConcurrentVector<T> &cVector) const;

      template <typename T>
      void getSelectVals(
          const Array2D<T> &morphedWaveforms,
          Array2D<bool> &selectVals,
          const T threshold) const;

      template <typename T>
      void ParallelDenoiseSimple(
          const Array2D<T> &inputArray2D,
          Array2D<T> &waveLessCoherent,
          const size_t grouping,
          const size_t groupingOffset) const;

      template <typename T>
      void ParallelDenoiseMorph2D(
          const Array2D<T> &inputArray2D,
          Array2D<bool> &selectVals,
          const char filterName,
          Array2D<T> &waveLessCoherent,
          const size_t fStructuringElementx,
          const size_t fStructuringElementy,
          const size_t grouping,
          const size_t groupingOffset,
          const T threshold) const;

      template <typename T>
      void ParallelDenoiseHough2D(
          const Array2D<T> &inputArray2D,
          const char filterName,
          Array2D<T> &waveLessCoherent) const;
    
  };
}

#endif
/** @} */ // end of doxygen group 

