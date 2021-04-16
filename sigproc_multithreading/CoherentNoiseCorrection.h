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
      
      /// Default constructor
      CoherentNoiseCorrection(){}

        float computeMedian(const ConcurrentVector<float>&) const;

        void getSelectVals(
          const ConcurrentArray2D<float>&,
          ConcurrentArray2D<bool>&,
          const unsigned int,
          const float) const;

        void denoiseMorph2D(
          ConcurrentArray2D<float>&,
          ConcurrentArray2D<float>&,
          const ConcurrentArray2D<float>&,
          ConcurrentArray2D<bool>&,
          const char,
          const unsigned int,
          const unsigned int,
          const unsigned int,
          const unsigned int,
          const unsigned int,
          const float) const;

        void denoiseHough2D(
          ConcurrentArray2D<float>& waveLessCoherent,
          ConcurrentArray2D<float>& morphedWaveforms,
          const ConcurrentArray2D<float>& fullEvent,
          ConcurrentArray2D<bool>& selectVals,
          ConcurrentArray2D<bool>& refinedSelectVals,
          const char filterName = 'g',
          const unsigned int grouping = 64,
          const unsigned int groupingOffset = 0,
          const unsigned int structuringElementx = 5,
          const unsigned int structuringElementy = 20,
          const unsigned int window = 0,
          const float thresholdFactor = 2.5,
          const size_t thetaSteps = 360,
          const unsigned int houghThreshold = 300,
          const unsigned int nmsWindowSize = 10,
          const unsigned int angleWindow = 50, 
          const unsigned int dilationX = 5,
          const unsigned int dilationY = 20,
          const unsigned int maxLines = 20,
          const float eps = 0.00001) const;
      
      /// Default destructor
      ~CoherentNoiseCorrection(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

