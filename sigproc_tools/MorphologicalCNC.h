/**
 * \file MorphologicalCNC.h
 *
 * \ingroup sigproc_tools
 *
 * \brief Class def header for a class MorphologicalCNC
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_MORPHOLOGICALCNC_H__
#define __SIGPROC_TOOLS_MORPHOLOGICALCNC_H__

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>

#include "Morph1DFast.h"
#include "Morph2DFast.h"
#include "MiscUtils.h"
#include "LineDetection.h"

namespace sigproc_tools {

  /**
     \class MorphologicalCNC
     User defined class MorphologicalCNC ... these comments are used to generate
     doxygen documentation!
  */

  template <class T> using Array2D = std::vector<std::vector<T>>;

  class MorphologicalCNC{

    public:

      /// Default constructor
      MorphologicalCNC(){}

      void getSelectVals(
        const Array2D<short>&,
        Array2D<bool>&,
        Array2D<bool>&,
        const unsigned int,
        const float) const;

      void getSelectVals(
        const Array2D<float>&,
        Array2D<bool>&,
        Array2D<bool>&,
        const unsigned int,
        const float) const;

      void getSelectVals(
        const Array2D<double>&,
        Array2D<bool>&,
        Array2D<bool>&,
        const unsigned int,
        const float) const;

      void denoiseMorph1D(
        Array2D<short>&,
        const Array2D<short>&,
        Array2D<bool>&,
        Array2D<bool>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float) const;

      void denoiseMorph1D(
        Array2D<float>&,
        const Array2D<float>&,
        Array2D<bool>&,
        Array2D<bool>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float) const;

      void denoiseMorph1D(
        Array2D<double>&,
        const Array2D<double>&,
        Array2D<bool>&,
        Array2D<bool>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float) const;


      void denoiseMorph2D(
        Array2D<short>&,
        const Array2D<short>&,
        Array2D<bool>&,
        Array2D<bool>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float) const;

      void denoiseMorph2D(
        Array2D<float>&,
        const Array2D<float>&,
        Array2D<bool>&,
        Array2D<bool>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float) const;

      void denoiseMorph2D(
        Array2D<double>&,
        const Array2D<double>&,
        Array2D<bool>&,
        Array2D<bool>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float) const;


      void denoiseHough2D(
        Array2D<short>& waveLessCoherent,
        const Array2D<short>& fullEvent,
        Array2D<bool>& selectVals,
        Array2D<bool>& refinedSelectVals,
        Array2D<bool>& roi,
        const char filterName = 'g',
        const unsigned int grouping = 64,
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

      void denoiseHough2D(
        Array2D<float>& waveLessCoherent,
        const Array2D<float>& fullEvent,
        Array2D<bool>& selectVals,
        Array2D<bool>& refinedSelectVals,
        Array2D<bool>& roi,
        const char filterName = 'g',
        const unsigned int grouping = 64,
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

      void denoiseHough2D(
        Array2D<double>& waveLessCoherent,
        const Array2D<double>& fullEvent,
        Array2D<bool>& selectVals,
        Array2D<bool>& refinedSelectVals,
        Array2D<bool>& roi,
        const char filterName = 'g',
        const unsigned int grouping = 64,
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


    private:

      template <typename T>
      void getSelectVals(
        const Array2D<T>& morphedWaveforms,
        Array2D<bool>& selectVals,
        Array2D<bool>& roi,
        const unsigned int window,
        const float thresholdFactor) const;

      template <typename T>
      void denoiseMorph1D(
        Array2D<T>& waveLessCoherent,
        const Array2D<T>& fullEvent,
        Array2D<bool>& selectVals,
        Array2D<bool>& roi,
        const char filterName='g',
        const unsigned int grouping=64,
        const unsigned int structuringElement=5,
        const unsigned int window=0,
        const float thresholdFactor=2.5) const;

      template <typename T>
      void denoiseMorph2D(
        Array2D<T>& waveLessCoherent,
        const Array2D<T>& fullEvent,
        Array2D<bool>& selectVals,
        Array2D<bool>& roi,
        const char filterName='g',
        const unsigned int grouping=64,
        const unsigned int structuringElementx=5,
        const unsigned int structuringElementy=20,
        const unsigned int window=0,
        const float thresholdFactor=2.5) const;

      template <typename T>
      void denoiseHough2D(
        Array2D<T>& waveLessCoherent,
        const Array2D<T>& fullEvent,
        Array2D<bool>& selectVals,
        Array2D<bool>& refinedSelectVals,
        Array2D<bool>& roi,
        const char filterName = 'g',
        const unsigned int grouping = 64,
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
        const float eps = 0.00001
      ) const;
      /// Default destructor
      ~MorphologicalCNC(){}

  };
}

#endif
/** @} */ // end of doxygen group
