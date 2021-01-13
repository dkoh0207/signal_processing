/**
 * \file FindROI2D.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class FindROI2D
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_FINDROI2D_H__
#define __SIGPROC_TOOLS_FINDROI2D_H__


#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <iterator>

#include "AdaptiveWiener.h"
#include "Thresholding.h"
#include "MiscUtils.h"
#include "Morph2DFast.h"
#include "LineDetection.h"
#include "FrequencyFilters1D.h"
#include "MorphologicalCNC.h"
#include "Morph2DFast.h"

namespace sigproc_tools {

  /**
     \class FindROI2D
     User defined class FindROI2D ... these comments are used to generate
     doxygen documentation!
  */
  template <class T> using Array2D = std::vector<std::vector<T>>;

  class FindROI2D{
    
  public:
    
    /// Default constructor
    FindROI2D(){}

    void applyChainFilter(
        const Array2D<float>& waveform2D,
        Array2D<float>& fullEvent,
        Array2D<bool>& outputROI,
        Array2D<float>& waveLessCoherent,
        Array2D<float>& morphedWaveform2D,
        // Default Parameters
        size_t FREQUENCY_THRESHOLD,
        size_t FREQUENCY_FILTER_SMOOTHNESS_ORDER,
        size_t FREQUENCY_FILTER_MODE,

        char MORPHOLOGICAL_FILTER_NAME,
        const unsigned int CHANNEL_GROUPING,
        const unsigned int STRUCTURING_ELEMENT_X,
        const unsigned int STRUCTURING_ELEMENT_Y,
        const unsigned int ROI_EXPAND_WINDOW_SIZE,
        const float MORPHOLOGICAL_THRESHOLD_FACTOR,

        const size_t THETASTEPS,
        const unsigned int HOUGH_THRESHOLD,
        const unsigned int NMS_WINDOW_SIZE,
        const unsigned int ANGLE_WINDOW,

        // float NOISE_VARIANCE = 20.0;
        const unsigned int ADFILTER_SX,
        const unsigned int ADFILTER_SY,

        const float GLOBAL_THRESHOLDING_FACTOR) const;
    
    /// Default destructor
    ~FindROI2D(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

