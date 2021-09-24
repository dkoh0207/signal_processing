/**
 * \file ROIFinder.h
 *
 * \ingroup sigproc_multithreading
 * 
 * \brief Class def header for a class ROIFinder
 *
 * @author koh0207
 */

/** \addtogroup sigproc_multithreading

    @{*/
#ifndef __SIGPROC_MULTITHREADING_ROIFINDER_H__
#define __SIGPROC_MULTITHREADING_ROIFINDER_H__

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <iterator>

#include "../sigproc_tools/Thresholding.h"
#include "../sigproc_tools/MiscUtils.h"
#include "../sigproc_tools/Morph2DFast.h"
#include "../sigproc_tools/LineDetection.h"
#include "../sigproc_tools/FrequencyFilters1D.h"
#include "../sigproc_tools/MorphologicalCNC.h"
#include "../sigproc_tools/EdgeDetection.h"
#include "../sigproc_tools/BilateralFilters.h"
#include "../sigproc_tools/Denoising.h"

#include "Morph2DFast.h"
#include "SigprocParallelDefs.h"

namespace sigproc_multithreading {

  /**
     \class ROIFinder
     User defined class ROIFinder ... these comments are used to generate
     doxygen documentation!
  */
  class ROIFinder1{
    
  public:

    ROIFinder1(){};
    
    /// Default constructor
    void applyChain(const Array2D<float>& input2D,
                    Array2D<bool>& output2D,
                    size_t FREQUENCY_THRESHOLD,
                    size_t FREQUENCY_FILTER_SMOOTHNESS_ORDER,
                    char MORPHOLOGICAL_FILTER_NAME,
                    const unsigned int CHANNEL_GROUPING,
                    const unsigned int CHANNEL_GROUPING_OFFSET,
                    const unsigned int STRUCTURING_ELEMENT_X,
                    const unsigned int STRUCTURING_ELEMENT_Y,
                    const float MORPHOLOGICAL_THRESHOLD_FACTOR,
                    
                    const size_t THETASTEPS,
                    const unsigned int HOUGH_THRESHOLD,
                    const unsigned int NMS_WINDOW_SIZE,
                    const unsigned int ANGLE_WINDOW,
                    
                    const unsigned int ADFILTER_SX,
                    const unsigned int ADFILTER_SY,
                    const float sigma_x,
                    const float sigma_y,
                    const float sigma_r,
                    const float lowThreshold,
                    const float highThreshold,
                    
                    const unsigned int BINARY_CLOSING_SX,
                    const unsigned int BINARY_CLOSING_SY) const;
    
    /// Default destructor
    ~ROIFinder1(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

