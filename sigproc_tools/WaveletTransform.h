/**
 * \file WaveletTransform.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class WaveletTransform
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_WAVELETTRANSFORM_H__
#define __SIGPROC_TOOLS_WAVELETTRANSFORM_H__

#include <vector>
#include <cmath>
#include "Wavelet.h"

namespace sigproc_tools {

  /**
     \class WaveletTransform
     User defined class WaveletTransform ... these comments are used to generate
     doxygen documentation!
  */
  class WaveletTransform{
    
  public:
    
    /// Default constructor
    WaveletTransform(){}

    // Adapted from Numerical Recipes in C++
    void wt1(
      std::vector<float>& wf,
      Daubechies4& wavelet,
      const int levels,
      const int isign) const;

    void wt2(
      const std::vector<std::vector<float>>& inputWaveform,
      std::vector<std::vector<float>>& transform,
      const size_t numChannels,
      const size_t nTicks,
      const int isign,
      Daubechies4& wavelet,
      const int levels) const;
    
    /// Default destructor
    ~WaveletTransform(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

