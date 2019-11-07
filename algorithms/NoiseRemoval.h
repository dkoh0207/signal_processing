/**
 * \file NoiseRemoval.h
 *
 * \ingroup algorithms
 * 
 * \brief Class def header for a class NoiseRemoval
 *
 * @author koh0207
 */

/** \addtogroup algorithms

    @{*/
#ifndef __ALGORITHMS_NOISEREMOVAL_H__
#define __ALGORITHMS_NOISEREMOVAL_H__

#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include "WaveformUtils.h"

namespace algorithms {

  /**
     \class NoiseRemoval
     User defined class NoiseRemoval ... these comments are used to generate
     doxygen documentation!
  */
  class NoiseRemoval{
    
  public:
    
    /// Default constructor
    NoiseRemoval(){}

    std::vector<std::vector<float>> removeCoherentNoise(
                          std::vector<std::vector<float>> &waveforms, 
                          const unsigned int grouping, 
                          const unsigned int nTicks);

    void filterWaveforms(const std::vector<std::vector<short>>& waveforms,
                                                  const unsigned int grouping,
                                                  const unsigned int nTicks,
                                                  std::vector<std::vector<float>>& noiseRemovedWfs,
                                                  std::vector<float>& means,
                                                  std::vector<float>& medians,
                                                  std::vector<float>& rmss);
    
    /// Default destructor
    ~NoiseRemoval(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

