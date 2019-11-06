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

    std::vector<std::vector<short>> removeCoherentNoise(
                          std::vector<std::vector<short>> &waveforms, 
                          const unsigned int grouping, 
                          const unsigned int nTicks);

    void getWaveformParams(const std::vector<short>& waveform,
                                                 short& mean,
                                                 short& median,
                                                 short& mode,
                                                 float& skewness,
                                                 float& rms);

    void filterWaveforms(const std::vector<std::vector<short>>& waveforms,
                                                  const unsigned int grouping,
                                                  const unsigned int nTicks,
                                                  std::vector<std::vector<short>>& noiseRemovedWfs,
                                                  std::vector<short>& means,
                                                  std::vector<short>& medians,
                                                  std::vector<float>& rmss);
    
    /// Default destructor
    ~NoiseRemoval(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

