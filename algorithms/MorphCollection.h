/**
 * \file MorphCollection.h
 *
 * \ingroup algorithms
 * 
 * \brief Class def header for a class MorphCollection
 *
 * @author koh0207
 */

/** \addtogroup algorithms

    @{*/
#ifndef __ALGORITHMS_MORPHCOLLECTION_H__
#define __ALGORITHMS_MORPHCOLLECTION_H__

#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include <string>
#include "WaveformUtils.h"

namespace algorithms {

  /**
     \class MorphCollection
     User defined class MorphCollection ... these comments are used to generate
     doxygen documentation!
  */
  class MorphCollection{
    
  public:
    
    /// Default constructor
    MorphCollection(){}
    
    std::vector<std::vector<float>> removeCoherentNoise(
                         std::vector<std::vector<float>>& filteredWaveforms, 
                         const unsigned int grouping, 
                         const unsigned int nTicks,
                         const unsigned int structuringElement,
                         std::vector<std::vector<float>>& intrinsicRMS,
                         std::vector<std::vector<bool>>& selectVals,
                         std::vector<std::vector<float>>& correctedMedians,
                         const float thresholdFactor=2.5);

    void filterWaveforms(const std::vector<std::vector<short>>& waveforms,
                                               const unsigned int grouping,
                                               const unsigned int nTicks,
                                               const unsigned int structuringElement,
                                               std::vector<std::vector<float>>& noiseRemovedWfs,
                                               std::vector<float>& means,
                                               std::vector<float>& medians,
                                               std::vector<float>& totalRMS,
                                               std::vector<std::vector<float>>& intrinsicRMS,
                                               std::vector<float>& cleanRMS);

    void getSelectVals(const std::vector<std::vector<float>>& waveforms,
                                                const unsigned int grouping,
                                                const unsigned int nTicks,
                                                const unsigned int structuringElement,
                                                std::vector<std::vector<bool>>& selectVals,
                                                const float thresholdFactor=2.5);

    void getSelectVals2D(const std::vector<std::vector<float> >& waveforms,
                                                const unsigned int grouping,
                                                const unsigned int nTicks,
                                                const unsigned int structuringElementx,
                                                const unsigned int structuringElementy,
                                                std::vector<std::vector<bool> >& selectVals,
                                                const unsigned int window,
                                                const float thresholdFactor=2.5);
  
    /// Default destructor
    ~MorphCollection(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

