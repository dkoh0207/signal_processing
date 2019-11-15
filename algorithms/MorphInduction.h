/**
 * \file MorphInduction.h
 *
 * \ingroup algorithms
 * 
 * \brief Class def header for a class MorphInduction
 *
 * @author koh0207
 */

/** \addtogroup algorithms

    @{*/
#ifndef __ALGORITHMS_MORPHINDUCTION_H__
#define __ALGORITHMS_MORPHINDUCTION_H__

#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "WaveformUtils.h"

namespace algorithms {

  /**
     \class MorphInduction
     User defined class MorphInduction ... these comments are used to generate
     doxygen documentation!
  */
  class MorphInduction{
    
  public:
    
    /// Default constructor
    MorphInduction(){}

    std::vector<std::vector<float>> removeCoherentNoise(
                         std::vector<std::vector<float>>& filteredWaveforms, 
                         const unsigned int grouping, 
                         const unsigned int nTicks,
                         const unsigned int structuringElement,
                         const unsigned int window,
                         std::vector<std::vector<float>>& intrinsicRMS,
                         std::vector<std::vector<bool>>& selectVals);

    std::vector<std::vector<float>> removeCoherentNoise2D(
                         std::vector<std::vector<float>>& filteredWaveforms, 
                         const unsigned int grouping, 
                         const unsigned int nTicks,
                         const unsigned int structuringElementx,
                         const unsigned int structuringElementy,
                         const unsigned int window,
                         std::vector<std::vector<float>>& intrinsicRMS,
                         std::vector<std::vector<bool>>& selectVals);

    void filterWaveforms(const std::vector<std::vector<short>>& waveforms,
                                               const unsigned int grouping,
                                               const unsigned int nTicks,
                                               const unsigned int structuringElement,
                                               const unsigned int window,
                                               std::vector<std::vector<float>>& noiseRemovedWfs,
                                               std::vector<float>& means,
                                               std::vector<float>& medians,
                                               std::vector<float>& rmss,
                                               std::vector<std::vector<float>>& intrinsicRMS);

    void getSelectVals(const std::vector<std::vector<float>>& waveforms,
                                                const unsigned int grouping,
                                                const unsigned int nTicks,
                                                const unsigned int structuringElement,
                                                std::vector<std::vector<bool>>& selectVals,
                                                const unsigned int window);

    void getSelectVals2D(const std::vector<std::vector<float>>& waveforms,
                                                const unsigned int grouping,
                                                const unsigned int nTicks,
                                                const unsigned int structuringElementx,
                                                const unsigned int structuringElementy,
                                                std::vector<std::vector<bool>>& selectVals,
                                                const unsigned int window);
    
    /// Default destructor
    ~MorphInduction(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

