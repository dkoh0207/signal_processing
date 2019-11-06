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

    std::vector<std::vector<double>> removeCoherentNoise(
                          std::vector<std::vector<double>> &waveforms, 
                          const int grouping, 
                          const int nTicks);
    
    /// Default destructor
    ~NoiseRemoval(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

