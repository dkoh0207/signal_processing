/**
 * \file MorphologicalFilter.h
 *
 * \ingroup algorithms
 * 
 * \brief Class def header for a class MorphologicalFilter
 *
 * @author koh0207
 */

/** \addtogroup algorithms

    @{*/
#ifndef __ALGORITHMS_MORPHOLOGICALFILTER_H__
#define __ALGORITHMS_MORPHOLOGICALFILTER_H__

#include <vector>
#include <algorithm>

namespace algorithms {

  /**
     \class MorphologicalFilter
     User defined class MorphologicalFilter ... these comments are used to generate
     doxygen documentation!
  */
  class MorphologicalFilter{
    
  public:
    
    /// Default constructor
    MorphologicalFilter(){}

    void FilterWaveform(RawDigitVector& waveform, size_t channel, size_t cnt, float pedestal) const {}

    void smoothInputWaveform(const RawDigitVector& inputWaveform, RawDigitVector& outputWaveform){}

    
    /// Default destructor
    ~MorphologicalFilter(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 

