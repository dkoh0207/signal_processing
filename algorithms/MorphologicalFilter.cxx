#ifndef __ALGORITHMS_MORPHOLOGICALFILTER_CXX__
#define __ALGORITHMS_MORPHOLOGICALFILTER_CXX__

#include "MorphologicalFilter.h"

void FilterWaveform(RawDigitVector& waveform, size_t channel, size_t cnt, float pedestal) const
    {
      std::vector<short> erosionVec;
      std::vector<short> dilationVec;
      std::vector<short> averageVec;
      std::vector<short> differenceVec;
    }

#endif
