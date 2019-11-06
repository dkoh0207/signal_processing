#ifndef __ALGORITHMS_NOISEREMOVAL_CXX__
#define __ALGORITHMS_NOISEREMOVAL_CXX__

#include "NoiseRemoval.h"

std::vector<std::vector<double>> algorithms::NoiseRemoval::removeCoherentNoise(std::vector<std::vector<double>> &waveforms, 
                         const int grouping, 
                         const int nTicks) 
{
  std::vector<std::vector<double>> waveLessCoherent(
    waveforms.size(), std::vector<double>(waveforms.at(0).size(), 0.0));

  auto numChannels = (int) waveforms.size();

  for (auto i=0; i<nTicks; ++i) {
    for (auto j=0; j<numChannels/grouping; ++j) {
      int group_start = j * grouping;
      int group_end = (j+1) * grouping;
      // Compute median.
      std::vector<double> v(grouping, 0.0);
      for (auto c=group_start; c<group_end; ++c) {
        v.push_back(waveforms[c][i]);
      }
      std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
      double median = v[v.size() / 2];
      for (auto c=group_start; c<group_end; ++c) {
        waveLessCoherent[c][i] = waveforms[c][i] - median;
      }
    }
  }
  return waveLessCoherent;
}

#endif
