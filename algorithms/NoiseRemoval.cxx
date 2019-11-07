#ifndef __ALGORITHMS_NOISEREMOVAL_CXX__
#define __ALGORITHMS_NOISEREMOVAL_CXX__

#include "NoiseRemoval.h"

std::vector<std::vector<float>> algorithms::NoiseRemoval::removeCoherentNoise(
                         std::vector<std::vector<float>>& filteredWaveforms, 
                         const unsigned int grouping, 
                         const unsigned int nTicks) 
{
  std::vector<std::vector<float>> waveLessCoherent(
    filteredWaveforms.size(), std::vector<float>(filteredWaveforms.at(0).size(), 0.0));

  auto numChannels = (int) filteredWaveforms.size();

  for (unsigned int i=0; i<nTicks; ++i) {
    for (unsigned int j=0; j<numChannels/grouping; ++j) {
      int group_start = j * grouping;
      int group_end = (j+1) * grouping;
      // Compute median.
      std::vector<float> v;
      v.resize(grouping);
      short counter = 0;
      for (auto c=group_start; c<group_end; ++c) {
        v[counter] = filteredWaveforms[c][i];
        ++counter;
      }
      float median = 0.0;
      if (v.size() % 2 == 0) {
        const auto m1 = v.begin() + v.size() / 2 - 1;
        const auto m2 = v.begin() + v.size() / 2;
        std::nth_element(v.begin(), m1, v.end());
        const auto e1 = *m1;
        std::nth_element(v.begin(), m2, v.end());
        const auto e2 = *m2;
        median = (e1 + e2) / 2.0;
      } else {
        const auto m = v.begin() + v.size() / 2;
        std::nth_element(v.begin(), m, v.end());
        median = *m;
      }
      for (auto c=group_start; c<group_end; ++c) {
        waveLessCoherent[c][i] = filteredWaveforms[c][i] - median;
      }
    }
  }
  return waveLessCoherent;
}

void algorithms::NoiseRemoval::filterWaveforms(const std::vector<std::vector<short>>& waveforms,
                                               const unsigned int grouping,
                                               const unsigned int nTicks,
                                               std::vector<std::vector<float>>& noiseRemovedWfs,
                                               std::vector<float>& means,
                                               std::vector<float>& medians,
                                               std::vector<float>& rmss)
{
  auto numChannels = waveforms.size();
  std::vector<std::vector<float>> filteredWaveforms;
  filteredWaveforms.resize(numChannels);
  means.resize(numChannels);
  medians.resize(numChannels);
  rmss.resize(numChannels);

  for (unsigned int i=0; i<numChannels; ++i) {
    float mean = 0;
    float median = 0;
    float mode = 0;
    float rms = 0.0;
    float skewness = 0.0;
    filteredWaveforms[i].resize(waveforms.at(0).size());
    WaveformUtils wUtils;
    wUtils.getWaveformParams(waveforms[i], mean, median, mode, skewness, rms);
    // Subtract Pedestals (Median of given waveform, one channel)
    std::transform(waveforms[i].begin(),waveforms[i].end(),filteredWaveforms[i].begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));
    means[i] = mean;
    medians[i] = median;
    rmss[i] = rms;
  }
  noiseRemovedWfs = removeCoherentNoise(filteredWaveforms, grouping, nTicks);
  return;
}



#endif
