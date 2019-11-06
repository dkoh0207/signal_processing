#ifndef __ALGORITHMS_NOISEREMOVAL_CXX__
#define __ALGORITHMS_NOISEREMOVAL_CXX__

#include "NoiseRemoval.h"

std::vector<std::vector<short>> algorithms::NoiseRemoval::removeCoherentNoise(
                         std::vector<std::vector<short>>& waveforms, 
                         const unsigned int grouping, 
                         const unsigned int nTicks) 
{
  std::vector<std::vector<short>> waveLessCoherent(
    waveforms.size(), std::vector<short>(waveforms.at(0).size(), 0.0));

  auto numChannels = (int) waveforms.size();

  for (unsigned int i=0; i<nTicks; ++i) {
    for (unsigned int j=0; j<numChannels/grouping; ++j) {
      int group_start = j * grouping;
      int group_end = (j+1) * grouping;
      // Compute median.
      std::vector<short> v(grouping, 0.0);
      for (auto c=group_start; c<group_end; ++c) {
        v.push_back(waveforms[c][i]);
      }
      std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
      short median = v[v.size() / 2];
      for (auto c=group_start; c<group_end; ++c) {
        waveLessCoherent[c][i] = waveforms[c][i] - median;
      }
    }
  }
  return waveLessCoherent;
}

void algorithms::NoiseRemoval::getWaveformParams(const std::vector<short>& waveform,
                                                 short& mean,
                                                 short& median,
                                                 short& mode,
                                                 float& skewness,
                                                 float& rms)
{
  /*
  Calculate waveform parameters for a given 1D waveform.

  INPUTS:
    - waveform: 1D RawDigit Waveform.
  */
  std::vector<short> localWaveform = waveform;
  std::sort(localWaveform.begin(),localWaveform.end(),[](const short& left, const short& right){return std::fabs(left) < std::fabs(right);});
  float realMean(float(std::accumulate(localWaveform.begin(),localWaveform.end(),0))/float(localWaveform.size()));
  median = localWaveform[localWaveform.size()/2];
  mean   = std::round(realMean);
  std::vector<float> adcLessPedVec;
  adcLessPedVec.resize(localWaveform.size());
  std::transform(localWaveform.begin(),localWaveform.end(),adcLessPedVec.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));
  rms      = std::sqrt(std::inner_product(adcLessPedVec.begin(), adcLessPedVec.end(), adcLessPedVec.begin(), 0.) / float(adcLessPedVec.size()));
  skewness = 3. * float(realMean - median) / rms;
  return;
}

void algorithms::NoiseRemoval::filterWaveforms(const std::vector<std::vector<short>>& waveforms,
                                               const unsigned int grouping,
                                               const unsigned int nTicks,
                                               std::vector<std::vector<short>>& noiseRemovedWfs,
                                               std::vector<short>& means,
                                               std::vector<short>& medians,
                                               std::vector<float>& rmss)
{
  auto numChannels = waveforms.size();
  std::vector<std::vector<short>> filteredWaveforms;
  filteredWaveforms.resize(numChannels);
  means.resize(numChannels);
  medians.resize(numChannels);
  rmss.resize(numChannels);

  for (unsigned int i=0; i<numChannels; ++i) {
    short mean = 0;
    short median = 0;
    short mode = 0;
    float rms = 0.0;
    float skewness = 0.0;
    filteredWaveforms[i].resize(waveforms.at(0).size());
    getWaveformParams(waveforms[i], mean, median, mode, skewness, rms);
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
