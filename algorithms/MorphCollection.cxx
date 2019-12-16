#ifndef __ALGORITHMS_MORPHCOLLECTION_CXX__
#define __ALGORITHMS_MORPHCOLLECTION_CXX__

#include "MorphCollection.h"


std::vector<std::vector<float>> algorithms::MorphCollection::removeCoherentNoise(
                         std::vector<std::vector<float>>& filteredWaveforms,
                         const unsigned int grouping, 
                         const unsigned int nTicks,
                         const unsigned int structuringElement,
                         std::vector<std::vector<float>>& intrinsicRMS,
                         std::vector<std::vector<bool>>& selectVals,
                         std::vector<std::vector<float>>& correctedMedians,
                         const float thresholdFactor)
{
  // Initialize output array that will hold coherent noise removed waveforms.
  std::vector<std::vector<float>> waveLessCoherent(
    filteredWaveforms.size(), std::vector<float>(filteredWaveforms.at(0).size(), 0.0));

  // Get signal protection regions "selectVals"
  selectVals.resize(filteredWaveforms.size());
  for (auto& v : selectVals) {
    v.resize(filteredWaveforms.at(0).size());
  }
  for (size_t i=0; i<selectVals.size(); ++i) {
    for (size_t j=0; j<selectVals.at(0).size(); ++j) {
      selectVals[i][j] = true;
    }
  }
  getSelectVals(filteredWaveforms, nTicks, structuringElement, selectVals, thresholdFactor);

  auto numChannels = (int) filteredWaveforms.size();
  auto nGroups = numChannels / grouping;

  correctedMedians.resize(nGroups);
  for (auto& v : correctedMedians) {
    v.resize(nTicks);
  }

  for (size_t i=0; i<nTicks; ++i) {
    for (size_t j=0; j<numChannels/grouping; ++j) {
      int group_start = j * grouping;
      int group_end = (j+1) * grouping;
      // Compute median.
      std::vector<float> v;
      for (auto c=group_start; c<group_end; ++c) {
        if (selectVals[c][i]) {
          v.push_back(filteredWaveforms[c][i]);
        }
      }
      float median = 0.0;
      if (v.size() > 0) {
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
      }
      correctedMedians[j][i] = median;
      for (auto c=group_start; c<group_end; ++c) {
        if (selectVals[c][i]) {
          waveLessCoherent[c][i] = filteredWaveforms[c][i] - median;
        } else {
          waveLessCoherent[c][i] = filteredWaveforms[c][i];
        }
      }
    }
  }
  intrinsicRMS.resize(nGroups);
  for (auto& v : intrinsicRMS) {
    v.resize(nTicks);
  }
  float rms = 0.0;
  for (size_t i=0; i<nGroups; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      std::vector<float> v;
      for (auto k=i*grouping; k<(i+1)*grouping; ++k) {
        v.push_back(waveLessCoherent[k][j]);
      }
      rms = std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.) / float(v.size()));
      intrinsicRMS[i][j] = rms;
    }
  }
  return waveLessCoherent;
}

void algorithms::MorphCollection::getSelectVals(const std::vector<std::vector<float>>& waveforms,
                                                const unsigned int nTicks,
                                                const unsigned int structuringElement,
                                                std::vector<std::vector<bool>>& selectVals,
                                                const float thresholdFactor)
{
  WaveformUtils wUtils;
  auto numChannels = waveforms.size();
  selectVals.resize(numChannels);

  for (size_t i=0; i<numChannels; ++i) {
    std::vector<float> dilation;
    wUtils.getDilation(waveforms[i], structuringElement, dilation);
    float dilationMed = 0.0;
    std::vector<float> localDilation;
    localDilation.resize(nTicks);
    for (size_t j=0; j<nTicks; ++j) {
      localDilation[j] = dilation[j];
    }
    if (localDilation.size() % 2 == 0) {
      const auto m1 = localDilation.begin() + localDilation.size() / 2 - 1;
      const auto m2 = localDilation.begin() + localDilation.size() / 2;
      std::nth_element(localDilation.begin(), m1, localDilation.end());
      const auto e1 = *m1;
      std::nth_element(localDilation.begin(), m2, localDilation.end());
      const auto e2 = *m2;
      dilationMed = (e1 + e2) / 2.0;
    } else {
      const auto m = localDilation.begin() + localDilation.size() / 2;
      std::nth_element(localDilation.begin(), m, localDilation.end());
      dilationMed = *m;
    }
    std::vector<float> dilationBase;
    dilationBase.resize(dilation.size());
    for (size_t i=0; i<dilationBase.size(); ++i) {
      dilationBase[i] = dilation[i] - dilationMed;
    }
    float dilationRMS;
    dilationRMS = std::sqrt(std::inner_product(dilationBase.begin(), 
      dilationBase.end(), dilationBase.begin(), 0.) / float(dilationBase.size()));
    float threshold;
    threshold = dilationRMS * thresholdFactor;
    for (size_t j=0; j<nTicks; ++j) {
      bool sVal = true;
      if (waveforms[i][j] >= threshold) {
        sVal = false;
      }
      selectVals[i][j] = sVal;
    }
  }
  return;
}


void algorithms::MorphCollection::filterWaveforms(const std::vector<std::vector<short>>& waveforms,
                                               const unsigned int grouping,
                                               const unsigned int nTicks,
                                               const unsigned int structuringElement,
                                               std::vector<std::vector<float>>& noiseRemovedWfs,
                                               std::vector<float>& means,
                                               std::vector<float>& medians,
                                               std::vector<float>& totalRMS,
                                               std::vector<std::vector<float>>& intrinsicRMS,
                                               std::vector<float>& cleanRMS)
{
  auto numChannels = waveforms.size();
  auto nGroups = (int) numChannels / grouping;
  std::vector<std::vector<float>> filteredWaveforms;
  std::vector<std::vector<float>> correctedMedians;
  std::vector<std::vector<bool>> selectVals;
  filteredWaveforms.resize(numChannels);
  means.resize(numChannels);
  medians.resize(numChannels);
  totalRMS.resize(numChannels);
  cleanRMS.resize(numChannels);
  intrinsicRMS.resize(nGroups);
  correctedMedians.resize(nGroups);

  for (size_t i=0; i<nGroups; ++i) {
    intrinsicRMS[i].resize(nTicks);
  }
  for (size_t i=0; i<nGroups; ++i) {
    correctedMedians[i].resize(nTicks);
  }

  float mean = 0;
  float median = 0;
  float rms = 0.0;
  float skewness = 0.0;

  for (size_t i=0; i<numChannels; ++i) {
    filteredWaveforms[i].resize(waveforms.at(0).size());
    WaveformUtils wUtils;
    wUtils.getWaveformParams(waveforms[i], mean, median, skewness, rms);
    // Subtract Pedestals (Median of given waveform, one channel)
    std::transform(waveforms[i].begin(),waveforms[i].end(),filteredWaveforms[i].begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));
    wUtils.getWaveformParams(filteredWaveforms[i], mean, median, skewness, rms);
    means[i] = mean;
    medians[i] = median;
    totalRMS[i] = rms;
  }

  noiseRemovedWfs = removeCoherentNoise(filteredWaveforms, grouping, nTicks, 
                                        structuringElement, intrinsicRMS, selectVals, correctedMedians);

  for (size_t i=0; i<numChannels; ++i) {
    WaveformUtils wUtils;
    wUtils.getWaveformParams(noiseRemovedWfs[i], mean, median, skewness, rms);
    cleanRMS[i] = rms;
  }

  return;
}

#endif
