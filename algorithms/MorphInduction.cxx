#ifndef __ALGORITHMS_MORPHINDUCTION_CXX__
#define __ALGORITHMS_MORPHINDUCTION_CXX__

#include "MorphInduction.h"

std::vector<std::vector<float>> algorithms::MorphInduction::removeCoherentNoise(
                         std::vector<std::vector<float>>& filteredWaveforms,
                         const unsigned int grouping, 
                         const unsigned int nTicks,
                         const unsigned int structuringElement,
                         const unsigned int window,
                         std::vector<std::vector<float>>& intrinsicRMS,
                         std::vector<std::vector<bool>>& selectVals)
{
  std::vector<std::vector<float>> waveLessCoherent(
    filteredWaveforms.size(), std::vector<float>(filteredWaveforms.at(0).size(), 0.0));

  selectVals.resize(filteredWaveforms.size());
  for (auto& v : selectVals) {
    v.resize(filteredWaveforms.at(0).size());
  }
  for (auto i=0; i<selectVals.size(); ++i) {
    for (auto j=0; j<selectVals.at(0).size(); ++j) {
      selectVals[i][j] = true;
    }
  }

  getSelectVals(filteredWaveforms, grouping, nTicks, structuringElement, selectVals, window);

  auto numChannels = (int) filteredWaveforms.size();

  for (auto i=0; i<nTicks; ++i) {
    for (auto j=0; j<numChannels/grouping; ++j) {
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
      for (auto c=group_start; c<group_end; ++c) {
        if (selectVals[c][i]) {
          waveLessCoherent[c][i] = filteredWaveforms[c][i] - median;
        } else {
          waveLessCoherent[c][i] = filteredWaveforms[c][i];
        }
      }
    }
  }
  auto nGroups = numChannels / grouping;
  intrinsicRMS.resize(nGroups);
  for (auto& v : intrinsicRMS) {
    v.resize(nTicks);
  }
  float rms = 0.0;
  for (auto i=0; i<nGroups; ++i) {
    for (auto j=0; j<nTicks; ++j) {
      std::vector<float> v(nTicks, 0.0);
      for (auto k=i*grouping; k<(i+1)*grouping; ++k) {
        v[k] = waveLessCoherent[k][j];
      }
      rms = std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.) / float(v.size()));
      intrinsicRMS[i][j] = rms;
    }
  }
  return waveLessCoherent;
}

void algorithms::MorphInduction::getSelectVals(const std::vector<std::vector<float>>& waveforms,
                                                const unsigned int grouping,
                                                const unsigned int nTicks,
                                                const unsigned int structuringElement,
                                                std::vector<std::vector<bool>>& selectVals,
                                                const unsigned int window)
{
  WaveformUtils wUtils;
  auto numChannels = waveforms.size();

  for (auto i=0; i<numChannels; ++i) {
    std::vector<float> dilation;
    std::vector<float> erosion;
    std::vector<float> average;
    std::vector<float> gradient;
    wUtils.getErosionDilationAverageDifference(waveforms[i], structuringElement, erosion, dilation, average, gradient);
    float gradientMed = 0.0;
    if (gradient.size() % 2 == 0) {
      const auto m1 = gradient.begin() + gradient.size() / 2 - 1;
      const auto m2 = gradient.begin() + gradient.size() / 2;
      std::nth_element(gradient.begin(), m1, gradient.end());
      const auto e1 = *m1;
      std::nth_element(gradient.begin(), m2, gradient.end());
      const auto e2 = *m2;
      gradientMed = (e1 + e2) / 2.0;
    } else {
      const auto m = gradient.begin() + gradient.size() / 2;
      std::nth_element(gradient.begin(), m, gradient.end());
      gradientMed = *m;
    }
    std::vector<float> gradientBase;
    gradientBase.resize(gradient.size());
    for (auto i=0; i<gradientBase.size(); ++i) {
      gradientBase[i] = gradient[i] - gradientMed;
    }
    float gradientRMS;
    gradientRMS = std::sqrt(std::inner_product(gradientBase.begin(), 
      gradientBase.end(), gradientBase.begin(), 0.) / float(gradientBase.size()));
    float threshold;
    threshold = gradientMed + gradientRMS * 2.5;
    // Should be better ways to do this.
    for (auto j=0; j<nTicks; ++j) {
      if (waveforms[i][j] >= threshold) {
        selectVals[i][j] = false;
        if (j > window-1) {
          for (auto k=0; k<window; ++k) {
            if (!selectVals[i][j-k]) continue;
            selectVals[i][j-k] = false;
          }
        }
        if (j < nTicks-window) {
          for (auto k=0; k<window; ++k) {
            if (!selectVals[i][j+k]) continue;
            selectVals[i][j+k] = false;
          }
        }
      }
    }
  }
  return;
}


void algorithms::MorphInduction::filterWaveforms(const std::vector<std::vector<short>>& waveforms,
                                               const unsigned int grouping,
                                               const unsigned int nTicks,
                                               const unsigned int structuringElement,
                                               const unsigned int window,
                                               std::vector<std::vector<float>>& noiseRemovedWfs,
                                               std::vector<float>& means,
                                               std::vector<float>& medians,
                                               std::vector<float>& rmss,
                                               std::vector<std::vector<float>>& intrinsicRMS)
{
  auto numChannels = waveforms.size();
  auto nGroups = numChannels / grouping;
  std::vector<std::vector<bool>> selectVals;
  std::vector<std::vector<float>> filteredWaveforms;
  filteredWaveforms.resize(numChannels);
  means.resize(numChannels);
  medians.resize(numChannels);
  rmss.resize(numChannels);

  for (auto i=0; i<numChannels; ++i) {
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
  noiseRemovedWfs = removeCoherentNoise(filteredWaveforms, grouping, nTicks, structuringElement, window, intrinsicRMS, selectVals);
  return;
}



#endif
