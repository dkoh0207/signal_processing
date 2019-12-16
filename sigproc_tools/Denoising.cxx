#ifndef __SIGPROC_TOOLS_DENOISING_CXX__
#define __SIGPROC_TOOLS_DENOISING_CXX__

#include "Denoising.h"


void sigproc_tools::Denoising::getSelectVals(
  const std::vector<std::vector<short>>& waveforms,
  const std::vector<std::vector<short>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  const unsigned int window,
  const float thresholdFactor)
{
  getSelectVals<short>(waveforms, morphedWaveforms,
    selectVals, window, thresholdFactor);
}

void sigproc_tools::Denoising::getSelectVals(
  const std::vector<std::vector<float>>& waveforms,
  const std::vector<std::vector<float>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  const unsigned int window,
  const float thresholdFactor)
{
  getSelectVals<float>(waveforms, morphedWaveforms,
    selectVals, window, thresholdFactor);
}

void sigproc_tools::Denoising::getSelectVals(
  const std::vector<std::vector<double>>& waveforms,
  const std::vector<std::vector<double>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  const unsigned int window,
  const float thresholdFactor)
{
  getSelectVals<double>(waveforms, morphedWaveforms,
    selectVals, window, thresholdFactor);
}

template <typename T>
void sigproc_tools::Denoising::getSelectVals(
  const std::vector<std::vector<T>>& waveforms,
  const std::vector<std::vector<T>>& morphedWaveforms,
  std::vector<std::vector<bool>>& selectVals,
  const unsigned int window,
  const float thresholdFactor)
{
  auto numChannels = waveforms.size();
  auto nTicks = waveforms.at(0).size();

  for (auto i=0; i<numChannels; ++i) {
    T median = 0.0;
    std::vector<T> localVec = morphedWaveforms[i];
    if (localVec.size() % 2 == 0) {
      const auto m1 = localVec.begin() + localVec.size() / 2 - 1;
      const auto m2 = localVec.begin() + localVec.size() / 2;
      std::nth_element(localVec.begin(), m1, localVec.end());
      const auto e1 = *m1;
      std::nth_element(localVec.begin(), m2, localVec.end());
      const auto e2 = *m2;
      median = (e1 + e2) / 2.0;
    } else {
      const auto m = localVec.begin() + localVec.size() / 2;
      std::nth_element(localVec.begin(), m, localVec.end());
      median = *m;
    }
    std::vector<T> baseVec;
    baseVec.resize(localVec.size());
    for (auto j=0; j<baseVec.size(); ++j) {
      baseVec[j] = morphedWaveforms[i][j] - median;
    }
    float rms;
    rms = std::sqrt(std::inner_product(baseVec.begin(), 
      baseVec.end(), baseVec.begin(), 0.) / float(baseVec.size()));
    float threshold;
    threshold = rms * thresholdFactor;

    for (int j=0; j<nTicks; ++j) {
      if (morphedWaveforms[i][j] > threshold) {
        // Check Bounds
        int lb = j - (int) window;
        int ub = j + (int) window + 1;
        int lowerBound = std::max(lb, 0);
        int upperBound = std::min(ub, (int) nTicks);
        for (auto k=lowerBound; k<upperBound; ++k) {
          selectVals[i][k] = true;
        }
      } else {
        selectVals[i][j] = false;
      }
    }
  }
  return;
}


void removeCoherentNoise1D(
  std::vector<std::vector<short>>& waveLessCoherent,
  const std::vector<std::vector<short>>& filteredWaveforms,
  std::vector<std::vector<short>>& morphedWaveforms,
  std::vector<std::vector<short>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<short>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise1D(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, correctedMedians,
    filterName, grouping, structuringElement, window, thresholdFactor);
  return;
}

void removeCoherentNoise1D(
  std::vector<std::vector<float>>& waveLessCoherent,
  const std::vector<std::vector<float>>& filteredWaveforms,
  std::vector<std::vector<float>>& morphedWaveforms,
  std::vector<std::vector<float>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<float>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise1D(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, correctedMedians,
    filterName, grouping, structuringElement, window, thresholdFactor);
  return;
}

void removeCoherentNoise1D(
  std::vector<std::vector<double>>& waveLessCoherent,
  const std::vector<std::vector<double>>& filteredWaveforms,
  std::vector<std::vector<double>>& morphedWaveforms,
  std::vector<std::vector<double>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<double>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise1D(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, correctedMedians,
    filterName, grouping, structuringElement, window, thresholdFactor);
  return;
}

template <typename T>
void removeCoherentNoise1D(
  std::vector<std::vector<T>>& waveLessCoherent,
  const std::vector<std::vector<T>>& filteredWaveforms,
  std::vector<std::vector<T>>& morphedWaveforms,
  std::vector<std::vector<T>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<T>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  auto numChannels = filteredWaveforms.size();
  auto nTicks = filteredWaveforms.at(0).size();
  auto nGroups = numChannels / grouping;

  // Coherent noise subtracted denoised waveforms
  waveLessCoherent.resize(filteredWaveforms.size());
  for (auto& v : waveLessCoherent) {
    v.resize(filteredWaveforms.at(0).size());
  }

  // Waveform with morphological filter applied
  morphedWaveforms.resize(filteredWaveforms.size());
  for (auto& v : morphedWaveforms) {
    v.resize(filteredWaveforms.at(0).size());
  }

  // Regions to protect waveform from coherent noise subtraction.
  selectVals.resize(filteredWaveforms.size());
  for (auto& v : selectVals) {
    v.resize(filteredWaveforms.at(0).size());
  }

  correctedMedians.resize(nGroups);
  for (auto& v : correctedMedians) {
    v.resize(nTicks);
  }

  intrinsicRMS.resize(nGroups);
  for (auto& v : intrinsicRMS) {
    v.resize(nTicks);
  }

  sigproc_tools::Morph1D denoiser;

  switch (filterName) {
    case 'd':
      for (auto i=0; i<numChannels; ++i) {
        denoiser.getDilation(filteredWaveforms[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    case 'e':
      for (auto i=0; i<numChannels; ++i) {
        denoiser.getErosion(filteredWaveforms[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    case 'a':
      for (auto i=0; i<numChannels; ++i) {
        denoiser.getAverage(filteredWaveforms[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    case 'g':
      for (auto i=0; i<numChannels; ++i) {
        denoiser.getGradient(filteredWaveforms[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
    default:
      for (auto i=0; i<numChannels; ++i) {
        denoiser.getDilation(filteredWaveforms[i],
          structuringElement, morphedWaveforms[i]);
      };
      break;
  }

  getSelectVals(filteredWaveforms, morphedWaveforms, 
    selectVals, window, thresholdFactor);

  for (auto i=0; i<nTicks; ++i) {
    for (auto j=0; j<nGroups; ++j) {
      int group_start = j * grouping;
      int group_end = (j+1) * grouping;
      // Compute median.
      std::vector<T> v;
      for (auto c=group_start; c<group_end; ++c) {
        if (!selectVals[c][i]) {
          v.push_back(filteredWaveforms[c][i]);
        }
      }
      T median = (T) 0;
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
      for (auto k=group_start; k<group_end; ++k) {
        if (!selectVals[k][i]) {
          waveLessCoherent[k][i] = filteredWaveforms[k][i] - median;
        } else {
          waveLessCoherent[k][i] = filteredWaveforms[k][i];
        }
      }
    }
  }

  T rms = (T) 0;
  for (auto i=0; i<nGroups; ++i) {
    for (auto j=0; j<nTicks; ++j) {
      std::vector<T> v;
      for (auto k=i*grouping; k<(i+1)*grouping; ++k) {
        v.push_back(waveLessCoherent[k][j]);
      }
      rms = std::sqrt(
        std::inner_product(
          v.begin(), v.end(), v.begin(), 0.) / T(v.size()));
      intrinsicRMS[i][j] = (T) rms;
    }
  }
  return;
}


void removeCoherentNoise2D(
  std::vector<std::vector<short>>& waveLessCoherent,
  const std::vector<std::vector<short>>& filteredWaveforms,
  std::vector<std::vector<short>>& morphedWaveforms,
  std::vector<std::vector<short>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<short>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise2D(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, correctedMedians,
    filterName, grouping, structuringElement, window, thresholdFactor);
  return;
}

void removeCoherentNoise2D(
  std::vector<std::vector<float>>& waveLessCoherent,
  const std::vector<std::vector<float>>& filteredWaveforms,
  std::vector<std::vector<float>>& morphedWaveforms,
  std::vector<std::vector<float>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<float>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise2D(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, correctedMedians,
    filterName, grouping, structuringElement, window, thresholdFactor);
  return;
}

void removeCoherentNoise2D(
  std::vector<std::vector<double>>& waveLessCoherent,
  const std::vector<std::vector<double>>& filteredWaveforms,
  std::vector<std::vector<double>>& morphedWaveforms,
  std::vector<std::vector<double>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<double>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElement,
  const unsigned int window,
  const float thresholdFactor)
{
  removeCoherentNoise2D(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, correctedMedians,
    filterName, grouping, structuringElement, window, thresholdFactor);
  return;
}

template <typename T>
void removeCoherentNoise2D(
  std::vector<std::vector<T>>& waveLessCoherent,
  const std::vector<std::vector<T>>& filteredWaveforms,
  std::vector<std::vector<T>>& morphedWaveforms,
  std::vector<std::vector<T>>& intrinsicRMS,
  std::vector<std::vector<bool>>& selectVals,
  std::vector<std::vector<T>>& correctedMedians,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window,
  const float thresholdFactor)
{
  auto numChannels = filteredWaveforms.size();
  auto nTicks = filteredWaveforms.at(0).size();
  auto nGroups = numChannels / grouping;

  // Coherent noise subtracted denoised waveforms
  waveLessCoherent.resize(filteredWaveforms.size());
  for (auto& v : waveLessCoherent) {
    v.resize(filteredWaveforms.at(0).size());
  }

  // Regions to protect waveform from coherent noise subtraction.
  selectVals.resize(filteredWaveforms.size());
  for (auto& v : selectVals) {
    v.resize(filteredWaveforms.at(0).size());
  }

  correctedMedians.resize(nGroups);
  for (auto& v : correctedMedians) {
    v.resize(nTicks);
  }

  intrinsicRMS.resize(nGroups);
  for (auto& v : intrinsicRMS) {
    v.resize(nTicks);
  }

  sigproc_tools::Morph2D denoiser;

  std::vector<std::vector<T>> dilation;
  std::vector<std::vector<T>> erosion;
  std::vector<std::vector<T>> average;
  std::vector<std::vector<T>> gradient;

  denoiser.getFilter2D(filteredWaveforms, structuringElementx,
    structuringElementy, dilation, erosion, average, gradient);

  switch (filterName) {
    case 'd':
      getSelectVals(filteredWaveforms, dilation, 
        selectVals, window, thresholdFactor);
      morphedWaveforms = dilation;
      break;
    case 'e':
      getSelectVals(filteredWaveforms, erosion, 
        selectVals, window, thresholdFactor);
      morphedWaveforms = erosion;
      break;
    case 'a':
      getSelectVals(filteredWaveforms, average, 
        selectVals, window, thresholdFactor);
      morphedWaveforms = average;
      break;
    case 'g':
      getSelectVals(filteredWaveforms, gradient, 
        selectVals, window, thresholdFactor);
      morphedWaveforms = gradient;
      break;
    default:
      getSelectVals(filteredWaveforms, gradient, 
        selectVals, window, thresholdFactor);
      morphedWaveforms = gradient;
      break;
  }

  for (auto i=0; i<nTicks; ++i) {
    for (auto j=0; j<nGroups; ++j) {
      int group_start = j * grouping;
      int group_end = (j+1) * grouping;
      // Compute median.
      std::vector<T> v;
      for (auto c=group_start; c<group_end; ++c) {
        if (!selectVals[c][i]) {
          v.push_back(filteredWaveforms[c][i]);
        }
      }
      T median = (T) 0;
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
      for (auto k=group_start; k<group_end; ++k) {
        if (!selectVals[k][i]) {
          waveLessCoherent[k][i] = filteredWaveforms[k][i] - median;
        } else {
          waveLessCoherent[k][i] = filteredWaveforms[k][i];
        }
      }
    }
  }

  float rms = 0.0;
  for (auto i=0; i<nGroups; ++i) {
    for (auto j=0; j<nTicks; ++j) {
      std::vector<T> v;
      for (auto k=i*grouping; k<(i+1)*grouping; ++k) {
        v.push_back(waveLessCoherent[k][j]);
      }
      rms = std::sqrt(
        std::inner_product(
          v.begin(), v.end(), v.begin(), 0.) / T(v.size()));
      intrinsicRMS[i][j] = (T) rms;
    }
  }
  return;
}


#endif
